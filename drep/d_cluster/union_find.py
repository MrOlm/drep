"""
Streaming, low-memory primary clustering via union-find (disjoint set).

Primary clustering in dRep is single-linkage hierarchical clustering at a fixed
distance cutoff. That is mathematically identical to finding the connected
components of the graph whose nodes are genomes and whose edges are the pairs
with ``distance <= cutoff``.

The classic dRep path materializes every pairwise MASH comparison into a long
``Mdb`` DataFrame (N^2 rows) and then pivots it into a dense N x N matrix before
handing it to scipy. For tens of thousands of genomes this is tens of GiB of RAM
and is the source of the crashes in issue #259 and the "Big dRep issue".

This module never builds the dense matrix and never needs to hold all N^2 pairs
in memory. It streams the MASH ``dist`` output, discards the ~99.9% of pairs that
are above the cutoff the instant they are read, and unions the survivors. Memory
is O(genomes + kept_edges) instead of O(genomes^2).
"""

import logging

import numpy as np
import pandas as pd

import drep.d_cluster.utils


class UnionFind:
    """
    Disjoint-set / union-find with path compression and union by rank.

    ``union`` and ``find`` are effectively O(alpha(N)) ~ O(1), so clustering the
    surviving edges is linear in the number of edges.
    """

    def __init__(self):
        self.parent = {}
        self.rank = {}

    def add(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0

    def find(self, x):
        # Find root
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        # Path compression (iterative, no recursion depth limit)
        while self.parent[x] != root:
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1

    def components(self):
        """
        Return {root: [members...]} for every set.
        """
        comps = {}
        for node in self.parent:
            comps.setdefault(self.find(node), []).append(node)
        return comps


def _components_to_cdb(uf):
    """
    Turn a populated UnionFind into a Cdb (columns: genome, primary_cluster).

    Clusters are numbered deterministically: largest first, ties broken by the
    alphabetically-smallest member. This makes runs reproducible regardless of
    the order edges happened to stream in.
    """
    comps = uf.components()

    ordered = sorted(
        comps.values(),
        key=lambda members: (-len(members), min(members)),
    )

    genomes = []
    clusters = []
    for cluster_id, members in enumerate(ordered, start=1):
        for genome in sorted(members):
            genomes.append(genome)
            clusters.append(cluster_id)

    Cdb = pd.DataFrame({'genome': genomes, 'primary_cluster': clusters})
    Cdb['primary_cluster'] = Cdb['primary_cluster'].astype(int)
    return Cdb


def cluster_long_df(db, cutoff, all_genomes=None):
    """
    Cluster an in-memory long-format MASH table with union-find (no pivot).

    This is the drop-in, single-linkage replacement for the pivot -> squareform
    -> scipy path in ``cluster_mash_database`` and for the ``low_ram`` path that
    used to ``stack()`` an already-dense matrix back into long format.

    Args:
        db: DataFrame with columns 'genome1', 'genome2', and either 'dist' or
            'similarity'.
        cutoff: distance cutoff (1 - P_ani). Pairs with dist <= cutoff are edges.
        all_genomes: optional iterable of every genome name, so singletons that
            never appear in an above-cutoff edge still get their own cluster. If
            not given, it is inferred from the genome1/genome2 columns.

    Returns:
        Cdb: DataFrame with columns 'genome', 'primary_cluster'.
    """
    if 'dist' in db.columns:
        dist = db['dist'].values
    else:
        dist = 1 - db['similarity'].values

    uf = UnionFind()

    # Seed every genome so singletons are represented
    if all_genomes is None:
        all_genomes = pd.unique(
            pd.concat([db['genome1'], db['genome2']], ignore_index=True)
        )
    for g in all_genomes:
        uf.add(g)

    mask = dist <= cutoff
    g1 = db['genome1'].values[mask]
    g2 = db['genome2'].values[mask]
    for a, b in zip(g1, g2):
        uf.add(a)
        uf.add(b)
        uf.union(a, b)

    return _components_to_cdb(uf)


def cluster_mash_files(dist_files, cutoff, all_genomes=None, chunksize=5_000_000,
                       name_from_fasta=True, progress=False):
    """
    Stream one or more MASH ``dist`` output files and cluster with union-find.

    Never builds the dense matrix and never holds all N^2 pairs in memory. Only
    genome names, the union-find bookkeeping, and one ``chunksize`` block of rows
    are resident at a time.

    Args:
        dist_files: path (str) or list of paths to MASH dist tsv output
            (columns: genome1, genome2, dist, p, kmers).
        cutoff: distance cutoff (1 - P_ani).
        all_genomes: optional iterable of every genome name to seed singletons.
        chunksize: rows per streamed block.
        name_from_fasta: if True, map file paths in the table to genome names via
            drep's basename logic (matches parse_mash_table behavior).
        progress: if True, show a tqdm progress bar over streamed rows.

    Returns:
        (Cdb, stats) where stats is a dict of counters for benchmarking/logging.
    """
    if isinstance(dist_files, (str, bytes)):
        dist_files = [dist_files]

    uf = UnionFind()
    if all_genomes is not None:
        for g in all_genomes:
            uf.add(g)

    if progress:
        try:
            from tqdm import tqdm
        except ImportError:
            logging.warning("tqdm not installed; primary-clustering progress bar disabled")
            progress = False

    total_rows = 0
    kept_edges = 0

    name_cache = {}

    def to_name(x):
        n = name_cache.get(x)
        if n is None:
            n = drep.d_cluster.utils._get_genome_name_from_fasta(x)
            name_cache[x] = n
        return n

    bar = tqdm(desc="  Primary clustering (streaming)", unit=" pairs") if progress else None

    for dist_file in dist_files:
        reader = pd.read_csv(
            dist_file,
            names=['genome1', 'genome2', 'dist', 'p', 'kmers'],
            usecols=['genome1', 'genome2', 'dist'],
            dtype={'genome1': str, 'genome2': str, 'dist': np.float32},
            sep='\t',
            chunksize=chunksize,
        )
        for chunk in reader:
            n = len(chunk)
            total_rows += n
            if bar is not None:
                bar.update(n)

            hits = chunk[chunk['dist'] <= cutoff]
            if len(hits) == 0:
                continue

            g1 = hits['genome1'].values
            g2 = hits['genome2'].values
            for a, b in zip(g1, g2):
                if name_from_fasta:
                    a = to_name(a)
                    b = to_name(b)
                uf.add(a)
                uf.add(b)
                uf.union(a, b)
                kept_edges += 1

    if bar is not None:
        bar.close()

    Cdb = _components_to_cdb(uf)

    stats = {
        'total_pairs_streamed': total_rows,
        'edges_kept': kept_edges,
        'genomes': len(uf.parent),
        'primary_clusters': Cdb['primary_cluster'].nunique() if len(Cdb) else 0,
    }
    return Cdb, stats


def build_ndb_from_edges(edges, Cdb):
    """
    Build a secondary-clustering Ndb out of primary's edge table, without running
    any new comparisons.

    Primary clustering (skani, sparse) already computed the exact ANI for every
    pair above skani's screening threshold. Secondary clustering only ever
    compares genomes *within* a primary cluster, and those pairs are a subset of
    what primary already has -- so re-running skani per primary cluster
    recomputes numbers that are already known, bit for bit.

    dRep's hierarchical secondary clustering needs a complete matrix per primary
    cluster, so pairs absent from the sparse edge list (i.e. below skani's
    screen, meaning no meaningful similarity) are filled in as ani=0 and
    coverage=0, and self-comparisons as 1.

    Args:
        edges: DataFrame from load_skani_sparse_edges.
        Cdb: primary clustering result with ['genome', 'primary_cluster'].

    Returns:
        Ndb: ['reference', 'querry', 'ani', 'alignment_coverage', 'primary_cluster']
    """
    g2p = Cdb.set_index('genome')['primary_cluster'].to_dict()

    e = edges.copy()
    e['pc'] = e['genome1'].map(g2p)
    # secondary only ever compares within a primary cluster
    e = e[e['pc'].notna() & (e['pc'] == e['genome2'].map(g2p))]

    have = set(zip(e['genome1'].values, e['genome2'].values))

    fill_r, fill_q, fill_pc, fill_ani, fill_cov = [], [], [], [], []
    for pc, sub in Cdb.groupby('primary_cluster'):
        gs = list(sub['genome'])
        for x in gs:
            for y in gs:
                if x == y:
                    fill_r.append(x); fill_q.append(y); fill_pc.append(pc)
                    fill_ani.append(1.0); fill_cov.append(1.0)
                elif (x, y) not in have:
                    fill_r.append(x); fill_q.append(y); fill_pc.append(pc)
                    fill_ani.append(0.0); fill_cov.append(0.0)

    kept = e.rename(columns={'genome1': 'reference', 'genome2': 'querry',
                             'pc': 'primary_cluster'})[
        ['reference', 'querry', 'ani', 'alignment_coverage', 'primary_cluster']]
    filled = pd.DataFrame({
        'reference': fill_r, 'querry': fill_q, 'ani': fill_ani,
        'alignment_coverage': fill_cov, 'primary_cluster': fill_pc,
    })

    Ndb = pd.concat([kept, filled], ignore_index=True)
    Ndb['primary_cluster'] = Ndb['primary_cluster'].astype(int)
    return Ndb


SKANI_SPARSE_COLUMNS = ['Ref_file', 'Query_file', 'ANI',
                        'Align_fraction_ref', 'Align_fraction_query']


def load_skani_sparse_edges(sparse_files, progress=False):
    """
    Load `skani triangle --sparse` output into a symmetric edge table.

    skani's sparse output is already only the pairs above its screening
    threshold, so unlike Mash's N^2 output it is small enough to hold in memory
    (~390k rows for 10,000 genomes) and can be reused rather than recomputed.

    Each input pair is emitted in both directions, because dRep treats
    alignment_coverage as the aligned fraction of the genome named in the first
    column, and the two directions have different coverages.

    Returns:
        DataFrame with ['genome1', 'genome2', 'ani', 'alignment_coverage'],
        where ani and alignment_coverage are 0-1 fractions.
    """
    if isinstance(sparse_files, (str, bytes)):
        sparse_files = [sparse_files]

    frames = []
    for f in sparse_files:
        d = pd.read_csv(f, sep='\t', usecols=SKANI_SPARSE_COLUMNS,
                        dtype={'Ref_file': str, 'Query_file': str,
                               'ANI': np.float32,
                               'Align_fraction_ref': np.float32,
                               'Align_fraction_query': np.float32})
        if len(d) == 0:
            continue
        a = np.array([drep.d_cluster.utils._get_genome_name_from_fasta(x)
                      for x in d['Ref_file'].values])
        b = np.array([drep.d_cluster.utils._get_genome_name_from_fasta(x)
                      for x in d['Query_file'].values])
        ani = (d['ANI'].values / 100).astype(np.float32)
        af_a = (d['Align_fraction_ref'].values / 100).astype(np.float32)
        af_b = (d['Align_fraction_query'].values / 100).astype(np.float32)

        frames.append(pd.DataFrame({
            'genome1': np.concatenate([a, b]),
            'genome2': np.concatenate([b, a]),
            'ani': np.concatenate([ani, ani]),
            # coverage is always the aligned fraction of the genome1 genome
            'alignment_coverage': np.concatenate([af_a, af_b]),
        }))

    if not frames:
        return pd.DataFrame(columns=['genome1', 'genome2', 'ani', 'alignment_coverage'])

    edges = pd.concat(frames, ignore_index=True)
    # drop self comparisons; they are added back explicitly where needed
    return edges[edges['genome1'] != edges['genome2']].reset_index(drop=True)


def cluster_edges(edges, ani_threshold, all_genomes, cov_threshold=0.0):
    """
    Union-find clustering of an in-memory edge table (see load_skani_sparse_edges).

    Args:
        edges: DataFrame with ['genome1', 'genome2', 'ani', 'alignment_coverage'].
        ani_threshold: minimum ANI (0-1) for a pair to be an edge.
        all_genomes: every genome name, so singletons get their own cluster.
        cov_threshold: minimum alignment_coverage (0-1) for a pair to be an edge.
            See run_skani_triangle_sparse for why this matters -- without it,
            genomes sharing a small conserved region chain together.

    Returns:
        (Cdb, stats)
    """
    uf = UnionFind()
    for g in all_genomes:
        uf.add(g)

    keep = edges['ani'].values >= ani_threshold
    if cov_threshold > 0:
        keep &= edges['alignment_coverage'].values >= cov_threshold

    g1 = edges['genome1'].values[keep]
    g2 = edges['genome2'].values[keep]
    for a, b in zip(g1, g2):
        uf.add(a)
        uf.add(b)
        uf.union(a, b)

    Cdb = _components_to_cdb(uf)
    stats = {
        'edges_total': len(edges),
        'edges_kept': int(keep.sum()),
        'genomes': len(uf.parent),
        'primary_clusters': Cdb['primary_cluster'].nunique() if len(Cdb) else 0,
    }
    return Cdb, stats


def cluster_skani_sparse_files(sparse_files, ani_threshold, all_genomes,
                               cov_threshold=0.0, chunksize=2_000_000, progress=False):
    """
    Stream `skani triangle --sparse` output and cluster with union-find.

    Unlike MASH ``dist``, the sparse skani output already contains *only* the
    above-screening-threshold pairs (an edge list), so there is no N^2 to stream
    at all -- just the surviving edges. This is the low-RAM, low-disk primary
    clustering path for very large genome sets.

    Expected columns (skani >= 0.2, with a header row):
        Ref_file, Query_file, ANI, Align_fraction_ref, Align_fraction_query,
        Ref_name, Query_name

    Args:
        sparse_files: path or list of paths to sparse skani output.
        ani_threshold: percent ANI (e.g. 90.0 for P_ani=0.9). Pairs at or above
            this are treated as edges.
        all_genomes: iterable of every genome name, so singletons that skani
            screened out still get their own primary cluster.
        cov_threshold: minimum aligned fraction (0-1) for a pair to count as an
            edge. skani reports two percentages; the larger is used (matching
            dRep's 'larger' coverage convention). 0 disables the filter.
        chunksize: rows per streamed block.
        progress: show a tqdm bar over streamed edges.

    Returns:
        (Cdb, Mdb, stats):
            Cdb: ['genome', 'primary_cluster']
            Mdb: reduced long-format table of the surviving edges only
                 (['genome1', 'genome2', 'similarity', 'dist']), for storage/plots.
            stats: dict of counters.
    """
    if isinstance(sparse_files, (str, bytes)):
        sparse_files = [sparse_files]

    uf = UnionFind()
    for g in all_genomes:
        uf.add(g)

    if progress:
        try:
            from tqdm import tqdm
        except ImportError:
            logging.warning("tqdm not installed; primary-clustering progress bar disabled")
            progress = False
    bar = tqdm(desc="  Primary clustering (sparse skani)", unit=" edges") if progress else None

    edges_seen = 0
    kept_edges = 0
    cov_pct = cov_threshold * 100.0
    mdb_g1, mdb_g2, mdb_sim = [], [], []

    name_cache = {}

    def to_name(x):
        n = name_cache.get(x)
        if n is None:
            n = drep.d_cluster.utils._get_genome_name_from_fasta(x)
            name_cache[x] = n
        return n

    for sparse_file in sparse_files:
        reader = pd.read_csv(
            sparse_file,
            sep='\t',
            usecols=['Ref_file', 'Query_file', 'ANI',
                     'Align_fraction_ref', 'Align_fraction_query'],
            dtype={'Ref_file': str, 'Query_file': str, 'ANI': np.float32,
                   'Align_fraction_ref': np.float32, 'Align_fraction_query': np.float32},
            chunksize=chunksize,
        )
        for chunk in reader:
            edges_seen += len(chunk)
            if bar is not None:
                bar.update(len(chunk))

            hits = chunk[chunk['ANI'] >= ani_threshold]
            if cov_pct > 0:
                larger_af = np.maximum(hits['Align_fraction_ref'].values,
                                       hits['Align_fraction_query'].values)
                hits = hits[larger_af >= cov_pct]
            if len(hits) == 0:
                continue

            for r, q, ani in zip(hits['Ref_file'].values,
                                 hits['Query_file'].values,
                                 hits['ANI'].values):
                a, b = to_name(r), to_name(q)
                if a == b:
                    continue
                uf.add(a)
                uf.add(b)
                uf.union(a, b)
                mdb_g1.append(a)
                mdb_g2.append(b)
                mdb_sim.append(ani / 100.0)
                kept_edges += 1

    if bar is not None:
        bar.close()

    Cdb = _components_to_cdb(uf)

    Mdb = pd.DataFrame({'genome1': mdb_g1, 'genome2': mdb_g2,
                        'similarity': np.array(mdb_sim, dtype=np.float32)})
    Mdb['dist'] = 1 - Mdb['similarity']

    stats = {
        'edges_seen': edges_seen,
        'edges_kept': kept_edges,
        'genomes': len(uf.parent),
        'primary_clusters': Cdb['primary_cluster'].nunique() if len(Cdb) else 0,
    }
    return Cdb, Mdb, stats
