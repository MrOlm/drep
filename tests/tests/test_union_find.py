"""
Unit tests for streaming union-find primary clustering (drep.d_cluster.union_find).
"""
import glob
import os
import shutil
import tempfile

import numpy as np
import pandas as pd
import pytest

import drep.d_cluster.union_find as uf
import drep.d_cluster.compare_utils as cu
import drep.d_cluster.utils


def _test_genomes():
    here = os.path.dirname(os.path.abspath(__file__))
    return [g for g in glob.glob(os.path.join(here, '../genomes/*'))
            if os.path.isfile(g)]


def test_union_find_basic():
    u = uf.UnionFind()
    for x in 'abcde':
        u.add(x)
    u.union('a', 'b')
    u.union('b', 'c')
    u.union('d', 'e')
    comps = {frozenset(v) for v in u.components().values()}
    assert comps == {frozenset('abc'), frozenset('de')}


def test_cluster_long_df_matches_expectation():
    # a-b close, c-d close, everything else far; e is a singleton
    rows = [
        ('a', 'b', 0.01), ('b', 'a', 0.01),
        ('c', 'd', 0.02), ('d', 'c', 0.02),
        ('a', 'c', 0.30), ('a', 'd', 0.30), ('a', 'e', 0.30),
        ('b', 'c', 0.30), ('b', 'e', 0.30), ('c', 'e', 0.30),
        ('d', 'e', 0.30),
    ]
    db = pd.DataFrame(rows, columns=['genome1', 'genome2', 'dist'])
    Cdb = uf.cluster_long_df(db, cutoff=0.1, all_genomes=list('abcde'))
    g2c = Cdb.set_index('genome')['primary_cluster'].to_dict()
    assert g2c['a'] == g2c['b']
    assert g2c['c'] == g2c['d']
    assert g2c['a'] != g2c['c']
    assert g2c['e'] != g2c['a'] and g2c['e'] != g2c['c']
    assert set(Cdb['genome']) == set('abcde')
    # deterministic numbering: largest clusters first
    assert Cdb['primary_cluster'].min() == 1


def test_cluster_mash_files_streaming(tmp_path):
    # write a small symmetric mash-style dist tsv
    f = tmp_path / "mash.tsv"
    names = [f"g{i}.fasta" for i in range(6)]
    block = [0, 0, 0, 1, 1, 1]  # two true clusters of 3
    with open(f, 'w') as o:
        for i, gi in enumerate(names):
            for j, gj in enumerate(names):
                d = 0.0 if i == j else (0.01 if block[i] == block[j] else 0.30)
                o.write(f"{gi}\t{gj}\t{d:.4f}\t0\t1000/1000\n")

    Cdb, stats = uf.cluster_mash_files(str(f), cutoff=0.1)
    assert stats['total_pairs_streamed'] == 36
    assert Cdb['primary_cluster'].nunique() == 2
    # genome names keep their basename (incl. extension), like parse_mash_table
    assert set(Cdb['genome']) == {f"g{i}.fasta" for i in range(6)}
    g2c = Cdb.set_index('genome')['primary_cluster'].to_dict()
    assert g2c['g0.fasta'] == g2c['g1.fasta'] == g2c['g2.fasta']
    assert g2c['g3.fasta'] == g2c['g4.fasta'] == g2c['g5.fasta']
    assert g2c['g0.fasta'] != g2c['g3.fasta']


def test_union_find_matches_scipy_membership():
    # Build a random symmetric similarity table and confirm union-find
    # and scipy single-linkage produce identical cluster membership.
    rng = np.random.default_rng(1)
    n = 40
    block = np.arange(n) // 4
    sim = np.where(block[:, None] == block[None, :], 0.99, 0.70)
    noise = np.triu(rng.normal(0, 0.01, (n, n)), 1)
    sim = np.clip(sim + noise + noise.T, 0, 1)
    np.fill_diagonal(sim, 1.0)

    rows = []
    for i in range(n):
        for j in range(n):
            rows.append((f"g{i}", f"g{j}", sim[i, j]))
    db = pd.DataFrame(rows, columns=['genome1', 'genome2', 'similarity'])

    scipy_Cdb, _ = cu.cluster_mash_database(db.copy(), P_ani=0.9,
                                            primary_clusterAlg='average',
                                            classic_primary_clustering=True)
    uf_Cdb, _ = cu.cluster_mash_database(db.copy(), P_ani=0.9,
                                         primary_clusterAlg='single')

    def membership(Cdb):
        return {frozenset(sub['genome']) for _, sub in Cdb.groupby('primary_cluster')}

    assert membership(scipy_Cdb) == membership(uf_Cdb)


def test_skani_sparse_min_af_filters_low_alignment_edges():
    """
    The aligned-fraction filter is what keeps skani ANI comparable to MASH's
    whole-genome similarity. skani reports ANI within aligned regions only, so a
    pair sharing one small conserved region looks like a high-ANI edge; under
    single linkage a few such bridges chain unrelated genomes into one giant
    primary cluster (measured on 10k UHGG genomes: min-af 0 collapsed 59% of the
    dataset into a single cluster).

    Here a and b are genuinely similar, and c is joined to each only by a
    high-ANI/low-alignment bridge. With the filter on, c must stay separate.
    """
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        f = os.path.join(td, 'sparse.tsv')
        rows = [
            # ref, query, ANI, af_ref, af_query
            ('a.fna', 'b.fna', 99.0, 90.0, 92.0),   # real relationship
            ('a.fna', 'c.fna', 98.0, 2.0, 3.0),     # spurious bridge (tiny overlap)
            ('b.fna', 'c.fna', 97.5, 2.5, 2.0),     # spurious bridge
        ]
        with open(f, 'w') as o:
            o.write("Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\n")
            for r in rows:
                o.write("\t".join(str(x) for x in r) + "\n")

        allg = ['a.fna', 'b.fna', 'c.fna']

        # No filter: the bridges chain a, b and c into one cluster
        C0, _, _ = uf.cluster_skani_sparse_files(f, 90.0, allg, cov_threshold=0.0)
        assert C0['primary_cluster'].nunique() == 1

        # With a 15% aligned-fraction floor, c is correctly left on its own
        C1, _, s1 = uf.cluster_skani_sparse_files(f, 90.0, allg, cov_threshold=0.15)
        g2c = C1.set_index('genome')['primary_cluster'].to_dict()
        assert g2c['a.fna'] == g2c['b.fna']
        assert g2c['c.fna'] != g2c['a.fna']
        assert s1['edges_kept'] == 1


def test_build_ndb_from_edges_fills_matrix():
    """
    Secondary clustering needs a complete matrix per primary cluster, so pairs
    absent from the sparse edge list must come back as ani=0 and self-pairs as 1.
    """
    edges = pd.DataFrame({
        'genome1': ['a', 'b'],
        'genome2': ['b', 'a'],
        'ani': [0.99, 0.99],
        'alignment_coverage': [0.9, 0.92],
    })
    Cdb = pd.DataFrame({'genome': ['a', 'b', 'c'], 'primary_cluster': [1, 1, 1]})
    Ndb = uf.build_ndb_from_edges(edges, Cdb)

    # complete 3x3 matrix for the one primary cluster
    assert len(Ndb) == 9
    g = Ndb.set_index(['reference', 'querry'])
    assert g.loc[('a', 'b'), 'ani'] == pytest.approx(0.99)
    assert g.loc[('a', 'a'), 'ani'] == 1.0
    assert g.loc[('a', 'a'), 'alignment_coverage'] == 1.0
    # c had no edges -> filled as no similarity
    assert g.loc[('a', 'c'), 'ani'] == 0.0
    assert g.loc[('c', 'a'), 'ani'] == 0.0
    # coverage is directional: fraction of the 'reference' genome
    assert g.loc[('a', 'b'), 'alignment_coverage'] == pytest.approx(0.9)
    assert g.loc[('b', 'a'), 'alignment_coverage'] == pytest.approx(0.92)


def test_build_ndb_from_edges_only_within_primary_clusters():
    """Secondary never compares across primary clusters."""
    edges = pd.DataFrame({
        'genome1': ['a', 'b'], 'genome2': ['b', 'a'],
        'ani': [0.99, 0.99], 'alignment_coverage': [0.9, 0.9],
    })
    Cdb = pd.DataFrame({'genome': ['a', 'b'], 'primary_cluster': [1, 2]})
    Ndb = uf.build_ndb_from_edges(edges, Cdb)
    # a and b are in different primary clusters: only self-comparisons survive
    assert set(zip(Ndb['reference'], Ndb['querry'])) == {('a', 'a'), ('b', 'b')}


@pytest.mark.skipif(shutil.which('skani') is None, reason="skani not installed")
def test_reused_edges_match_rerunning_skani():
    """
    The one-pass path must give exactly what re-running skani per primary cluster
    gives -- that is the whole premise for skipping the second pass.
    """
    import drep.d_cluster.compare_utils as cu
    import drep.d_cluster.utils

    genomes = _test_genomes()
    Bdb = drep.d_cluster.utils.load_genomes(genomes)
    workdir = tempfile.mkdtemp()
    try:
        Mdb, Cdb, _ = cu.primary_cluster_skani_sparse(
            Bdb, os.path.join(workdir, 'p'), P_ani=0.9, S_ani=0.99,
            cov_thresh=0.1, processors=4, primary_progress=False)

        # reuse primary's edges
        Ndb_r, Cdb_r, _ = cu.secondary_clustering_from_primary_edges(
            Bdb, Cdb, Mdb, S_ani=0.99, cov_thresh=0.1, clusterAlg='average')

        # re-run skani per primary cluster (the classic path)
        Ndb_c, Cdb_c, _ = cu.secondary_clustering(
            Bdb, Cdb, 'skani', os.path.join(workdir, 's'),
            S_ani=0.99, cov_thresh=0.1, clusterAlg='average', processors=4)

        def part(C):
            return {frozenset(s['genome']) for _, s in C.groupby('secondary_cluster')}

        assert part(Cdb_r) == part(Cdb_c)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


@pytest.mark.skipif(shutil.which('skani') is None, reason="skani not installed")
def test_sparse_skani_primary_matches_mash():
    """
    Sparse-skani primary clustering should recover the same primary partition as
    the classic MASH path on the bundled test genomes.
    """
    genomes = _test_genomes()
    Bdb = drep.d_cluster.utils.load_genomes(genomes)
    workdir = tempfile.mkdtemp()
    try:
        _, Cdb_sk, cret = cu.primary_cluster_skani_sparse(
            Bdb, os.path.join(workdir, 'sk'), P_ani=0.9, processors=4,
            primary_progress=False)
        _, Cdb_mash, _ = cu.all_vs_all_MASH(
            Bdb, os.path.join(workdir, 'mash'), P_ani=0.9, processors=4)

        # Every input genome is represented (including singletons skani screened out)
        assert set(Cdb_sk['genome']) == set(Bdb['genome'])
        # Small genome sets still get a scipy linkage so the primary dendrogram
        # can be drawn, and it must cover every genome -- not just the ones that
        # appear in the sparse edge list
        assert not isinstance(cret[0], str), "expected a real linkage matrix for a small set"
        assert list(cret[1].columns) == sorted(Bdb['genome'])

        def part(Cdb):
            return {frozenset(sub['genome']) for _, sub in Cdb.groupby('primary_cluster')}

        assert part(Cdb_sk) == part(Cdb_mash)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def test_classic_primary_clustering_uses_dense_path():
    """--classic_primary_clustering must produce a real scipy linkage matrix."""
    rng = np.random.default_rng(2)
    n = 20
    block = np.arange(n) // 4
    sim = np.where(block[:, None] == block[None, :], 0.99, 0.70)
    noise = np.triu(rng.normal(0, 0.01, (n, n)), 1)
    sim = np.clip(sim + noise + noise.T, 0, 1)
    np.fill_diagonal(sim, 1.0)
    rows = [(f"g{i:02d}", f"g{j:02d}", sim[i, j]) for i in range(n) for j in range(n)]
    db = pd.DataFrame(rows, columns=['genome1', 'genome2', 'similarity'])

    _, cret = cu.cluster_mash_database(db.copy(), P_ani=0.9,
                                       primary_clusterAlg='average',
                                       classic_primary_clustering=True)
    # cret[0] is a real linkage matrix (ndarray), not the streaming marker
    assert not isinstance(cret[0], str)
    assert cret[2]['linkage_method'] == 'average'
