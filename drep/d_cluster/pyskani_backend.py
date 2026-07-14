"""
In-process skani comparisons via pyskani (https://github.com/althonos/pyskani).

The subprocess-based comparison algorithms re-sketch genomes on every
invocation. That is especially wasteful during greedy secondary clustering,
where each new genome is compared against the growing set of cluster
representatives by spawning a fresh subprocess that re-sketches *every*
representative -- O(N * R) sketching work and N subprocess spawns for N genomes.

pyskani lets us sketch each genome exactly once, keep the representatives in an
in-memory database, and query it directly. Sketching becomes O(N), there are no
subprocesses, and no temporary files are written.

This module is imported lazily: pyskani is an optional dependency, and dRep
falls back to the subprocess skani/fastANI implementations without it.
"""

import functools
import logging
import os

import pandas as pd

import drep.d_cluster.utils

# Column layout every secondary-clustering comparison must return
NDB_COLUMNS = ['reference', 'querry', 'ani', 'alignment_coverage']


def import_pyskani():
    """
    Import pyskani, raising an actionable error if it isn't installed.
    """
    try:
        import pyskani
    except ImportError:
        raise ImportError(
            "The 'pyskani' S_algorithm requires the pyskani package, which is not "
            "installed. Install it with `pip install pyskani` (or `pip install "
            "drep[pyskani]`), or choose a different --S_algorithm (e.g. skani, "
            "which uses the skani executable instead)."
        )
    return pyskani


def load_contigs(location):
    """
    Read a FASTA file into a list of contig sequences (as bytes), which is what
    pyskani's sketch/query expect.
    """
    from Bio import SeqIO
    return [bytes(record.seq) for record in SeqIO.parse(location, 'fasta')]


# During greedy clustering a genome is queried and then, if it founds a new
# cluster, immediately sketched as a representative. A tiny cache avoids parsing
# the same FASTA twice in a row. Deliberately kept at maxsize=2 -- caching every
# genome's contigs would hold the entire input set in memory.
@functools.lru_cache(maxsize=2)
def _load_contigs_cached(location):
    return load_contigs(location)


class PyskaniDatabase:
    """
    A pyskani database that sketches each genome exactly once.

    Sketching is the expensive part of an ANI comparison, so genomes are sketched
    on the way in and the resulting database is queried many times. Contigs are
    cached only while needed to sketch/query, not retained for the lifetime of
    the object.
    """

    def __init__(self, **kwargs):
        pyskani = import_pyskani()
        self.db = pyskani.Database()
        self.names = []
        # Screening cutoff. Mirrors skani's -s flag; the subprocess pairwise path
        # uses -s 1 (i.e. compare essentially everything), so default low here and
        # let callers raise it when they only care about close relatives.
        self.cutoff = kwargs.get('pyskani_cutoff', 0.01)
        self.learned_ani = kwargs.get('pyskani_learned_ani', None)

    def add(self, name, contigs):
        """Sketch a genome once and store it as a reference."""
        self.db.sketch(name, *contigs)
        self.names.append(name)

    def add_genome(self, location):
        """Sketch a genome from a FASTA path; returns its dRep genome name."""
        name = drep.d_cluster.utils._get_genome_name_from_fasta(location)
        self.add(name, _load_contigs_cached(location))
        return name

    def query_hits(self, name, contigs):
        """
        Query the database with a genome, returning the raw pyskani hits.
        """
        kw = {'cutoff': self.cutoff}
        if self.learned_ani is not None:
            kw['learned_ani'] = self.learned_ani
        return self.db.query(name, *contigs, **kw)

    def query(self, name, contigs, query_as_reference=False):
        """
        Query the database, returning (reference, querry, ani, alignment_coverage)
        tuples.

        Everywhere in dRep, alignment_coverage is the aligned fraction of the
        genome named in the 'reference' column (skani's .af matrix cell [A][B] is
        the aligned fraction of A; fastANI's matched/total is the fraction of the
        genome that load_fastani puts in 'reference'). Both orientations below
        respect that rule.

        Args:
            query_as_reference: if False (default, pairwise use), emit
                reference=the database genome and coverage=its aligned fraction.
                If True (greedy use), emit reference=the queried genome and
                coverage=the queried genome's aligned fraction -- which is what
                get_cluster_rep expects, since it reads the representative out of
                the 'querry' column.
        """
        rows = []
        for hit in self.query_hits(name, contigs):
            if query_as_reference:
                rows.append((hit.query_name, hit.reference_name,
                             hit.identity, hit.query_fraction))
            else:
                rows.append((hit.reference_name, hit.query_name,
                             hit.identity, hit.reference_fraction))
        return rows

    def query_genome(self, location, query_as_reference=False):
        name = drep.d_cluster.utils._get_genome_name_from_fasta(location)
        return self.query(name, _load_contigs_cached(location),
                          query_as_reference=query_as_reference)


def _fill_missing_pairs(rows, names):
    """
    skani only reports pairs that share enough k-mers. dRep's hierarchical
    secondary clustering needs a complete matrix, so absent pairs are filled in
    as ani=0 / coverage=0 (no detectable relatedness), matching what the
    subprocess `skani triangle --min-af 0` path yields for unrelated genomes.
    """
    have = {(r[0], r[1]) for r in rows}
    filled = list(rows)
    for a in names:
        for b in names:
            if (a, b) not in have:
                filled.append((a, b, 0.0, 0.0))
    return filled


def run_pairwise_pyskani(genome_list, **kwargs):
    """
    All-vs-all ANI within a set of genomes, in-process.

    Each genome is sketched exactly once and then queried against the database,
    so the sketching cost is linear in the number of genomes rather than
    quadratic.

    Args:
        genome_list: list of genome file locations.

    Returns:
        Ndb: DataFrame with ['reference', 'querry', 'ani', 'alignment_coverage'].
    """
    db = PyskaniDatabase(**kwargs)

    contigs = {}
    for location in genome_list:
        name = drep.d_cluster.utils._get_genome_name_from_fasta(location)
        cs = load_contigs(location)
        contigs[name] = cs
        db.add(name, cs)

    logging.debug(f"pyskani: sketched {len(contigs)} genomes once; querying")

    rows = []
    for name, cs in contigs.items():
        rows.extend(db.query(name, cs))

    rows = _fill_missing_pairs(rows, list(contigs.keys()))
    Ndb = pd.DataFrame(rows, columns=NDB_COLUMNS)

    # A genome can hit itself with identity slightly below 1 depending on
    # sketching; force exact self-identity like the other algorithms do.
    self_mask = Ndb['reference'] == Ndb['querry']
    Ndb.loc[self_mask, 'ani'] = 1.0
    Ndb.loc[self_mask, 'alignment_coverage'] = 1.0

    # Keep one row per ordered pair (query returns each direction once)
    Ndb = Ndb.drop_duplicates(subset=['reference', 'querry'], keep='first')
    return Ndb.reset_index(drop=True)


def pyskani_one_vs_many(location, db, **kwargs):
    """
    Compare one genome against an existing PyskaniDatabase of representatives.

    Used by greedy secondary clustering: the database holds every current cluster
    representative, already sketched, so this is a single in-process query rather
    than a subprocess that re-sketches all representatives.

    Emits rows with the representative in the 'querry' column (mirroring the
    fastANI greedy path), because get_cluster_rep reads the winning
    representative from there.

    Returns an Ndb-shaped DataFrame (possibly empty if nothing is similar).
    """
    rows = db.query_genome(location, query_as_reference=True)
    if len(rows) == 0:
        return pd.DataFrame(columns=NDB_COLUMNS)
    return pd.DataFrame(rows, columns=NDB_COLUMNS)
