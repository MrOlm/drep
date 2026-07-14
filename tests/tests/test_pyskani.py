"""
Tests for the in-process pyskani backend (drep.d_cluster.pyskani_backend).

pyskani is an optional dependency, so every test here skips cleanly when it
isn't installed.
"""
import glob
import os
import shutil
import tempfile

import pandas as pd
import pytest

import drep.d_cluster.compare_utils as cu
import drep.d_cluster.cluster_utils as clu
import drep.d_cluster.external as ext
import drep.d_cluster.greedy_clustering as gc
import drep.d_cluster.utils
import drep.d_filter

def _has_pyskani():
    try:
        import pyskani  # noqa: F401
        return True
    except ImportError:
        return False


requires_pyskani = pytest.mark.skipif(not _has_pyskani(), reason="pyskani not installed")
requires_skani = pytest.mark.skipif(shutil.which('skani') is None, reason="skani not installed")
requires_fastani = pytest.mark.skipif(shutil.which('fastANI') is None, reason="fastANI not installed")


def _test_genomes():
    here = os.path.dirname(os.path.abspath(__file__))
    return sorted([g for g in glob.glob(os.path.join(here, '../genomes/*'))
                   if os.path.isfile(g)])


def _partition(Cdb, col='secondary_cluster'):
    return {frozenset(sub['genome']) for _, sub in Cdb.groupby(col)}


@requires_pyskani
def test_pyskani_ndb_shape():
    """run_pairwise_pyskani returns a complete, well-formed Ndb."""
    import drep.d_cluster.pyskani_backend as pb
    genomes = _test_genomes()
    Ndb = pb.run_pairwise_pyskani(genomes)

    assert list(Ndb.columns) == ['reference', 'querry', 'ani', 'alignment_coverage']
    # Every ordered pair present (dRep's hierarchical clustering needs a full matrix)
    assert len(Ndb) == len(genomes) ** 2
    # ANI and coverage are on a 0-1 scale
    assert Ndb['ani'].between(0, 1).all()
    assert Ndb['alignment_coverage'].between(0, 1).all()
    # Self comparisons are exactly 1
    selfs = Ndb[Ndb['reference'] == Ndb['querry']]
    assert (selfs['ani'] == 1).all()
    assert (selfs['alignment_coverage'] == 1).all()


@requires_pyskani
@requires_skani
def test_pyskani_agrees_with_subprocess_skani():
    """
    pyskani and the skani executable should report the same ANI/coverage, and
    produce the same secondary clusters.
    """
    import drep.d_cluster.pyskani_backend as pb
    genomes = _test_genomes()
    workdir = tempfile.mkdtemp()
    try:
        sub = ext.run_pairwise_skani(genomes, os.path.join(workdir, 'skani/'), processors=4)
        py = pb.run_pairwise_pyskani(genomes)

        m = pd.merge(sub, py, on=['reference', 'querry'], suffixes=('_sub', '_py'))
        # Only close pairs matter for clustering; skani's min-af filter drops
        # distant pairs from pyskani's output (they're filled in as ani=0).
        close = m[m['ani_sub'] >= 0.95]
        assert len(close) > 0
        assert (close['ani_sub'] - close['ani_py']).abs().max() < 0.001
        assert (close['alignment_coverage_sub'] - close['alignment_coverage_py']).abs().max() < 0.001

        for sa, nc in [(0.99, 0.1), (0.95, 0.1)]:
            c1, _ = clu.genome_hierarchical_clustering(sub, S_ani=sa, cov_thresh=nc,
                                                       comp_method='skani', cluster='X')
            c2, _ = clu.genome_hierarchical_clustering(py, S_ani=sa, cov_thresh=nc,
                                                       comp_method='pyskani', cluster='X')
            assert _partition(c1) == _partition(c2), f"clusters differ at S_ani={sa}"
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


@requires_pyskani
@requires_fastani
def test_greedy_pyskani_matches_greedy_fastani():
    """Greedy clustering should give the same answer via pyskani as via fastANI."""
    genomes = _test_genomes()
    bdb = drep.d_cluster.utils.load_genomes(genomes)
    bdb = drep.d_filter._add_lengthN50(bdb, bdb)

    workdir = tempfile.mkdtemp()
    try:
        parts = {}
        for alg in ['fastANI', 'pyskani']:
            d = os.path.join(workdir, alg + '/')
            os.makedirs(d, exist_ok=True)
            Ndb, Cdb, _ = gc.compare_genomes_greedy(
                bdb, alg, d, S_ani=0.95, cov_thresh=0.1, cluster='P1', processors=4)
            parts[alg] = _partition(Cdb)
        assert parts['fastANI'] == parts['pyskani']
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


@requires_pyskani
def test_greedy_pyskani_sketches_each_genome_once():
    """
    The whole point of the pyskani greedy path: every genome is sketched exactly
    once, no matter how many representatives accumulate.
    """
    import drep.d_cluster.pyskani_backend as pb
    genomes = _test_genomes()
    bdb = drep.d_cluster.utils.load_genomes(genomes)
    bdb = drep.d_filter._add_lengthN50(bdb, bdb)

    sketch_calls = []
    orig = pb.PyskaniDatabase.add

    def counting_add(self, name, contigs):
        sketch_calls.append(name)
        return orig(self, name, contigs)

    pb.PyskaniDatabase.add = counting_add
    workdir = tempfile.mkdtemp()
    try:
        gc.compare_genomes_greedy(bdb, 'pyskani', os.path.join(workdir, 'g/'),
                                  S_ani=0.95, cov_thresh=0.1, cluster='P1')
        # One sketch per representative, and never the same genome twice
        assert len(sketch_calls) == len(set(sketch_calls)), "a genome was sketched more than once"
        assert len(sketch_calls) <= len(genomes)
    finally:
        pb.PyskaniDatabase.add = orig
        shutil.rmtree(workdir, ignore_errors=True)


def test_estimate_time_handles_every_S_algorithm():
    """
    estimate_time only drives a log line, but it used to raise UnboundLocalError
    for any algorithm it didn't know about -- which killed a 10k-genome run at
    the start of secondary clustering. Every --S_algorithm choice must work, and
    unknown ones must not raise.
    """
    from drep.d_cluster.utils import estimate_time

    for alg in ['ANIn', 'gANI', 'goANI', 'ANImf', 'fastANI', 'skani', 'pyskani']:
        t = estimate_time(100, alg)
        assert t > 0, f"{alg} gave {t!r}"

    # An unrecognized algorithm must degrade gracefully, not raise
    assert estimate_time(100, 'some_future_algorithm') > 0


@requires_pyskani
def test_compare_genomes_dispatches_pyskani():
    """--S_algorithm pyskani is reachable through the normal dispatch path."""
    genomes = _test_genomes()
    bdb = drep.d_cluster.utils.load_genomes(genomes)
    workdir = tempfile.mkdtemp()
    try:
        Ndb = cu.compare_genomes(bdb, 'pyskani', workdir)
        assert len(Ndb) == len(genomes) ** 2
        assert 'ani' in Ndb.columns
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
