import os
import tempfile
import pandas as pd
from collections import defaultdict

import drep
import drep.d_cluster.external

import tests.test_utils as test_utils

import pytest
class Empty():
    pass

@pytest.fixture()
def self():
    # Set up
    self = Empty()

    yield self

    # self.genomes = test_utils.load_test_genomes()
    # self.wd_loc = test_utils.load_test_wd_loc()
    # self.s_wd_loc = test_utils.load_solutions_wd()
    # self.testdir = test_utils.load_random_test_dir()
    #
    # importlib.reload(logging)
    # if os.path.isdir(self.wd_loc):
    #     shutil.rmtree(self.wd_loc)
    #
    # yield self
    #
    # # Teardown
    # logging.shutdown()
    # if os.path.isdir(self.wd_loc):
    #     shutil.rmtree(self.wd_loc)

def test_compare_dfs(self):
    '''
    test compare_dfs
    '''
    df1 = pd.DataFrame({'genome':['a', 'b', 'c'], 'value':[0, 0, 1]})
    df2 = pd.DataFrame({'genome':['c', 'b', 'a'], 'value':[1, 0, 0]})
    df3 = pd.DataFrame({'genome':['a', 'b', 'c'], 'value':[0, 0, 0]})
    df4 = pd.DataFrame({'value':[0, 0, 1], 'genome':['a', 'b', 'c']})

    assert not df1.equals(df2)
    assert not df1.equals(df3)
    assert df1.sort_index(axis=1).equals(df4.sort_index(axis=1))

    assert test_utils.compare_dfs(df1, df2, verbose=True)
    assert test_utils.compare_dfs(df1, df4, verbose=True)

    assert not test_utils.compare_dfs(df1, df3, verbose=True)

def test_scoring(self):
    '''
    test scoring calculation
    '''
    table = defaultdict(list)
    table['testnumber'].append("1")
    table['completeness'].append(100)
    table['contamination'].append(0)
    table['N50'].append(100)
    table['length'].append(100000)
    table['strain_heterogeneity'].append(0)

    table['testnumber'].append("2")
    table['completeness'].append(100)
    table['contamination'].append(0)
    table['N50'].append(100)
    table['length'].append(100000)
    table['strain_heterogeneity'].append(50)

    table['testnumber'].append("3")
    table['completeness'].append(88.2)
    table['contamination'].append(10)
    table['N50'].append(100)
    table['length'].append(100000)
    table['strain_heterogeneity'].append(50)
    df = pd.DataFrame(table)

    for i, row in df.groupby('testnumber'):
        score = drep.d_choose.score_row(row)
        if i == "1":
            assert score == 107.0
        elif i == "2":
            assert score == 107.0
        elif i == "3":
            assert score == 90.2

def test_n50(self):
    '''
    test N50 calculation
    '''
    import drep.d_filter
    genomes = test_utils.load_test_genomes()

    # test a genome with a single scaffold
    genome = [x for x in genomes if 'EC20' in x][0]
    n50 = drep.d_filter.calc_n50(genome)
    assert n50 == 3427276

    # test a real genome
    genome = [x for x in genomes if 'T2' in x][0]
    n50 = drep.d_filter.calc_n50(genome)
    assert n50 == 774663, n50

def test_load_fastani_pandas3_compat(self):
    '''
    Regression test for GitHub issue #299: load_fastani crashes with pandas>=3.0
    because delim_whitespace=True was removed. Use sep=r'\s+' instead.
    '''
    # Minimal fastANI output format: reference query ANI j1 j2 (whitespace-delimited)
    fastani_output = (
        "/path/to/genome_A.fasta\t/path/to/genome_B.fasta\t98.5\t900\t1000\n"
        "/path/to/genome_B.fasta\t/path/to/genome_A.fasta\t98.5\t850\t1000\n"
    )
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write(fastani_output)
        tmp_path = f.name

    try:
        fdb = drep.d_cluster.external.load_fastani(tmp_path)
        assert len(fdb) == 2
        assert set(fdb.columns) == {'reference', 'querry', 'ani', 'alignment_coverage'}
        assert abs(fdb['ani'].iloc[0] - 0.985) < 0.001
    finally:
        os.unlink(tmp_path)
def test_load_genomeInfo_delimiters(self):
    '''
    Regression test for GitHub issue #305: a tab-delimited genomeInfo file was
    parsed as a .csv, mangling the column names, because pandas doesn't raise
    when reading with the wrong delimiter
    '''
    header = ['genome', 'completeness', 'contamination']
    rows = [
        ['Enterococcus_casseliflavus_EC20.fasta', '98.28', '0.0'],
        ['Escherichia_coli_Sakai.fna', '100.0', '1.5'],
    ]

    for sep, suffix in [(',', '.csv'), ('\t', '.tsv'), (';', '.csv')]:
        table = '\n'.join([sep.join(r) for r in [header] + rows]) + '\n'
        with tempfile.NamedTemporaryFile(mode='w', suffix=suffix, delete=False) as f:
            f.write(table)
            tmp_path = f.name

        try:
            Idb = drep.d_filter.load_genomeInfo(tmp_path)
            assert set(Idb.columns) == set(header), (sep, list(Idb.columns))
            assert len(Idb) == 2
            assert Idb['completeness'].tolist() == [98.28, 100.0]
            assert Idb['genome'].tolist()[0] == 'Enterococcus_casseliflavus_EC20.fasta'
        finally:
            os.unlink(tmp_path)

def test_load_genomeInfo_missing_columns(self):
    '''
    A genomeInfo file without the required columns (e.g. raw CheckM2 output,
    which calls them "Name" and "Completeness") should still load, and then be
    rejected by _validate_genomeInfo with a message listing what was found
    '''
    table = "Name\tCompleteness\tContamination\ngenome_A\t98.28\t0.0\n"
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        f.write(table)
        tmp_path = f.name

    try:
        Idb = drep.d_filter.load_genomeInfo(tmp_path)
        assert list(Idb.columns) == ['Name', 'Completeness', 'Contamination']

        bdb = pd.DataFrame({'genome': ['genome_A'], 'location': ['/tmp/genome_A']})
        with pytest.raises(KeyError) as e:
            drep.d_filter._validate_genomeInfo(Idb, bdb)
        assert 'Completeness' in str(e.value)
    finally:
        os.unlink(tmp_path)
