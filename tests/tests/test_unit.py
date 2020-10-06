import pandas as pd
from collections import defaultdict

import drep

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