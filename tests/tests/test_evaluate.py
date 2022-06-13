import os
import glob
import shutil
import pandas as pd
import importlib
import logging
import pytest

import tests.test_utils as test_utils

import drep
from drep import argumentParser
from drep.controller import Controller
from drep.WorkDirectory import WorkDirectory

@pytest.fixture()
def self():
    self = test_utils.load_common_self()
    yield self
    self.teardown()

def test_tertiary_clustering_1(self):
    '''
    Test --run_tertiary_clustering fully
    '''
    test_dir = self.test_dir

    # Check that wont run without dereplicate
    args = drep.argumentParser.parse_args(
        ['compare', self.wd_loc, '--run_tertiary_clustering', '-g'] + self.genomes)
    try:
        drep.controller.Controller().parseArguments(args)
        assert False
    except ValueError:
        pass

    args = drep.argumentParser.parse_args(
        ['dereplicate', self.wd_loc, '--run_tertiary_clustering', '--ignoreGenomeQuality', '-g'] + self.genomes)
    drep.controller.Controller().parseArguments(args)

    # Load test results
    wd = drep.WorkDirectory.WorkDirectory(self.wd_loc)
    Cdb = wd.get_db('Cdb').sort_values('genome').reset_index(drop=True)

    # Load solutions
    wdS = drep.WorkDirectory.WorkDirectory(self.s_wd_loc)
    CdbS = wdS.get_db('Cdb').sort_values('genome').reset_index(drop=True)

    assert 'original_secondary_cluster' not in CdbS.columns
    assert 'original_secondary_cluster' in Cdb.columns

def test_tertiary_clustering_2(self):
    '''
    Quick tests for --run_tertiary_clustering fully
    '''
    test_dir = self.test_dir

    # Edit Cdb and Wdb
    wd = drep.WorkDirectory.WorkDirectory(self.working_wd_loc)
    Cdb = wd.get_db('Cdb')
    Cdb['secondary_cluster'] = [c if g != 'Enterococcus_faecalis_T2.fna' else '1_3' for c, g in zip(Cdb['secondary_cluster'], Cdb['genome'])]

    Wdb = wd.get_db('Wdb')
    db = pd.DataFrame({'genome':['Enterococcus_faecalis_T2.fna'], 'cluster':['1_3'], 'score':[50]})
    Wdb = pd.concat([Wdb, db])
    assert len(Wdb) == 5

    wd.store_db(Wdb, 'Wdb')
    wd.store_db(Cdb, 'Cdb')

    # Run tertiary clustering
    args = drep.argumentParser.parse_args(
        ['dereplicate', self.working_wd_loc, '--run_tertiary_clustering', '--S_algorithm', 'ANImf', '-sa', '0.99', '-g'] + self.genomes)
    drep.d_evaluate.d_evaluate_wrapper(args.work_directory, evaluate=['2'], **vars(args))

    wd = drep.WorkDirectory.WorkDirectory(self.working_wd_loc)
    Cdb = wd.get_db('Cdb').sort_values('genome').reset_index(drop=True)
    Wdb = wd.get_db('Wdb').sort_values('genome').reset_index(drop=True)
    assert len(Cdb['secondary_cluster'].unique()) == 4
    assert len(Cdb['original_secondary_cluster'].unique()) == 5
    assert len(Wdb) == 4
    assert '1_1.3' in Cdb['secondary_cluster'].tolist()