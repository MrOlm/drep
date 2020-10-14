import os
import glob
import shutil
import pandas as pd
import importlib
import logging

import tests.test_utils as test_utils

import pytest

import drep
from drep import argumentParser
from drep.controller import Controller
from drep.WorkDirectory import WorkDirectory

# class Empty():
#     pass
#
# @pytest.fixture()
# def self():
#     # Set up
#     self = Empty()
#     self.s_wd_loc = test_utils.load_solutions_wd()
#     self.working_wd_loc = test_utils.load_test_wd_loc()
#     logging.shutdown()
#     if os.path.isdir(self.working_wd_loc):
#         shutil.rmtree(self.working_wd_loc)
#
#     # copy over the data from solutions directory
#     os.mkdir(self.working_wd_loc)
#     shutil.copytree(os.path.join(self.s_wd_loc, 'data'), \
#                     os.path.join(self.working_wd_loc, 'data'))
#     shutil.copytree(os.path.join(self.s_wd_loc, 'data_tables'), \
#                     os.path.join(self.working_wd_loc, 'data_tables'))
#     shutil.copytree(os.path.join(self.s_wd_loc, 'log'), \
#                     os.path.join(self.working_wd_loc, 'log'))
#     importlib.reload(logging)
#
#     yield self
#
#     # Teardown
#     logging.shutdown()
#     if os.path.isdir(self.working_wd_loc):
#         shutil.rmtree(self.working_wd_loc)

@pytest.fixture()
def self():
    self = test_utils.load_common_self()
    yield self
    self.teardown()

def test_choose_1(self):
    '''
    Ensure choose can handle when Chdb is not present, running checkM automatically
    '''
    # Delete Chdb
    wd_loc = self.working_wd_loc
    os.remove(wd_loc + '/data_tables/Chdb.csv')

    # Modify Bdb so the genome locations are right
    genomes = test_utils.load_test_genomes()
    g2l = {os.path.basename(g): g for g in genomes}

    Bdb = pd.read_csv(wd_loc + '/data_tables/Bdb.csv')
    Bdb['location'] = Bdb['genome'].map(g2l)
    Bdb.to_csv(wd_loc + '/data_tables/Bdb.csv', index=False)

    # Run choose - this should re-run checkM and re-generate chdb
    drep.d_choose.d_choose_wrapper(wd_loc, checkM_method='taxonomy_wf')
    # args = argumentParser.parse_args(['choose', wd_loc, '--checkM_method', \
    #                                   'taxonomy_wf'])
    # controller = Controller()
    # controller.parseArguments(args)

    Swd = WorkDirectory(self.s_wd_loc)
    wd = WorkDirectory(self.working_wd_loc)
    for db in ['Chdb', 'genomeInformation']:
        db1 = Swd.get_db(db)
        db2 = wd.get_db(db)
        assert test_utils.compare_dfs(db1, db2), "{0} is not the same!".format(db)

def test_choose_2(self):
    '''
    Try out the --skipCheckM argument for choose
    '''
    # Delete Chdb
    wd_loc = self.working_wd_loc
    os.remove(wd_loc + '/data_tables/Chdb.csv')
    os.remove(wd_loc + '/data_tables/Sdb.csv')
    os.remove(wd_loc + '/data_tables/Wdb.csv')

    # Run choose with --skipCheckM
    args = argumentParser.parse_args(['dereplicate', wd_loc, '--ignoreGenomeQuality'])
    kwargs = vars(args)
    del kwargs['genomes']
    drep.d_choose.d_choose_wrapper(wd_loc, **kwargs)
    #
    # controller = Controller()
    # controller.parseArguments(args)

    Swd  = WorkDirectory(self.s_wd_loc)
    wd   = WorkDirectory(self.working_wd_loc)
    for db in ['Sdb', 'Wdb', 'genomeInformation']:
        db1 = Swd.get_db(db)
        db2 =  wd.get_db(db)
        assert not test_utils.compare_dfs(db1, db2), "{0} is the same!".format(db)

    sdb = wd.get_db('Sdb')
    Swd.get_db(db)
    for s in sdb['score'].tolist():
        assert (s > 0) & (s < 5)

    gdb = wd.get_db('genomeInformation')
    assert 'centrality' in gdb.columns

def test_centrality_1(self):
    """
    Test the methods drep.d_choose.add_centrality and "choose_winners" on a small set of genomes
    """
    wd = drep.WorkDirectory.WorkDirectory(self.working_wd_loc)
    kwargs = vars(argumentParser.parse_args(['dereplicate', self.working_wd_loc, '--ignoreGenomeQuality']))
    del kwargs['genomes']

    # Modify Cdb
    cdb = wd.get_db('Cdb')
    cdb['secondary_cluster'] = [x.replace('1_2', '1_1') for x in cdb['secondary_cluster']]
    wd.store_db(cdb, 'Cdb')

    # Run calculation
    bdb = wd.get_db('Bdb')
    Gdb = drep.d_filter.calc_genome_info(bdb['location'].tolist())
    Gdb = drep.d_choose.add_centrality(wd, Gdb, **kwargs)

    # Test result of add_centrality
    assert 'centrality' in list(Gdb.columns)
    assert len(Gdb[Gdb['centrality'] > 0]) > 0
    assert len(Gdb[Gdb['centrality'] > 1]) == 0
    assert len(Gdb[Gdb['centrality'].isna()]) == 0

    # Run choose winners
    Sdb, Wdb = drep.d_choose.choose_winners(cdb, Gdb, **kwargs)

    # Compare against choose winners with no centrality weight
    kwargs = vars(argumentParser.parse_args(['dereplicate', self.working_wd_loc, '--ignoreGenomeQuality', '-centW', '0']))
    del kwargs['genomes']
    Sdb2, Wdb2 = drep.d_choose.choose_winners(cdb, Gdb, **kwargs)

    # Make sure you get different values, and make sure they're not too different
    assert not test_utils.compare_dfs2(Sdb, Sdb2)
    assert abs(Sdb['score'].mean() - Sdb2['score'].mean()) < 1

    # Make sure S_ani is being loaded properly
    kwargs = vars(
        argumentParser.parse_args(['dereplicate', self.working_wd_loc, '--ignoreGenomeQuality', '-sa', '0.95']))
    del kwargs['genomes']
    Sdb2, Wdb2 = drep.d_choose.choose_winners(cdb, Gdb, **kwargs)
    assert not test_utils.compare_dfs2(Sdb, Sdb2)
    assert abs(Sdb['score'].mean()) < Sdb2['score'].mean()

def test_centrality_2(self):
    """
    Test calculating centrality from scratch
    """
    wd = drep.WorkDirectory.WorkDirectory(self.working_wd_loc)
    bdb = wd.get_db('Bdb')

    # Modify Cdb
    cdb = wd.get_db('Cdb')
    cdb['secondary_cluster'] = [x.replace('1_2', '1_1') for x in cdb['secondary_cluster']]

    # Run calculation
    Ndb = drep.d_choose.calc_centrality_from_scratch(bdb, cdb, os.path.join(wd.get_dir('MASH'), 'centrality_calculations/'))
    assert len(Ndb) == 9
    assert set(Ndb.columns.tolist()) == set(['reference', 'querry', 'ani', 'cluster'])