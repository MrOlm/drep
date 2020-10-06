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

import pytest

class Empty():
    pass

@pytest.fixture()
def self():
    # Set up
    self = Empty()
    self.s_wd_loc = test_utils.load_solutions_wd()
    self.working_wd_loc = test_utils.load_test_wd_loc()
    self.genomes = test_utils.load_test_genomes()

    logging.shutdown()
    if os.path.isdir(self.working_wd_loc):
        shutil.rmtree(self.working_wd_loc, ignore_errors=True)

    os.mkdir(self.working_wd_loc)
    shutil.copytree(os.path.join(self.s_wd_loc, 'data'), \
                    os.path.join(self.working_wd_loc, 'data'))
    shutil.rmtree(os.path.join(self.working_wd_loc, 'data', 'Clustering_files'))

    importlib.reload(logging)

    yield self

    # Teardown
    logging.shutdown()
    if os.path.isdir(self.working_wd_loc):
        shutil.rmtree(self.working_wd_loc, ignore_errors=True)
#
# class test_rerun():
#     '''
#     These tests require no calling of external programs, so are faster
#     '''
#
#     def setUp(self):
#         '''
#         copy the solutions data to the new work directory
#         '''
#         self.s_wd_loc = test_utils.load_solutions_wd()
#         self.working_wd_loc = test_utils.load_test_wd_loc()
#         self.genomes = test_utils.load_test_genomes()
#
#         self.tearDown()
#
#         # copy over the data from solutions directory
#         os.mkdir(self.working_wd_loc)
#         shutil.copytree(os.path.join(self.s_wd_loc, 'data'), \
#                         os.path.join(self.working_wd_loc, 'data'))
#         shutil.rmtree(os.path.join(self.working_wd_loc, 'data', 'Clustering_files'))
#
#         importlib.reload(logging)
#         if os.path.isdir(self.working_wd_loc):
#             shutil.rmtree(self.working_wd_loc)
#
#     def tearDown(self):
#         logging.shutdown()
#         if os.path.isdir(self.working_wd_loc):
#             shutil.rmtree(self.working_wd_loc)
#
#     def run(self):
#         # self.setUp()
#         # self.unit_tests_1()
#         # self.tearDown()
#
#         self.setUp()
#         self.unit_tests_2()
#         self.tearDown()
#
#         self.setUp()
#         self.unit_tests_3()
#         self.tearDown()
#
#         self.setUp()
#         self.unit_tests_4()
#         self.tearDown()
#
#         # self.setUp()
#         # self.unit_tests_5()
#         # self.tearDown()
#
#         self.setUp()
#         self.unit_tests_6()
#         self.tearDown()

def test_unit_1(self):
    '''
    Test a normal run of cluster
    '''
    # normal complete run
    args = argumentParser.parse_args(['dereplicate', self.working_wd_loc, '--S_algorithm', 'ANImf', '-g'] + \
                                     self.genomes)
    kwargs = vars(args)
    drep.d_cluster.controller.d_cluster_wrapper(self.working_wd_loc, **kwargs)

    # Verify
    Swd = WorkDirectory(self.s_wd_loc)
    wd = WorkDirectory(self.working_wd_loc)

    # Confirm the following are correct:
    # for db in ['Cdb', 'Mdb', 'Ndb']:
    for db in ['Cdb', 'Ndb']:
        db1 = Swd.get_db(db)
        db2 = wd.get_db(db)

        # get rid of some precision on the ANI; you are comparing fastANI with ANImf
        if db == 'Ndb':
            db1['ani'] = [round(x, 3) for x in db1['ani']]
            db2['ani'] = [round(x, 3) for x in db2['ani']]
            db1['alignment_length'] = [round(x, -6) for x in db1['alignment_length']]
            db2['alignment_length'] = [round(x, -6) for x in db2['alignment_length']]

            #db1 = db1[db2.columns]
            db1 = db1[['ani', 'alignment_length', 'querry', 'reference']]
            db2 = db2[['ani', 'alignment_length', 'querry', 'reference']]

            db1 = db1.sort_values(['querry', 'reference']).reset_index(drop=True)
            db2 = db2.sort_values(['querry', 'reference']).reset_index(drop=True)

        if db == 'Cdb':
            db1 = db1[['genome', 'secondary_cluster']].sort_values('genome').reset_index(drop=True)
            db2 = db2[['genome', 'secondary_cluster']].sort_values('genome').reset_index(drop=True)

        assert test_utils.compare_dfs2(db1, db2, verbose=True), "{0} is not the same!".format(db)

def test_unit_2(self):
    '''
    Test cluster with --SkipSecondary
    '''
    # run
    args = argumentParser.parse_args(['dereplicate', self.working_wd_loc, '-g'] + \
                                     self.genomes + ['--SkipSecondary'])
    kwargs = vars(args)
    drep.d_cluster.controller.d_cluster_wrapper(self.working_wd_loc, **kwargs)

    # Verify
    Swd = WorkDirectory(self.s_wd_loc)
    wd = WorkDirectory(self.working_wd_loc)

    # Confirm the following are the same:
    # for db in ['Mdb']:
    #     db1 = Swd.get_db(db)
    #     db2 =  wd.get_db(db)
    #     assert test_utils.compare_dfs(db1, db2), "{0} is not the same!".format(db)

    # Confirm the following are not the same:
    for db in ['Cdb', 'Ndb']:
        db1 = Swd.get_db(db)
        db2 = wd.get_db(db)
        assert not test_utils.compare_dfs(db1, db2), "{0} is the same! (and shouldn't be)".format(db)

def test_unit_3(self):
    '''
    Test cluster with --skipMash
    '''

    # normal complete run
    args = argumentParser.parse_args(['dereplicate', self.working_wd_loc, '-g'] + \
                                     self.genomes + ['--SkipMash'])
    kwargs = vars(args)
    drep.d_cluster.controller.d_cluster_wrapper(self.working_wd_loc, **kwargs)

    # Verify
    Swd = WorkDirectory(self.s_wd_loc)
    wd = WorkDirectory(self.working_wd_loc)

    # Confirm the following are not the same:
    for db in ['Cdb', 'Ndb']:  # , 'Mdb']:
        db1 = Swd.get_db(db)
        db2 = wd.get_db(db)

        assert not test_utils.compare_dfs(db1, db2), "{0} is the same! (and shouldn't be)".format(db)

def test_unit_4(self):
    '''
    Test changing cluster -pa
    '''
    # normal complete run
    args = argumentParser.parse_args(['dereplicate', self.working_wd_loc, '-g'] + \
                                     self.genomes + ['-pa', '0.10'])
    kwargs = vars(args)
    drep.d_cluster.controller.d_cluster_wrapper(self.working_wd_loc, **kwargs)

    # Verify
    Swd = WorkDirectory(self.s_wd_loc)
    wd = WorkDirectory(self.working_wd_loc)

    # Confirm the following are correct:
    # for db in ['Mdb']:
    #     db1 = Swd.get_db(db)
    #     db2 =  wd.get_db(db)
    #     assert test_utils.compare_dfs(db1, db2), "{0} is not the same!".format(db)

    # Confirm the following are not the same:
    for db in ['Ndb', 'Cdb']:
        db1 = Swd.get_db(db)
        db2 = wd.get_db(db)
        assert not test_utils.compare_dfs(db1, db2), "{0} is the same! (and shouldn't be)".format(db)

def test_unit_5(self):
    '''
    Test changing cluster --S_algorithm gANI
    '''
    loc, works = drep.d_bonus.find_program('ANIcalculator')
    if not works:
        return


    # normal complete run
    args = argumentParser.parse_args(['dereplicate', self.working_wd_loc, '-g'] + \
                                     self.genomes + ['--S_algorithm', 'gANI'])
    kwargs = vars(args)
    drep.d_cluster.controller.d_cluster_wrapper(self.working_wd_loc, **kwargs)

    # Verify
    Swd = WorkDirectory(self.s_wd_loc)
    wd = WorkDirectory(self.working_wd_loc)

    # Confirm the following are correct:
    for db in ['Cdb', 'Mdb']:
        db1 = Swd.get_db(db)
        db2 = wd.get_db(db)
        assert test_utils.compare_dfs(db1, db2), "{0} is not the same!".format(db)

def test_unit_6(self):
    '''
    Test drep call commands
    '''
    # try on single mash command

    wd = WorkDirectory(self.working_wd_loc)
    MASH_folder = wd.get_dir('MASH')
    log_folder = wd.get_dir('cmd_logs')

    mash_exe = 'mash'
    all_file = MASH_folder + 'ALL.msh'

    cmd = [mash_exe, 'dist', all_file, all_file, '>', MASH_folder
           + 'MASH_table.tsv']
    cmd = ' '.join(cmd)
    drep.run_cmd(cmd, shell=True, logdir=log_folder)

    assert len(glob.glob(log_folder + '*')) == 3


def test_unit_7(self):
    '''
    Test cluster with --SkipSecondary
    '''
    # run
    args = argumentParser.parse_args(['dereplicate', self.working_wd_loc, '-g'] + \
                                     self.genomes + ['--SkipSecondary'])
    kwargs = vars(args)
    drep.d_cluster.controller.d_cluster_wrapper(self.working_wd_loc, **kwargs)

    # Verify
    Swd = WorkDirectory(self.s_wd_loc)
    wd = WorkDirectory(self.working_wd_loc)

    # Confirm the following are the same:
    # for db in ['Mdb']:
    #     db1 = Swd.get_db(db)
    #     db2 =  wd.get_db(db)
    #     assert test_utils.compare_dfs(db1, db2), "{0} is not the same!".format(db)

    # Confirm the following are not the same:
    for db in ['Cdb', 'Ndb']:
        db1 = Swd.get_db(db)
        db2 = wd.get_db(db)
        assert not test_utils.compare_dfs(db1, db2), "{0} is the same! (and shouldn't be)".format(db)