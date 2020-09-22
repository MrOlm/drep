import os
import glob
import shutil
import pandas as pd
import importlib
import logging

import tests.test_utils as test_utils

import drep
from drep import argumentParser
from drep.controller import Controller
from drep.WorkDirectory import WorkDirectory

class test_choose():
    def __init__(self):
        pass

    def setUp(self):
        self.s_wd_loc = test_utils.load_solutions_wd()
        self.working_wd_loc = test_utils.load_test_wd_loc()

        self.tearDown()

        # copy over the data from solutions directory
        os.mkdir(self.working_wd_loc)
        shutil.copytree(os.path.join(self.s_wd_loc, 'data'), \
            os.path.join(self.working_wd_loc, 'data'))
        shutil.copytree(os.path.join(self.s_wd_loc, 'data_tables'), \
            os.path.join(self.working_wd_loc, 'data_tables'))
        shutil.copytree(os.path.join(self.s_wd_loc, 'log'), \
            os.path.join(self.working_wd_loc, 'log'))

        importlib.reload(logging)

    def tearDown(self):
        logging.shutdown()
        if os.path.isdir(self.working_wd_loc):
            shutil.rmtree(self.working_wd_loc)

    def run(self):
        self.setUp()
        self.unit_test_1()
        self.tearDown()

        self.setUp()
        self.unit_test_2()
        self.tearDown()

    def unit_test_1(self):
        '''
        Ensure choose can handle when Chdb is not present, running checkM automatically
        '''
        # Delete Chdb
        wd_loc = self.working_wd_loc
        os.remove(wd_loc + '/data_tables/Chdb.csv')

        # Modify Bdb so the genome locations are right
        genomes = test_utils.load_test_genomes()
        g2l = {os.path.basename(g):g for g in genomes}

        Bdb = pd.read_csv(wd_loc + '/data_tables/Bdb.csv')
        Bdb['location'] = Bdb['genome'].map(g2l)
        Bdb.to_csv(wd_loc + '/data_tables/Bdb.csv', index=False)

        # Run choose - this should re-run checkM and re-generate chdb
        args = argumentParser.parse_args(['choose', wd_loc, '--checkM_method',\
            'taxonomy_wf'])
        controller = Controller()
        controller.parseArguments(args)

        Swd  = WorkDirectory(self.s_wd_loc)
        wd   = WorkDirectory(self.working_wd_loc)
        for db in ['Chdb', 'genomeInformation']:
            db1 = Swd.get_db(db)
            db2 =  wd.get_db(db)
            assert test_utils.compare_dfs(db1, db2), "{0} is not the same!".format(db)

    def unit_test_2(self):
        '''
        Try out the --skipCheckM argument for choose
        '''
        # Delete Chdb
        wd_loc = self.working_wd_loc
        os.remove(wd_loc + '/data_tables/Chdb.csv')
        os.remove(wd_loc + '/data_tables/Sdb.csv')
        os.remove(wd_loc + '/data_tables/Wdb.csv')

        # Run choose with --skipCheckM
        args = argumentParser.parse_args(['choose', wd_loc, '--ignoreGenomeQuality'])
        controller = Controller()
        controller.parseArguments(args)

        Swd  = WorkDirectory(self.s_wd_loc)
        wd   = WorkDirectory(self.working_wd_loc)
        for db in ['Sdb', 'Wdb', 'genomeInformation']:
            db1 = Swd.get_db(db)
            db2 =  wd.get_db(db)
            assert not test_utils.compare_dfs(db1, db2), "{0} is the same!".format(db)

        sdb = wd.get_db('Sdb')
        for s in sdb['score'].tolist():
            assert (s > 0) & (s < 5)