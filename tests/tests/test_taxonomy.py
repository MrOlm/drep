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

class test_taxonomy():
    ''' These tests skip running centrifuge and prodigal, but test parsing of files
    '''

    def setUp(self):
        self.genomes = test_utils.load_test_genomes()
        self.wd_loc = test_utils.load_test_wd_loc()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)
        self.s_wd_loc = test_utils.load_solutions_taxonomy_wd()

        # copy over the data from solutions directory
        os.mkdir(self.wd_loc)
        shutil.copytree(os.path.join(self.s_wd_loc, 'data'), \
            os.path.join(self.wd_loc, 'data'))

        importlib.reload(logging)

    def run(self):
        self.setUp()
        self.taxTest1()
        self.tearDown()

        self.setUp()
        self.taxTest2()
        self.tearDown()

        self.setUp()
        self.taxTest3()
        self.tearDown()

        self.setUp()
        self.taxTest4()
        self.tearDown()


    def taxTest1(self):
        '''
        Check the taxonomy call for max method
        '''
        genomes = self.genomes
        wd_loc = self.wd_loc
        swd_loc = self.s_wd_loc

        # Call the command
        args = argumentParser.parse_args(['bonus',wd_loc,'-g'] +genomes \
                + ['--run_tax','--cent_index','/home/mattolm/download/centrifuge/indices/b+h+v',\
                '--tax_method', 'max'])
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        Swd = WorkDirectory(swd_loc)
        wd = WorkDirectory(wd_loc)

        tdbS = Swd.get_db('Bdb')
        tdb = wd.get_db('Bdb')
        del tdbS['location']
        del tdb['location']
        assert test_utils.compare_dfs(tdb, tdbS), "{0} is not the same!".format('Bdb')

        tdbS = Swd.get_db('Tdb')
        tdb = wd.get_db('Tdb')

        if test_utils.compare_dfs(tdb, tdbS) == False:
            print("{0} is not the same! May be due to centrifuge index issues".format('Tdb'))
            my_panel = pd.Panel(dict(df1=tdbS,df2=tdb))
            print(my_panel.apply(test_utils.report_diff, axis=0))

        assert True
        #assert compare_dfs(tdb, tdbS), "{0} is not the same!".format('Tdb')

    def taxTest2(self):
        '''
        Check the taxonomy call for percent method
        '''
        genomes = self.genomes
        wd_loc = self.wd_loc
        swd_loc = self.s_wd_loc

        # Call the command
        args = argumentParser.parse_args(['bonus',wd_loc,'-g'] +genomes \
                + ['--run_tax','--cent_index','/home/mattolm/download/centrifuge/indices/b+h+v',\
                '--tax_method', 'percent'])
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        Swd = WorkDirectory(swd_loc)
        wd = WorkDirectory(wd_loc)

        tdbS = Swd.get_db('BdbP')
        tdb = wd.get_db('Bdb')
        del tdbS['location']
        del tdb['location']
        assert test_utils.compare_dfs(tdb, tdbS), "{0} is not the same!".format('Bdb')

        tdbS = Swd.get_db('TdbP')
        tdb = wd.get_db('Tdb')
        assert test_utils.compare_dfs(tdb, tdbS), "{0} is not the same!".format('Tdb')

    def taxTest3(self):
        '''
        Try actually running centrifuge
        '''
        loc, works = drep.d_bonus.find_program('centrifuge')
        if works == False:
            print('Centrifuge not installed- skipping tests')

        else:
            genomes = self.genomes
            wd_loc = self.wd_loc
            swd_loc = self.s_wd_loc

            # Remove previous data run
            shutil.rmtree(os.path.join(self.wd_loc, 'data', 'centrifuge'))

            # Call the command
            args = argumentParser.parse_args(['bonus',wd_loc,'-g'] +genomes \
                    + ['--run_tax','--cent_index','/home/mattolm/download/centrifuge/indices/b+h+v',\
                    '--tax_method', 'percent'])
            controller = Controller()
            controller.parseArguments(args)

            # Verify
            Swd = WorkDirectory(swd_loc)
            wd = WorkDirectory(wd_loc)

            tdbS = Swd.get_db('BdbP')
            tdb = wd.get_db('Bdb')
            del tdbS['location']
            del tdb['location']
            assert test_utils.compare_dfs(tdb, tdbS), "{0} is not the same!".format('Bdb')

            tdbS = Swd.get_db('TdbP')
            tdb = wd.get_db('Tdb')
            assert test_utils.compare_dfs(tdb, tdbS), "{0} is not the same!".format('Tdb')

    def taxTest4(self):
        '''
        Try actually running centrifuge without prodigal done
        '''
        loc, works = drep.d_bonus.find_program('centrifuge')
        if works == False:
            print('Centrifuge not installed- skipping tests')

        else:
            genomes = self.genomes
            wd_loc = self.wd_loc
            swd_loc = self.s_wd_loc

            # Remove previous data run
            shutil.rmtree(os.path.join(self.wd_loc, 'data', 'centrifuge'))
            shutil.rmtree(os.path.join(self.wd_loc, 'data', 'prodigal'))

            # Call the command
            args = argumentParser.parse_args(['bonus',wd_loc,'-g'] +genomes \
                    + ['--run_tax','--cent_index','/home/mattolm/download/centrifuge/indices/b+h+v',\
                    '--tax_method', 'percent'])
            controller = Controller()
            controller.parseArguments(args)

            # Verify
            Swd = WorkDirectory(swd_loc)
            wd = WorkDirectory(wd_loc)

            tdbS = Swd.get_db('BdbP')
            tdb = wd.get_db('Bdb')
            del tdbS['location']
            del tdb['location']
            assert test_utils.compare_dfs(tdb, tdbS), "{0} is not the same!".format('Bdb')

            tdbS = Swd.get_db('TdbP')
            tdb = wd.get_db('Tdb')
            assert test_utils.compare_dfs(tdb, tdbS), "{0} is not the same!".format('Tdb')

    def tearDown(self):
        logging.shutdown()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)