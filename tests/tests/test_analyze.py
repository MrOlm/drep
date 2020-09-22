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

class test_analyze():
    def __init__(self):
        pass

    def setUp(self):
        self.s_wd_loc = test_utils.load_solutions_wd()
        self.working_wd_loc = test_utils.load_test_wd_loc()
        self.test_dir = test_utils.load_random_test_dir()

        self.tearDown()

        # copy over the data from solutions directory
        os.mkdir(self.working_wd_loc)
        os.mkdir(self.test_dir)
        shutil.copytree(os.path.join(self.s_wd_loc, 'data'), \
            os.path.join(self.working_wd_loc, 'data'))
        shutil.copytree(os.path.join(self.s_wd_loc, 'data_tables'), \
            os.path.join(self.working_wd_loc, 'data_tables'))
        shutil.copytree(os.path.join(self.s_wd_loc, 'log'), \
            os.path.join(self.working_wd_loc, 'log'))

        # edit Bbd to point to the right genomes
        bdb = drep.d_cluster.utils.load_genomes(test_utils.load_test_genomes())
        bdb.to_csv(os.path.join(self.working_wd_loc, 'data_tables', 'Bdb.csv'), \
            index=False)

        importlib.reload(logging)

    def tearDown(self):
        logging.shutdown()
        if os.path.isdir(self.working_wd_loc):
            shutil.rmtree(self.working_wd_loc)
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run(self):
        self.setUp()
        self.plot_1_test_1()
        self.tearDown()

        self.setUp()
        self.plot_5_test_1()
        self.tearDown()

        self.setUp()
        self.plot_5_test_2()
        self.tearDown()

        self.setUp()
        self.plot_6_test_1()
        self.tearDown()

        self.setUp()
        self.functional_test_1()
        self.tearDown()

        self.setUp()
        self.functional_test_2()
        self.tearDown()

    def plot_1_test_1(self):
        '''
        Test drep.d_analyze.mash_dendrogram_from_wd
        '''
        wd_loc = self.working_wd_loc
        wd = drep.WorkDirectory.WorkDirectory(wd_loc)
        test_dir = self.test_dir

        # Make sure it works
        assert len(glob.glob(test_dir + '/*')) == 0
        drep.d_analyze.mash_dendrogram_from_wd(wd, plot_dir=test_dir)
        assert len(glob.glob(test_dir + '/*')) == 1
        for f in glob.glob(test_dir + '/*'):
            assert os.path.getsize(f) > 0

        # Make sure it crashes gracefully if can't make plot
        os.remove(os.path.join(wd.get_dir('data_tables'), 'Mdb.csv'))
        wd = drep.WorkDirectory.WorkDirectory(wd_loc)
        drep.d_analyze.mash_dendrogram_from_wd(wd, plot_dir=test_dir)

    def plot_5_test_1(self):
        '''
        Test drep.d_analyze.plot_binscoring_from_wd

        Make sure it works without any genomeInfo
        '''
        wd_loc = self.working_wd_loc
        wd = drep.WorkDirectory.WorkDirectory(wd_loc)
        test_dir = self.test_dir

        assert len(glob.glob(test_dir + '/*')) == 0
        drep.d_analyze.plot_binscoring_from_wd(wd, plot_dir=test_dir)
        assert len(glob.glob(test_dir + '/*')) == 1
        for f in glob.glob(test_dir + '/*'):
            assert os.path.getsize(f) > 0

    def plot_5_test_2(self):
        '''
        Test drep.d_analyze.plot_binscoring_from_wd

        Make sure it works without any genomeInfo
        '''
        wd_loc = self.working_wd_loc
        wd = drep.WorkDirectory.WorkDirectory(wd_loc)
        test_dir = self.test_dir

        # Make sure it works whithout any genome info
        os.remove(os.path.join(wd.get_dir('data_tables'), 'Chdb.csv'))
        wd = drep.WorkDirectory.WorkDirectory(wd_loc)

        assert len(glob.glob(test_dir + '/*')) == 0
        drep.d_analyze.plot_binscoring_from_wd(wd, plot_dir=test_dir)
        assert len(glob.glob(test_dir + '/*')) == 1
        for f in glob.glob(test_dir + '/*'):
            assert os.path.getsize(f) > 0

    def plot_6_test_1(self):
        '''
        Test plot 6 with different things missing
        '''
        # Test with everything there
        args = argumentParser.parse_args(['analyze',self.working_wd_loc,'-pl', '6'])
        controller = Controller()
        controller.parseArguments(args)
        fig_dir = os.path.join(self.working_wd_loc, 'figures', '')

        figs = [os.path.basename(f) for f in glob.glob(fig_dir + '*')]
        FIGS = ['Winning_genomes.pdf']

        assert sorted(figs) == sorted(FIGS)
        for fig in glob.glob(fig_dir + '*'):
            assert os.path.getsize(fig) > 0

        # Test with removing Widb
        db_loc = os.path.join(self.working_wd_loc, 'data_tables', 'Widb.csv')
        os.remove(db_loc)
        for f in glob.glob(fig_dir + '*'):
            os.remove(f)

        args = argumentParser.parse_args(['analyze',self.working_wd_loc,'-pl', '6'])
        controller = Controller()
        controller.parseArguments(args)
        fig_dir = os.path.join(self.working_wd_loc, 'figures', '')

        figs = [os.path.basename(f) for f in glob.glob(fig_dir + '*')]
        FIGS = ['Winning_genomes.pdf']

        assert sorted(figs) == sorted(FIGS)
        for fig in glob.glob(fig_dir + '*'):
            assert os.path.getsize(fig) > 0

    def functional_test_1(self):
        '''
        Ensure analyze produces all plots
        '''
        args = argumentParser.parse_args(['analyze',self.working_wd_loc,'-pl'] + \
            ['a'])
        controller = Controller()
        controller.parseArguments(args)

        FIGS = ['Cluster_scoring.pdf', 'Clustering_scatterplots.pdf', \
            'Primary_clustering_dendrogram.pdf', 'Secondary_clustering_dendrograms.pdf', \
            'Winning_genomes.pdf', 'Secondary_clustering_MDS.pdf']

        fig_dir = os.path.join(self.working_wd_loc, 'figures', '')
        figs = [os.path.basename(f) for f in glob.glob(fig_dir + '*')]

        assert sorted(figs) == sorted(FIGS)
        for fig in glob.glob(fig_dir + '*'):
            assert os.path.getsize(fig) > 0

    def functional_test_2(self):
        '''
        Ensure analyze crashes gracefully
        '''
        wd_loc = self.working_wd_loc
        wd = drep.WorkDirectory.WorkDirectory(wd_loc)
        os.remove(os.path.join(wd.get_dir('data_tables'), 'Mdb.csv'))
        os.remove(os.path.join(wd.get_dir('data_tables'), 'Cdb.csv'))
        os.remove(os.path.join(wd.get_dir('data_tables'), 'Bdb.csv'))

        args = argumentParser.parse_args(['analyze',self.working_wd_loc,'-pl'] + \
            ['a'])
        controller = Controller()
        controller.parseArguments(args)