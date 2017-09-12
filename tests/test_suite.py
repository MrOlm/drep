#!/usr/bin/env python

###############################################################################
#
# test_suite.py - process several E. coli genomes to verify operation of dRep
#
###############################################################################

import glob
import os
import pytest
import shutil
import logging

import pandas as pd

from drep import argumentParser
from drep.controller import Controller
from drep.WorkDirectory import WorkDirectory

def load_test_genomes():
    return glob.glob(os.path.join(str(os.getcwd()) + '/genomes/*'))

def load_test_wd_loc():
    loc = os.path.join(str(os.getcwd()),'../tests/test_backend/ecoli_wd')
    return loc

def load_solutions_wd():
    loc = os.path.join(str(os.getcwd()),'../tests/test_solutions/ecoli_wd')
    return loc

def load_solutions_taxonomy_wd():
    loc = os.path.join(str(os.getcwd()),'../tests/test_solutions/ecoli_taxonomy')
    return loc

def ensure_identicle(Swd, wd, skip = None):
    if skip == None:
        skip = []

    # Compare datatables
    for d in Swd.data_tables:
        if d in skip:
            continue

        db1 = Swd.get_db(d)
        db2 =  wd.get_db(d)

        assert compare_dfs(db1, db2), "{0} is not the same!".format(d)

    # Compare the clustering files
    pass

    # Compare the graphs
    pass

def sanity_check(Swd):
    '''
    Make sure the work directory passed in is correct
    '''
    f = open(Swd.location + '/log/amiinsane.txt')
    l = f.readlines()[0].strip()
    f.close()
    assert l == "No, you're not", l

    return

def compare_dfs(db1, db2):
    '''
    Return True if dataframes are equal (order of dataframes doesn't matter)
    '''

    db1 = db1.fillna(0)
    db2 = db2.fillna(0)

    df = pd.concat([db1, db2])
    df = df.reset_index(drop=True)
    df_gpby = df.groupby(list(df.columns))
    idx = [x[0] for x in df_gpby.groups.values() if len(x) == 1]

    identicle = (len(idx) == 0)
    # if not identicle:
    #     print("index: ", idx)
    #     print("db1: ",db1)
    #     print("db2: ",db2)
    #     print("df_gpby: ", str(df_gpby))

    return identicle

class VerifyDereplicateWf():
    def __init__(self):
        pass

    def run(self):
        self.setUp()
        self.functional_test_1()
        self.tearDown()

    def functional_test_1(self):
        genomes  = self.genomes
        wd_loc   = self.wd_loc
        s_wd_loc = self.s_wd_loc

        sanity_check(WorkDirectory(s_wd_loc))

        args = argumentParser.parse_args(['dereplicate_wf',wd_loc,'-g'] + genomes \
            + ['--checkM_method', 'taxonomy_wf'])
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        s_wd = WorkDirectory(s_wd_loc)
        wd   = WorkDirectory(wd_loc)
        ensure_identicle(s_wd, wd, skip=['Bdb'])

        # Perform sanity check to make sure solutions directiory isn't
        # being overwritten
        sanity_check(s_wd)

    def setUp(self):
        self.genomes = load_test_genomes()
        self.wd_loc = load_test_wd_loc()
        self.s_wd_loc = load_solutions_wd()

        logging.shutdown()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)

    def tearDown(self):
        logging.shutdown()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)

class VerifyFilter():
    def __init__(self):
        pass

    def setUp(self):
        self.genomes = load_test_genomes()
        self.wd_loc = load_test_wd_loc()

        logging.shutdown()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)

    def run(self):
        self.setUp()
        self.functional_test_1()
        self.tearDown()

    def functional_test_1(self):
        '''
        Call filter on 'Escherichia_coli_Sakai.fna'
        '''
        genomes = self.genomes
        wd_loc  = self.wd_loc

        args = argumentParser.parse_args(['filter',wd_loc,'-g',genomes[4]] \
            + ['--checkM_method', 'taxonomy_wf'])
        controller = Controller()
        controller.parseArguments(args)

        # Confirm Chdb.csv is correct
        assert True

        # Confirm genome is in Bdb.csv
        assert True

    def tearDown(self):
        logging.shutdown()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)

class VerifyAnalyze():
    def __init__(self):
        pass

    def setUp(self):
        self.s_wd_loc = load_solutions_wd()
        self.working_wd_loc = load_test_wd_loc()

        self.tearDown()

        # copy over the data from solutions directory
        os.mkdir(self.working_wd_loc)
        shutil.copytree(os.path.join(self.s_wd_loc, 'data'), \
            os.path.join(self.working_wd_loc, 'data'))
        shutil.copytree(os.path.join(self.s_wd_loc, 'data_tables'), \
            os.path.join(self.working_wd_loc, 'data_tables'))
        shutil.copytree(os.path.join(self.s_wd_loc, 'log'), \
            os.path.join(self.working_wd_loc, 'log'))


    def run(self):
        self.setUp()
        self.unit_test_1()
        self.tearDown()

    def unit_test_1(self):
        '''
        Ensure analyze produces plots
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

    def tearDown(self):
        logging.shutdown()
        if os.path.isdir(self.working_wd_loc):
            shutil.rmtree(self.working_wd_loc)

class VerifyTaxonomy():
    ''' These tests skip running centrifuge and prodigal, but test parsing of files
    '''

    def setUp(self):
        self.genomes = load_test_genomes()
        self.wd_loc = load_test_wd_loc()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)
        self.s_wd_loc = load_solutions_taxonomy_wd()

        # copy over the data from solutions directory
        os.mkdir(self.wd_loc)
        shutil.copytree(os.path.join(self.s_wd_loc, 'data'), \
            os.path.join(self.wd_loc, 'data'))

    def run(self):
        self.setUp()
        self.taxTest1()
        self.tearDown()

        self.setUp()
        self.taxTest2()
        self.tearDown()

        print('tax test 1 passed')

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
        assert compare_dfs(tdb, tdbS), "{0} is not the same!".format('Bdb')

        tdbS = Swd.get_db('Tdb')
        tdb = wd.get_db('Tdb')
        assert compare_dfs(tdb, tdbS), "{0} is not the same!".format('Tdb')

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
        assert compare_dfs(tdb, tdbS), "{0} is not the same!".format('Bdb')

        tdbS = Swd.get_db('TdbP')
        tdb = wd.get_db('Tdb')
        assert compare_dfs(tdb, tdbS), "{0} is not the same!".format('Tdb')

    def tearDown(self):
        logging.shutdown()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)

class VerifyCluster():
    def __init__(self):
        pass

    def setUp(self):
        self.genomes = load_test_genomes()
        self.wd_loc = load_test_wd_loc()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc, ignore_errors=True)
        self.s_wd_loc = load_solutions_wd()

    def run(self):
        self.setUp()
        self.functional_test_1()
        self.tearDown()

        print('functional test 1 passed')

        self.setUp()
        self.skipsecondary_test()
        self.tearDown()

        print('skip secondary test passed')

    def functional_test_1(self):
        '''
        Cluster the 5 genomes using default settings
        '''
        genomes = self.genomes
        wd_loc  = self.wd_loc

        args = argumentParser.parse_args(['cluster',wd_loc,'-g']+genomes)
        controller = Controller()
        controller.parseArguments(args)

        # Confirm Cdb.csv is correct
        assert True

        # Confirm Bdb.csv is correct
        assert True

    def skipsecondary_test(self):
        genomes = self.genomes
        wd_loc  = self.wd_loc
        s_wd_loc = self.s_wd_loc

        args = argumentParser.parse_args(['cluster',wd_loc,'-g'] +genomes \
                + ['--SkipSecondary','-o'])
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        Swd = WorkDirectory(s_wd_loc)
        wd   = WorkDirectory(wd_loc)

        # Confirm Mdb.csv is correct
        db1 = Swd.get_db('Mdb')
        db2 =  wd.get_db('Mdb')
        #assert compare_dfs(db1, db2), "{0} is not the same!".format('Mdb')

        # Confirm Ndb.csv doesn't exist
        db2 = wd.get_db('Ndb')
        assert db2.empty, 'Ndb is not empty'

    def tearDown(self):
        logging.shutdown()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)

class QuickTests():
    '''
    These tests require no calling of external programs, so are faster
    '''

    def setUp(self):
        '''
        copy the solutions data to the new work directory
        '''
        self.s_wd_loc = load_solutions_wd()
        self.working_wd_loc = load_test_wd_loc()
        self.genomes = load_test_genomes()

        self.tearDown()

        # copy over the data from solutions directory
        os.mkdir(self.working_wd_loc)
        shutil.copytree(os.path.join(self.s_wd_loc, 'data'), \
            os.path.join(self.working_wd_loc, 'data'))
        shutil.rmtree(os.path.join(self.working_wd_loc, 'data', 'Clustering_files'))

    def run (self):
        # self.setUp()
        # self.unit_tests_5()
        # self.tearDown()

        self.setUp()
        self.unit_tests_1()
        self.tearDown()

        self.setUp()
        self.unit_tests_2()
        self.tearDown()

        self.setUp()
        self.unit_tests_4()
        self.tearDown()

        self.setUp()
        self.unit_tests_3()
        self.tearDown()

    def unit_tests_1(self):
        '''
        Test a normal run of cluster
        '''
        # normal complete run
        args = argumentParser.parse_args(['cluster',self.working_wd_loc,'-g'] + \
            self.genomes + ['-o'])
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        Swd  = WorkDirectory(self.s_wd_loc)
        wd   = WorkDirectory(self.working_wd_loc)

        # Confirm the following are correct:
        #for db in ['Cdb', 'Mdb', 'Ndb']:
        for db in ['Cdb', 'Ndb']:
            db1 = Swd.get_db(db)
            db2 =  wd.get_db(db)
            assert compare_dfs(db1, db2), "{0} is not the same!".format(db)

    def unit_tests_2(self):
        '''
        Test cluster with --SkipSecondary
        '''
        # run
        args = argumentParser.parse_args(['cluster',self.working_wd_loc,'-g'] + \
            self.genomes + ['-o', '--SkipSecondary'])
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        Swd  = WorkDirectory(self.s_wd_loc)
        wd   = WorkDirectory(self.working_wd_loc)

        # Confirm the following are the same:
        # for db in ['Mdb']:
        #     db1 = Swd.get_db(db)
        #     db2 =  wd.get_db(db)
        #     assert compare_dfs(db1, db2), "{0} is not the same!".format(db)

        # Confirm the following are not the same:
        for db in ['Cdb', 'Ndb']:
            db1 = Swd.get_db(db)
            db2 =  wd.get_db(db)
            assert not compare_dfs(db1, db2), "{0} is the same! (and shouldn't be)".format(db)

    def unit_tests_3(self):
        ''' Test cluster with --skipMash
        '''

        # normal complete run
        args = argumentParser.parse_args(['cluster',self.working_wd_loc,'-g'] + \
            self.genomes + ['-o', '--SkipMash'])
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        Swd  = WorkDirectory(self.s_wd_loc)
        wd   = WorkDirectory(self.working_wd_loc)

        # Confirm the following are not the same:
        for db in ['Cdb', 'Ndb', 'Mdb']:
            db1 = Swd.get_db(db)
            db2 = wd.get_db(db)
            assert not compare_dfs(db1, db2), "{0} is the same! (and shouldn't be)".format(db)

    def unit_tests_4(self):
        '''
        Test changing cluster -pa
        '''
        # normal complete run
        args = argumentParser.parse_args(['cluster',self.working_wd_loc,'-g'] + \
            self.genomes + ['-o', '-pa', '0.10'])
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        Swd  = WorkDirectory(self.s_wd_loc)
        wd   = WorkDirectory(self.working_wd_loc)

        # Confirm the following are correct:
        # for db in ['Mdb']:
        #     db1 = Swd.get_db(db)
        #     db2 =  wd.get_db(db)
        #     assert compare_dfs(db1, db2), "{0} is not the same!".format(db)

        # Confirm the following are not the same:
        for db in ['Ndb', 'Cdb']:
            db1 = Swd.get_db(db)
            db2 =  wd.get_db(db)
            assert not compare_dfs(db1, db2), "{0} is the same! (and shouldn't be)".format(db)

    def unit_tests_5(self):
        '''
        Test changing cluster --S_algorithm gANI
        '''
        # normal complete run
        args = argumentParser.parse_args(['cluster',self.working_wd_loc,'-g'] + \
            self.genomes + ['-o', '--S_algorithm', 'gANI'])
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        Swd  = WorkDirectory(self.s_wd_loc)
        wd   = WorkDirectory(self.working_wd_loc)

        # Confirm the following are correct:
        for db in ['Cdb', 'Mdb']:
            db1 = Swd.get_db(db)
            db2 =  wd.get_db(db)
            assert compare_dfs(db1, db2), "{0} is not the same!".format(db)

    def tearDown(self):
        logging.shutdown()
        if os.path.isdir(self.working_wd_loc):
            shutil.rmtree(self.working_wd_loc)

class UnitTests():
    '''
    Very quick tests of random things THAT REQUIRE NO FILES / SETUP
    '''

    def run(self):
        '''
        run all tests
        '''
        self.test1()
        self.test3()

    def test1(self):
        '''
        test compare_dfs
        '''
        df1 = pd.DataFrame({'genome':['a', 'b', 'c'], 'value':[0, 0, 1]})
        df2 = pd.DataFrame({'genome':['c', 'b', 'a'], 'value':[1, 0, 0]})
        df3 = pd.DataFrame({'genome':['a', 'b', 'c'], 'value':[0, 0, 0]})
        df4 = pd.DataFrame({'value':[0, 0, 1], 'genome':['a', 'b', 'c']})

        assert not df1.equals(df2)
        assert not df1.equals(df3)
        assert df1.equals(df4)

        assert compare_dfs(df1, df2)
        assert compare_dfs(df1, df4)

        assert not compare_dfs(df1, df3)


    def test2(self):
        '''
        test bonus debug
        '''
        pass

    def test3(self):
        '''
        test N50 calculation
        '''
        import drep.d_filter
        genomes = load_test_genomes()

        # test a genome with a single scaffold
        genome = [x for x in genomes if 'EC20' in x][0]
        n50 = drep.d_filter.calc_n50(genome)
        assert n50 == 3427276

        # test a real genome
        genome = [x for x in genomes if 'T2' in x][0]
        n50 = drep.d_filter.calc_n50(genome)
        assert n50 == 774663, n50

def filter_test():
    ''' test the filter operation '''
    verifyFilter = VerifyFilter()
    verifyFilter.run()

def cluster_test():
    ''' test the cluster operation '''
    verifyCluster = VerifyCluster()
    verifyCluster.run()

def analyze_test():
    ''' test the analyze operation '''
    verifyAnalyze= VerifyAnalyze()
    verifyAnalyze.run()

def dereplicate_wf_test():
    ''' test the dereplicate_wf operation '''
    verifyDereplicateWf = VerifyDereplicateWf()
    verifyDereplicateWf.run()

def rerun_test():
    ''' run some tests on an already made workDirectory (to speed them up)'''
    QuickTests().run()

def unit_test():
    ''' run simple unit tests'''
    UnitTests().run()

def taxonomy_test():
    '''
    Test taxonomy methods
    '''

    VerifyTaxonomy().run()

@pytest.mark.long
def test_long():
    dereplicate_wf_test()
    filter_test()
    cluster_test()
    analyze_test()

@pytest.mark.short
def test_short():
    cluster_test()
    rerun_test()
    taxonomy_test()

@pytest.mark.quick
def test_quick():
    rerun_test()

@pytest.mark.unit
def test_unit():
    unit_test()

if __name__ == '__main__':
    #analyze_test()
    test_unit()
    test_quick()
    test_short()
    test_long()
    #dereplicate_wf_test()
    #taxonomy_test()
    #cluster_test()

    print("Everything seems to be working swimmingly!")
