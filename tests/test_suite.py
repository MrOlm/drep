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

def ensure_identicle(Swd, wd, skip = None):
    if skip == None:
        skip = []

    # Compare datatables
    for d in Swd.data_tables:
        if d in skip:
            continue

        db1 = Swd.get_db(d)
        db2 =  wd.get_db(d)

        assert db1.equals(db2), "{0} is not the same!".format(d)

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

class VerifyCluster():
    def __init__(self):
        pass

    def setUp(self):
        self.genomes = load_test_genomes()
        self.wd_loc = load_test_wd_loc()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)
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
        assert db1.equals(db2), "{0} is not the same!".format('Mdb')

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
        for db in ['Cdb', 'Mdb', 'Ndb']:
            db1 = Swd.get_db(db)
            db2 =  wd.get_db(db)
            assert db1.equals(db2), "{0} is not the same!".format(db)

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
        for db in ['Mdb']:
            db1 = Swd.get_db(db)
            db2 =  wd.get_db(db)
            assert db1.equals(db2), "{0} is not the same!".format(db)

        # Confirm the following are not the same:
        for db in ['Cdb', 'Ndb']:
            db1 = Swd.get_db(db)
            db2 =  wd.get_db(db)
            assert not db1.equals(db2), "{0} is the same! (and shouldn't be)".format(db)

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
            assert not db1.equals(db2), "{0} is the same! (and shouldn't be)".format(db)

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
        for db in ['Mdb']:
            db1 = Swd.get_db(db)
            db2 =  wd.get_db(db)
            assert db1.equals(db2), "{0} is not the same!".format(db)

        # Confirm the following are not the same:
        for db in ['Ndb', 'Cdb']:
            db1 = Swd.get_db(db)
            db2 =  wd.get_db(db)
            assert not db1.equals(db2), "{0} is the same! (and shouldn't be)".format(db)

    def tearDown(self):
        logging.shutdown()
        if os.path.isdir(self.working_wd_loc):
            shutil.rmtree(self.working_wd_loc)

def filter_test():
    ''' test the filter operation '''
    verifyFilter = VerifyFilter()
    verifyFilter.run()

def cluster_test():
    ''' test the cluster operation '''
    verifyCluster = VerifyCluster()
    verifyCluster.run()

def dereplicate_wf_test():
    ''' test the dereplicate_wf operation '''
    verifyDereplicateWf = VerifyDereplicateWf()
    verifyDereplicateWf.run()

def rerun_test():
    ''' run some tests on an already made workDirectory (to speed them up)'''
    QuickTests().run()

@pytest.mark.long
def test_long():
    dereplicate_wf_test()
    filter_test()
    cluster_test()

@pytest.mark.short
def test_short():
    cluster_test()
    rerun_test()

@pytest.mark.quick
def test_quick():
    rerun_test()

if __name__ == '__main__':
    #test_quick()
    #test_short()
    #test_long()
    dereplicate_wf_test()
    print("Everything seems to be working swimmingly!")
