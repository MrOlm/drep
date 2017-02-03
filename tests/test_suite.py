#!/usr/bin/env python

###############################################################################
#
# test_suite.py - process several E. coli genomes to verify operation of dRep
#
###############################################################################

import glob
import os
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

def ensure_identicle(Swd, wd):
    # Compare datatables
    for d in Swd.data_tables:
        db1 = Swd.get_db(d)
        db2 =  wd.get_db(d)

        assert db1.equals(db2), "{0} is not the same!".format(d)

        # TO DO: http://stackoverflow.com/questions/17095101/outputting-difference-in-two-pandas-dataframes-side-by-side-highlighting-the-d

    # Compare the clustering files
    pass

    # Compare the graphs
    pass

def sanity_check(Swd):
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

        args = argumentParser.parse_args(['dereplicate_wf',wd_loc,'-g'] + genomes)
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        s_wd = WorkDirectory(s_wd_loc)
        wd   = WorkDirectory(wd_loc)
        ensure_identicle(s_wd, wd)

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

        args = argumentParser.parse_args(['filter',wd_loc,'-g',genomes[4]])
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

    def run(self):
        self.setUp()
        self.functional_test_1()
        self.tearDown()

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

    def tearDown(self):
        logging.shutdown()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)

def test_filter():
    ''' test the filter operation '''
    verifyFilter = VerifyFilter()
    verifyFilter.run()

def test_cluster():
    ''' test the cluster operation '''
    verifyCluster = VerifyCluster()
    verifyCluster.run()

def test_dereplicate_wf():
    ''' test the dereplicate_wf operation '''
    verifyDereplicateWf = VerifyDereplicateWf()
    verifyDereplicateWf.run()

if __name__ == '__main__':
    test_dereplicate_wf()
    print("Everything seems to be working swimmingly!")
    #print("Run with py.test")
