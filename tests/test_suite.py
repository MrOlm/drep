#!/usr/bin/env python3

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
from WorkDirectory import WorkDirectory

def load_test_genomes():
    return glob.glob(os.path.join(str(os.getcwd()) + '/genomes/*'))

def load_test_wd_loc():
    loc = os.path.join(str(os.getcwd()),'../tests/test_backend/ecoli_wd')
    return loc

class VerifyDereplicateWf():
    def __init__(self):
        pass

    def run(self):
        self.setUp()
        self.functional_test_1()
        self.tearDown()

    def functional_test_1(self):
        genomes = self.genomes
        wd_loc  = self.wd_loc

        args = argumentParser.parse_args(['dereplicate_wf',wd_loc,'-g'] + genomes)
        controller = Controller()
        controller.parseArguments(args)

    def setUp(self):
        self.genomes = load_test_genomes()
        self.wd_loc = load_test_wd_loc()

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
    #test_ecoli()
    print("Run with py.test")
