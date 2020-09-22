import importlib
import logging
import os
import shutil

import tests.test_utils as test_utils

from drep import argumentParser
from drep.controller import Controller
from drep.WorkDirectory import WorkDirectory

class test_dereplicate():
    def __init__(self):
        pass

    def setUp(self):
        self.genomes = test_utils.load_test_genomes()
        self.wd_loc = test_utils.load_test_wd_loc()
        self.s_wd_loc = test_utils.load_solutions_wd()

        importlib.reload(logging)
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)

    def tearDown(self):
        logging.shutdown()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)

    def run(self):
        # self.setUp()
        # self.functional_test_1()
        # self.tearDown()
        #
        # self.setUp()
        # self.functional_test_2()
        # self.tearDown()

        self.setUp()
        self.functional_test_3()
        self.tearDown()

    def functional_test_1(self):
        genomes  = self.genomes
        wd_loc   = self.wd_loc
        s_wd_loc = self.s_wd_loc

        test_utils.sanity_check(WorkDirectory(s_wd_loc))

        args = argumentParser.parse_args(['dereplicate',wd_loc,'-g'] + genomes \
            + ['--checkM_method', 'taxonomy_wf', '--debug', '--S_algorithm',
                        'ANImf'])
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        s_wd = WorkDirectory(s_wd_loc)
        wd   = WorkDirectory(wd_loc)
        test_utils.ensure_identicle(s_wd, wd, skip=['Bdb', 'Mdb'])

        # Perform sanity check to make sure solutions directiory isn't
        # being overwritten
        test_utils.sanity_check(s_wd)

    def functional_test_2(self):
        genomes  = self.genomes
        wd_loc   = self.wd_loc
        s_wd_loc = self.s_wd_loc

        test_utils.sanity_check(WorkDirectory(s_wd_loc))

        args = argumentParser.parse_args(['compare',wd_loc,'-g'] + genomes)
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        s_wd = WorkDirectory(s_wd_loc)
        wd   = WorkDirectory(wd_loc)
        test_utils.ensure_identicle(s_wd, wd, skip=['Bdb', 'Chdb', 'Sdb', 'Wdb', 'Widb',\
            'genomeInformation', 'Mdb'])

        # Perform sanity check to make sure solutions directiory isn't
        # being overwritten
        test_utils.sanity_check(s_wd)

    def functional_test_3(self):
        '''
        Use goANI
        '''
        genomes  = self.genomes
        wd_loc   = self.wd_loc
        s_wd_loc = self.s_wd_loc

        test_utils.sanity_check(WorkDirectory(s_wd_loc))

        args = argumentParser.parse_args(['compare',wd_loc,'--S_algorithm',
                    'goANI','-g'] + genomes)
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        s_wd = WorkDirectory(s_wd_loc)
        wd   = WorkDirectory(wd_loc)
        Ndb = wd.get_db('Ndb')
        assert len(Ndb) > 0

        # Perform sanity check to make sure solutions directiory isn't
        # being overwritten
        test_utils.sanity_check(s_wd)