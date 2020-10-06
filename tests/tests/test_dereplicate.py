import importlib
import logging
import os
import shutil
import pytest

import tests.test_utils as test_utils

from drep import argumentParser
from drep.controller import Controller
from drep.WorkDirectory import WorkDirectory


class Empty():
    pass

@pytest.fixture()
def self():
    # Set up
    self = Empty()
    self.genomes = test_utils.load_test_genomes()
    self.zipped_genomes = test_utils.load_zipped_genomes()
    self.large_genome_set = test_utils.load_large_genome_set()
    self.wd_loc = test_utils.load_test_wd_loc()
    self.wd_loc2 = test_utils.load_test_wd_loc_2()
    self.s_wd_loc = test_utils.load_solutions_wd()

    importlib.reload(logging)
    if os.path.isdir(self.wd_loc):
        shutil.rmtree(self.wd_loc)
    if os.path.isdir(self.wd_loc2):
        shutil.rmtree(self.wd_loc2)

    yield self

    # Teardown
    logging.shutdown()
    if os.path.isdir(self.wd_loc):
        shutil.rmtree(self.wd_loc)
    if os.path.isdir(self.wd_loc2):
        shutil.rmtree(self.wd_loc2)


# class test_dereplicate():
#     def __init__(self):
#         pass
#
#     def setUp(self):
#         self.genomes = test_utils.load_test_genomes()
#         self.zipped_genomes = test_utils.load_zipped_genomes()
#         self.large_genome_set = test_utils.load_large_genome_set()
#         self.wd_loc = test_utils.load_test_wd_loc()
#         self.wd_loc2 = test_utils.load_test_wd_loc_2()
#         self.s_wd_loc = test_utils.load_solutions_wd()
#
#         importlib.reload(logging)
#         if os.path.isdir(self.wd_loc):
#             shutil.rmtree(self.wd_loc)
#         if os.path.isdir(self.wd_loc2):
#             shutil.rmtree(self.wd_loc2)
#
#     def tearDown(self):
#         logging.shutdown()
#         if os.path.isdir(self.wd_loc):
#             shutil.rmtree(self.wd_loc)
#         if os.path.isdir(self.wd_loc2):
#             shutil.rmtree(self.wd_loc2)
#
#     def run(self):
#         self.setUp()
#         self.functional_test_1()
#         self.tearDown()
#
#         self.setUp()
#         self.functional_test_2()
#         self.tearDown()
#
#         self.setUp()
#         self.functional_test_3()
#         self.tearDown()
#
#         self.setUp()
#         self.functional_test_4()
#         self.tearDown()
#
#         self.setUp()
#         self.functional_test_5()
#         self.tearDown()
#
#         self.setUp()
#         self.functional_test_6()
#         self.tearDown()
#
#         self.setUp()
#         self.functional_test_7()
#         self.tearDown()
#
#         self.setUp()
#         self.functional_test_8()
#         self.tearDown()

def test_dereplicate_1(self):
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

def test_dereplicate_2(self):
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

def test_dereplicate_3(self):
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

def test_dereplicate_4(self):
    '''
    Test the ability of primary clustering to take a large genome set and break it into chunks
    '''
    genomes  = self.large_genome_set
    wd_loc   = self.wd_loc
    wd_loc2  = self.wd_loc2

    if len(genomes) == 0:
        print("*** THIS TEST ONLY WORKS ON MO'S DEVELOPMENT MACHINE ***")
        return

    # Get normal results
    args = argumentParser.parse_args(['compare', wd_loc2, '--S_algorithm',
                                      'fastANI', '--SkipSecondary',
                                      '--primary_chunksize', '50', '-g'] + genomes)
    Controller().parseArguments(args)
    wd = WorkDirectory(wd_loc2)
    CSdb = wd.get_db('Cdb')

    # Run with chunking
    args = argumentParser.parse_args(['compare',wd_loc,'--S_algorithm',
                'fastANI','--SkipSecondary', '--multiround_primary_clustering',
                '--primary_chunksize', '50', '-g'] + genomes)
    Controller().parseArguments(args)

    # Verify they're the same
    wd = WorkDirectory(wd_loc)
    Cdb = wd.get_db('Cdb')

    assert len(CSdb) == len(Cdb)
    for c in ['primary_cluster', 'secondary_cluster']:
        assert set(CSdb[c].value_counts().to_dict().values()) == set(Cdb[c].value_counts().to_dict().values())
        assert set(CSdb[c].value_counts().to_dict().keys()) == set(Cdb[c].value_counts().to_dict().keys())
    assert set(CSdb['genome'].tolist()) == set(Cdb['genome'].tolist())
    assert set(Cdb.columns) - set(CSdb.columns) == set(['length', 'subcluster', 'primary_representitive'])

def test_dereplicate_5(self):
    '''
    Test greedy clustering
    '''
    genomes = self.large_genome_set[:10]
    wd_loc = self.wd_loc
    wd_loc2 = self.wd_loc2

    if len(genomes) == 0:
        print("*** THIS TEST ONLY WORKS ON MO'S DEVELOPMENT MACHINE ***")
        return

    # Get greedy results
    args = argumentParser.parse_args(['compare', wd_loc2, '--S_algorithm',
                                      'fastANI', '--SkipMash',
                                      '--clusterAlg', 'greedy', '-sa', '0.95', '-g'] + genomes)
    Controller().parseArguments(args)
    wd = WorkDirectory(wd_loc2)
    CSdb = wd.get_db('Cdb')

    # Run normal
    args = argumentParser.parse_args(['compare', wd_loc, '--S_algorithm',
                                      'fastANI', '--SkipMash',
                                      '-sa', '0.95', '-g'] + genomes)
    Controller().parseArguments(args)

    # Verify they're the same
    wd = WorkDirectory(wd_loc)
    Cdb = wd.get_db('Cdb')

    assert len(CSdb) == len(Cdb)
    for c in ['primary_cluster', 'secondary_cluster']:
        assert set(CSdb[c].value_counts().to_dict().values()) == set(Cdb[c].value_counts().to_dict().values()), c
        assert set(CSdb[c].value_counts().to_dict().keys()) == set(Cdb[c].value_counts().to_dict().keys()), c
    assert set(CSdb['genome'].tolist()) == set(Cdb['genome'].tolist())
    assert set(CSdb.columns) - set(Cdb.columns) == set(['greedy_representative'])

def test_dereplicate_6(self):
    '''
    Test zipped genomes
    '''
    genomes = self.zipped_genomes
    wd_loc = self.wd_loc
    s_wd_loc = self.s_wd_loc

    test_utils.sanity_check(WorkDirectory(s_wd_loc))

    args = argumentParser.parse_args(['compare', wd_loc, '--S_algorithm',
                                      'fastANI', '-g'] + genomes)
    controller = Controller()
    controller.parseArguments(args)

    # Verify
    wd = WorkDirectory(wd_loc)
    anis = wd.get_db('Ndb')['ani'].tolist()
    assert max(anis) <= 1
    assert min(anis) >= 0
    assert len(set(anis)) > 1

def test_dereplicate_7(self):
    '''
    Test greedy clustering fully
    '''
    genomes = self.large_genome_set
    wd_loc = self.wd_loc
    wd_loc2 = self.wd_loc2

    if len(genomes) == 0:
        print("*** THIS TEST ONLY WORKS ON MO'S DEVELOPMENT MACHINE ***")
        return

    # Get greedy results
    args = argumentParser.parse_args(['compare', wd_loc2, '--S_algorithm',
                                      'fastANI', '--multiround_primary_clustering', '--primary_chunksize', '50',
                                      '--clusterAlg', 'greedy', '-sa', '0.95', '-g'] + genomes)
    Controller().parseArguments(args)
    wd = WorkDirectory(wd_loc2)
    CSdb = wd.get_db('Cdb')

    # Run normal
    args = argumentParser.parse_args(['compare', wd_loc, '--S_algorithm',
                                      'fastANI', '--multiround_primary_clustering', '--primary_chunksize', '50',
                                      '-sa', '0.95', '-g'] + genomes)
    Controller().parseArguments(args)

    # Verify they're the same
    wd = WorkDirectory(wd_loc)
    Cdb = wd.get_db('Cdb')

    assert len(CSdb) == len(Cdb)
    for c in ['primary_cluster', 'secondary_cluster']:
        assert set(CSdb[c].value_counts().to_dict().values()) == set(Cdb[c].value_counts().to_dict().values()), c
        if c != 'secondary_cluster':
            assert set(CSdb[c].value_counts().to_dict().keys()) == set(Cdb[c].value_counts().to_dict().keys()), c
    assert set(CSdb['genome'].tolist()) == set(Cdb['genome'].tolist())
    assert set(CSdb.columns) - set(Cdb.columns) == set(['greedy_representative'])

def test_dereplicate_8(self):
    '''
    Test greedy clustering with some primary clusters only having a single member
    '''
    if len(self.large_genome_set) == 0:
        print("*** THIS TEST ONLY WORKS ON MO'S DEVELOPMENT MACHINE ***")
        return


    genomes = [self.large_genome_set[0], self.large_genome_set[20]]
    wd_loc = self.wd_loc
    wd_loc2 = self.wd_loc2

    # Get greedy results
    args = argumentParser.parse_args(['compare', wd_loc2, '--S_algorithm',
                                      'fastANI', '--multiround_primary_clustering', '--primary_chunksize', '50',
                                      '--clusterAlg', 'greedy', '-sa', '0.95', '-pa', '0.99', '-g'] + genomes)
    Controller().parseArguments(args)
    wd = WorkDirectory(wd_loc2)
    CSdb = wd.get_db('Cdb')

    # Run normal
    args = argumentParser.parse_args(['compare', wd_loc, '--S_algorithm',
                                      'fastANI', '--multiround_primary_clustering', '--primary_chunksize', '50',
                                      '-sa', '0.95', '-pa', '0.99', '-g'] + genomes)
    Controller().parseArguments(args)

    # Verify they're the same
    wd = WorkDirectory(wd_loc)
    Cdb = wd.get_db('Cdb')

    assert len(CSdb) == len(Cdb)
    for c in ['primary_cluster', 'secondary_cluster']:
        assert set(CSdb[c].value_counts().to_dict().values()) == set(Cdb[c].value_counts().to_dict().values()), c
        if c != 'secondary_cluster':
            assert set(CSdb[c].value_counts().to_dict().keys()) == set(Cdb[c].value_counts().to_dict().keys())#, [set(CSdb[c].value_counts().to_dict().keys()), set(Cdb[c].value_counts().to_dict().keys())]
    assert set(CSdb['genome'].tolist()) == set(Cdb['genome'].tolist())
    assert set(CSdb.columns) - set(Cdb.columns) == set(['greedy_representative'])
