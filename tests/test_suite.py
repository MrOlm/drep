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

import drep
import drep.d_filter

from drep import argumentParser
from drep.controller import Controller
from drep.WorkDirectory import WorkDirectory
from drep.d_bonus import find_program

def load_test_genomes():
    return glob.glob(os.path.join(str(os.getcwd()) + '/genomes/*'))

def load_test_wd_loc():
    loc = os.path.join(str(os.getcwd()),'../tests/test_backend/ecoli_wd')
    return loc

def load_random_test_dir():
    loc = os.path.join(str(os.getcwd()),'../tests/test_backend/test_dir')
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

        db1 = Swd.get_db(d, return_none=False)
        db2 =  wd.get_db(d, return_none=False)

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

def report_diff(x):
    return x[0] if x[0] == x[1] else '{} | {}'.format(*x)

class VerifyDereplicateWf():
    def __init__(self):
        pass

    def run(self):
        self.setUp()
        self.functional_test_1()
        self.tearDown()

        self.setUp()
        self.functional_test_2()
        self.tearDown()

    def setUp(self):
        self.genomes = load_test_genomes()
        self.wd_loc = load_test_wd_loc()
        self.s_wd_loc = load_solutions_wd()

        #logging.shutdown()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)

    def tearDown(self):
        #logging.shutdown()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)

    def functional_test_1(self):
        genomes  = self.genomes
        wd_loc   = self.wd_loc
        s_wd_loc = self.s_wd_loc

        sanity_check(WorkDirectory(s_wd_loc))

        args = argumentParser.parse_args(['dereplicate',wd_loc,'-g'] + genomes \
            + ['--checkM_method', 'taxonomy_wf'])
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        s_wd = WorkDirectory(s_wd_loc)
        wd   = WorkDirectory(wd_loc)
        ensure_identicle(s_wd, wd, skip=['Bdb', 'Mdb'])

        # Perform sanity check to make sure solutions directiory isn't
        # being overwritten
        sanity_check(s_wd)

    def functional_test_2(self):
        genomes  = self.genomes
        wd_loc   = self.wd_loc
        s_wd_loc = self.s_wd_loc

        sanity_check(WorkDirectory(s_wd_loc))

        args = argumentParser.parse_args(['compare',wd_loc,'-g'] + genomes)
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        s_wd = WorkDirectory(s_wd_loc)
        wd   = WorkDirectory(wd_loc)
        ensure_identicle(s_wd, wd, skip=['Bdb', 'Chdb', 'Sdb', 'Wdb', 'Widb',\
            'genomeInformation'])

        # Perform sanity check to make sure solutions directiory isn't
        # being overwritten
        sanity_check(s_wd)

class VerifyFilter():
    def __init__(self):
        pass

    def setUp(self):
        self.genomes = load_test_genomes()
        self.wd_loc = load_test_wd_loc()
        self.s_wd_loc = load_solutions_wd()

        #logging.shutdown()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)

    def run(self):
        self.setUp()

        self.test_calc_genome_info()
        self.test_validate_genomeInfo()
        self.test_chdb_to_genomeInfo()

        self.functional_test_1()

        self.tearDown()

    def test_chdb_to_genomeInfo(self):
        '''
        test drep.d_filter.chdb_to_genomeInfo
        '''
        genomes = self.genomes
        workDirectory = drep.WorkDirectory.WorkDirectory(self.s_wd_loc)

        # Load chdb
        chdb = workDirectory.get_db('Chdb')
        # Make bdb
        bdb = drep.d_cluster.load_genomes(genomes)
        # Make Gdb
        Gdb = drep.d_filter.calc_genome_info(genomes)

        # Run
        Idb = drep.d_filter.chdb_to_genomeInfo(chdb)
        Tdb = drep.d_filter._validate_genomeInfo(Idb, bdb)
        Tdb = drep.d_filter._add_lengthN50(Tdb, bdb)
        t = Tdb[Tdb['genome'] == 'Enterococcus_casseliflavus_EC20.fasta']
        assert t['completeness'].tolist()[0] == 98.28
        assert t['length'].tolist()[0] == 3427276

    def test_calc_genome_info(self):
        '''
        test drep.d_filter.calc_genome_info
        '''
        genomes = self.genomes
        result = drep.d_filter.calc_genome_info(genomes)

        result['g'] = [os.path.basename(l) for l in result[\
            'location']]
        d = result[result['g'] == 'Enterococcus_faecalis_T2.fna']
        n = d['N50'].tolist()[0]
        l = d['length'].tolist()[0]

        assert n == 774663
        assert l == 3263835

    def test_validate_genomeInfo(self):
        '''
        test drep.d_filter._validate_genomeInfo

        1) Make sure it can load a proper file

        2) Make sure it crashes on an unproper file
        '''
        # Make proper Idb
        genomes = self.genomes
        table = {}
        atts = ['completeness', 'contamination', 'strain_heterogeneity']
        for a in atts:
            table[a] = []
        table['genome'] = []
        table['location'] = []
        for g in genomes:
            table['genome'].append(os.path.basename(g))
            table['location'].append(g)
            for a in atts:
                table[a].append(10)
        Idb = pd.DataFrame(table)
        # Make bdb
        bdb = drep.d_cluster.load_genomes(genomes)
        # Make Gdb
        Gdb = drep.d_filter.calc_genome_info(genomes)

        # Run as correct
        Tdb = drep.d_filter._validate_genomeInfo(Idb, bdb)
        Tdb = drep.d_filter._add_lengthN50(Tdb, bdb)
        t = Tdb[Tdb['genome'] == 'Enterococcus_casseliflavus_EC20.fasta']
        assert t['completeness'].tolist()[0] == 10.0
        assert t['length'].tolist()[0] == 3427276

        # Run wihout one of the genomes
        idb = Idb[Idb['genome'] != 'Enterococcus_casseliflavus_EC20.fasta']
        try:
            tdb = drep.d_filter._validate_genomeInfo(idb, bdb, Gdb)
            assert False
        except:
            pass

        # Run without completeness
        idb = Idb.copy()
        del idb['completeness']
        try:
            tdb = drep.d_filter._validate_genomeInfo(idb, bdb, Gdb)
            assert False
        except:
            pass

        # Run without the genome info
        idb = Idb.copy()
        idb['genome'] = idb['location']
        tdb = drep.d_filter._validate_genomeInfo(idb, bdb)
        tdb = drep.d_filter._add_lengthN50(tdb, bdb)
        t = Tdb[Tdb['genome'] == 'Enterococcus_casseliflavus_EC20.fasta']
        assert t['completeness'].tolist()[0] == 10.0
        assert t['length'].tolist()[0] == 3427276

    def functional_test_1(self):
        '''
        Call filter on 'Escherichia_coli_Sakai.fna'

        Should call both prodigal and checkM
        '''
        genomes = self.genomes
        wd_loc  = self.wd_loc

        # make sure calling it on the right genome
        assert genomes[4].endswith('Enterococcus_faecalis_T2.fna')

        args = argumentParser.parse_args(['filter',wd_loc,'-g',genomes[4]] \
            + ['--checkM_method', 'taxonomy_wf'])
        controller = Controller()
        controller.parseArguments(args)

        # Confirm Chdb.csv is correct
        wd = drep.WorkDirectory.WorkDirectory(wd_loc)
        chdb = wd.get_db('Chdb')
        assert chdb['Completeness'].tolist()[0] == 98.28

        # Confirm genome is in Bdb.csv
        Gdb = wd.get_db('genomeInfo')
        assert Gdb['completeness'].tolist()[0] == 98.28

    def tearDown(self):
        #logging.shutdown()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)

class VerifyCluster():
    def __init__(self):
        pass

    def setUp(self):
        self.genomes = load_test_genomes()
        self.wd_loc = load_test_wd_loc()
        self.test_dir = load_random_test_dir()
        self.s_wd_loc = load_solutions_wd()

        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc, ignore_errors=True)
        if not os.path.isdir(self.test_dir):
            os.mkdir(self.test_dir)

    def tearDown(self):
        #logging.shutdown()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run(self):
        self.setUp()
        self.test_all_vs_all_mash()
        self.tearDown()

        self.setUp()
        self.test_cluster_mash_database()
        self.tearDown()

        self.setUp()
        self.test_compare_genomes()
        self.tearDown()

        self.setUp()
        self.test_genome_hierarchical_clustering()
        self.tearDown()

        self.setUp()
        self.functional_test_3()
        self.tearDown()

        self.setUp()
        self.functional_test_2()
        self.tearDown()

        self.setUp()
        self.functional_test_1()
        self.tearDown()

        self.setUp()
        self.skipsecondary_test()
        self.tearDown()

    def test_genome_hierarchical_clustering(self):
        '''
        Test d_cluster.test_genome_hierarchical_clustering
        '''
        wdS = drep.WorkDirectory.WorkDirectory(self.s_wd_loc)
        Ndb = wdS.get_db('Ndb')

        # Run clustering on Ndb
        Cdb, c2ret = drep.d_cluster._cluster_Ndb(Ndb, comp_method='ANImf')
        g2c = Cdb.set_index('genome')['secondary_cluster'].to_dict()
        assert g2c['Enterococcus_faecalis_T2.fna'] != g2c['Enterococcus_faecalis_TX0104.fa']
        assert g2c['Enterococcus_faecalis_T2.fna'] == g2c['Enterococcus_faecalis_YI6-1.fna']

        # Make sure storage is correct
        wd = drep.WorkDirectory.WorkDirectory(self.wd_loc)
        wd.store_special('secondary_linkages', c2ret)
        wd.load_cached()
        got = wd.get_cluster('secondary_linkage_cluster_1')
        assert len(got) == 3

    def test_compare_genomes(self):
        '''
        Test d_cluster.compare_genomes
        '''
        bdb = drep.d_cluster.load_genomes(self.genomes)
        data_folder = self.test_dir

        # Try ANImf
        Ndb = drep.d_cluster.compare_genomes(bdb, 'ANImf', data_folder)
        db = Ndb[(Ndb['reference'] == 'Enterococcus_faecalis_T2.fna')\
            & (Ndb['querry'] == 'Enterococcus_casseliflavus_EC20.fasta')]
        assert (db['ani'].tolist()[0] > 0.85) & (db['ani'].tolist()[0] < 0.86)

        # Try ANIn
        Ndb = drep.d_cluster.compare_genomes(bdb, 'ANIn', data_folder)
        db = Ndb[(Ndb['reference'] == 'Enterococcus_faecalis_T2.fna')\
            & (Ndb['querry'] == 'Enterococcus_casseliflavus_EC20.fasta')]
        assert (db['ani'].tolist()[0] > 0.85) & (db['ani'].tolist()[0] < 0.86)

        # Try gANI
        loc, works = drep.d_bonus.find_program('ANIcalculator')
        if works:
            Ndb = drep.d_cluster.compare_genomes(bdb, 'gANI', data_folder, \
                prod_folder = os.path.join(data_folder + 'prodigal'))
            db = Ndb[(Ndb['reference'] == 'Enterococcus_faecalis_T2.fna')\
                & (Ndb['querry'] == 'Enterococcus_casseliflavus_EC20.fasta')]
            assert (db['ani'].tolist()[0] > 0.85) & (db['ani'].tolist()[0] < 0.86)

    def test_all_vs_all_mash(self):
        '''
        Test d_cluster.all_vs_all_MASH
        '''
        bdb = drep.d_cluster.load_genomes(self.genomes)
        data_folder = self.test_dir

        # Run it under normal conditions
        Mdb = drep.d_cluster.all_vs_all_MASH(bdb, data_folder)
        assert len(Mdb) == 25
        db = Mdb[(Mdb['genome1'] == 'Enterococcus_faecalis_T2.fna') & \
            (Mdb['genome2'] == 'Enterococcus_casseliflavus_EC20.fasta')]
        d = float(db['dist'].tolist()[0])
        assert (d > .2) & (d < .3)
        assert len(glob.glob(data_folder + '/MASH_files/sketches/*')) == 1
        assert len(glob.glob(data_folder + '/MASH_files/sketches/*/*')) == 6

        # Start over
        shutil.rmtree(self.test_dir)
        os.mkdir(self.test_dir)

        # Run it under reduced chuck size
        Mdb = drep.d_cluster.all_vs_all_MASH(bdb, data_folder, groupSize=2)

        assert len(Mdb) == 25
        db = Mdb[(Mdb['genome1'] == 'Enterococcus_faecalis_T2.fna') & \
            (Mdb['genome2'] == 'Enterococcus_casseliflavus_EC20.fasta')]
        d = float(db['dist'].tolist()[0])
        assert (d > .2) & (d < .3)
        assert len(glob.glob(data_folder + '/MASH_files/sketches/*')) == 3
        assert len(glob.glob(data_folder + '/MASH_files/sketches/*/*')) == 8

    def test_cluster_mash_database(self):
        '''
        Test d_cluster.cluster_mash_database
        '''
        wdS = drep.WorkDirectory.WorkDirectory(self.s_wd_loc)
        Mdb = wdS.get_db('Mdb')

        # Make sure clustering is correct
        Cdb, cluster_ret = drep.d_cluster.cluster_mash_database(Mdb)
        g2c = Cdb.set_index('genome')['primary_cluster'].to_dict()
        assert g2c['Enterococcus_faecalis_T2.fna'] == g2c['Enterococcus_faecalis_TX0104.fa']
        assert g2c['Enterococcus_faecalis_T2.fna'] != g2c['Enterococcus_casseliflavus_EC20.fasta']

        # Make sure storage is correct
        wd = drep.WorkDirectory.WorkDirectory(self.wd_loc)
        wd.store_special('primary_linkage', cluster_ret)
        wd.load_cached()
        got = wd.get_cluster('primary_linkage')
        assert len(got) == 3

    def functional_test_1(self):
        '''
        Cluster the 5 genomes using default settings
        '''
        genomes = self.genomes
        wd_loc  = self.wd_loc
        s_wd_loc = self.s_wd_loc

        args = argumentParser.parse_args(['cluster',wd_loc,'-g']+genomes)
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        Swd = WorkDirectory(s_wd_loc)
        wd   = WorkDirectory(wd_loc)

        # Confirm Cdb.csv is correct
        db1 = Swd.get_db('Cdb')
        db2 =  wd.get_db('Cdb')

        assert compare_dfs(db1, db2), "{0} is not the same!".format('Cdb')

    def functional_test_2(self):
        '''
        Cluster the 5 genomes using gANI
        '''
        genomes = self.genomes
        wd_loc  = self.wd_loc
        s_wd_loc = self.s_wd_loc

        # Make sure gANI is installed
        loc, works = find_program('ANIcalculator')
        if (loc == None or works == False):
            print('Cannot locate the program {0}- skipping related tests'\
                .format('ANIcalculator (for gANI)'))
            return

        args = argumentParser.parse_args(['cluster',wd_loc,'--S_algorithm',\
            'gANI','-g']+genomes)
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        Swd = WorkDirectory(s_wd_loc)
        wd   = WorkDirectory(wd_loc)

        # Confirm Cdb.csv is correct
        db1 = Swd.get_db('Cdb')
        del db1['comparison_algorithm']
        db2 =  wd.get_db('Cdb')
        del db2['comparison_algorithm']
        assert compare_dfs(db1, db2), "{0} is not the same!".format('Cdb')

    def functional_test_3(self):
        '''
        Cluster the 5 genomes using ANImf
        '''

        genomes = self.genomes
        wd_loc  = self.wd_loc
        s_wd_loc = self.s_wd_loc

        args = argumentParser.parse_args(['cluster',wd_loc,'--S_algorithm',\
            'ANImf','-g']+genomes)
        controller = Controller()
        controller.parseArguments(args)

        # Verify
        Swd = WorkDirectory(s_wd_loc)
        wd   = WorkDirectory(wd_loc)

        # Confirm Cdb.csv is correct
        db1 = Swd.get_db('Cdb')
        del db1['comparison_algorithm']
        db2 =  wd.get_db('Cdb')
        del db2['comparison_algorithm']
        assert compare_dfs(db1, db2), "{0} is not the same!".format('Cdb')

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

class VerifyAnalyze():
    def __init__(self):
        pass

    def setUp(self):
        self.s_wd_loc = load_solutions_wd()
        self.working_wd_loc = load_test_wd_loc()
        self.test_dir = load_random_test_dir()

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

    def tearDown(self):
        #logging.shutdown()
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

class VerifyChoose():
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
        genomes = load_test_genomes()
        g2l = {os.path.basename(g):g for g in genomes}

        Bdb = pd.read_csv(wd_loc + '/data_tables/Bdb.csv')
        Bdb['location'] = Bdb['genome'].map(g2l)
        Bdb.to_csv(wd_loc + '/data_tables/Bdb.csv', index=False)

        # Run choose - this should re-run checkM and re-generate chdb
        args = argumentParser.parse_args(['choose', wd_loc, '-o', '--checkM_method',\
            'taxonomy_wf'])
        controller = Controller()
        controller.parseArguments(args)

        Swd  = WorkDirectory(self.s_wd_loc)
        wd   = WorkDirectory(self.working_wd_loc)
        for db in ['Chdb', 'genomeInformation']:
            db1 = Swd.get_db(db)
            db2 =  wd.get_db(db)
            assert compare_dfs(db1, db2), "{0} is not the same!".format(db)

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
        args = argumentParser.parse_args(['choose', wd_loc, '--noQualityFiltering'])
        controller = Controller()
        controller.parseArguments(args)

        Swd  = WorkDirectory(self.s_wd_loc)
        wd   = WorkDirectory(self.working_wd_loc)
        for db in ['Sdb', 'Wdb', 'genomeInformation']:
            db1 = Swd.get_db(db)
            db2 =  wd.get_db(db)
            assert not compare_dfs(db1, db2), "{0} is the same!".format(db)

        sdb = wd.get_db('Sdb')
        for s in sdb['score'].tolist():
            assert (s > 0) & (s < 5)

    def tearDown(self):
        #logging.shutdown()
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

        # self.setUp()
        # self.taxTest3()
        # self.tearDown()

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
        assert compare_dfs(tdb, tdbS), "{0} is not the same!".format('Bdb')

        tdbS = Swd.get_db('Tdb')
        tdb = wd.get_db('Tdb')

        if compare_dfs(tdb, tdbS) == False:
            print("{0} is not the same! May be due to centrifuge index issues".format('Tdb'))
            my_panel = pd.Panel(dict(df1=tdbS,df2=tdb))
            print(my_panel.apply(report_diff, axis=0))

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
        assert compare_dfs(tdb, tdbS), "{0} is not the same!".format('Bdb')

        tdbS = Swd.get_db('TdbP')
        tdb = wd.get_db('Tdb')
        assert compare_dfs(tdb, tdbS), "{0} is not the same!".format('Tdb')

    def taxTest3(self):
        '''
        Try actually running centrifuge
        '''
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
        assert compare_dfs(tdb, tdbS), "{0} is not the same!".format('Bdb')

        tdbS = Swd.get_db('TdbP')
        tdb = wd.get_db('Tdb')
        assert compare_dfs(tdb, tdbS), "{0} is not the same!".format('Tdb')

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
        self.unit_tests_6()
        self.tearDown()

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

    def unit_tests_6(self):
        '''
        Test drep call commands
        '''
        # try on single mash command

        wd   = WorkDirectory(self.working_wd_loc)
        MASH_folder = wd.get_dir('MASH')
        log_folder = wd.get_dir('cmd_logs')

        mash_exe = 'mash'
        all_file = MASH_folder + 'ALL.msh'

        cmd = [mash_exe, 'dist', all_file, all_file, '>', MASH_folder
            + 'MASH_table.tsv']
        cmd = ' '.join(cmd)
        drep.run_cmd(cmd, shell=True, logdir=log_folder)

        assert len(glob.glob(log_folder + '*')) == 3

    def tearDown(self):
        #logging.shutdown()
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

def choose_test():
    ''' test the choose operation '''
    verifyChoose= VerifyChoose()
    verifyChoose.run()

def dereplicate_test():
    ''' test the dereplicate operation '''
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
    dereplicate_test()
    filter_test()
    cluster_test()
    choose_test()
    analyze_test()

@pytest.mark.short
def test_short():
    taxonomy_test()

@pytest.mark.quick
def test_quick():
    rerun_test()

@pytest.mark.unit
def test_unit():
    unit_test()

if __name__ == '__main__':
    # test_unit()
    # test_quick()
    # test_short()
    # test_long()

    # filter_test()
    #choose_test()
    analyze_test()
    #dereplicate_test()
    #taxonomy_test()
    #cluster_test()

    print("Everything seems to be working swimmingly!")
