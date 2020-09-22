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
from drep.d_bonus import find_program

class test_cluster():
    def __init__(self):
        pass

    def setUp(self):
        self.genomes = test_utils.load_test_genomes()
        self.broken_genome = test_utils.load_broken_genome()
        self.wd_loc = test_utils.load_test_wd_loc()
        self.test_dir = test_utils.load_random_test_dir()
        self.s_wd_loc = test_utils.load_solutions_wd()

        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc, ignore_errors=True)
        if not os.path.isdir(self.test_dir):
            os.mkdir(self.test_dir)

        importlib.reload(logging)

    def tearDown(self):
        logging.shutdown()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run(self):
        # self.setUp()
        # self.test_all_vs_all_mash()
        # self.tearDown()
        #
        # self.setUp()
        # self.test_list_genome_load()
        # self.tearDown()
        #
        # self.setUp()
        # self.test_all_vs_all_mash()
        # self.tearDown()
        #
        # self.setUp()
        # self.test_cluster_mash_database()
        # self.tearDown()
        #
        # # self.setUp()
        # # self.time_compare_genomes()
        # # self.tearDown()
        #
        # self.setUp()
        # self.test_goANI()
        # self.tearDown()
        #
        # self.setUp()
        # self.test_goANI2()
        # self.tearDown()
        #
        # self.setUp()
        # self.test_fastANI()
        # self.tearDown()
        #
        # self.setUp()
        # self.test_compare_genomes()
        # self.tearDown()

        self.setUp()
        self.test_genome_hierarchical_clustering()
        self.tearDown()

        self.setUp()
        self.functional_test_4()
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

    def test_list_genome_load(self):
        '''
        Test inputing a list of genomes via a text file
        '''
        bdb = drep.d_cluster.utils.load_genomes(self.genomes)
        data_folder = self.test_dir

        # Make the list of genomes
        if not os.path.exists(data_folder):
            os.mkdir(data_folder)
        genome_loc = os.path.join(data_folder, 'genomes.txt')
        with open(genome_loc, 'w') as o:
            for i, row in bdb.iterrows():
                o.write(row['location'] + '\n')

        # Test it out
        wd_loc  = self.wd_loc
        s_wd_loc = self.s_wd_loc

        args = argumentParser.parse_args(['cluster',wd_loc,'--S_algorithm',\
            'fastANI','-g',genome_loc])
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
        assert test_utils.compare_dfs(db1, db2), "{0} is not the same!".format('Cdb')

        Ndb = drep.d_cluster.compare_utils.compare_genomes(bdb, 'fastANI', data_folder)
        db = Ndb[(Ndb['reference'] == 'Enterococcus_faecalis_T2.fna')\
            & (Ndb['querry'] == 'Enterococcus_casseliflavus_EC20.fasta')]

        assert (db['ani'].tolist()[0] > 0.7) & (db['ani'].tolist()[0] < 0.8)

    def test_genome_hierarchical_clustering(self):
        '''
        Test d_cluster.test_genome_hierarchical_clustering
        '''
        wdS = drep.WorkDirectory.WorkDirectory(self.s_wd_loc)
        Ndb = wdS.get_db('Ndb')

        # Run clustering on Ndb
        Cdb, c2ret = drep.d_cluster.utils._cluster_Ndb(Ndb, comp_method='ANImf')
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
        bdb = drep.d_cluster.utils.load_genomes(self.genomes)
        data_folder = self.test_dir

        # Try gANI
        loc, works = drep.d_bonus.find_program('ANIcalculator')
        if works:
            p_folder = os.path.join(data_folder, 'prodigal')
            #print(p_folder)
            Ndb = drep.d_cluster.compare_utils.compare_genomes(bdb, 'gANI', data_folder, \
                                                               prod_folder = p_folder)
            db = Ndb[(Ndb['reference'] == 'Enterococcus_faecalis_T2.fna')\
                & (Ndb['querry'] == 'Enterococcus_casseliflavus_EC20.fasta')]

            assert (db['ani'].tolist()[0] > 0.7) & (db['ani'].tolist()[0] < 0.75)

        # Try ANImf
        Ndb = drep.d_cluster.compare_utils.compare_genomes(bdb, 'ANImf', data_folder)
        db = Ndb[(Ndb['reference'] == 'Enterococcus_faecalis_T2.fna')\
            & (Ndb['querry'] == 'Enterococcus_casseliflavus_EC20.fasta')]
        assert (db['ani'].tolist()[0] > 0.85) & (db['ani'].tolist()[0] < 0.86)

        # Try ANIn
        Ndb = drep.d_cluster.compare_utils.compare_genomes(bdb, 'ANIn', data_folder)
        db = Ndb[(Ndb['reference'] == 'Enterococcus_faecalis_T2.fna')\
            & (Ndb['querry'] == 'Enterococcus_casseliflavus_EC20.fasta')]
        assert (db['ani'].tolist()[0] > 0.85) & (db['ani'].tolist()[0] < 0.86)

    def test_goANI(self):
        '''
        Test goANI
        '''
        import time

        bdb = drep.d_cluster.utils.load_genomes(self.genomes)
        data_folder = self.test_dir

        # Copy over prodigal
        self.s_wd_loc = test_utils.load_solutions_wd()
        p_folder = os.path.join(data_folder, 'data/prodigal/')
        shutil.copytree(os.path.join(self.s_wd_loc, 'data/prodigal'), \
            p_folder)

        # Try goANI
        p_folder = os.path.join(data_folder, 'data/prodigal/')
        Ndb = drep.d_cluster.compare_utils.compare_genomes(bdb, 'goANI', data_folder, \
                                                           prod_folder = p_folder)
        db = Ndb[(Ndb['reference'] == 'Enterococcus_faecalis_T2.fna')\
            & (Ndb['querry'] == 'Enterococcus_casseliflavus_EC20.fasta')]

        assert (db['ani'].tolist()[0] > 0.7) & (db['ani'].tolist()[0] < 0.8)

    def test_goANI2(self):
        '''
        Test goANI in the case where the genomes share no genes
        '''
        import time

        bdb = drep.d_cluster.utils.load_genomes(self.genomes)
        data_folder = self.test_dir

        # Copy over prodigal
        self.s_wd_loc = test_utils.load_solutions_wd()
        p_folder = os.path.join(data_folder, 'data/prodigal/')
        shutil.copytree(os.path.join(self.s_wd_loc, 'data/prodigal'), \
            p_folder)

        # Remove all but one gene in one of the prodigal files
        p_folder = os.path.join(data_folder, 'data/prodigal/')
        for f in glob.glob(p_folder + '*'):
            if 'Escherichia_coli_Sakai.fna.fna' in f:
                new_file = open(f + '.2', 'w')
                old_file = open(f, 'r')
                j = 0
                for line in old_file.readlines():
                    if ((line[0] == '>') & (j != 0)):
                        break
                    j += 1
                    new_file.write(line.strip() + '/')
                new_file.close()
                old_file.close()
                os.remove(f)
                shutil.copy(f + '.2', f)

        # Try goANI
        p_folder = os.path.join(data_folder, 'data/prodigal/')
        Ndb = drep.d_cluster.compare_utils.compare_genomes(bdb, 'goANI', data_folder, \
                                                           prod_folder = p_folder)
        db = Ndb[(Ndb['reference'] == 'Enterococcus_faecalis_T2.fna')\
            & (Ndb['querry'] == 'Enterococcus_casseliflavus_EC20.fasta')]

        assert (db['ani'].tolist()[0] > 0.7) & (db['ani'].tolist()[0] < 0.8)

    def test_fastANI(self):
        '''
        Test fastANI
        '''
        bdb = drep.d_cluster.utils.load_genomes(self.genomes)
        data_folder = self.test_dir

        Ndb = drep.d_cluster.compare_utils.compare_genomes(bdb, 'fastANI', data_folder)
        db = Ndb[(Ndb['reference'] == 'Enterococcus_faecalis_T2.fna')\
            & (Ndb['querry'] == 'Enterococcus_casseliflavus_EC20.fasta')]

        assert (db['ani'].tolist()[0] > 0.7) & (db['ani'].tolist()[0] < 0.8)

    def time_compare_genomes(self):
        '''
        Time d_cluster.compare_genomes
        '''
        import time

        bdb = drep.d_cluster.utils.load_genomes(self.genomes)
        data_folder = self.test_dir

        for method in ['fastANI', 'ANIn', 'ANImf']:
            # Try ANImf
            start = time.time()
            Ndb = drep.d_cluster.compare_utils.compare_genomes(bdb, method, data_folder, processors=1)
            db = Ndb[(Ndb['reference'] == 'Enterococcus_faecalis_T2.fna')\
                & (Ndb['querry'] == 'Enterococcus_casseliflavus_EC20.fasta')]
            assert (db['ani'].tolist()[0] > 0.7) & (db['ani'].tolist()[0] < 0.9)
            end = time.time()

            comps = len(bdb) * len(bdb)
            print("{1} time: {0:.2f} seconds for {2} comparisons ({3:.2f} seconds per comparison)".format(end-start, method, comps, (end-start)/comps))

    def test_all_vs_all_mash(self):
        '''
        Test d_cluster.all_vs_all_MASH
        '''
        bdb = drep.d_cluster.utils.load_genomes(self.genomes)
        bdb = drep.d_filter._add_lengthN50(bdb, bdb)
        data_folder = self.test_dir

        # Run it under normal conditions
        Mdb, Cdb, cluster_ret = drep.d_cluster.compare_utils.all_vs_all_MASH(bdb, data_folder)
        assert len(Mdb) == 25
        db = Mdb[(Mdb['genome1'] == 'Enterococcus_faecalis_YI6-1.fna') & \
            (Mdb['genome2'] == 'Enterococcus_faecalis_TX0104.fa')]

        d = float(db['dist'].tolist()[0])
        assert (d > .01) & (d < .02)
        assert len(glob.glob(data_folder + '/MASH_files/sketches/*')) == 1
        assert len(glob.glob(data_folder + '/MASH_files/sketches/*/*')) == 6

        # Start over
        shutil.rmtree(self.test_dir)
        os.mkdir(self.test_dir)

        # Run it under reduced chuck size
        Mdb, Cdb, cluster_ret = drep.d_cluster.compare_utils.all_vs_all_MASH(bdb, data_folder, primary_chunksize=2, multiround_primary_clustering=True)

        assert len(Mdb) == 25
        db = Mdb[(Mdb['genome1'] == 'Enterococcus_faecalis_YI6-1.fna') & \
            (Mdb['genome2'] == 'Enterococcus_faecalis_TX0104.fa')]
        d = float(db['dist'].tolist()[0])
        assert (d > .01) & (d < .02)
        assert len(glob.glob(data_folder + '/MASH_files/sketches/*')) == 3
        assert len(glob.glob(data_folder + '/MASH_files/sketches/*/*')) == 8

    def test_cluster_mash_database(self):
        '''
        Test d_cluster.cluster_mash_database
        '''
        wdS = drep.WorkDirectory.WorkDirectory(self.s_wd_loc)
        Mdb = wdS.get_db('Mdb')

        # Make sure clustering is correct
        Cdb, cluster_ret = drep.d_cluster.compare_utils.cluster_mash_database(Mdb)
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

        assert test_utils.compare_dfs(db1, db2), "{0} is not the same!".format('Cdb')

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
        assert test_utils.compare_dfs(db1, db2), "{0} is not the same!".format('Cdb')

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
        assert test_utils.compare_dfs(db1, db2), "{0} is not the same!".format('Cdb')

    def functional_test_4(self):
        '''
        Cluster the 5 genomes using fastANI
        '''

        genomes = self.genomes
        wd_loc  = self.wd_loc
        s_wd_loc = self.s_wd_loc

        args = argumentParser.parse_args(['cluster',wd_loc,'--S_algorithm',\
            'fastANI','-g']+genomes)
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
        assert test_utils.compare_dfs(db1, db2), "{0} is not the same!".format('Cdb')

    def skipsecondary_test(self):
        genomes = self.genomes
        wd_loc  = self.wd_loc
        s_wd_loc = self.s_wd_loc

        args = argumentParser.parse_args(['cluster',wd_loc,'-g'] +genomes \
                + ['--SkipSecondary'])
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