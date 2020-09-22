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

class test_filter():
    def __init__(self):
        pass

    def setUp(self):
        self.genomes = test_utils.load_test_genomes()
        self.wd_loc = test_utils.load_test_wd_loc()
        self.s_wd_loc = test_utils.load_solutions_wd()
        self.testdir = test_utils.load_random_test_dir()

        importlib.reload(logging)
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)

    def tearDown(self):
        logging.shutdown()
        if os.path.isdir(self.wd_loc):
            shutil.rmtree(self.wd_loc)

    def run(self):
        self.setUp()
        self.test_calc_genome_info()
        self.tearDown()

        self.setUp()
        self.test_validate_genomeInfo()
        self.tearDown()

        self.setUp()
        self.test_chdb_to_genomeInfo()
        self.tearDown()

        self.setUp()
        self.functional_test_1()
        self.tearDown()

        self.setUp()
        self.functional_test_2()
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
        bdb = drep.d_cluster.utils.load_genomes(genomes)
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
        bdb = drep.d_cluster.utils.load_genomes(genomes)
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
        genome = [g for g in genomes if g.endswith('Enterococcus_faecalis_T2.fna')]
        assert len(genome) == 1
        genome = genome[0]

        args = argumentParser.parse_args(['filter',wd_loc,'-g',genome] \
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

    def functional_test_2(self):
        '''
        Call filter on 'Escherichia_coli_Sakai.fna' with GenomeInfo provivded
        '''
        genomes = self.genomes
        wd_loc  = self.wd_loc

        # make sure calling it on the right genome
        genome = [g for g in genomes if g.endswith('Enterococcus_faecalis_T2.fna')]
        assert len(genome) == 1
        genome = genome[0]

        table = {}
        atts = ['completeness', 'contamination', 'strain_heterogeneity']
        for a in atts:
            table[a] = []
        table['genome'] = []
        table['location'] = []
        for g in [genome]:
            table['genome'].append(os.path.basename(g))
            table['location'].append(g)
            for a in atts:
                table[a].append(10)
        Idb = pd.DataFrame(table)

        if not os.path.isdir(self.testdir):
            os.mkdir(self.testdir)
        GI_loc = os.path.join(self.testdir, 'genomeInfo.csv')
        Idb.to_csv(GI_loc, index=False)

        args = argumentParser.parse_args(['filter',wd_loc,'-g',genome] \
            + ['--genomeInfo', GI_loc])
        controller = Controller()
        controller.parseArguments(args)

        # Confirm Chdb.csv is correct
        wd = drep.WorkDirectory.WorkDirectory(wd_loc)

        # Confirm genome is in Bdb.csv
        Gdb = wd.get_db('genomeInfo')
        assert Gdb['completeness'].tolist()[0] == 10