#!/usr/bin/env python3

import logging
import glob
import pandas as pd
import os
import sys
import shutil
import multiprocessing

import drep
import drep.WorkDirectory
import drep.d_cluster
import drep.d_choose

def d_filter_wrapper(wd, **kwargs):
    '''
    Controller for the dRep filter operation

    Args:
        wd (WorkDirectory): The current workDirectory
        **kwargs: Command line arguments

    Keyword Args:
        genomes: genomes to filter in .fasta format
        genomeInfo: location of .csv file with the columns: ["genome"(basename of .fasta file of that genome), "completeness"(0-100 value for completeness of the genome), "contamination"(0-100 value of the contamination of the genome)]

        processors: Threads to use with checkM / prodigal
        overwrite: Overwrite existing data in the work folder

        length: minimum genome length when filtering
        completeness: minimum genome completeness when filtering
        contamination: maximum genome contamination when filtering

        skipCheckM: run checkM on the genome list to calculate comp and con
        checkM_method: Either lineage_wf (more accurate) or taxonomy_wf (faster)

    Returns:
        Bdb.csv: A dataframe of filtered genomes in the workDirectory
        Chdb.csv: A dataframe of raw checkM results in the workDirectory
        GenomeInfo.csv: A dataframe of genome information in the workDirectory

    '''
    # Load the WorkDirectory.
    logging.debug("Loading work directory in filter")
    workDirectory = drep.WorkDirectory.WorkDirectory(wd)
    logging.debug(str(workDirectory))

    # Validate arguments
    logging.debug("Validating filter arguments")
    _validate_arguments(workDirectory, **kwargs)

    # Get the thing to filter - bdb (columns = genome, location)
    if wd.hasDb('Bdb'):
        logging.info("Will filter Bdb")
        bdb = wd.get_db('Bdb')
    else:
        logging.info("Will filter the genome list")
        bdb = drep.d_cluster.load_genomes(kwargs['genomes'])

    # Calculate the length and N50 of all genomes
    GenomeInfo = calc_genome_info(bdb)

    # Filter by size if need be
    if kwargs.get('length', 0) > 1:
        logging.debug("Filtering genomes by size")
        bdb = filter_bdb_length(bdb, kwargs['length'], verbose=True)

    # Run checkM if need be
    #if (not kwargs.get('skipCheckM',False)):
    if (not kwargs.get('skipCheckM',False)) & (kwargs.get('')):

        # Run checkM
        if kwargs.get('Chdb',None) != None:
            logging.debug("Loading provided CheckM data...")
            Chdb = kwargs.get('Chdb')
            Chdb = pd.read_csv(Chdb)
            validate_chdb(Chdb, bdb)

        elif workDirectory.hasDb('Chdb'):
            logging.debug("Loading CheckM data from work directory...")
            Chdb = workDirectory.get_db('Chdb')
            validate_chdb(Chdb, bdb)

        else:
            logging.debug("Running CheckM")
            Chdb = run_checkM_wrapper(bdb, workDirectory, **kwargs)


        # Save genome information db
        pass

    # Filter
        # Filter bdb
        bdb = filter_bdb(bdb, Chdb, **kwargs)

    # Save Bdb
    loc = workDirectory.location + '/data_tables/'
    if saveAs == 'Bdb':
        bdb.to_csv(loc + 'Bdb.csv',index=False)
        workDirectory.data_tables['Bdb'] = bdb

def calc_genome_info(genomes: list):
    '''
    Calculate the length and N50 of a list of genome locations

    Args:
        genomes: list of locations of genomes

    Returns:
        DataFrame: pandas dataframe with ["location", "length", "N50"]
    '''
    table = pd.DataFrame()
    for att in ['location', 'length', 'N50']:
        table[att] = []

    for loc in genomes:
        table['location'].append(loc)
        table['length'].append()

    return 'poopy'


def filter_bdb(bdb, chdb, **kwargs):
    min_comp = kwargs.get('completeness',False)
    max_con = kwargs.get('contamination',False)
    # min_strain_htr = kwargs.get('strain_htr',False)
    start_genomes = list(bdb['genome'].unique())
    assert len(start_genomes) > 0

    db = chdb.copy()

    if min_comp != False:
        db = db[(db['Completeness'] >= min_comp)]
    if max_con != False:
        db = db[db['Contamination'] <= max_con]
    # if min_strain_htr != False:
    #     db = db[db['Strain heterogeneity'] <= min_strain_htr]
    keep_genomes = list(db['Bin Id'].unique())
    bdb = bdb[bdb['genome'].isin(keep_genomes)]

    logging.info("{0:.2f}% of genomes passed checkM filtering".format((len(keep_genomes)/len(start_genomes))*100))

    return bdb

def filter_bdb_length(bdb, min_length, verbose=False):
    start = len(bdb['location'].unique())
    bdb['length'] =  bdb['location'].map(drep.fasta_length)
    x = bdb[bdb['length'] >= min_length]
    end = len(x['location'].unique())

    logging.info("{0:.2f}% of genomes passed length filtering".format((end/start)*100))

    return x

def validate_chdb(Chdb, bdb):
    quit = False
    b_genomes = bdb['genome'].tolist()
    for genome in b_genomes:
        if genome not in Chdb['Bin Id'].tolist():
            logging.error("{0} is not in checkM db".format(genome))
            quit = True
    if quit:
        logging.error("New checkM db needs to be made")
        sys.exit()
    return

def _validate_arguments(wd,**kwargs):
    '''
    Validate arguments to the filter wrapper

    1) Make sure you either have a genome list or Bdb to filter. Bdb is the rarer
    case, but can be done. Throw warnings if needed.

    Args:
        wd (WorkDirectory): The current workDirectory
        **kwargs: Command line arguments

    Returns:
        Nothing
    '''

    # Figure out what you're going to filter
    if wd.hasDb('Wdb'):
        logging.warning("NOTE: Wdb already exists! This will not be filtered! Be sure you know what you're doing")
    if wd.hasDb('Bdb'):
        if kwargs.get('genomes',None) != None:
            logging.error("Both Bdb and a genome list are found- either don't include "\
                    + "a genome list or start a new work directory!")
            sys.exit()
        if wd.hasDb('Cdb'):
            logging.warning("NOTE: Clustering already exists! This will not be filtered! Be sure you know what you're doing")
    else:
        if kwargs.get('genomes',None) == None:
            logging.error("I don't have anything to filter! Give me a genome list")
            sys.exit()

    # # Validate checkM stuff
    # if kwargs.get('genomeInfo', None) != None:
    #

    return

def run_checkM_wrapper(bdb, workDirectory, **kwargs):
    # Make the folder to house the checkM info
    checkM_loc = workDirectory.location + '/data/checkM/'
    drep.make_dir(checkM_loc,dry=kwargs.get('dry',False),\
                overwrite=kwargs.get('overwrite',False))

    # Run prodigal
    prod_folder = workDirectory.location + '/data/prodigal/'
    if not os.path.exists(prod_folder):
        os.makedirs(prod_folder)
    logging.info("Running prodigal")
    run_prodigal(bdb, prod_folder, wd=workDirectory, **kwargs)

    # Run checkM
    checkM_outfolder = checkM_loc + 'checkM_outdir/'
    Chdb = run_checkM(prod_folder, checkM_outfolder, wd=workDirectory, **kwargs)
    validate_chdb(Chdb, bdb)

    # Fix genome size and N50 of Chdb
    Chdb = fix_chdb(Chdb,bdb)

    # Save checkM run
    workDirectory.store_db(Chdb,'Chdb')

    return Chdb

def fix_chdb(Chdb, Bdb):
    g2s = {}
    g2n = {}
    for genome in Chdb['Bin Id'].unique():
        loc = Bdb['location'][Bdb['genome'] == genome].tolist()[0]
        g2s[genome] = drep.fasta_length(loc)
        g2n[genome] = calc_n50(loc)
    Chdb['Genome size (bp)'] = Chdb['Bin Id'].map(g2s)
    Chdb['N50 (scaffolds)'] = Chdb['Bin Id'].map(g2n)

    return Chdb

def calc_n50(loc):
    from Bio import SeqIO

    lengths = []
    for seq_record in SeqIO.parse(loc, "fasta"):
        lengths.append(len(seq_record))

    half = sum(lengths)/2
    tally = 0
    for l in sorted(lengths, reverse=True):
        tally += l
        if tally > half:
            return l

def run_prodigal(bdb, out_dir, **kwargs):
    t = kwargs.get('processors','6')
    loc = shutil.which('prodigal')
    if loc == None:
        logging.error('Cannot locate the program {0}- make sure its in the system path'\
            .format('prodigal'))
        sys.exit()

    cmds = []
    for genome in bdb['location'].unique():
        fna = "{0}{1}{2}".format(out_dir,os.path.basename(genome),'.fna')
        faa = "{0}{1}{2}".format(out_dir,os.path.basename(genome),'.faa')
        if os.path.exists(fna) and os.path.exists(faa):
            pass
        else:
            cmds.append(['prodigal','-i',genome,'-d',fna,'-a',faa,'-m','-p','meta'])

    if len(cmds) > 0:
        if 'wd' in kwargs:
            logdir = kwargs.get('wd').get_dir('cmd_logs')
        else:
            logdir = False
        drep.thread_cmds(cmds, shell=False, logdir=logdir, t=int(t))

    else:
        logging.info("Past prodigal runs found- will not re-run")

def run_checkM(genome_folder,checkm_outf,**kwargs):
    import drep.d_bonus as dBonus
    t = str(kwargs.get('processors','6'))
    loc, works = dBonus.find_program('checkm')
    if loc == None:
        logging.error('Cannot locate the program {0}- make sure its in the system path'\
            .format('checkm'))
        sys.exit()
    if works == False:
        logging.error('Program {0} is not working!! Im going to crash now'\
            .format('checkm'))
        sys.exit()
    check_exe = loc

    checkm_method = kwargs.get('checkM_method','lineage_wf')

    # Run checkM initial
    if checkm_method == 'taxonomy_wf':
         cmd = [check_exe,checkm_method,'domain','Bacteria',genome_folder,checkm_outf,'-f',\
            checkm_outf + '/results.tsv','--tab_table','-t',str(t),'-g','-x','faa']
    else:
         cmd = [check_exe,checkm_method,genome_folder,checkm_outf,'-f',\
            checkm_outf + '/results.tsv','--tab_table','-t',str(t),'--pplacer_threads',\
            str(t),'-g','-x','faa']

    logging.debug("Running CheckM with command: {0}".format(cmd))

    if 'wd' in kwargs:
        logdir = kwargs.get('wd').get_dir('cmd_logs')
    else:
        logdir = False
    drep.run_cmd(cmd, shell=False, logdir=logdir)

    # Run checkM again for the better table
    if checkm_method == 'taxonomy_wf':
        lineage = checkm_outf + 'Bacteria.ms'
    else:
        lineage = checkm_outf + 'lineage.ms'
    desired_file = checkm_outf + 'Chdb.tsv'
    cmd = [check_exe,'qa', lineage, checkm_outf, '-f', desired_file, '-t',\
            str(t), '--tab_table','-o', '2']
    logging.debug("Running CheckM with command: {0}".format(cmd))

    if 'wd' in kwargs:
        logdir = kwargs.get('wd').get_dir('cmd_logs')
    else:
        logdir = False
    drep.run_cmd(cmd, shell=False, logdir=logdir)

    # Load table and return it
    try:
        chdb = pd.read_table(desired_file,sep='\t')
    except:
        logging.error("!!! checkM failed !!!\nIf using pyenv, make sure both python2 and " +\
            "python3 are available (for example: pyenv global 3.5.1 2.7.9)")
        sys.exit()
    return chdbfix_
