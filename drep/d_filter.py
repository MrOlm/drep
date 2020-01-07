#!/usr/bin/env python3
'''
d_filter - a subset of drep

Filter genomes based on genome length or quality. Also can run prodigal and checkM
'''

import logging
import glob
import pandas as pd
import os
import sys
import shutil
import multiprocessing
from Bio import SeqIO

import drep
import drep.WorkDirectory
import drep.d_cluster
import drep.d_choose
import drep.d_bonus

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
        debug: If True, make extra output when running external scripts

        length: minimum genome length when filtering
        completeness: minimum genome completeness when filtering
        contamination: maximum genome contamination when filtering

        ignoreGenomeQuality: Don't run checkM or do any quality-based filtering (not recommended)
        checkM_method: Either lineage_wf (more accurate) or taxonomy_wf (faster)

    Returns:
        Nothing: stores Bdb.csv, Chdb.csv, and GenomeInfo.csv in the work directory
    '''
    # Load the WorkDirectory.
    logging.debug("Loading work directory in filter")
    workDirectory = drep.WorkDirectory.WorkDirectory(wd)
    logging.debug(str(workDirectory))
    wd = workDirectory

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
    logging.info("Calculating genome info of genomes")
    Gdb = calc_genome_info(bdb['location'].tolist())

    # Filter by size if need be
    if kwargs.get('length', 0) > 1:
        logging.debug("Filtering genomes by size")
        bdb = _filter_bdb_length(bdb, Gdb, kwargs['length'])

    # Get comp/con information
    if kwargs.get('ignoreGenomeQuality', False):
        logging.debug("Skipping all quality-based filtering")

    else:
        Gdb = _get_run_genomeInfo(wd, bdb, **kwargs)

    # Filter
    if not kwargs.get('ignoreGenomeQuality', False):
        logging.debug("Filtering genomes")
        bdb = filter_bdb(bdb, Gdb, **kwargs)

    # Save Bdb and genomeInfo
    logging.debug("Storing resulting files")

    workDirectory.store_db(bdb, 'Bdb')
    if not kwargs.get('ignoreGenomeQuality', False):
        workDirectory.store_db(Gdb, 'genomeInfo')

def _get_run_genomeInfo(workDirectory, bdb, **kwargs):
    '''
    Through kwargs and the wd, get genomeInfo

    Args:
        wd: workDirectory
        bdb: current bdb
        kwrags: keyword arguments

    Keyword arguments:
        no_run: if True, don't run any programs and return None if genomeInfo can't be made

    Returns:
        DataFrame: genomeInfo
    '''
    if kwargs.get('genomeInfo', None) != None:
        logging.debug("Loading provided genome quality information")
        try:
            Idb = pd.read_csv(kwargs.get('genomeInfo'))
            Tdb = _validate_genomeInfo(Idb, bdb)
            Gdb = _add_lengthN50(Tdb, bdb)
        except:
            Idb = pd.read_table(kwargs.get('genomeInfo'))
            Tdb = _validate_genomeInfo(Idb, bdb)
            Gdb = _add_lengthN50(Tdb, bdb)

    elif workDirectory.hasDb('genomeInfo'):
        logging.debug("Loading genomeInfo.csv from work directory")
        Idb = workDirectory.get_db('genomeInfo')
        Tdb = _validate_genomeInfo(Idb, bdb)
        Gdb = _add_lengthN50(Tdb, bdb)

    elif workDirectory.hasDb('Chdb'):
        logging.debug("Loading Chdb.csv from work directory")
        Chdb = workDirectory.get_db('Chdb')
        Idb = chdb_to_genomeInfo(Chdb)
        Tdb = _validate_genomeInfo(Idb, bdb)
        Gdb = _add_lengthN50(Tdb, bdb)

    elif kwargs.get('no_run', False):
        logging.debug("Making basic genomeInfo")
        Gdb = calc_genome_info(bdb['location'].tolist())
        del Gdb['location']

    else:
        logging.debug("Running CheckM")
        Chdb = _run_checkM_wrapper(bdb, workDirectory, **kwargs)
        Idb = chdb_to_genomeInfo(Chdb)
        Tdb = _validate_genomeInfo(Idb, bdb)
        Gdb = _add_lengthN50(Tdb, bdb)

    return Gdb

def _validate_genomeInfo(Idb, bdb):
    '''
    Validate genomeinfo file; merge in Gdb if needed

    Args:
        genomeInfo: DataFrame of genome info
        bdb: DataFrame with genomes to confirm

    Returns:
        DataFrame: Validated genomeInfo
    '''
    # Make sure it has required columns
    for r in ['completeness', 'contamination', 'genome']:
        assert r in Idb.columns, "{0} missing from GenomeInfo".format(r)

    # Make sure correct datatypes
    for r in ['completeness', 'contamination', 'strain_heterogeneity']:
        if r in Idb:
            Idb[r] = Idb[r].astype(float)

    # See if full path was used
    for genome in list(bdb['genome'].unique()):
        if genome not in Idb['genome'].tolist():
            if os.path.basename(genome) in [os.path.basename(x) for x in Idb['genome']]:
                logging.warning("Provided genome info has full genome path- correcting")
                Idb['location'] = Idb['genome']
                Idb['genome'] = [os.path.basename(x) for x in Idb['location']]
                break

    # Make sure it matchs up with bdb
    for genome in list(bdb['genome'].unique()):
        if genome not in Idb['genome'].tolist():
            assert genome in Idb['genome'].tolist(), "{0} missing from GenomeInfo".format(genome)

    # Throw warnings if you think this is weird
    for r in ['completeness', 'contamination', 'strain_heterogeneity']:
        if r in Idb:
            if (Idb[r].max() <= 1) & (Idb[r].min() > 0):
                logging.warning("GenomeInfo has no values over 1 for {0}".format(r) + \
                    "- these should be 0-100, not 0-1!")

    return Idb

def _add_lengthN50(Idb, bdb):
    '''
    If Gdb doesn't have length or N50, add it
    '''
    # Figure out if you need to add them
    add = False
    for r in ['length', 'N50']:
        if r not in Idb:
            add = True

    if not add:
        return Idb

    Gdb = calc_genome_info(bdb['location'].tolist())
    for r in ['length', 'N50']:
        if r not in Idb:
            Idb[r] = Idb['genome'].map(Gdb.set_index('genome')[r].to_dict())

    return Idb

def calc_genome_info(genomes: list):
    '''
    Calculate the length and N50 of a list of genome locations

    Args:
        genomes: list of locations of genomes

    Returns:
        DataFrame: pandas dataframe with ["location", "length", "N50", "genome"]
    '''
    table = {}
    for att in ['location', 'length', 'N50', 'genome']:
        table[att] = []

    for loc in genomes:
        table['location'].append(loc)
        table['genome'].append(os.path.basename(loc))
        table['length'].append(calc_fasta_length(loc))
        table['N50'].append(calc_n50(loc))

    return pd.DataFrame(table)


def filter_bdb(bdb, Gdb, **kwargs):
    '''
    Filter bdb based on Gdb

    Args:
        bdb: DataFrame with ["genome"]
        Gdb: DataFrame with ["genome", "completeness", "contamination"]

    Keyword args:
        min_comp: Minimum genome completeness (%)
        max_con: Maximum genome contamination (%)

    Returns:
        DataFrame: bdb filtered based on completeness and contamination
    '''
    # Get kwargs set up
    min_comp = kwargs.get('completeness',False)
    max_con = kwargs.get('contamination',False)

    # log the number of staring genomes
    start_genomes = list(bdb['genome'].unique())
    assert len(start_genomes) > 0

    # Filter Gdb based on the required metrics
    db = Gdb.copy()
    if min_comp != False:
        db = db[(db['completeness'] >= min_comp)]
    if max_con != False:
        db = db[db['contamination'] <= max_con]

    # Make a list of the remaining genomes
    keep_genomes = list(db['genome'].unique())

    # filter bdb according to these genomes
    bdb = bdb[bdb['genome'].isin(keep_genomes)]

    logging.info("{0:.2f}% of genomes passed checkM filtering".format(\
            (len(keep_genomes)/len(start_genomes))*100))
    return bdb

def _filter_bdb_length(bdb, Gdb, min_length):
    '''
    Filter bdb (['genome', 'location']) based on the minimum length in Gdb
    (['location', 'length', 'N50'])

    Args:
        bdb: DataFrame
        Gdb: DataFrame
        min_length: length to filter on

    Returns:
        DataFrame: bdb filtered by length
    '''

    start = len(bdb['location'].unique())
    bdb['length'] =  bdb['location'].map(Gdb.set_index('location')['length'].to_dict())
    x = bdb[bdb['length'] >= min_length]
    end = len(x['location'].unique())

    logging.info("{0:.2f}% of genomes passed length filtering".format((end/start)*100))

    return x

def validate_chdb(Chdb, bdb):
    '''
    Make sure all genomes in bdb are in Chdb

    Args:
        Chdb: dataframe of checkM information
        bdb: dataframe with ['genome']
    '''
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

    return

def _run_checkM_wrapper(bdb, workDirectory, **kwargs):
    '''
    Controller for running checkM through a workDirectory
    '''
    # Run prodigal
    prod_folder = workDirectory.get_dir('prodigal')
    logging.info("Running prodigal")
    run_prodigal(bdb['location'].tolist(), prod_folder, wd=workDirectory, **kwargs)

    # Make the folder to house the checkM info
    checkM_loc = workDirectory.get_dir('checkM')

    # Run checkM
    logging.info("Running checkM")
    checkM_outfolder = checkM_loc + 'checkM_outdir/'
    Chdb = run_checkM(prod_folder, checkM_outfolder, wd=workDirectory, **kwargs)
    validate_chdb(Chdb, bdb)

    # Fix genome size and N50 of Chdb
    Chdb = _fix_chdb(Chdb,bdb)

    # Save checkM run
    workDirectory.store_db(Chdb,'Chdb')

    return Chdb

def _fix_chdb(Chdb, Bdb):
    '''
    Add corrent genome length and N50 info when running checkM with proteins
    '''
    g2s = {}
    g2n = {}
    for genome in Chdb['Bin Id'].unique():
        loc = Bdb['location'][Bdb['genome'] == genome].tolist()[0]
        g2s[genome] = calc_fasta_length(loc)
        g2n[genome] = calc_n50(loc)
    Chdb['Genome size (bp)'] = Chdb['Bin Id'].map(g2s)
    Chdb['N50 (scaffolds)'] = Chdb['Bin Id'].map(g2n)

    return Chdb

def calc_n50(loc):
    '''
    Calculate the N50 of a .fasta file

    Args:
        fasta_loc: location of .fasta file.

    Returns:
        int: N50 of .fasta file.
    '''
    lengths = []
    for seq_record in SeqIO.parse(loc, "fasta"):
        lengths.append(len(seq_record))

    half = sum(lengths)/2
    tally = 0
    for l in sorted(lengths, reverse=True):
        tally += l
        if tally > half:
            return l

def run_prodigal(genome_list, out_dir, **kwargs):
    '''
    Run prodigal on a set of genomes, store the output in the out_dir

    Args:
        genome_list: list of genomes to run prodigal on
        out_dir: output directory to store prodigal output

    Keyword Args:
        processors: number of processors to multithread with
        exe_loc: location of prodigal excutible (will try and find with shutil if not provided)
        debug: log all of the commands
        wd: if you want to log commands, you also need the wd

    '''
    # Get set up
    t = kwargs.get('processors','6')
    loc = kwargs.get('exe_loc', None)
    if loc == None:
        loc = drep.get_exe('prodigal')

    # Make sure it's a list
    assert type(genome_list) == list

    # Make list of commands
    cmds = []
    for genome in genome_list:
        fna = "{0}{1}".format(os.path.join(out_dir,os.path.basename(genome)),'.fna')
        faa = "{0}{1}".format(os.path.join(out_dir,os.path.basename(genome)),'.faa')
        if os.path.exists(fna) and os.path.exists(faa):
            pass
        else:
            cmds.append(['prodigal','-i',genome,'-d',fna,'-a',faa,'-m','-p','meta'])

    # Run commands
    if len(cmds) > 0:
        if ('wd' in kwargs) and (kwargs.get('debug', False) == True):
            logdir = kwargs.get('wd').get_dir('cmd_logs')
        else:
            logdir = False
            #logdir = "/home/mattolm/Programs/drep/tests/test_backend/logs/"

        drep.thread_cmds(cmds, shell=False, logdir=logdir, t=int(t))

    else:
        logging.info("Past prodigal runs found- will not re-run")

def run_checkM(genome_folder,checkm_outf,**kwargs):
    '''
    Run checkM

    WARNING- this will result in wrong genome lenth and genome N50 estimate, due to
    it being run on prodigal output

    Args:
        genome_folder: location of folder to run checkM on - should be full of files ending in .faa (result of prodigal)
        checkm_outf: location of folder to store checkM output

    Keyword args:
        processors: number of threads
        checkm_method: either lineage_wf or taxonomy_wf
        debug: log all of the commands
        wd: if you want to log commands, you also need the wd
        set_recursion: if not 0, set the python recursion
    '''
    # Find checkm exe
    loc, works = drep.d_bonus.find_program('checkm')
    if loc == None:
        logging.error('Cannot locate the program {0}- make sure its in the system path'\
            .format('checkm'))
        sys.exit()
    if works == False:
        logging.error('Program {0} is not working!! Im going to crash now'\
            .format('checkm'))
        sys.exit()
    check_exe = loc

    # Get set up
    t = str(kwargs.get('processors','6'))
    checkm_method = kwargs.get('checkM_method','lineage_wf')

    # Set recursion
    R = kwargs.get('set_recursion', '0')
    if R != '0':
        logging.warning('Setting Maximum Recursion depth to {0}'.format(R))
        sys.setrecursionlimit(int(R))

    # Run checkM initial
    if checkm_method == 'taxonomy_wf':
         cmd = [check_exe,checkm_method,'domain','Bacteria',genome_folder,checkm_outf,'-f',\
            checkm_outf + '/results.tsv','--tab_table','-t',str(t),'-g','-x','faa']
    else:
         cmd = [check_exe,checkm_method,genome_folder,checkm_outf,'-f',\
            checkm_outf + '/results.tsv','--tab_table','-t',str(t),'--pplacer_threads',\
            str(t),'-g','-x','faa']

    logging.debug("Running CheckM with command: {0}".format(cmd))

    if ('wd' in kwargs) & (kwargs.get('debug', False) == True):
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

    if ('wd' in kwargs) & (kwargs.get('debug', False) == True):
        logdir = kwargs.get('wd').get_dir('cmd_logs')
    else:
        logdir = False
    drep.run_cmd(cmd, shell=False, logdir=logdir)

    # Load table
    try:
        chdb = pd.read_table(desired_file,sep='\t')
    except:
        logging.error("!!! checkM failed !!!\nYou can run again with the --debug option to see what went wrong (command logs will be created in the log folder)")
        sys.exit()

    # Return table
    return chdb

def calc_fasta_length(fasta_loc):
    '''
    Calculate the length of the .fasta file and retun length

    Args:
        fasta_loc: location of .fasta file

    Returns:
        int: total length of all .fasta files
    '''
    total = 0
    for seq_record in SeqIO.parse(fasta_loc, "fasta"):
        total += len(seq_record)
    return total

def chdb_to_genomeInfo(chdb):
    '''
    Convert the output of checkM (chdb) into genomeInfo

    Args:
        chdb: dataframe of checkM

    Returns:
        DataFrame: genomeInfo
    '''
    Gdb = chdb.copy()
    Gdb.rename(columns={'Bin Id': 'genome', 'Completeness': 'completeness',\
        'Contamination':'contamination', \
        'Strain heterogeneity':'strain_heterogeneity'}, inplace=True)
    Gdb = Gdb[['genome', 'completeness', 'contamination', 'strain_heterogeneity']]
    return Gdb
