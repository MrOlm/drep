#!/usr/bin/env python3
'''
d_cluster - a subset of dRep

Clusters a list of genomes with both primary and secondary clustering
'''

import pandas as pd
import os
import glob
import shutil
import logging
import random
import string
import sys
import json
import scipy.cluster.hierarchy
import scipy.spatial.distance as ssd
import numpy as np
import pickle
import time
import glob

import drep
import drep.d_filter
import drep.d_bonus

# from memory_profiler import profile
# import psutil
# process = psutil.Process(os.getpid())

# This is to make pandas shut up with it's warnings
pd.options.mode.chained_assignment = None

def d_cluster_wrapper(workDirectory, **kwargs):
    '''
    Controller for the dRep cluster operation

    Validates arguments, calls cluster_genomes, then stores the output

    Args:
        wd (WorkDirectory): The current workDirectory
        **kwargs: Command line arguments

    Keyword Args:
        processors: Threads to use with checkM / prodigal
        overwrite: Overwrite existing data in the work folder
        debug: If True, make extra output when running external scripts

        MASH_sketch: MASH sketch size
        S_algorithm: Algorithm for secondary clustering comaprisons {ANImf,gANI,ANIn}
        n_PRESET: Presets to pass to nucmer {normal, tight}

        P_ANI: ANI threshold to form primary (MASH) clusters (default: 0.9)
        S_ANI: ANI threshold to form secondary clusters (default: 0.99)
        SkipMash: Skip MASH clustering, just do secondary clustering on all genomes
        SkipSecondary: Skip secondary clustering, just perform MASH clustering
        COV_THRESH: Minmum level of overlap between genomes when doing secondary comparisons (default: 0.1)
        coverage_method: Method to calculate coverage of an alignment {total,larger}
        CLUSTERALG: Algorithm used to cluster genomes (passed to scipy.cluster.hierarchy.linkage (default: average)

        genomes: genomes to cluster in .fasta format. Not necessary if already loaded sequences with the "filter" operation

    Returns:
        Stores Cdb, Mdb, Ndb, Bdb
    '''
    # Load the WorkDirectory.
    logging.debug("Loading work directory")
    workDirectory = drep.WorkDirectory.WorkDirectory(workDirectory)
    logging.debug(str(workDirectory))
    wd = workDirectory

    # Parse arguments
    Bdb = _parse_cluster_arguments(workDirectory, **kwargs)
    data_folder = wd.get_dir('data')

    # Run the main program
    logging.debug("Calling cluster genomes")
    Cdb, Mdb, Ndb = cluster_genomes(list(Bdb['location'].tolist()), \
        data_folder, wd=workDirectory, **kwargs)

    # Save the output
    logging.debug("Main program run complete- saving output to {0}".format(data_folder))
    wd.store_db(Cdb, 'Cdb')
    wd.store_db(Mdb, 'Mdb')
    wd.store_db(Ndb, 'Ndb')
    if not wd.hasDb('Bdb'):
        wd.store_db(Bdb, 'Bdb')

    # Log arguments
    wd.store_special('cluster_log', kwargs)

def cluster_genomes(genome_list, data_folder, **kwargs):
    """
    Clusters a set of genomes using the dRep primary and secondary clustering method

    Takes a number of command line arguments and returns a couple pandas
    dataframes. Done in a number of steps

    Args:
        genomes: list of genomes to be clustered
        data_folder: location where MASH and ANIn data will be stored

    Keyword Args:
        processors: Threads to use with checkM / prodigal
        overwrite: Overwrite existing data in the work folder
        debug: If True, make extra output when running external scripts

        MASH_sketch: MASH sketch size
        S_algorithm: Algorithm for secondary clustering comaprisons {ANImf,gANI,ANIn}
        n_PRESET: Presets to pass to nucmer {normal, tight}

        P_ANI: ANI threshold to form primary (MASH) clusters (default: 0.9)
        S_ANI: ANI threshold to form secondary clusters (default: 0.99)
        SkipMash: Skip MASH clustering, just do secondary clustering on all genomes
        SkipSecondary: Skip secondary clustering, just perform MASH clustering
        COV_THRESH: Minmum level of overlap between genomes when doing secondary comparisons (default: 0.1)
        coverage_method: Method to calculate coverage of an alignment {total,larger}
        CLUSTERALG: Algorithm used to cluster genomes (passed to scipy.cluster.hierarchy.linkage (default: average)

        n_c: nucmer argument c
        n_maxgap: nucmer argument maxgap
        n_noextend: nucmer argument noextend
        n_method: nucmer argument method
        n_preset: preset nucmer arrangements: {tight,normal}

        wd: workDirectory (needed to store clustering results and run prodigal)

    Returns:
        list: [Mdb(db of primary clustering), Ndb(db of secondary clustering, Cdb(clustering information))]
    """
    # Handle debug mode
    if kwargs.get('debug', False):
        debug = True
        wd = kwargs.get('wd')
    else:
        debug = False

    logging.info("Clustering Step 1. Parse Arguments")

    Bdb = load_genomes(genome_list)
    algorithm = kwargs.get('S_algorithm', 'ANImf')

    # Deal with nucmer presets
    if kwargs.get('n_preset', None) != None:
        kwargs['n_c'], kwargs['n_maxgap'], kwargs['n_noextend'], kwargs['n_method'] \
        = _nucmer_preset(kwargs['n_PRESET'])
    logging.debug("kwargs to cluster: {0}".format(kwargs))

    logging.info("Clustering Step 2. Perform MASH (primary) clustering")
    #logging.debug('total memory - {0:.2f} Mbp'.format(int(process.memory_info().rss)/1000000))
    # figure out if you have cached Mdb / CdbF
    cached = (debug and wd.hasDb('Mdb') and wd.hasDb('CdbF'))

    if kwargs.get('SkipMash', False):
        logging.info("2. Nevermind! Skipping Mash")
        # Make a "Cdb" where all genomes are in the same cluster
        Cdb = _gen_nomash_cdb(Bdb)
        # Make a blank "Mdb" for storage anyways
        Mdb = pd.DataFrame({'Blank':[]})

    elif cached:
        logging.info('2. Nevermind! Loading cached primary clustering')
        Mdb = wd.get_db('Mdb')
        Cdb = wd.get_db('CdbF')
        logging.info('2. Primary clustering cache loaded')

    else:
        logging.info("2a. Run pair-wise MASH clustering")
        Mdb = all_vs_all_MASH(Bdb, data_folder, **kwargs)

        if debug:
            logging.debug("Debug mode on - saving Mdb ASAP")
            wd.store_db(Mdb, 'Mdb')

        logging.info("2b. Cluster pair-wise MASH clustering")
        Cdb, cluster_ret = cluster_mash_database(Mdb, **kwargs)

        if debug:
            logging.debug("Debug mode on - saving CdbF ASAP")
            wd.store_db(Cdb, 'CdbF')

        # Store the primary clustering results
        if kwargs.get('wd', None) != None:
            kwargs.get('wd').store_special('primary_linkage', cluster_ret)

    logging.info("{0} primary clusters made".format(len(Cdb['primary_cluster'].unique())))
    logging.info("Step 3. Perform secondary clustering")
    #logging.debug('total memory - {0:.2f} Mbp'.format(int(process.memory_info().rss)/1000000))

    # Wipe any old secondary clusters
    if kwargs.get('wd', None) != None:
        kwargs.get('wd')._wipe_secondary_clusters()

    if not kwargs.get('SkipSecondary', False):
        # See if cached
        cached = (debug and wd.hasDb('Ndb'))

        if cached:
            logging.info('3. Loading cached secondary clustering')
            Ndb = wd.get_db('Ndb')

            # Get rid of broken ones
            Ndb = Ndb.dropna(subset=['reference'])

            logging.info('3. Secondary clustering cache loaded')

        # Run comparisons, make Ndb
        else:
            _print_time_estimate(Bdb, Cdb, algorithm, kwargs.get('processors', 6))
            Ndb = pd.DataFrame()
            for bdb, name in iteratre_clusters(Bdb,Cdb):
                logging.debug('running cluster {0}'.format(name))
                #logging.debug('total memory - {0:.2f} Mbp'.format(int(process.memory_info().rss)/1000000))
                ndb = compare_genomes(bdb, algorithm, data_folder, **kwargs)

                if len(ndb) == 0:
                    logging.error("CRITICAL ERROR WITH PRIMARY CLUSTER {0}; TRYING AGAIN".format(name))
                    ndb = compare_genomes(bdb, algorithm, data_folder, **kwargs)

                if len(ndb) > 0:
                    ndb['primary_cluster'] = name
                    Ndb = Ndb.append(ndb)
                else:
                    logging.error("DOUBLE CRITICAL ERROR AGAIN WITH PRIMARY CLUSTER {0}; SKIPPING".format(name))

            if debug:
                logging.debug("Debug mode on - saving Ndb ASAP")
                wd.store_db(Ndb, 'Ndb')

        # Run clustering on Ndb
        Cdb, c2ret = _cluster_Ndb(Ndb, comp_method=algorithm, **kwargs)

        # Store the secondary clustering results
        if kwargs.get('wd', None) != None:
            kwargs.get('wd').store_special('secondary_linkages', c2ret)

    else:
        logging.info("3. Nevermind! Skipping secondary clustering")
        Cdb = _gen_nomani_cdb(Cdb, data_folder = data_folder, **kwargs)
        Ndb = pd.DataFrame({'Blank':[]})

    logging.info(
    "Step 4. Return output")

    return Cdb, Mdb, Ndb

def _cluster_Ndb(Ndb, id='primary_cluster', **kwargs):
    '''
    Cluster Ndb on id

    Args:
        Ndb: obvious
        id: thing that defines subclusters in Ndb (defauls = 'primary_cluster')

    Keyword arguments:
        all: passed on to genome_hierarchical_clustering

    Returns:
        list: [Cdb, c2ret(cluster name -> clusering files)]
    '''
    Cdb = pd.DataFrame()
    c2ret = {}
    for name, ndb in Ndb.groupby(id):
        cdb, cluster_ret = genome_hierarchical_clustering(ndb, cluster=name, **kwargs)
        cdb[id] = name
        Cdb = pd.concat([Cdb,cdb], ignore_index=True)
        c2ret[name] = cluster_ret

    return Cdb, c2ret

def genome_hierarchical_clustering(Ndb, **kwargs):
    '''
    Cluster ANI database

    Args:
        Ndb: result of secondary clustering

    Keyword arguments:
        clusterAlg: how to cluster the database (default = single)
        S_ani: thershold to cluster at (default = .99)
        cov_thresh: minumum coverage to be included in clustering (default = .5)
        cluster: name of the cluster
        comp_method: comparison algorithm used

    Returns:
        list: [Cdb, {cluster:[linkage, linkage_db, arguments]}]
    '''
    logging.debug('Clustering ANIn database')

    S_Lmethod = kwargs.get('clusterAlg', 'single')
    S_Lcutoff = 1 - kwargs.get('S_ani', .99)
    cov_thresh = float(kwargs.get('cov_thresh',0.5))
    cluster = kwargs.get('cluster','')
    comp_method = kwargs.get('comp_method', 'unk')

    Table = {'genome':[],'secondary_cluster':[]}

    # Handle the case where there's only one genome
    if len(Ndb['reference'].unique()) == 1:
        Table['genome'].append(os.path.basename(Ndb['reference'].unique().tolist()[0]))
        Table['secondary_cluster'].append("{0}_0".format(cluster))
        cluster_ret = []

    else:
        # Make linkage Ndb
        Ldb = make_linkage_Ndb(Ndb, **kwargs)

        # 3) Cluster the linkagedb
        Gdb, linkage = cluster_hierarchical(Ldb, linkage_method= S_Lmethod, \
                                    linkage_cutoff= S_Lcutoff)

        # 4) Extract secondary clusters
        for clust, d in Gdb.groupby('cluster'):
            for genome in d['genome'].tolist():
                Table['genome'].append(genome)
                Table['secondary_cluster'].append("{0}_{1}".format(cluster,clust))

        # 5) Save the linkage
        arguments = {'linkage_method':S_Lmethod,'linkage_cutoff':S_Lcutoff,\
                    'comparison_algorithm':comp_method,'minimum_coverage':cov_thresh}
        cluster_ret = [linkage, Ldb, arguments]

    # Return the database
    Gdb = pd.DataFrame(Table)
    Gdb['threshold'] = S_Lcutoff
    Gdb['cluster_method'] = S_Lmethod
    Gdb['comparison_algorithm'] = comp_method

    return Gdb, cluster_ret


def make_linkage_Ndb(Ndb, **kwargs):
    '''
    Filter the Ndb in accordance with the kwargs. Average reciprical ANI values

    Args:
        Ndb: result of secondary clutsering

    Keyword arguments:
        cov_thresh: minimum coverage threshold (default = 0.5)

    Return:
        DataFrame: pivoted DataFrame ready for clustering
    '''
    d = Ndb.copy()
    cov_thresh = float(kwargs.get('cov_thresh',0.5))

    # Remove values without enough coverage
    d.loc[d['alignment_coverage'] <= cov_thresh, 'ani'] = 0

    # Make a linkagedb by averaging values and setting self-compare to 1
    add_avani(d)
    d['dist'] = 1 - d['av_ani']
    db = d.pivot("reference", "querry", "dist")

    return db

def iteratre_clusters(Bdb, Cdb, id='primary_cluster'):
    '''
    An iterator: Given Bdb and Cdb, yeild smaller Bdb's in the same cluster

    Args:
        Bdb: [genome, location]
        Cdb: [genome, id]
        id: what to iterate on (default = 'primary_cluster')

    Returns:
        list: [d(subset of b), cluster(name of cluster)]
    '''
    Bdb = pd.merge(Bdb,Cdb)
    for cluster, d in Bdb.groupby(id):
        yield d, cluster

def estimate_time(comps, alg):
    '''
    Estimate time, in minutes, based on comparison algorithm and number of comparisons

    Args:
        comps: number of genomes comparisons to perform
        alg: algorthm used

    Return:
        float: time to perfom comparison (in minutes)
    '''
    if alg == 'ANIn':
        time = comps * .33
    elif alg == 'gANI':
        time = comps * .1
    elif alg == 'goANI':
        time = comps * .1
    elif alg == 'ANImf':
        time = comps * .5
    elif alg == 'fastANI':
        time = comps * 0.00667
    return time

def _parse_cluster_arguments(workDirectory, **kwargs):
    '''
    Parse and validate clustering arguments

    Figure out what genomes you're going to be clustering, make sure there's no
    conflicts with the workDirectory

    Args:
        workDirectory: self explainatory

    Returns:
        DataFrame: Bdb
    '''
    # Make sure you have the required program installed
    loc = shutil.which('mash')
    if loc == None:
        logging.error('Cannot locate the program {0}- make sure its in the system path'\
            .format('mash'))

    # If genomes are provided, load them
    if kwargs.get('genomes',None) != None:
        assert workDirectory.hasDb("Bdb") == False, \
            "Don't provide new genomes- you already have them in the work directory"
        Bdb = load_genomes(kwargs['genomes'])

    # If genomes are not provided, don't load them
    if kwargs.get('genomes',None) == None:
        assert workDirectory.hasDb("Bdb") != False, \
        "Must either provide a genome list, or run the 'filter' operation with the same work directory"
        Bdb = workDirectory.get_db('Bdb')

    # Make sure this isn't going to overwrite old data
    # overwrite = kwargs.get('overwrite', True)
    # for db in ['Cdb','Mdb','Ndb']:
    #     if workDirectory.hasDb(db):
    #         if not overwrite:
    #             logging.error("clustering already exists; run with --overwrite to continue")
    #             sys.exit()
    #         logging.debug("THIS WILL OVERWRITE {0}".format(db))

    # Make sure people weren't dumb with their cutoffs
    for v in ['P_ani', 'S_ani']:
        if kwargs.get(v) > 1:
            logging.warning("{0} is set to {1}- this should be \
                between 0-1, not 1-100".format(v, kwargs.get(v)))

    return Bdb

def cluster_hierarchical(db, linkage_method= 'single', linkage_cutoff= 0.10):
    '''
    Perform hierarchical clustering on a symmetrical distiance matrix

    Args:
        db: result of db.pivot usually
        linkage_method: passed to scipy.cluster.hierarchy.fcluster
        linkage_cutoff: distance to draw the clustering line (default = .1)

    Returns:
        list: [Cdb, linkage]
    '''
    # Save names
    names = list(db.columns)

    # Generate linkage dataframe
    arr =  np.asarray(db)
    try:
        arr = ssd.squareform(arr)
    except:
        logging.error("The database passed in is not symmetrical!")
        logging.error(arr)
        logging.error(names)
        sys.exit()
    linkage = scipy.cluster.hierarchy.linkage(arr, method= linkage_method)

    # Form clusters
    fclust = scipy.cluster.hierarchy.fcluster(linkage,linkage_cutoff, \
                    criterion='distance')
    # Make Cdb
    Cdb = _gen_cdb_from_fclust(fclust,names)

    return Cdb, linkage

def _gen_cdb_from_fclust(fclust,names):
    '''
    Make Cdb from the result of scipy.cluster.hierarchy.fcluster

    Args:
        fclust: result of scipy.cluster.hierarchy.fcluster
        names: list(db.columns) of the input dataframe

    Returns:
        DataFrame: Cdb
    '''
    Table={'cluster':[],'genome':[]}
    for i, c in enumerate(fclust):
        Table['cluster'].append(c)
        Table['genome'].append(names[i])

    return pd.DataFrame(Table)

def all_vs_all_MASH(Bdb, data_folder, **kwargs):
    """
    Run MASH pairwise within all samples in Bdb

    Args:
        Bdb: dataframe with genome, location
        data_folder: location to store output files

    Keyword Args:
        MASH_sketch: size of mash sketches
        dry: dont actually run anything
        processors: number of processors to multithread with
        mash_exe: location of mash excutible (will try and find with shutil if not provided)
        groupSize: max number of mash sketches to hold in each folder
        debug: if True, log all of the commands
        wd: if you want to log commands, you also need the wd
    """

    MASH_s = kwargs.get('MASH_sketch',1000)
    dry = kwargs.get('dry',False)
    # overwrite = kwargs.get('overwrite', False)
    mash_exe = kwargs.get('mash_exe', None)
    p = kwargs.get('processors',6)
    groupSize = kwargs.get('groupSize', 1000)

    # set up logdir
    if ('wd' in kwargs) and (kwargs.get('debug', False) == True):
        logdir = kwargs.get('wd').get_dir('cmd_logs')
    else:
        logdir = False

    # Find mash
    mash_exe = kwargs.get('exe_loc', None)
    if mash_exe == None:
        mash_exe = drep.get_exe('mash')

    # Set up folders
    MASH_folder = os.path.join(data_folder, 'MASH_files/')
    if not os.path.exists(MASH_folder):
        os.makedirs(MASH_folder)

    sketch_folder = os.path.join(MASH_folder, 'sketches/')
    if not os.path.exists(sketch_folder):
        os.makedirs(sketch_folder)

    # Make chunks
    l2g = Bdb.set_index('location')['genome'].to_dict()
    locations = list(Bdb['location'].unique())
    chunks = [locations[x:x+groupSize] for x in range(0, len(locations), groupSize)]

    # Make the MASH sketches
    cmds = []
    chunk_folders = []
    for i, chunk in enumerate(chunks):
        chunk_folder = os.path.join(sketch_folder, "chunk_{0}".format(i))
        chunk_folders.append(chunk_folder)
        if not os.path.exists(chunk_folder):
            os.makedirs(chunk_folder)
        for fasta in chunk:
            genome = l2g[fasta]
            file = os.path.join(chunk_folder, genome)
            if not os.path.isfile(file + '.msh'):
                cmd = [mash_exe, 'sketch', fasta, '-s', str(MASH_s), '-o',
                    file]
                cmds.append(cmd)

    if not dry:
        if len(cmds) > 0:
            drep.thread_cmds(cmds, logdir=logdir, t=int(p))

    # Combine MASH sketches within chunk
    cmds = []
    alls = []
    for chunk_folder in chunk_folders:
        all_file = os.path.join(chunk_folder, 'chunk_all.msh')
        cmd = [mash_exe, 'paste', all_file] \
                + glob.glob(os.path.join(chunk_folder, '*'))
        cmds.append(cmd)
        alls.append(all_file)
    if not dry:
        if len(cmds) > 0:
            drep.thread_cmds(cmds, logdir=logdir, t=int(p))

    # Combine MASH sketches of all chunks
    all_file = os.path.join(MASH_folder, 'ALL.msh')
    cmd = [mash_exe, 'paste', all_file] + alls
    drep.run_cmd(cmd, dry, shell=False, logdir=logdir)

    # Calculate distances
    cmd = [mash_exe, 'dist','-p', str(p), all_file, all_file, '>', MASH_folder
            + 'MASH_table.tsv']
    cmd = ' '.join(cmd)
    drep.run_cmd(cmd, dry, shell=True, logdir=logdir)

    # Make Mdb based on all genomes in the MASH folder
    file = MASH_folder + 'MASH_table.tsv'

    iniCols = ['genome1','genome2','dist','p','kmers']
    uCols = ['genome1','genome2','dist']
    dTypes = {'genome1':'category', 'genome2':'category', 'dist':np.float32}
    Mdb = pd.read_csv(file, names=iniCols, usecols=uCols, dtype=dTypes, sep='\t')
    Mdb['genome1'] = Mdb['genome1'].apply(_get_genome_name_from_fasta)
    Mdb['genome2'] = Mdb['genome2'].apply(_get_genome_name_from_fasta)
    Mdb['similarity'] = 1 - Mdb['dist']

    # Filter out those genomes that are in the MASH folder but shouldn't be in Mdb
    genomes = Bdb['genome'].unique()
    Mdb = Mdb[Mdb['genome1'].isin(genomes)]
    Mdb = Mdb[Mdb['genome2'].isin(genomes)]

    # Reorder categories to be correct
    for g in ['genome1', 'genome2']:
        Mdb[g] = Mdb[g].cat.remove_unused_categories()
        Mdb[g] = Mdb[g].cat.reorder_categories(sorted((Mdb[g].unique())), ordered=True)

    return Mdb

def cluster_mash_database(db, **kwargs):
    '''
    From a Mash database, cluster and return Cdb

    Args:
        db: Mdb (all_vs_all Mash results)

    Keyword arguments:
        clusterAlg: how to cluster database (default = single)
        P_ani: threshold to cluster at (default = 0.9)

    Returns:
        list: [Cdb, [linkage, linkage_db, arguments]]
    '''
    logging.debug('Clustering MASH database')

    # Load key words
    P_Lmethod = kwargs.get('clusterAlg','single')
    P_Lcutoff = 1 - kwargs.get('P_ani',.9)

    # Do the actual clustering
    db['dist'] = 1 - db['similarity']
    linkage_db = db.pivot("genome1","genome2","dist")
    Cdb, linkage = cluster_hierarchical(linkage_db, linkage_method= P_Lmethod, \
                                linkage_cutoff= P_Lcutoff)
    Cdb = Cdb.rename(columns={'cluster':'primary_cluster'})

    # Preparing clustering for return
    arguments = {'linkage_method':P_Lmethod,'linkage_cutoff':P_Lcutoff,\
                    'comparison_algorithm':'MASH'}
    cluster_ret = [linkage, linkage_db, arguments]

    return Cdb, cluster_ret

def _get_genome_name_from_fasta(fasta):
    '''
    Just do a os.path.basename command

    Args:
        fasta: location of .fasta file

    Returns:
        string: basename of .fasta files
    '''
    return str(os.path.basename(fasta))

def parse_delta(filename):
    '''
    Parse a .delta file from nucmer

    Args:
        filename: location of .delta file

    Returns:
        list: [alignment_length, similarity_errors]
    '''

    aln_length, sim_errors = 0, 0
    for line in [l.strip().split() for l in open(filename, 'r').readlines()]:
        if line[0] == 'NUCMER' or line[0].startswith('>'):  # Skip headers
            continue
        # We only process lines with seven columns:
        if len(line) == 7:
            aln_length += abs(int(line[1]) - int(line[0]))
            sim_errors += int(line[4])
    return aln_length, sim_errors

def gen_gANI_cmd(file, g1, g2, exe):
    '''
    Generate a command to run gANI (ANIcalculator)

    Will return [] if g1 == g2

    Args:
        file: output file name
        g1: location of the genes of genome1
        g2: location of the genes of genome2
        exe: location of gANI executible

    Returns:
        list: command to run
    '''
    # Don't do self-comparions
    if g1 == g2:
        return []

    else:
        cmd = [exe,'-genome1fna',g1,'-genome2fna',g2,'-outfile',file,'-outdir',file + 'TEMP']
    return cmd

def gen_goANI_cmd(file, g1, g2, exe):
    '''
    Generate a command to run goANI

    Will return [] if g1 == g2

    Args:
        file: output file name
        g1: location of the genes of genome1
        g2: location of the genes of genome2
        exe: location of nsimscan executible

    Returns:
        list: command to run
    '''
    # Don't do self-comparions
    if g1 == g2:
        return []

    else:
        cmd = [exe,'--om','TABX', g1, g2, file]
    return cmd

def process_deltafiles(deltafiles, org_lengths, logger=None, **kwargs):
    '''
    Parse a list of delta files into a pandas dataframe

    Files NEED to be named in the format genome1_vs_genome2.delta, or genome1_vs_genome2.filtered.delta

    Args:
        deltafiles: list of .delta files
        org_lengths: dictionary of genome to genome length

    Keyword Arguments:
        coverage_method: default is larger
        logger: if not None, will log 0 errors

    Returns:
        DataFrame: information about the alignments
    '''

    Table = {'querry':[],'reference':[],'alignment_length':[],'similarity_errors':[],
            'ref_coverage':[],'querry_coverage':[],'ani':[], 'reference_length':[],
            'querry_length':[],'alignment_coverage':[]}

    coverage_method = kwargs.get('coverage_method', 'larger')

    for deltafile in deltafiles:
        qname, sname = os.path.splitext(os.path.split(deltafile)[-1])[0].split('_vs_')
        # Handle naming of "filtered"
        if sname.endswith('.delta'):
            sname = sname[:-6]

        tot_length, tot_sim_error = parse_delta(deltafile)
        if tot_length == 0 and logger is not None:
            logging.info("Total alignment length reported in " +
                               "%s is zero!" % deltafile)
        query_cover = float(tot_length) / org_lengths[qname]
        sbjct_cover = float(tot_length) / org_lengths[sname]

        try:
            perc_id = 1 - float(tot_sim_error) / tot_length
        except ZeroDivisionError:
            #print("Alignment between {0} and {1} has 0 alignment!".format(qname,sname))
            perc_id = 0  # set arbitrary value of zero identity


        Table['querry'].append(qname)
        Table['querry_length'].append(org_lengths[qname])
        Table['reference'].append(sname)
        Table['reference_length'].append(org_lengths[sname])
        Table['alignment_length'].append(tot_length)
        Table['similarity_errors'].append(tot_sim_error)
        Table['ani'].append(perc_id)
        Table['ref_coverage'].append(sbjct_cover)
        Table['querry_coverage'].append(query_cover)

        if coverage_method == 'total':
            Table['alignment_coverage'].append((tot_length * 2)/(org_lengths[qname]\
                                                             + org_lengths[sname]))
        elif coverage_method == 'larger':
            Table['alignment_coverage'].append(max((tot_length/org_lengths[qname]),\
                (tot_length/org_lengths[sname])))

    dTypes = {'querry':'category','querry_length':np.int64, 'reference':'category',\
            'reference_length':np.int64, 'alignment_length':np.int64, 'similarity_errors':np.int64,\
            'ani':np.float32, 'ref_coverage':np.float32, 'querry_coverage':np.float32}
    df = pd.DataFrame(Table)
    for c, t in dTypes.items():
        df[c] = df[c].astype(t)
    return df

def process_goani_files(files):
    '''
    From a list of nsimscan files, return a parsed DataFrame like normal,
    and a special one for calculating dn/ds

    Args:
        list: files

    Returns:
        DataFrame: Ndb
    '''
    Table = {'querry':[],'reference':[],'ani':[],'alignment_coverage':[]}
    for file in files:
        results = parse_nsim_file(file)
        results['reference'] = os.path.basename(file).split('_vs_')[0]
        results['querry'] = os.path.basename(file).split('_vs_')[1][:-4]

        # Add forward results
        Table['reference'].append(results['reference'])
        Table['querry'].append(results['querry'])
        Table['ani'].append(float(results['ani'])/100)
        Table['alignment_coverage'].append(results['af'])
        #Table['alignment_coverage'].append(results['rq_coverage'])

        # Add reverse results [they're the same]
        # Table['reference'].append(results['querry'])
        # Table['querry'].append(results['reference'])
        # Table['ani'].append(float(results['ani'])/100)
        # Table['alignment_coverage'].append(results['af'])
        #Table['alignment_coverage'].append(results['qr_coverage'])

    # Add self comparisons
    for g in set(Table['reference']):
        Table['reference'].append(g)
        Table['querry'].append(g)
        Table['ani'].append(1)
        Table['alignment_coverage'].append(1)

    Gdb = pd.DataFrame(Table)
    return Gdb

def process_gani_files(files):
    '''
    From a list of gANI output files, return a parsed DataFrame

    Args:
        list: files

    Returns:
        DataFrame: Ndb
    '''
    Table = {'querry':[],'reference':[],'ani':[],'alignment_coverage':[]}
    for file in files:
        results = parse_gani_file(file)

        # calculate coverage using the "larger" method
        cov = max(results['rq_coverage'], results['qr_coverage'])

        # Add forward results
        Table['reference'].append(results['reference'])
        Table['querry'].append(results['querry'])
        Table['ani'].append(float(results['rq_ani'])/100)
        Table['alignment_coverage'].append(cov)
        #Table['alignment_coverage'].append(results['rq_coverage'])

        # Add reverse results [they're the same]
        Table['reference'].append(results['querry'])
        Table['querry'].append(results['reference'])
        Table['ani'].append(float(results['qr_ani'])/100)
        Table['alignment_coverage'].append(cov)
        #Table['alignment_coverage'].append(results['qr_coverage'])

    # Add self comparisons
    for g in set(Table['reference']):
        Table['reference'].append(g)
        Table['querry'].append(g)
        Table['ani'].append(1)
        Table['alignment_coverage'].append(1)

    Gdb = pd.DataFrame(Table)
    return Gdb


def compare_genomes(bdb, algorithm, data_folder, **kwargs):
    '''
    Compare a list of genomes using the algorithm specified

    This method takes in bdb (a table with the columns location and genome), runs
    pair-wise comparisons between all genomes in the sample, and returns a table
    with at least the columns 'reference', 'querry', 'ani','coverage', depending
    on what algorithm is called

    Args:
        bdb: DataFrame with ['genome', 'location'] (drep.d_filter.load_genomes)
        algorithm: options are ANImf, ANIn, gANI
        data_folder: location to store output files

    Keyword Arguments:
        wd: either this or prod_folder needed for gANI
        prod_folder: either this or wd needed for gANI

    Return:
        DataFrame: Ndb (['reference', 'querry', 'ani','coverage'])
    '''
    # To handle other versions of this method which passed in a WorkDirectory
    # instead of data_folder string
    if isinstance(data_folder, drep.WorkDirectory.WorkDirectory):
        data_folder = data_folder.get_dir('data')

    if algorithm == 'ANImf':
        genome_list = bdb['location'].tolist()
        working_data_folder = os.path.join(data_folder, 'ANImf_files/')
        df = run_pairwise_ANImf(genome_list, working_data_folder, **kwargs)
        return df

    elif algorithm == 'ANIn':
        genome_list = bdb['location'].tolist()
        working_data_folder = os.path.join(data_folder, 'ANIn_files/')
        df = run_pairwise_ANIn(genome_list, working_data_folder, **kwargs)
        return df

    elif algorithm == 'fastANI':
        genome_list = bdb['location'].tolist()
        working_data_folder = os.path.join(data_folder, 'fastANI_files/')
        df = run_pairwise_fastANI(genome_list, working_data_folder, **kwargs)
        return df

    elif algorithm == 'gANI':
        # Figure out prodigal folder
        wd = kwargs.get('wd', False)
        if not wd:
            prod_folder = kwargs.pop('prod_folder', False)
            assert prod_folder != False
        else:
            prod_folder = wd.get_dir('prodigal')

        working_data_folder = os.path.join(data_folder, 'gANI_files/')
        df = run_pairwise_gANI(bdb, working_data_folder, \
                prod_folder=prod_folder, **kwargs)
        return df

    elif algorithm == 'goANI':
        # Figure out prodigal folder
        wd = kwargs.get('wd', False)
        if not wd:
            prod_folder = kwargs.pop('prod_folder', False)
            assert prod_folder != False
        else:
            prod_folder = wd.get_dir('prodigal')

        working_data_folder = os.path.join(data_folder, 'goANI_files/')
        df = run_pairwise_goANI(bdb, working_data_folder, \
                prod_folder=prod_folder, **kwargs)
        return df

    else:
        logging.error("{0} not supportedd".format(algorithm))
        sys.exit()

def run_pairwise_ANIn(genome_list, ANIn_folder, **kwargs):
    '''
    Given a list of genomes and an output folder, compare all genomes using ANImf

    Args:
        genome_list: list of locations of genome files
        ANIn_folder: folder to store the output of comparison

    Keyword arguments:
        processors: threads to use
        debug: if true save extra output
        wd: needed if debug is True
    '''
    p = kwargs.get('processors',6)
    genomes = genome_list

    # Make folder
    if not os.path.exists(ANIn_folder):
        os.makedirs(ANIn_folder)

    # Gen commands
    cmds = []
    files = []
    for g1 in genomes:
        # Make it so each reference is it's own folder, to spread out .delta files
        cur_folder = os.path.join(ANIn_folder, _get_genome_name_from_fasta(g1))
        if not os.path.exists(cur_folder):
            os.makedirs(cur_folder)

        for g2 in genomes:
            file_name = "{0}{1}_vs_{2}".format(ANIn_folder, \
                        _get_genome_name_from_fasta(g1),\
                        _get_genome_name_from_fasta(g2))
            files.append(file_name)

            # If the file doesn't already exist, add it to what needs to be run
            if not os.path.isfile(file_name + '.delta'):
                cmds.append(gen_nucmer_cmd(file_name,g1,g2))

    # Run commands
    if len(cmds) > 0:
        for c in cmds:
            logging.debug(' '.join(c))

        if ('wd' in kwargs) and (kwargs.get('debug', False)):
            logdir = kwargs.get('wd').get_dir('cmd_logs')
        else:
            logdir = False
        drep.thread_cmds(cmds, shell=False, logdir=logdir, t=int(p))

    # Make dictionary of genome lengths
    org_lengths = {}
    for genome in genomes:
        org_lengths[_get_genome_name_from_fasta(genome)] = \
            drep.d_filter.calc_fasta_length(genome)



    deltafiles = ["{0}.delta".format(file) for file in files]
    df = process_deltafiles(deltafiles, org_lengths, **kwargs)

    return df

def run_pairwise_ANImf(genome_list, ANIn_folder, **kwargs):
    '''
    Given a list of genomes and an output folder, compare all genomes using ANImf

    Args:
        genome_list: list of locations of genome files
        ANIn_folder: folder to store the output of comparison

    Keyword arguments:
        processors: threads to use
        debug: if true save extra output
        wd: needed if debug is True
    '''
    p = kwargs.get('processors',6)
    genomes = genome_list

    # Make folder if doesnt exist
    if not os.path.exists(ANIn_folder):
        os.makedirs(ANIn_folder)

    # Gen commands
    cmds = []
    files = []
    for g1 in genomes:
        # Make it so each reference is it's own folder, to spread out .delta files
        cur_folder = os.path.join(ANIn_folder, _get_genome_name_from_fasta(g1))
        if not os.path.exists(cur_folder):
            os.makedirs(cur_folder)

        for g2 in genomes:
            file_name = "{0}/{1}_vs_{2}".format(cur_folder, \
                        _get_genome_name_from_fasta(g1),\
                        _get_genome_name_from_fasta(g2))
            files.append(file_name)

            # If the file doesn't already exist, add it to what needs to be run
            if not os.path.isfile(file_name + '.delta.filtered'):
                cmds.append(gen_animf_cmd(file_name,g1,g2))

    # Run commands
    if len(cmds) > 0:
        for c in cmds:
            logging.debug(c)

        if ('wd' in kwargs) and (kwargs.get('debug', False)):
            logdir = kwargs.get('wd').get_dir('cmd_logs')
        else:
            logdir = False
        drep.thread_cmds(cmds, shell=True, logdir=logdir, t=int(p))

    # Make dictionary of genome lengths
    org_lengths = {}
    for genome in genomes:
        org_lengths[_get_genome_name_from_fasta(genome)] = drep.d_filter.calc_fasta_length(genome)

    # Parse output
    deltafiles = ["{0}.delta.filtered".format(f) for f in files]
    df = process_deltafiles(deltafiles, org_lengths, **kwargs)

    return df

def run_pairwise_fastANI(genome_list, outdir, **kwargs):
    p = kwargs.get('processors',6)
    code = _randomString(stringLength=10)

    # Make folders
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    tmp_dir = os.path.join(outdir, 'tmp/')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # Make genome list
    glist = os.path.join(tmp_dir, 'genomeList')
    glist = _make_glist(genome_list, glist)

    # Gen command
    exe_loc = drep.get_exe('fastANI')
    out_base = os.path.join(outdir, 'fastANI_out_{0}'.format(code))
    cmd = [exe_loc, '--ql', glist, '--rl', glist, '-o', out_base, '--matrix', '-t', str(p), '--minFraction', str(0)]
    logging.debug(' '.join(cmd) + ' ' + code)

    # Run command
    if ('wd' in kwargs) and (kwargs.get('debug', False)):
        logdir = kwargs.get('wd').get_dir('cmd_logs')
    else:
        logdir = False
    drep.thread_cmds([cmd], shell=False, logdir=logdir, t=1)

    # Load results
    fdb = load_fastani(out_base)

    # fix missing ones
    try:
        fdb = _fix_fastani(fdb)
        return fdb

    # handle broken self
    except:
        logging.error("CRITICAL ERROR WITH SECONDARY CLUSTERING CODE {0}; SKIPPING".format(code))
        return pd.DataFrame()

def load_fastani(file):
    fdb = pd.read_csv(file, names=['reference', 'querry', 'ANI', 'j1', 'j2'], delim_whitespace=True)
    for c in ['reference', 'querry']:
        fdb[c] = [_get_genome_name_from_fasta(x) for x in fdb[c]]
    fdb = fdb.rename(columns={'ANI':'ani'})
    fdb['alignment_coverage'] = [(j1/j2) for j1, j2 in zip(fdb['j1'], fdb['j2'])]
    fdb = fdb[['reference', 'querry', 'ani', 'alignment_coverage']]
    fdb['ani'] = [x/100 for x in fdb['ani']]

    return fdb

def _fix_fastani(odb):
    # Add back missing genomes
    fdb = odb.pivot('reference', 'querry', 'ani')
    fdb.reset_index(level=0, inplace=True)
    fdb.fillna(0, inplace=True)
    fdb = fdb.melt(id_vars=['reference']).rename(
            columns={'value':'ani'})

    # Add back alignment coverage
    fdb = pd.merge(fdb, odb[['reference', 'querry', 'alignment_coverage']], on=['reference', 'querry'], how='outer')
    fdb['alignment_coverage'] = fdb['alignment_coverage'].fillna(0)

    assert len(fdb['reference'].unique()) == len(fdb['querry'].unique())
    assert len(fdb) == (len(fdb['reference'].unique()) * len(fdb['querry'].unique()))

    return fdb

def _make_glist(genomes, floc):
    o = open(floc, 'w')
    for g in genomes:
        loc = os.path.abspath(g)
        assert os.path.isfile(g)
        o.write(loc + '\n')
    o.close()
    return floc

def run_pairwise_gANI(bdb, gANI_folder, prod_folder, **kwargs):
    '''
    Run pairwise gANI on a list of Genomes

    Args:
        bdb: DataFrame with ['genome', 'location']
        gANI_folder: folder to store gANI output
        prod_folder: folder containing prodigal output from genomes (will run if needed)

    Keyword arguments:
        debug: log all of the commands
        wd: if you want to log commands, you also need the wd
        processors: threads to use

    Returns:
        DataFrame: Ndb for gANI
    '''
    p = kwargs.get('processors',6)
    gANI_exe = drep.get_exe('ANIcalculator')
    genomes = bdb['location'].tolist()

    # Make folders
    if not os.path.exists(gANI_folder):
        os.makedirs(gANI_folder)
    if not os.path.exists(prod_folder):
        os.makedirs(prod_folder)

    # Remove crap folders- they shouldn't exist and if they do it messes things up
    crap_folders = glob.glob(gANI_folder + '*.gANITEMP')
    for crap_folder in crap_folders:
        if os.path.exists(crap_folder):
            shutil.rmtree(crap_folder)

    # Run prodigal
    logging.debug("Running prodigal...")
    drep.d_filter.run_prodigal(bdb['location'].tolist(), prod_folder, **kwargs)

    # Gen gANI commands
    logging.debug("Running gANI...")
    cmds = []
    files = []
    for i, g1 in enumerate(genomes):
        # Make it so each reference is it's own folder, to spread out .delta files
        cur_folder = os.path.join(gANI_folder, _get_genome_name_from_fasta(g1))
        if not os.path.exists(cur_folder):
            os.makedirs(cur_folder)

        for j, g2 in enumerate(genomes):
            if i > j:
                name1= _get_genome_name_from_fasta(g1)
                name2= _get_genome_name_from_fasta(g2)
                file_name = "{0}/{1}_vs_{2}.gANI".format(cur_folder, \
                            name1, name2)
                files.append(file_name)

                # If the file doesn't already exist, add it to what needs to be run
                if not os.path.isfile(file_name):
                    fna1 = "{0}.fna".format(os.path.join(prod_folder,name1))
                    fna2 = "{0}.fna".format(os.path.join(prod_folder,name2))
                    cmds.append(gen_gANI_cmd(file_name,fna1,fna2,gANI_exe))

    # Run commands
    if len(cmds) > 0:
        logging.debug('Running gANI commands: {0}'.format('\n'.join([' '.join(x) for x in cmds])))
        if ('wd' in kwargs) and (kwargs.get('debug', False) == True):
            logdir = kwargs.get('wd').get_dir('cmd_logs')
        else:
            logdir = False
            #logdir = "/home/mattolm/Programs/drep/tests/test_backend/logs/"
        drep.thread_cmds(cmds, logdir=logdir, t=int(p))

    else:
        logging.debug("gANI already run- will not re-run")

    # Parse output
    df = process_gani_files(files)

    # Add self-comparisons if there is only one genome
    if len(genomes) == 1:
        Table = {'querry':[],'reference':[],'ani':[],'alignment_coverage':[]}
        for g in genomes:
            Table['reference'].append(_get_genome_name_from_fasta(g))
            Table['querry'].append(_get_genome_name_from_fasta(g))
            Table['ani'].append(1)
            Table['alignment_coverage'].append(1)
        d = pd.DataFrame(Table)
        df = pd.concat([df,d],ignore_index=True)

    return df

def run_pairwise_goANI(bdb, goANI_folder, prod_folder, **kwargs):
    '''
    Run pairwise goANI on a list of Genomes

    Args:
        bdb: DataFrame with ['genome', 'location']
        goANI_folder: folder to store gANI output
        prod_folder: folder containing prodigal output from genomes (will run if needed)

    Keyword arguments:
        debug: log all of the commands
        wd: if you want to log commands, you also need the wd
        processors: threads to use

    Returns:
        DataFrame: Ndb for gANI
    '''
    p = kwargs.get('processors',6)
    nsimscan_exe = drep.get_exe('nsimscan')
    genomes = bdb['location'].tolist()

    # Make folders
    if not os.path.exists(goANI_folder):
        os.makedirs(goANI_folder)
    if not os.path.exists(prod_folder):
        os.makedirs(prod_folder)

    # Run prodigal
    logging.debug("Running prodigal...")
    drep.d_filter.run_prodigal(bdb['location'].tolist(), prod_folder, **kwargs)

    # Gen gANI commands
    logging.debug("Running goANI...")
    cmds = []
    files = []
    for i, g1 in enumerate(genomes):
        # Make it so each reference is it's own folder, to spread out .delta files
        cur_folder = os.path.join(goANI_folder, _get_genome_name_from_fasta(g1))
        if not os.path.exists(cur_folder):
            os.makedirs(cur_folder)

        for j, g2 in enumerate(genomes):
            if i != j:
                name1= _get_genome_name_from_fasta(g1)
                name2= _get_genome_name_from_fasta(g2)
                file_name = "{0}/{1}_vs_{2}.sim".format(cur_folder, \
                            name1, name2)
                files.append(file_name)

                # If the file doesn't already exist, add it to what needs to be run
                if not os.path.isfile(file_name):
                    fna1 = "{0}.fna".format(os.path.join(prod_folder,name1))
                    fna2 = "{0}.fna".format(os.path.join(prod_folder,name2))
                    cmds.append(gen_goANI_cmd(file_name,fna1,fna2,nsimscan_exe))

    # Run commands
    if len(cmds) > 0:
        logging.debug('Running goANI commands: {0}'.format('\n'.join([' '.join(x) for x in cmds])))
        if ('wd' in kwargs) and (kwargs.get('debug', False) == True):
            logdir = kwargs.get('wd').get_dir('cmd_logs')
        else:
            logdir = False
            #logdir = "/home/mattolm/Programs/drep/tests/test_backend/logs/"
        drep.thread_cmds(cmds, logdir=logdir, t=int(p))

    else:
        logging.debug("goANI already run- will not re-run")

    # Parse output
    df = process_goani_files(files)

    # Add self-comparisons if there is only one genome
    if len(genomes) == 1:
        Table = {'querry':[],'reference':[],'ani':[],'alignment_coverage':[]}
        for g in genomes:
            Table['reference'].append(_get_genome_name_from_fasta(g))
            Table['querry'].append(_get_genome_name_from_fasta(g))
            Table['ani'].append(1)
            Table['alignment_coverage'].append(1)
        d = pd.DataFrame(Table)
        df = pd.concat([df,d],ignore_index=True)

    return df

def parse_gani_file(file):
    '''
    Parse gANI file, return dictionary of results

    Args:
        file: location of gANI file

    Returns:
        dict: results in the gANI file
    '''
    try:
        x = pd.read_table(file)
    except:
        logging.error('gANI file {0} does not exist. The most likely reason is that '.format(file)\
            + 'one of the genomes has a .fasta header than gANI doesnt like.'\
            + ' Known issues include having a header over 160 characaters long, or any '\
            + 'kind of special character besides "_" (including ".", ":", and "-").'\
            + ' To fix this error either use ANIn (which is not as picky with fasta headers)'\
            + ', or fix the .fasta headers to conform to those rules.')
        sys.exit()
    x = x.rename(columns={'GENOME1':'reference','GENOME2':'querry','AF(1->2)':'rq_coverage',\
                        'AF(2->1)':'qr_coverage','ANI(1->2)':'rq_ani','ANI(2->1)':'qr_ani'})
    dict = x.to_dict(orient='list')
    dict = {k:v[0] for k,v in dict.items()}
    dict['reference'] = dict['reference'][:-4]
    dict['querry'] = dict['querry'][:-4]
    return dict

def parse_nsim_file(file):
    '''
    Parse nsim file, return dictionary of results and gene datatable

    Args:
        file: location of nsimscan file

    Returns:
        dict: results in the gANI file
        db: gene-based alignment results
    '''
    # Load
    db = pd.read_csv(file, sep='\t')

    if len(db) == 0:
        logging.warning("File {0} is empty, indicating a nsimscan failure! Run with --debug and check the log folder for details".format(file))

        return _summarize_nsimsan(db)

    db = db.rename(columns={'#qry_id':'qry_id', '#Q_id':'qry_id', 'S_id':'sbj_id', 'trg_len':'sbj_len'})

    # Filter
    try:
        db = _filter_nsimscan(db)
    except:
        print("ERROR! FILE {0}".format(file))
        print(db)
        return _summarize_nsimsan(pd.DataFrame())

    # Make summary results
    x = _summarize_nsimsan(db)

    return x

def _filter_nsimscan(db1, minAF=0.7, minANI=70):
    '''
    Filter a single nsim scan file
    '''
    db1 = db1.sort_values(['al_len','qry_id', 'sbj_id'], ascending=False)

    # Filter like gANI
    db1['af'] = [a/min(o,t) for a,o,t in zip(db1['al_len'],
                                             db1['qry_len'],
                                             db1['sbj_len'])]
    db1 = db1[(db1['af'] >= minAF) & (db1['p_inden'] >= minANI)]

    return db1

def _summarize_nsimsan(db):
    '''
    Take all those aligned genes and return a summary ANI
    '''
    table = {}
    if len(db) > 0:
        table['ani'] = sum([p * l for p, l in zip(db['p_inden'], db['al_len'])]) \
                        / db['al_len'].sum()
        table['af'] = db['al_len'].sum() / db['qry_len'].sum()
    else:
        table['ani'] = 0
        table['af'] = 0

    return table


def gen_nucmer_cmd(prefix,ref,querry,c='65',noextend=False,maxgap='90',method='mum'):
    '''
    Generate command to run with nucmer

    Args:
        prefix: desired output name (will have .delta appended)
        ref: location of reference genome
        querry: location of querry genomes

    Keyword args:
        c: c value
        noextend: either True or False
        maxgap: maxgap
        method: detault is 'mum'

    Returns:
        list: command to run number
    '''
    cmd = ['nucmer','--' + method,'-p',prefix,'-c',str(c),'-g',str(maxgap)]
    if noextend: cmd.append('--noextend')
    cmd += [ref,querry]
    return cmd

def gen_animf_cmd(prefix,ref,querry, **kwargs):
    '''
    return animf command. It will be a single string, with the format:
    "nucmer cmd ; filter cmd"

    Args:
        prefix: desired file output name
        ref: location of reference genome
        querry: location of querry genome

    Keyword args:
        all: passed on to gen_nucmer_cmd

    Returns:
        string: format is "nucmer cmd ; filter cmd"
    '''

    nucmer_cmd = gen_nucmer_cmd(prefix,ref,querry, **kwargs)

    delta = prefix + '.delta'
    out = delta + '.filtered'

    filter_cmd = gen_filter_cmd(prefix + '.delta', out)
    return ' '.join(nucmer_cmd) + '; ' + ' '.join(filter_cmd)

def gen_filter_cmd(delta, out):
    '''
    return delta-filter command

    Args:
        delta: desired .delta file to filter
        out: desired output file

    Returns:
        list: cmd to run
    '''
    cmd = ["delta-filter", '-r', '-q', delta, '>', out]
    return cmd

def _gen_nomash_cdb(Bdb):
    '''
    From Bdb, just add a column of 'primary_cluster' = 0

    Args:
        Bdb: dataframe with [genome, location]

    Returns:
        DataFrame: Cdb
    '''
    Cdb = Bdb.copy()
    Cdb['primary_cluster'] = 0
    return Cdb

def load_genomes(genome_list):
    '''
    Takes a list of genome locations, returns a pandas dataframe with the
    locations and the basename of the genomes

    Args:
        genome_list: list of genome locations

    Returns:
        DataFrame: pandas dataframe with columns ['genome', 'location']
    '''
    assert type(genome_list) == type(list())

    # Load from a text file
    if len(genome_list) == 1:
        logging.info('Loading genomes from a list')
        try:
            Table = {'genome':[],'location':[]}
            with open(genome_list[0], 'r') as o:
                for line in o.readlines():
                    genome = line.strip()
                    assert os.path.isfile(genome), "{0} is not a file".format(genome)
                    Table['genome'].append(os.path.basename(genome))
                    Table['location'].append(os.path.abspath(genome))
            Bdb = pd.DataFrame(Table)
            return Bdb
        except:
            logging.info('Nevermind! Ill try loading as a genome now')

    Table = {'genome':[], 'location':[]}
    for genome in genome_list:
        assert os.path.isfile(genome), "{0} is not a file".format(genome)
        Table['genome'].append(os.path.basename(genome))
        Table['location'].append(os.path.abspath(genome))
    Bdb = pd.DataFrame(Table)
    return Bdb



def add_avani(db):
    '''
    add a column titled 'av_ani' to the passed in dataframe

    dataframe must have rows reference, querey, and ani

    Args:
        db: dataframe
    '''

    logging.debug('making dictionary for average_ani')
    combo2value = {}
    for i, row in db.iterrows():
        combo2value["{0}-vs-{1}".format(row['querry'], row['reference'])] \
            = row['ani']

    logging.debug('list comprehension for average_ani')
    db['av_ani'] = [np.mean([combo2value["{0}-vs-{1}".format(q,r)], \
                        combo2value["{0}-vs-{1}".format(r,q)]]) if r != q else 1\
                        for q, r in zip(db['querry'].tolist(), \
                        db['reference'].tolist())]

    logging.debug('averageing done')

def _nucmer_preset(preset):
    '''
    Return the values from a nucmer 'preset'

    Args:
        preset: either "tight" or "normal"

    Returns:
        list: [c, n_maxgap, n_noextend, n_method]
    '''
    assert preset in ['tight','normal']

    if preset == 'tight':
        return 65, 1, True, 'mum'

    elif preset == 'normal':
        return 65, 90, False, 'mum'

def _gen_nomani_cdb(Cdb, **kwargs):
    '''
    Make a Cdb that looks like the Cdb from secondary clustering

    Args:
        Cdb: Result from Mash clustering

    Returns:
        DataFrame: Cdb
    '''
    c_method = kwargs.get('clusterAlg','single')
    threshold = 1 - kwargs.get('P_ani',.9)

    # Make Cdb look like Cdb
    cdb = Cdb.copy()
    cdb['secondary_cluster'] = ["{0}_0".format(i) for i in cdb['primary_cluster']]
    cdb.rename(columns={'primary_cluster': 'primary_cluster'}, inplace=True)
    cdb['threshold'] = threshold
    cdb['cluster_method'] = c_method
    cdb['comparison_algorithm'] = 'MASH'

    return cdb

def _print_time_estimate(Bdb, Cdb, algorithm, cores):
    '''
    Print an estimate of how long genome comaprisons will take

    Args:
        Bdb: DataFrame with [genome, location]
        Cdb: Clustering DataFrame
        algorthm: algorithm to estimate time with
    '''
    comps = 0
    for bdb, name in iteratre_clusters(Bdb,Cdb):
        g = len(bdb['genome'].unique())
        comps += (g * g)

    time = estimate_time(comps, algorithm)
    time = time / int(cores)
    logging.info("Running {0} {1} comparisons- should take ~ {2:.1f} min".format(\
            comps, algorithm, time))

def _randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))
