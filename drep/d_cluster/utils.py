#!/usr/bin/env python3
'''
d_cluster - a subset of dRep

Clusters a list of genomes with both primary and secondary clustering
'''

import pandas as pd
import os
import logging
import random
import string
import numpy as np

import drep
import drep.d_cluster.cluster_utils
import drep.d_cluster.compare_utils
import drep.d_filter
import drep.d_bonus

# from drep.d_cluster import genomeChunk
import drep.d_cluster.controller
#from drep.d_cluster.controller import genomeChunk, genome_hierarchical_clustering, iteratre_clusters, cluster_hierarchical

# This is to make pandas shut up with it's warnings
pd.options.mode.chained_assignment = None

import drep.d_cluster.external
#from drep.d_cluster.external import parse_gani_file, parse_nsim_file, gen_nucmer_cmd, add_avani



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
        cdb, cluster_ret = drep.d_cluster.cluster_utils.genome_hierarchical_clustering(ndb, cluster=name, **kwargs)
        cdb[id] = name
        Cdb = pd.concat([Cdb,cdb], ignore_index=True)
        c2ret[name] = cluster_ret

    return Cdb, c2ret


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
    drep.d_cluster.external.add_avani(d)
    d['dist'] = 1 - d['av_ani']
    db = d.pivot(index="reference", columns="querry", values="dist")

    return db


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


def parse_mash_table(location):
    try:
        iniCols = ['genome1', 'genome2', 'dist', 'p', 'kmers']
        uCols = ['genome1', 'genome2', 'dist']
        dTypes = {'genome1': 'category', 'genome2': 'category', 'dist': np.float32}
        Mdb = pd.read_csv(location, names=iniCols, usecols=uCols, dtype=dTypes, sep='\t')
    except ValueError:
        logging.error("You're using an incompatible version of Mash! I'll try and fix it, but I recommend upgrading Mash to v2.2 or downgrading to v1 if problems occur")
        Mdb = pd.read_csv(location, sep='\t')

    Mdb['genome1'] = Mdb['genome1'].apply(_get_genome_name_from_fasta).astype('category')
    Mdb['genome2'] = Mdb['genome2'].apply(_get_genome_name_from_fasta).astype('category')
    Mdb['similarity'] = 1 - Mdb['dist']

    return Mdb


def merge_genome_chunks(mash_exe, genome_chunks, sketch_folder, MASH_folder):
    all_file = os.path.join(MASH_folder, 'ALL.msh')
    cmd = [mash_exe, 'paste', all_file] + [gc.all_file for gc in genome_chunks]

    new_locs = []
    new_names = []
    for gc in genome_chunks:
        new_locs += gc.genome_locations
        new_names += gc.genome_names

    new_gc = drep.d_cluster.compare_utils.genomeChunk(new_locs, 'all', sketch_folder, new_names, no_create=True)
    new_gc.all_file = all_file

    return cmd, new_gc


# def all_vs_all_MASH(Bdb, data_folder, **kwargs):
#     """
#     Run MASH pairwise within all samples in Bdb
#
#     Args:
#         Bdb: dataframe with genome, location
#         data_folder: location to store output files
#
#     Keyword Args:
#         MASH_sketch: size of mash sketches
#         dry: dont actually run anything
#         processors: number of processors to multithread with
#         mash_exe: location of mash excutible (will try and find with shutil if not provided)
#         groupSize: max number of mash sketches to hold in each folder
#         debug: if True, log all of the commands
#         wd: if you want to log commands, you also need the wd
#     """
#
#     MASH_s = kwargs.get('MASH_sketch',1000)
#     dry = kwargs.get('dry',False)
#     # overwrite = kwargs.get('overwrite', False)
#     p = kwargs.get('processors',6)
#     groupSize = kwargs.get('groupSize', 1000)
#
#     # set up logdir
#     if ('wd' in kwargs) and (kwargs.get('debug', False) == True):
#         logdir = kwargs.get('wd').get_dir('cmd_logs')
#     else:
#         logdir = False
#
#     # Find mash
#     mash_exe = kwargs.get('exe_loc', None)
#     if mash_exe == None:
#         mash_exe = drep.get_exe('mash')
#
#     # Set up folders
#     MASH_folder = os.path.join(data_folder, 'MASH_files/')
#     if not os.path.exists(MASH_folder):
#         os.makedirs(MASH_folder)
#
#     sketch_folder = os.path.join(MASH_folder, 'sketches/')
#     if not os.path.exists(sketch_folder):
#         os.makedirs(sketch_folder)
#
#     # Make chunks
#     l2g = Bdb.set_index('location')['genome'].to_dict()
#     locations = list(Bdb['location'].unique())
#     chunks = [locations[x:x+groupSize] for x in range(0, len(locations), groupSize)]
#
#     # Make the MASH sketches
#     cmds = []
#     chunk_folders = []
#     for i, chunk in enumerate(chunks):
#         chunk_folder = os.path.join(sketch_folder, "chunk_{0}".format(i))
#         chunk_folders.append(chunk_folder)
#         if not os.path.exists(chunk_folder):
#             os.makedirs(chunk_folder)
#         for fasta in chunk:
#             genome = l2g[fasta]
#             file = os.path.join(chunk_folder, genome)
#             if not os.path.isfile(file + '.msh'):
#                 cmd = [mash_exe, 'sketch', fasta, '-s', str(MASH_s), '-o',
#                     file]
#                 cmds.append(cmd)
#
#     if not dry:
#         if len(cmds) > 0:
#             drep.thread_cmds(cmds, logdir=logdir, t=int(p))
#
#     # Combine MASH sketches within chunk
#     cmds = []
#     alls = []
#     for chunk_folder in chunk_folders:
#         all_file = os.path.join(chunk_folder, 'chunk_all.msh')
#         cmd = [mash_exe, 'paste', all_file] \
#                 + glob.glob(os.path.join(chunk_folder, '*'))
#         cmds.append(cmd)
#         alls.append(all_file)
#     if not dry:
#         if len(cmds) > 0:
#             drep.thread_cmds(cmds, logdir=logdir, t=int(p))
#
#     # Combine MASH sketches of all chunks
#     all_file = os.path.join(MASH_folder, 'ALL.msh')
#     cmd = [mash_exe, 'paste', all_file] + alls
#     drep.run_cmd(cmd, dry, shell=False, logdir=logdir)
#
#     # Calculate distances
#     cmd = [mash_exe, 'dist','-p', str(p), all_file, all_file, '>', MASH_folder
#             + 'MASH_table.tsv']
#     cmd = ' '.join(cmd)
#     drep.run_cmd(cmd, dry, shell=True, logdir=logdir)
#
#     # Make Mdb based on all genomes in the MASH folder
#     file = MASH_folder + 'MASH_table.tsv'
#
#     iniCols = ['genome1','genome2','dist','p','kmers']
#     uCols = ['genome1','genome2','dist']
#     dTypes = {'genome1':'category', 'genome2':'category', 'dist':np.float32}
#     Mdb = pd.read_csv(file, names=iniCols, usecols=uCols, dtype=dTypes, sep='\t')
#     Mdb['genome1'] = Mdb['genome1'].apply(_get_genome_name_from_fasta).astype('category')
#     Mdb['genome2'] = Mdb['genome2'].apply(_get_genome_name_from_fasta).astype('category')
#     Mdb['similarity'] = 1 - Mdb['dist']
#
#     # Filter out those genomes that are in the MASH folder but shouldn't be in Mdb
#     genomes = Bdb['genome'].unique()
#     Mdb = Mdb[Mdb['genome1'].isin(genomes)]
#     Mdb = Mdb[Mdb['genome2'].isin(genomes)]
#
#     # Reorder categories to be correct
#     for g in ['genome1', 'genome2']:
#         Mdb[g] = Mdb[g].cat.remove_unused_categories()
#         Mdb[g] = Mdb[g].cat.reorder_categories(sorted((Mdb[g].unique())), ordered=True)
#
#     return Mdb

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
        results = drep.d_cluster.external.parse_nsim_file(file)
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
        results = drep.d_cluster.external.parse_gani_file(file)

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
                cmds.append(drep.d_cluster.external.gen_nucmer_cmd(file_name, g1, g2))

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
    if algorithm != 'fastANI':
        comps = 0
        for bdb, name in drep.d_cluster.cluster_utils.iteratre_clusters(Bdb, Cdb):
            g = len(bdb['genome'].unique())
            comps += (g * g)

        time = estimate_time(comps, algorithm)
        time = time / int(cores)

    else:
        time = 0
        comps = 0
        for bdb, name in drep.d_cluster.cluster_utils.iteratre_clusters(Bdb, Cdb):
            g = len(bdb['genome'].unique())
            time += predict_fastani_time(g, cores)
            comps += (g * g)
        # Convert to minutes
        time = time / 60

    logging.info("Running {0} {1} comparisons- should take ~ {2:.1f} min".format(\
            comps, algorithm, time))


def predict_fastani_time(genomes, cores=6):
    '''
    Forumula is f(x) = a*x**2 + b*x

    http://localhost:8888/doc/tree/Other/OneOffs/dRep_performance_fastANI_timing_1.ipynb
    '''
    a = 0.018969442838018713
    b = 14.604906958328073

    genomes = int(genomes)
    time_6cores = a * genomes ** 2 + b * genomes

    time = time_6cores * (6 / cores)

    return time

def _randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))
