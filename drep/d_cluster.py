#!/usr/bin/env python3

import pandas as pd
import os
import glob
import shutil
import multiprocessing
import logging
from subprocess import call
import sys
import json
import scipy.cluster.hierarchy
import scipy.spatial.distance as ssd
import numpy as np
import pickle
import time
import glob

import drep as dm
import drep
import drep.d_filter as dFilter

# This is to make pandas shut up with it's warnings
pd.options.mode.chained_assignment = None

"""
Bdb = pandas DataFrame with the columns genome and location
Mdb = pandas containing raw MASH information
Ndb = pandas containing raw ANIn information
Cdb = pandas containing clustering information (both MASH and ANIn)
"""

""" Program architecture:

*   The main method is cluster_genomes(). This method will take all necessary arguments
    for a normal clustering pipeline, and return all normal outputs

*   The method called by the main program is d_cluster_wrapper(). This method will take
    the raw arguments and workDirectory, parse and extract the needed info, and then call
    cluster_genomes(). It will then take the output information and save it to the
    WorkDirectory object.
"""

def cluster_genomes(Bdb, data_folder, **kwargs):

    """
    Takes a number of command line arguments and returns a couple pandas dataframes

    Required Input:
    * Bdb           - pandas dataframe with the columns "genome" and "location"
    * data_folder   - location where MASH and ANIn data will be stored

    Optional Input:

    ***** Clustering Arguments *****
    * skipMash      - If true, skip MASH altogether (all_vs_all for ANIn)
    * P_Lmethod      - MASH method of determining linkage
                    (to be passed to scipy.cluster.hierarchy.linkage)
    * P_Lcutoff      - MASH linkage cutoff
    * MASH_ANI      - ANI threshold for simple MASH clustering

    * SkipSecondary - If true, skip secondary clustering altogether
    * ANIn_cov      - Coverage threshold for secondary clustering
    * S_Lmethod     - ANIn method of determining linkage
                    (to be passed to scipy.cluster.hierarchy.linkage)
    * S_Lcutoff     - ANIn linkage cutoff


    ***** Comparison Algorithm Arguments *****
    * MASH_s        - MASH sketch size

    * n_c           - nucmer argument c
    * n_maxgap      - nucmer argument maxgap
    * n_noextend    - nucmer argument noextend
    * n_method      - nucmer argument method

    * n_preset      - preset nucmer arrangements: {tight,normal}

    ***** Other Arguments *****
    * dry           - don't actually do anything, just print commands
    * threads       - number of threads to use for multithreading steps
    * overwrite     - overwrite existing data

    Returned Output:
    * Mdb = MASH comparison specifics
    * Ndb = ANIn comparison specifics
    * Cdb = Clustering information
    """

    logging.info("Step 1. Parse Arguments")

    # Deal with nucmer presets
    if kwargs.get('n_preset', None) != None:
        kwargs['n_c'], kwargs['n_maxgap'], kwargs['n_noextend'], kwargs['n_method'] \
        = nucmer_preset(kwargs['n_PRESET'])

    logging.debug("kwargs: {0}".format(kwargs))

    logging.info("Step 2. Perform MASH (primary) clustering")
    if not kwargs.get('SkipMash', False):

        logging.info("2a. Run pair-wise MASH clustering")
        Mdb = all_vs_all_MASH(Bdb, data_folder, **kwargs)

        logging.info("2b. Cluster pair-wise MASH clustering")
        Cdb = cluster_mash_database(Mdb, data_folder= data_folder, **kwargs)

    else:
        logging.info("2. Nevermind! Skipping Mash")
        Cdb = gen_nomash_cdb(Bdb)
        Mdb = pd.DataFrame({'Blank':[]})

    logging.info("{0} primary clusters made".format(len(Cdb['MASH_cluster'].unique())))

    logging.info("Step 3. Perform secondary clustering")

    Ndb, Cdb = run_secondary_clustering(Bdb, Cdb, data_folder, Mdb = Mdb, \
                **kwargs)

    logging.info(
    "Step 4. Return output")

    return Cdb, Mdb, Ndb

def run_secondary_clustering(Bdb, Cdb, data_folder, **kwargs):
    SkipSecondary = kwargs.get('SkipSecondary', False)
    algorithm = kwargs.get('S_algorithm', 'ANIn')

    # This should only happen when called from the wrapper, and Mdb need to be provided
    if SkipSecondary:
        Mdb = kwargs.pop('Mdb')

        Cdb = gen_nomani_cdb(Cdb, Mdb, data_folder = data_folder, **kwargs)
        Ndb = pd.DataFrame({'Blank':[]})

        return Ndb, Cdb

    # Estimate time
    comps = 0
    for bdb, name in iteratre_clusters(Bdb,Cdb):
        g = len(bdb['genome'].unique())
        comps += (g * g)
    time = estimate_time(comps, algorithm)
    time = time / int(kwargs.get('processors'))
    logging.info("Running {0} {1} comparisons- should take ~ {2:.1f} min".format(\
            comps, algorithm, time))

    # Run comparisons
    Ndb = pd.DataFrame()
    for bdb, name in iteratre_clusters(Bdb,Cdb):
        ndb = compare_genomes(bdb, algorithm, data_folder, **kwargs)
        ndb['MASH_cluster'] = name
        Ndb = pd.concat([Ndb,ndb], ignore_index= True)

    # Clear out clustering folder
    c_folder = data_folder + 'Clustering_files/'
    if not os.path.exists(data_folder):
        os.makedirs(data_folder)
    logging.debug('clobbering {0}'.format(data_folder))
    if kwargs.get('overwrite',False):
        for fn in glob.glob(data_folder + 'secondary_linkage_cluster*'):
            os.remove(fn)

    # Run clustering
    Cdb = pd.DataFrame()
    for ndb, name in iteratre_clusters(Ndb,Ndb):
        cdb = genome_hierarchical_clustering(ndb, c_folder, algorithm,\
                cluster=name, **kwargs)
        cdb['primary_cluster'] = name
        Cdb = pd.concat([Cdb,cdb], ignore_index=True)

    return Ndb, Cdb

'''
Description
'''
def genome_hierarchical_clustering(Ndb, data_folder, comp_method, **kwargs):
    logging.debug('Clustering ANIn database')
    S_Lmethod = kwargs.get('clusterAlg', 'single')
    S_Lcutoff = 1 - kwargs.get('S_ani', .99)
    cov_thresh = float(kwargs.get('cov_thresh',0.5))

    cluster = kwargs.get('cluster','')

    Table = {'genome':[],'secondary_cluster':[]}

    # Handle the case where there's only one genome
    if len(Ndb['reference'].unique()) == 1:
        Table['genome'].append(os.path.basename(Ndb['reference'].unique().tolist()[0]))
        Table['secondary_cluster'].append("{0}_0".format(cluster))

    else:
        # 2) Make a linkage-db in algorithm-specific manner
        Ldb = make_linkage_Ndb(Ndb,comp_method,**kwargs)

        # 3) Cluster the linkagedb
        Gdb, linkage = cluster_hierarchical(Ldb, linkage_method= S_Lmethod, \
                                    linkage_cutoff= S_Lcutoff)

        # 4) Extract secondary clusters
        for clust in Gdb['cluster'].unique():
            d = Gdb[Gdb['cluster'] == clust]
            for genome in d['genome'].tolist():
                Table['genome'].append(genome)
                Table['secondary_cluster'].append("{0}_{1}".format(cluster,clust))

        # 5) Save the linkage
        arguments = {'linkage_method':S_Lmethod,'linkage_cutoff':S_Lcutoff,\
                    'comparison_algorithm':comp_method,'minimum_coverage':cov_thresh}
        pickle_name = "secondary_linkage_cluster_{0}.pickle".format(cluster)
        logging.debug('Saving secondary_linkage pickle {1} to {0}'.format(data_folder,\
                                                            pickle_name))
        with open(data_folder + pickle_name, 'wb') as handle:
            pickle.dump(linkage, handle)
            pickle.dump(Ldb,handle)
            pickle.dump(arguments,handle)

    # Return the database
    Gdb = pd.DataFrame(Table)
    Gdb['threshold'] = S_Lcutoff
    Gdb['cluster_method'] = S_Lmethod
    Gdb['comparison_algorithm'] = comp_method

    return Gdb

'''
Filter the Ndb folder in accordance with the kwargs, depending on the algorithm
'''
def make_linkage_Ndb(Ndb,algorithm,**kwargs):
    d = Ndb.copy()

    if algorithm == 'ANIn':
        cov_thresh = float(kwargs.get('cov_thresh',0.5))

        # Remove values without enough coverage
        d.loc[d['alignment_coverage'] <= cov_thresh, 'ani'] = 0

        # Make a linkagedb by averaging values and setting self-compare to 1
        d['av_ani'] = d.apply(lambda row: average_ani (row,d),axis=1)
        d['dist'] = 1 - d['av_ani']
        db = d.pivot("reference", "querry", "dist")

    if algorithm == 'gANI':
        cov_thresh = float(kwargs.get('cov_thresh',0.5))

        # Remove values without enough coverage
        d.loc[d['alignment_coverage'] <= cov_thresh, 'ani'] = 0

        # Make a linkagedb by averaging values and setting self-compare to 1
        d['av_ani'] = d.apply(lambda row: average_ani (row,d),axis=1)
        d['dist'] = 1 - d['av_ani']
        #print(d[d.duplicated()])

        db = d.pivot("reference", "querry", "dist")

    return db

def iteratre_clusters(Bdb, Cdb, id='MASH_cluster'):
    Bdb = pd.merge(Bdb,Cdb)
    for cluster in Bdb[id].unique():
        d = Bdb[Bdb[id] == cluster]
        yield d, cluster

def estimate_time(comps, alg):
    if alg == 'ANIn':
        time = comps * .33
    if alg == 'gANI':
        time = comps * .1
    return time

def d_cluster_wrapper(workDirectory, **kwargs):

    # Load the WorkDirectory.
    logging.debug("Loading work directory")
    workDirectory = drep.WorkDirectory.WorkDirectory(workDirectory)
    logging.debug(str(workDirectory))

    # Parse arguments
    Bdb, data_folder, kwargs = parse_arguments(workDirectory, **kwargs)
    if kwargs['n_PRESET'] != None:
        kwargs['n_c'], kwargs['n_maxgap'], kwargs['n_noextend'], kwargs['n_method'] \
        = nucmer_preset(kwargs['n_PRESET'])

    # Run the main program
    Cdb, Mdb, Ndb = cluster_genomes(Bdb, data_folder, **kwargs)

    # Save the output
    data_dir = workDirectory.location + '/data_tables/'
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    logging.debug("Main program run complete- saving output to {0}".format(data_dir))

    Cdb.to_csv(os.path.join(data_dir,'Cdb.csv'),index=False)
    Mdb.to_csv(os.path.join(data_dir,'Mdb.csv'),index=False)
    Ndb.to_csv(os.path.join(data_dir,'Ndb.csv'),index=False)
    Bdb.to_csv(os.path.join(data_dir,'Bdb.csv'),index=False)

    # Log arguments
    cluster_log = workDirectory.location + '/log/cluster_arguments.json'
    with open(cluster_log, 'w') as fp:
        json.dump(kwargs, fp)
    fp.close()

def parse_arguments(workDirectory, **kwargs):

    # Make sure you have the required program installed
    loc = shutil.which('mash')
    if loc == None:
        logging.error('Cannot locate the program {0}- make sure its in the system path'\
            .format('mash'))
    kwargs['mash_exe'] = loc

    # If genomes are provided, load them
    if kwargs.get('genomes',None) != None:
        assert workDirectory.hasDb("Bdb") == False, \
        "Don't provide new genomes- you already have them in the work directory"
        Bdb = load_genomes(kwargs['genomes'])
    # If genomes are not provided, don't load them
    if kwargs.get('genomes',None) == None:
        assert workDirectory.hasDb("Bdb") != False, \
        "Must either provide a genome list, or run the 'filter' operation with the same work directory"
        Bdb = workDirectory.data_tables['Bdb']

    # Make sure this isn't going to overwrite old data
    overwrite = kwargs.get('overwrite',False)
    for db in ['Cdb','Mdb','Ndb']:
        if workDirectory.hasDb(db):
            if not overwrite:
                logging.error("clustering already exists; run with --overwrite to continue")
                sys.exit()
            logging.debug("THIS WILL OVERWRITE {0}".format(db))


    # Make data_folder
    data_folder = os.path.join(workDirectory.location, 'data/')

    return Bdb, data_folder, kwargs

def cluster_hierarchical(db, linkage_method= 'single', linkage_cutoff= 0.10):

    # Save names
    names = list(db.columns)

    # Generate linkage dataframe
    arr =  np.asarray(db)
    try:
        arr = ssd.squareform(arr)
    except:
        logging.error("The database passed in is not symmetrical!")
        logging.error(arr)
        sys.exit()
    linkage = scipy.cluster.hierarchy.linkage(arr, method= linkage_method)

    # Form clusters
    fclust = scipy.cluster.hierarchy.fcluster(linkage,linkage_cutoff, \
                                                                    criterion='distance')

    # Make Cdb
    Cdb = gen_cdb_from_fclust(fclust,names)

    return Cdb, linkage

def gen_cdb_from_fclust(fclust,names):

    Table={'cluster':[],'genome':[]}
    for i, c in enumerate(fclust):
        Table['cluster'].append(c)
        Table['genome'].append(names[i])

    return pd.DataFrame(Table)


def cluster_anin_database(Cdb, Ndb, data_folder = False, **kwargs):

    logging.debug('Clustering ANIn database')

    cov_thresh = float(kwargs.get('cov_thresh',0.5))
    S_Lmethod = kwargs.get('clusterAlg', 'single')
    S_Lcutoff = 1 - kwargs.get('S_ani', .99)
    overwrite = kwargs.get('overwrite', False)


    if (data_folder != False):
        data_folder = data_folder + 'Clustering_files/'
        if not os.path.exists(data_folder):
            os.makedirs(data_folder)

        # Delete all existing pickles
        logging.debug('clobbering {0}'.format(data_folder))
        if kwargs.get('overwrite',False):
            for fn in glob.glob(data_folder + 'secondary_linkage_cluster*'):
                os.remove(fn)

    Table = {'genome':[],'ANIn_cluster':[]}

    # For every MASH cluster-
    for cluster in Cdb['MASH_cluster'].unique():
        # Filter the database to this cluster
        d = Ndb[Ndb['reference'].isin(Cdb['genome'][Cdb['MASH_cluster'] == cluster].tolist())]

        # Remove values without enough coverage
        d.loc[d['alignment_coverage'] <= cov_thresh, 'ani'] = 0

        # Handle case where cluster has one member
        if len(d['reference'].unique()) == 1:
            Table['genome'].append(d['reference'].unique().tolist()[0])
            Table['ANIn_cluster'].append("{0}_0".format(cluster))
            continue

        # Make a linkagedb
        d['av_ani'] = d.apply(lambda row: average_ani (row,d),axis=1)
        d['dist'] = 1 - d['av_ani']
        db = d.pivot("reference", "querry", "dist")

        Gdb, linkage = cluster_hierarchical(db, linkage_method= S_Lmethod, \
                                    linkage_cutoff= S_Lcutoff)

        # For every ANIn cluster
        for clust in Gdb['cluster'].unique():

            # Filter the database to this cluster
            d = Gdb[Gdb['cluster'] == clust]

            # For every genome in this cluster
            for genome in d['genome'].tolist():

                # Save cluster information
                Table['genome'].append(genome)
                Table['ANIn_cluster'].append("{0}_{1}".format(cluster,clust))

        # Save the linkage
        if (data_folder != False):
            arguments = {'linkage_method':S_Lmethod,'linkage_cutoff':S_Lcutoff,\
                        'comparison_algorithm':'ANIn','minimum_coverage':cov_thresh}
            pickle_name = "secondary_linkage_cluster_{0}.pickle".format(cluster)
            logging.debug('Saving secondary_linkage pickle {1} to {0}'.format(data_folder,\
                                                                pickle_name))
            with open(data_folder + pickle_name, 'wb') as handle:
                pickle.dump(linkage, handle)
                pickle.dump(db,handle)
                pickle.dump(arguments,handle)

    Gdb = pd.DataFrame(Table)
    Cdb = pd.merge(Gdb, Cdb)
    Cdb['threshold'] = S_Lcutoff
    Cdb['cluster_method'] = S_Lmethod
    Cdb['comparison_algorithm'] = 'ANIn'
    Cdb = Cdb.rename(columns={'MASH_cluster':'primary_cluster',\
                        'ANIn_cluster':'secondary_cluster'})

    return Cdb



def cluster_anin_simple(Cdb, Ndb, ANIn=.99, cov_thresh=0.5):

    Table = {'genome':[],'ANIn_cluster':[]}

    # For every MASH cluster-
    for cluster in Cdb['MASH_cluster'].unique():
        # Filter the database to this cluster
        d = Ndb[Ndb['reference'].isin(Cdb['genome'][Cdb['MASH_cluster'] == cluster].tolist())]

        # Make a graph of the genomes in this cluster based on Ndb
        g = make_graph_anin(d,cov_thresh=cov_thresh,anin_thresh=ANIn)
        df = cluster_graph(g)

        # For every ANIn cluster in this graph -
        for clust in df['cluster'].unique():
            # Filter the database to this cluster
            d = df[df['cluster'] == clust]

            # For every genome in this cluster-
            for genome in d['genome'].tolist():

                # Save the cluster information
                Table['genome'].append(genome)
                Table['ANIn_cluster'].append("{0}_{1}".format(cluster,clust))

    Gdb = pd.DataFrame(Table)

    return pd.merge(Gdb,Cdb)

def run_anin_on_clusters(Bdb, Cdb, data_folder, **kwargs):

    """
    For each cluster in Cdb, run pairwise ANIn
    """

    n_c = kwargs.get('n_c', 65)
    n_maxgap = kwargs.get('n_maxgap', 90)
    n_noextend = kwargs.get('n_noextend', False)
    n_method = kwargs.get('method', 'mum')
    p = kwargs.get('processors', 6)
    dry = kwargs.get('dry',False)
    overwrite = kwargs.get('overwrite', False)


    # Set up folders
    ANIn_folder = data_folder + 'ANIn_files/'
    if not os.path.exists(ANIn_folder):
        os.makedirs(ANIn_folder)

    # Add cluster information to Cdb
    Bdb = pd.merge(Bdb,Cdb)
    logging.info("{0} MASH clusters were made".format(len(Bdb['MASH_cluster']\
                .unique().tolist())))


    # Step 1. Make the directories and generate the list of commands to be run
    cmds = []
    files = []
    for cluster in Bdb['MASH_cluster'].unique():
        d = Bdb[Bdb['MASH_cluster'] == cluster]
        genomes = d['location'].tolist()
        for g1 in genomes:
            for g2 in genomes:
                file_name = "{0}{1}_vs_{2}".format(ANIn_folder, \
                            get_genome_name_from_fasta(g1),\
                            get_genome_name_from_fasta(g2))
                files.append(file_name + '.delta')

                # If the file doesn't already exist, add it to what needs to be run
                if not os.path.isfile(file_name + '.delta'):
                    cmds.append(gen_nucmer_cmd(file_name,g1,g2,c=n_c,noextend=n_noextend,\
                                maxgap=n_maxgap,method=n_method))

    # Step 2. Run the nucmer commands

    if not dry:
        if len(cmds) > 0:
            thread_nucmer_cmds_status(cmds,p)

    # Step 3. Parse the nucmer output

    org_lengths = {y:dm.fasta_length(x) for x,y in zip(Bdb['location'].tolist(),Bdb['genome'].tolist())}
    Ndb = process_deltadir(files, org_lengths)
    Ndb['MASH_cluster'] = None
    for cluster in Bdb['MASH_cluster'].unique():
        d = Bdb[Bdb['MASH_cluster'] == cluster]
        Ndb['MASH_cluster'][Ndb['reference'].isin(d['genome'].tolist())] = cluster

    return Ndb

def gen_nucmer_commands(genomes,outf,c=65,maxgap=90,noextend=False,method='mum'):
    cmds = []
    for g1 in genomes:
        for g2 in genomes:
            out = "{0}{1}_vs_{2}".format(outf,get_genome_name_from_fasta(g1),get_genome_name_from_fasta(g2))
            cmds.append(gen_nucmer_cmd(out,g1,g2,c=c,noextend=noextend,maxgap=maxgap,method=method))

    return cmds

def run_nucmer_genomeList(genomes,outf,b2s,c=65,maxgap=90,noextend=False,method='mum',dry=False):
    """
    genomes is a list of locations of genomes in .fasta file. This will do pair-wise
    comparisons of those genomes using the nucmer settings given
    """

    # Run commands on biotite
    cmds = []
    for g1 in genomes:
        for g2 in genomes:
            out = "{0}{1}_vs_{2}".format(outf,get_genome_name_from_fasta(g1),get_genome_name_from_fasta(g2))
            cmds.append(gen_nucmer_cmd(out,g1,g2,c=c,noextend=noextend,maxgap=maxgap,method=method))
    if not dry:
        thread_nucmer_cmds_status(cmds, t=p)

    # Parse resulting folder
    data = process_deltadir(outf, b2s)
    data['c'] = c
    data['maxgap'] = maxgap
    data['noextend'] = noextend
    data['method'] = method

    return data

def all_vs_all_MASH(Bdb, data_folder, **kwargs):
    """
    Run MASH pairwise within all samples in Bdb
    """

    MASH_s = kwargs.get('MASH_sketch',1000)
    dry = kwargs.get('dry',False)
    overwrite = kwargs.get('overwrite', False)
    mash_exe = kwargs.get('mash_exe', 'mash')

    # Set up folders
    MASH_folder = data_folder + 'MASH_files/'
    if not os.path.exists(MASH_folder):
        os.makedirs(MASH_folder)

    sketch_folder = MASH_folder + 'sketches/'
    if not os.path.exists(sketch_folder):
        os.makedirs(sketch_folder)

    # Make the MASH sketches
    cmds = []
    for fasta in Bdb['location'].unique():
        genome = Bdb['genome'][Bdb['location'] == fasta].tolist()[0]
        file = sketch_folder + genome
        if not os.path.isfile(file + '.msh'):
            cmd = [mash_exe, 'sketch', fasta, '-s', str(MASH_s), '-o',
                file]
            cmds.append(cmd)

    if not dry:
        if len(cmds) > 0:
            thread_mash_cmds_status(cmds)

    # Combine MASH sketches
    cmd = [mash_exe, 'paste', MASH_folder + 'ALL.msh', sketch_folder+ '*']
    cmd = ' '.join(cmd)
    dm.run_cmd(cmd,dry,True)

    # Calculate distances
    all_file = MASH_folder + 'ALL.msh'
    cmd = [mash_exe, 'dist', all_file, all_file, '>', MASH_folder
            + 'MASH_table.tsv']
    cmd = ' '.join(cmd)
    dm.run_cmd(cmd,dry,True)

    # Make Mdb based on all genomes in the MASH folder
    Mdb = pd.DataFrame()
    file = MASH_folder + 'MASH_table.tsv'

    table = pd.read_csv(file,sep='\t',header = None)
    table.columns = ['genome1','genome2','dist','p','kmers']
    table['genome1'] = table['genome1'].apply(get_genome_name_from_fasta)
    table['genome2'] = table['genome2'].apply(get_genome_name_from_fasta)
    Mdb = pd.concat([Mdb,table],ignore_index=True)
    Mdb['similarity'] = 1 - Mdb['dist'].astype(float)

    # Filter out those genomes that are in the MASH folder but shouldn't be in Mdb
    genomes = Bdb['genome'].unique()
    Mdb = Mdb[Mdb['genome1'].isin(genomes)]
    Mdb = Mdb[Mdb['genome2'].isin(genomes)]

    return Mdb

def cluster_mash_database(db, data_folder= False, **kwargs):

    logging.debug('Clustering MASH database')

    P_Lmethod = kwargs.get('clusterAlg','single')
    P_Lcutoff = 1 - kwargs.get('P_ani',.9)
    dry = kwargs.get('dry',False)
    overwrite = kwargs.get('overwrite', False)

    if data_folder != False:
        data_folder = data_folder + 'Clustering_files/'
        if not os.path.exists(data_folder):
            os.makedirs(data_folder)

    db['dist'] = 1 - db['similarity']
    linkage_db = db.pivot("genome1","genome2","dist")
    Cdb, linkage = cluster_hierarchical(linkage_db, linkage_method= P_Lmethod, \
                                linkage_cutoff= P_Lcutoff)
    Cdb = Cdb.rename(columns={'cluster':'MASH_cluster'})

    if (data_folder != False):
        arguments = {'linkage_method':P_Lmethod,'linkage_cutoff':P_Lcutoff,\
                        'comparison_algorithm':'MASH'}
        logging.debug('Saving primary_linkage pickle to {0}'.format(data_folder))
        with open(data_folder + 'primary_linkage.pickle', 'wb') as handle:
            pickle.dump(linkage, handle)
            pickle.dump(linkage_db, handle)
            pickle.dump(arguments, handle)

    return Cdb

def get_genome_name_from_fasta(fasta):
    return str(os.path.basename(fasta))

# Parse NUCmer delta file to get total alignment length and total sim_errors
def parse_delta(filename):
    """Returns (alignment length, similarity errors) tuple from passed .delta.
    - filename - path to the input .delta file
    Extracts the aligned length and number of similarity errors for each
    aligned uniquely-matched region, and returns the cumulative total for
    each as a tuple.
    """
    aln_length, sim_errors = 0, 0
    for line in [l.strip().split() for l in open(filename, 'rU').readlines()]:
        if line[0] == 'NUCMER' or line[0].startswith('>'):  # Skip headers
            continue
        # We only process lines with seven columns:
        if len(line) == 7:
            aln_length += abs(int(line[1]) - int(line[0]))
            sim_errors += int(line[4])
    return aln_length, sim_errors

def gen_gANI_cmd(file, g1, g2, dir, exe):
    # Handle self comparison
    # Did a pretty exhaustive test- same comaprisons always give 100% ANI
    if g1 == g2:
        # make a copy of g1
        #g1T = g1 + '.GANI_TEMP'
        #shutil.copyfile(g1,g1T)

        #cmd = [exe,'-genome1fna',g1T,'-genome2fna',g2,'-outfile',file,'-outdir',file + 'TEMP']
        return []
    else:
        cmd = [exe,'-genome1fna',g1,'-genome2fna',g2,'-outfile',file,'-outdir',file + 'TEMP']
    #cmd = [exe,'-genome1fna',g1,'-genome2fna',g2,'-outfile',file,'-outdir',dir]
    return cmd


def process_deltadir(deltafiles, org_lengths, logger=None):
    """Returns a tuple of ANIm results for .deltas in passed directory.
    - delta_dir - path to the directory containing .delta files
    - org_lengths - dictionary of total sequence lengths, keyed by sequence
    Returns the following pandas dataframes in a tuple; query sequences are
    rows, subject sequences are columns:
    - alignment_lengths - symmetrical: total length of alignment
    - percentage_identity - symmetrical: percentage identity of alignment
    - alignment_coverage - non-symmetrical: coverage of query and subject
    - similarity_errors - symmetrical: count of similarity errors
    May throw a ZeroDivisionError if one or more NUCmer runs failed, or a
    very distant sequence was included in the analysis.
    """
    # Process directory to identify input files
    #deltafiles = glob.glob(delta_dir + '*.delta')

    Table = {'querry':[],'reference':[],'alignment_length':[],'similarity_errors':[],
            'ref_coverage':[],'querry_coverage':[],'ani':[], 'reference_length':[],
            'querry_length':[],'alignment_coverage':[]}

    # Process .delta files assuming that the filename format holds:
    # org1_vs_org2.delta
    zero_error = False  # flag to register a divide-by-zero error
    for deltafile in deltafiles:
        qname, sname = os.path.splitext(os.path.split(deltafile)[-1])[0].split('_vs_')
        tot_length, tot_sim_error = parse_delta(deltafile)
        if tot_length == 0 and logger is not None:
            logging.info("Total alignment length reported in " +
                               "%s is zero!" % deltafile)
        query_cover = float(tot_length) / org_lengths[qname]
        sbjct_cover = float(tot_length) / org_lengths[sname]
        # Calculate percentage ID of aligned length. This may fail if
        # total length is zero.
        # The ZeroDivisionError that would arise should be handled
        # Common causes are that a NUCmer run failed, or that a very
        # distant sequence was included in the analysis.
        try:
            perc_id = 1 - float(tot_sim_error) / tot_length
        except ZeroDivisionError:
            #print("Alignment between {0} and {1} has 0 alignment!".format(qname,sname))
            perc_id = 0  # set arbitrary value of zero identity
            zero_error = True

        Table['querry'].append(qname)
        Table['querry_length'].append(org_lengths[qname])
        Table['reference'].append(sname)
        Table['reference_length'].append(org_lengths[sname])
        Table['alignment_length'].append(tot_length)
        Table['similarity_errors'].append(tot_sim_error)
        Table['ani'].append(perc_id)
        Table['ref_coverage'].append(sbjct_cover)
        Table['querry_coverage'].append(query_cover)
        Table['alignment_coverage'].append((tot_length * 2)/(org_lengths[qname]\
                                                             + org_lengths[sname]))

    df = pd.DataFrame(Table)
    return df

def process_deltafiles(deltafiles, org_lengths, logger=None, **kwargs):

    Table = {'querry':[],'reference':[],'alignment_length':[],'similarity_errors':[],
            'ref_coverage':[],'querry_coverage':[],'ani':[], 'reference_length':[],
            'querry_length':[],'alignment_coverage':[]}

    # Process .delta files assuming that the filename format holds:
    # org1_vs_org2.delta
    coverage_method = kwargs.get('coverage_method')
    logging.debug('coverage_method is {0}'.format(coverage_method))

    for deltafile in deltafiles:
        qname, sname = os.path.splitext(os.path.split(deltafile)[-1])[0].split('_vs_')
        tot_length, tot_sim_error = parse_delta(deltafile)
        if tot_length == 0 and logger is not None:
            logging.info("Total alignment length reported in " +
                               "%s is zero!" % deltafile)
        query_cover = float(tot_length) / org_lengths[qname]
        sbjct_cover = float(tot_length) / org_lengths[sname]
        # Calculate percentage ID of aligned length. This may fail if
        # total length is zero.
        # The ZeroDivisionError that would arise should be handled
        # Common causes are that a NUCmer run failed, or that a very
        # distant sequence was included in the analysis.
        try:
            perc_id = 1 - float(tot_sim_error) / tot_length
        except ZeroDivisionError:
            #print("Alignment between {0} and {1} has 0 alignment!".format(qname,sname))
            perc_id = 0  # set arbitrary value of zero identity
            zero_error = True

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

    df = pd.DataFrame(Table)
    return df

### MAKE IT SO THAT YOU REMOVE THE _TEMP MARKER FROM THE GENOMES, AND DELTE THE
### FILE WHILE YOU'RE AT IT
def process_gani_files(files):
    Table = {'querry':[],'reference':[],'ani':[],'alignment_coverage':[]}
    for file in files:
        results = parse_gani_file(file)

        # Add forward results
        Table['reference'].append(results['reference'])
        Table['querry'].append(results['querry'])
        Table['ani'].append(float(results['rq_ani'])/100)
        Table['alignment_coverage'].append(results['rq_coverage'])

        # Add reverse results [they're the same]
        Table['reference'].append(results['querry'])
        Table['querry'].append(results['reference'])
        Table['ani'].append(float(results['qr_ani'])/100)
        Table['alignment_coverage'].append(results['qr_coverage'])


    # Add self comparisons
    for g in set(Table['reference']):
        Table['reference'].append(g)
        Table['querry'].append(g)
        Table['ani'].append(1)
        Table['alignment_coverage'].append(1)

    Gdb = pd.DataFrame(Table)
    return Gdb

'''
This method takes in bdb (a table with the columns location and genome), runs
pair-wise comparisons between all genomes in the sample, and returns a table
with at least the columns 'reference', 'querry', 'ani','coverage', depending
on what algorithm is called
'''
def compare_genomes(bdb, algorithm, data_folder, **kwargs):
    # To handle other versions of this method which passed in a WorkDirectory
    # instead of data_folder string
    if isinstance(data_folder,drep.WorkDirectory.WorkDirectory):
        data_folder = data_folder.location + '/data/'

    if algorithm == 'ANIn':
        genome_list = bdb['location'].tolist()
        working_data_folder = data_folder + 'ANIn_files/'
        df = run_pairwise_ANIn(genome_list, working_data_folder, **kwargs)
        return df

    elif algorithm == 'gANI':
        working_data_folder = data_folder + 'gANI_files/'
        prod_folder = data_folder + 'prodigal/'
        df = run_pairwise_gANI(bdb, working_data_folder, \
                prod_folder = prod_folder, **kwargs)
        return df

    elif algorithm == 'mauve':
        working_data_folder = data_folder + 'mauve_files/'
        df = run_pairwise_mauve(bdb, working_data_folder, **kwargs)
        return df

    else:
        logging.error("{0} not supported".format(algorithm))
        sys.exit()

def run_pairwise_ANIn(genome_list, ANIn_folder, **kwargs):
    p = kwargs.get('processors',6)
    genomes = genome_list

    # Make folder
    if not os.path.exists(ANIn_folder):
        os.makedirs(ANIn_folder)

    # Gen commands
    cmds = []
    files = []
    for g1 in genomes:
        for g2 in genomes:
            file_name = "{0}{1}_vs_{2}".format(ANIn_folder, \
                        get_genome_name_from_fasta(g1),\
                        get_genome_name_from_fasta(g2))
            files.append(file_name)

            # If the file doesn't already exist, add it to what needs to be run
            if not os.path.isfile(file_name + '.delta'):
                cmds.append(gen_nucmer_cmd(file_name,g1,g2))

    # Run commands
    if len(cmds) > 0:
        thread_nucmer_cmds_status(cmds,p,verbose=False)

    # Parse output
    org_lengths = {get_genome_name_from_fasta(y):dm.fasta_length(x) \
                    for x,y in zip(genomes,genomes)}
    deltafiles = ["{0}.delta".format(file) for file in files]
    df = process_deltafiles(deltafiles, org_lengths, **kwargs)

    return df

def run_pairwise_mauve(bdb, data_folder, **kwargs):
    p = kwargs.get('processors',6)

    # Make folder
    if not os.path.exists(data_folder):
        os.makedirs(data_folder)

    # Gen commands
    cmds = []
    files = []
    for g1 in genomes:
        for g2 in genomes:
            '''
            file_name = "{0}{1}_vs_{2}".format(ANIn_folder, \
                        get_genome_name_from_fasta(g1),\
                        get_genome_name_from_fasta(g2))
            files.append(file_name)

            # If the file doesn't already exist, add it to what needs to be run
            if not os.path.isfile(file_name + '.delta'):
                cmds.append(gen_nucmer_cmd(file_name,g1,g2))
            '''

    # Run commands
    if len(cmds) > 0:
        t#hread_nucmer_cmds_status(cmds,p,verbose=False)

    # Parse output
    '''
    org_lengths = {get_genome_name_from_fasta(y):dm.fasta_length(x) \
                    for x,y in zip(genomes,genomes)}
    deltafiles = ["{0}.delta".format(file) for file in files]
    df = process_deltafiles(deltafiles, org_lengths)
    '''

    return df

def run_pairwise_gANI(bdb, gANI_folder, verbose = False, **kwargs):
    gANI_exe = kwargs.get('gANI_exe','ANIcalculator')
    p = kwargs.get('processors',6)
    prod_folder = kwargs.get('prod_folder')
    genomes = bdb['location'].tolist()

    # Make folder
    if not os.path.exists(gANI_folder):
        os.makedirs(gANI_folder)

    # Remove crap folders- they shouldn't exist and if they do it messes things up
    crap_folders = glob.glob(gANI_folder + '*.gANITEMP')
    for crap_folder in crap_folders:
        if os.path.exists(crap_folder):
            logging.debug("CRAP FOLDER EXISTS FOR gANI- removing {0}".format(crap_folder))
            shutil.rmtree(crap_folder)

    if not os.path.exists(prod_folder):
        os.makedirs(prod_folder)

    # Run prodigal
    if verbose:
        logging.info("Running prodigal...")
    dFilter.run_prodigal(bdb, prod_folder, verbose=verbose)

    # Gen gANI commands
    if verbose:
        logging.info("Running gANI...")
    cmds = []
    files = []
    for i, g1 in enumerate(genomes):
        for j, g2 in enumerate(genomes):
            if i > j:
                name1= get_genome_name_from_fasta(g1)
                name2= get_genome_name_from_fasta(g2)
                file_name = "{0}{1}_vs_{2}.gANI".format(gANI_folder, \
                            name1, name2)
                files.append(file_name)

                # If the file doesn't already exist, add it to what needs to be run
                if not os.path.isfile(file_name):
                    fna1 = "{0}{1}.fna".format(prod_folder,name1)
                    fna2 = "{0}{1}.fna".format(prod_folder,name2)
                    cmds.append(gen_gANI_cmd(file_name,fna1,fna2,gANI_folder,gANI_exe))

    # Run commands
    if len(cmds) > 0:
        logging.debug('Running gANI commands: {0}'.format('\n'.join([' '.join(x) for x in cmds])))
        thread_mash_cmds_status(cmds,p)
    else:
        if verbose:
            logging.info("gANI already run- will not re-run")

    # Parse output
    df = process_gani_files(files)

    # Handle self-comparisons
    #df['reference'] = [i.replace('.fna.GANI_','') for i in df['reference']]
    #df['querry'] = [i.replace('.fna.GANI_','') for i in df['querry']]
    #df = df.drop_duplicates()

    # Add self-comparisons if these is only one genome
    if len(genomes) == 1:
        Table = {'querry':[],'reference':[],'ani':[],'alignment_coverage':[]}
        for g in genomes:
            Table['reference'].append(g)
            Table['querry'].append(g)
            Table['ani'].append(1)
            Table['alignment_coverage'].append(1)
        d = pd.DataFrame(Table)
        df = pd.concat([df,d],ignore_index=True)

    return df

def parse_gani_file(file):
    x = pd.read_table(file)
    x = x.rename(columns={'GENOME1':'reference','GENOME2':'querry','AF(1->2)':'rq_coverage',\
                        'AF(2->1)':'qr_coverage','ANI(1->2)':'rq_ani','ANI(2->1)':'qr_ani'})
    dict = x.to_dict(orient='list')
    dict = {k:v[0] for k,v in dict.items()}
    dict['reference'] = dict['reference'][:-4]
    dict['querry'] = dict['querry'][:-4]
    return dict

def gen_nucmer_cmd(prefix,ref,querry,c='65',noextend=False,maxgap='90',method='mum'):
    cmd = ['nucmer','--' + method,'-p',prefix,'-c',str(c),'-g',str(maxgap)]
    if noextend: cmd.append('--noextend')
    cmd += [ref,querry]
    return cmd

def thread_nucmer_cmds(cmds,t=10):
    pool = multiprocessing.Pool(processes=t)
    pool.map(run_nucmer_cmd,cmds)
    pool.close()
    pool.join()
    return

def thread_nucmer_cmds_status(cmds,t=10,verbose=True):
    total = len(cmds)
    if verbose:
        minutes = ((float(total) * (float(0.33))) /float(t))
        logging.info("Running {0} mummer comparisons: should take ~ {1:.1f} min".format(total,minutes))
    pool = multiprocessing.Pool(processes=int(t))
    rs = pool.map_async(run_nucmer_cmd,cmds)
    pool.close()

    pool.join()
    return

def thread_mash_cmds_status(cmds,t=10):
    total = len(cmds)
    pool = multiprocessing.Pool(processes=t)
    rs = pool.map_async(run_nucmer_cmd,cmds)
    pool.close()

    while(True):
        done = total - (rs._number_left * rs._chunksize)
        percR = (done/total) * 100
        sys.stdout.write('\r')
        sys.stdout.write("[{0:20}] {1:3.2f}%".format('='*int(percR/5), percR))
        sys.stdout.flush()
        if (rs.ready()):
            sys.stdout.write('\n')
            break
        time.sleep(0.5)

    pool.join()
    return

def gen_nomash_cdb(Bdb):
    Cdb = Bdb.copy()
    Cdb['MASH_cluster'] = 0
    return Cdb

def run_nucmer_cmd(cmd,dry=False,shell=False):
    devnull = open(os.devnull, 'w')
    if shell:
        if not dry: call(cmd,shell=True, stderr=devnull,stdout=devnull)
        else: print(cmd)
    else:
        if not dry: call(cmd, stderr=devnull,stdout=devnull)
        else: print(' '.join(cmd))
    return

def run_mash_cmd(cmds):
    cmd = ' '.join(cmd)
    dm.run_cmd(cmd,dry=False,shell=True)

def load_genomes(genome_list):
    Table = {'genome':[],'location':[]}

    for genome in genome_list:
        assert os.path.isfile(genome), "{0} is not a file".format(genome)
        Table['genome'].append(os.path.basename(genome))
        Table['location'].append(os.path.abspath(genome))

    Bdb = pd.DataFrame(Table)
    return Bdb

def average_ani(row,db):
    g1 = row['querry']
    g2 = row['reference']
    if g1 == g2:
        return 1
    else:
        ani1 = float(row['ani'])
        ani2 = float(db['ani'][(db['reference'] == g1) & (db['querry'] == g2)].tolist()[0])
        avg = np.mean([ani1,ani2])
        return avg

def nucmer_preset(preset):
   #nucmer argument c, n_maxgap, n_noextend, n_method

    assert preset in ['tight','normal']

    if preset == 'tight':
        return 65, 1, True, 'mum'

    elif preset == 'normal':
        return 65, 90, False, 'mum'

def gen_nomani_cdb(Cdb, Mdb, **kwargs):
    c_method = kwargs.get('clusterAlg','single')
    threshold = 1 - kwargs.get('P_ani',.9)
    data_folder = kwargs.get('data_folder') + 'Clustering_files/'

    # Make Cdb look like Cdb
    cdb = Cdb.copy()
    cdb['secondary_cluster'] = ["{0}_0".format(i) for i in cdb['MASH_cluster']]
    cdb.rename(columns={'MASH_cluster': 'primary_cluster'}, inplace=True)
    cdb['threshold'] = threshold
    cdb['cluster_method'] = c_method
    cdb['comparison_algorithm'] = 'MASH'

    # Delete any only secondary clusters
    if kwargs.get('overwrite',False):
        for fn in glob.glob(data_folder + 'secondary_linkage_cluster*'):
            os.remove(fn)

    return cdb

def test_clustering():
    # Get test genomes
    test_genomes = glob.glob(str(os.getcwd()) + '/../test/genomes/*')
    names = [os.path.basename(g) for g in test_genomes]
    assert names == ['Enterococcus_faecalis_T2.fna',
	                 'Escherichia_coli_Sakai.fna',
	                 'Enterococcus_casseliflavus_EC20.fasta',
	                 'Enterococcus_faecalis_TX0104.fa',
	                 'Enterococcus_faecalis_YI6-1.fna']
    Bdb = load_genomes(test_genomes)

    # Set test directory
    test_directory = str(os.getcwd()) + '/../test/test_backend/'
    dm.make_dir(test_directory,overwrite=True)

    # Perform functional test
    Cdb, Mdb, Ndb = cluster_genomes(Bdb,test_directory)

    # Confirm it's right by showing there is one group of 3 and two groups of one
    group_lengths = sorted([len(Cdb['genome'][Cdb['ANIn_cluster'] == x].tolist()) for x in Cdb['ANIn_cluster'].unique()])
    assert group_lengths == [1, 1, 3]

    print("Functional test success!")

if __name__ == '__main__':

    test_clustering()
