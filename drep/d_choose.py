#!/usr/bin/env python3
'''
d_choose - a subset of drep

Choose best genome from each cluster
'''

import logging
import glob
import pandas as pd
import os
import sys
import shutil
import numpy as np

import drep
import drep.d_cluster
import drep.d_filter
import drep.WorkDirectory

def d_choose_wrapper(wd, **kwargs):
    '''
    Controller for the dRep choose operation

    Based off of the formula:
    A*Completeness - B*Contamination + C*(Contamination *
    (strain_heterogeneity/100)) + D*log(N50) + E*log(size)

    A = completeness_weight;
    B = contamination_weight;
    C = strain_heterogeneity_weight;
    D = N50_weight;
    E = size_weight

    Args:
        wd (WorkDirectory): The current workDirectory
        **kwargs: Command line arguments

    Keyword args:
        genomeInfo: .csv genomeInfo file
        ignoreGenomeQuality: Don't run checkM or do any quality-based filtering (not recommended)
        checkM_method: Either lineage_wf (more accurate) or taxonomy_wf (faster)

        completeness_weight: see formula
        contamination_weight: see formula
        strain_heterogeneity_weight: see formula
        N50_weight: see formula
        size_weight: see formula

    Returns:
        Makes Sdb (scoreDb) in the workDirectory
    '''
    # Load the WorkDirectory.
    logging.info("Loading work directory")
    workDirectory = drep.WorkDirectory.WorkDirectory(wd)
    logging.debug(str(workDirectory))
    wd = workDirectory

    # Make sure you have enough
    kwargs = _validate_choose_arguments(workDirectory, kwargs)
    Cdb = wd.get_db('Cdb')
    bdb = wd.get_db('Bdb')

    # Get comp/con information
    if kwargs.get('ignoreGenomeQuality', False):
        logging.debug("Skipping all quality-based filtering")
        Gdb = drep.d_filter.calc_genome_info(bdb['location'].tolist())
    else:
        Gdb = drep.d_filter._get_run_genomeInfo(wd, bdb, **kwargs)

    if kwargs.get('centrality_weight', 0) > 0:
        Gdb = add_centrality(wd, Gdb, **kwargs)

    wd.store_db(Gdb, 'genomeInformation', overwrite=True)

    # Call a method with Cdb and Chdb, returning Sdb (scored db) and Wdb (winner db)

    Sdb, Wdb = choose_winners(Cdb, Gdb, **kwargs)

    # Save Sdb and Wdb
    workDirectory.store_db(Sdb,'Sdb',overwrite=True)
    workDirectory.store_db(Wdb,'Wdb',overwrite=True)

    # Make a "winning genomes" folder and populate it
    logging.debug('saving dereplicated genomes')
    g2l = bdb.set_index('genome')['location'].to_dict()
    derep_genomes = [g2l[g] for g in list(Wdb['genome'].unique())]
    wd.store_special('dereplicated_genomes', derep_genomes)

def choose_winners(Cdb, Gdb, **kwargs):
    '''
    Make a scoring database and pick the winner of each cluster

    Args:
        Cdb: clustering database
        Gdb: genome information database

    Keyword Args:
        See wrapper

    Returns:
        List: [Sdb (scoring database), Wdb (winner database)]
    '''
    # If you have extra genome weights, load them
    ew_loc = kwargs.get('extra_weight_table', None)
    if ew_loc is not None:
        Edb = load_extra_weight_table(ew_loc, list(Cdb['genome'].unique()), **kwargs)
    else:
        Edb = None

    # Generate Sdb
    genomes = list(Cdb['genome'].unique())
    Sdb = score_genomes(genomes, Gdb, Edb=Edb, **kwargs)
    logging.debug("Sdb finished")

    # Generate Wdb
    Wdb = pick_winners(Sdb,Cdb)
    logging.debug("Wdb finished")

    return Sdb, Wdb

def load_extra_weight_table(loc, genomes, **kwargs):
    """

    Args:
        loc: location of extra weight table
        genomes: list of genomes you have in Cdb
        **kwargs: nothing realld

    Returns:
        dataframe with columns "genome" and "extra_weight"

    """
    if not os.path.join(loc):
        logging.error(f"COULD NOT FIND FILE {loc}; WILL NOT PROCESS EXTRA WEIGHTS")
        return None
    else:
        try:
            db = pd.read_csv(loc, sep='\t', names=['genome', 'extra_weight'])
            db['extra_weight'] = db['extra_weight'].astype(float)
        except:
            f = ''
            with open(loc, 'r') as o:
                x = 0
                for line in o.readlines():
                    f += line + '\n'
                    x += 1

                    if x > 10:
                        break

            logging.error(f"COULD NOT LOAD FILE {loc}; WILL NOT PROCESS EXTRA WEIGHTS. FILE LOOKS LIKE:\n{f}")
            return None

        wg = set(db['genome'].tolist())
        gg = set(genomes)

        missing = "\n".join(list(wg-gg)[:10])

        logging.info(f'Loaded {len(db)} extra weights. {len(gg.intersection(wg))} of {len(gg)} genomes have an extra weight. {len(wg - gg)} genomes have a weight but ARE NOT KNOWN BY DREP; here are some examples:\n{missing}')

        return db[db['genome'].isin(gg)]

def pick_winners(Sdb, Cdb):
    '''
    Based on clustering and scores, pick the best genome from every cluster

    Args:
        Sdb: score of every genome
        Cdb: clustering

    Returns:
        DataFrame: Wdb (winner database)
    '''
    Table = {'genome':[],'cluster':[],'score':[]}
    for cluster in Cdb['secondary_cluster'].unique():
        d = Sdb[Sdb['genome'].isin(Cdb['genome'][Cdb['secondary_cluster'] == cluster].tolist())]
        score = d['score'].max()
        genome = d['genome'][d['score'] == score].tolist()[0]
        Table['genome'].append(genome)
        Table['score'].append(score)
        Table['cluster'].append(cluster)

    Wdb = pd.DataFrame(Table)
    return Wdb

def score_genomes(genomes, Gdb, Edb=None, **kwargs):
    '''
    Calculate the scores for a list of genomes

    Args:
        genomes: list of genomes
        Gdb: genome information database

    Keyword Args:
        See wrapper

    Returns:
        DataFrame: Sdb (scoring database)
    '''

    if Edb is not None:
        g2e = Edb.set_index('genome')['extra_weight'].to_dict()
    else:
        g2e = {}

    Table = {'genome':[],'score':[]}
    for genome in genomes:
        if genome in g2e:
            extra = g2e[genome]
        else:
            extra = 0
        row = Gdb[Gdb['genome'] == genome]

        score = score_row(row, extra=extra, **kwargs)
        Table['genome'].append(genome)
        Table['score'].append(score)

    Sdb = pd.DataFrame(Table)
    return Sdb

def score_row(row, extra=0, **kwargs):
    '''
    Perform the scoring of a row based on kwargs

    Args:
        row: row of genome information

    Keyword Args:
        ignoreGenomeQuality: Don't run checkM or do any quality-based filtering (not recommended)

        completeness_weight: see formula
        contamination_weight: see formula
        strain_heterogeneity_weight: see formula
        N50_weight: see formula
        size_weight: see formula
        extra = extra weight to apply

    Returns:
        float: score
    '''
    comW =  kwargs.get('completeness_weight',1)
    conW =  kwargs.get('contamination_weight',1)
    n50W =  kwargs.get('N50_weight',1)
    sizeW = kwargs.get('size_weight',1)
    strW =  kwargs.get('strain_heterogeneity_weight',1)
    centW = kwargs.get('centrality_weight', 0)

    # For centrality calculations
    S_ani = kwargs.get('S_ani', 0.99)
    if centW > 0:
        cent = row['centrality'].tolist()[0]
    else:
        cent = 0

    n50 = float(row['N50'].tolist()[0])
    size = float(row['length'].tolist()[0])

    if kwargs.get('ignoreGenomeQuality', False):
        score = (np.log10(n50) * n50W) + (np.log10(size) * sizeW) + ((cent - S_ani) * centW) + float(extra)
        return score

    com = float(row['completeness'].tolist()[0])
    con = float(row['contamination'].tolist()[0])

    if 'strain_heterogeneity' in row:
        strh = float(row['strain_heterogeneity'].tolist()[0])
    else:
        strh = 0

    score = (com * comW) - (con * conW) + (strW * (con * (strh/100))) \
        + (np.log10(n50) * n50W) + (np.log10(size) * sizeW) + ((cent - S_ani) * centW) + float(extra)
    return score

def _validate_choose_arguments(wd, kwargs):
    '''
    Validate choose arguments

    Make sure you have a Cdb

    Args:
        wd: WorkDirectory
        kwargs: keyword arguments
    '''
    if not wd.hasDb('Cdb'):
        logging.error("Can't find Cdb- quitting")
        logging.error("Cdb is not found in the work directory- you must run cluster before you ",
                + "can choose")
        sys.exit()

    # Validate centrality arguments
    if (kwargs.get('SkipSecondary', True)) & (kwargs.get('centrality_weight', 0) > 0):
        logging.error(
            "You skipped secondary clustering but have centrality weight above 0. You cant do that. I will now set the centrality weight to 0 to avoid a crash")
        kwargs['centrality_weight'] = 0

    return kwargs

def add_centrality(wd, Gdb, **kwargs):
    """
    Add a columns named "centrality" to genome info
    """
    Ndb = wd.get_db('Ndb')
    Cdb = wd.get_db('Cdb')

    if Cdb['cluster_method'].iloc[0] == 'greedy':
        Ndb = calc_centrality_from_scratch(wd.get_db('Bdb'), Cdb, os.path.join(wd.get_dir('MASH'), 'centrality_calculations/'))

    g2c = Cdb.set_index('genome')['secondary_cluster'].to_dict()
    c2s = Cdb['secondary_cluster'].value_counts().to_dict()

    Ndb['cluster_1'] = Ndb['reference'].map(g2c)
    Ndb['cluster_2'] = Ndb['querry'].map(g2c)
    Ndb = Ndb[Ndb['cluster_1'] == Ndb['cluster_2']]
    Ndb = Ndb[Ndb['reference'] != Ndb['querry']]

    genome2centrality = {}
    for cluster, ndb in Ndb.groupby('cluster_1'):
        #print(f"Cluster {cluster} has {c2s[cluster]} members and {len(ndb)} comps")

        mlen = c2s[cluster]
        assert len(ndb) == (mlen * mlen) - mlen
        for genome, db in ndb.groupby('reference'):
            genome2centrality[genome] = db['ani'].mean()

    Gdb['centrality'] = Gdb['genome'].map(genome2centrality).fillna(0)
    return Gdb

def calc_centrality_from_scratch(Bdb, Cdb, data_folder):
    """
    Calculate centrality from scratch using Mash
    """
    logging.info("Calculating centrality using Mash")

    # 1) Run calculations
    dbs = []
    Xdb = pd.merge(Cdb, Bdb, on='genome', how='left')
    for cluster, bdb in Xdb.groupby('secondary_cluster'):
        if len(bdb) <= 1:
            continue

        df = os.path.join(data_folder, cluster + '/')
        mdb, cdb, cluster_ret = drep.d_cluster.compare_utils.all_vs_all_MASH(bdb, df, MASH_sketch=10000)
        mdb['ani'] = 1 - mdb['dist']
        mdb['cluster'] = cluster
        mdb = mdb.rename(columns={'genome1':'reference', 'genome2':'querry'})
        dbs.append(mdb[['reference', 'querry', 'ani', 'cluster']])

    Mdb = pd.concat(dbs).reset_index(drop=True)
    return Mdb
