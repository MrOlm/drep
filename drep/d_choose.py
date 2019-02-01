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
    _validate_choose_arguments(workDirectory, **kwargs)
    Cdb = wd.get_db('Cdb')
    bdb = wd.get_db('Bdb')

    # Get comp/con information
    if kwargs.get('ignoreGenomeQuality', False):
        logging.debug("Skipping all quality-based filtering")
        Gdb = drep.d_filter.calc_genome_info(bdb['location'].tolist())
    else:
        Gdb = drep.d_filter._get_run_genomeInfo(wd, bdb, **kwargs)
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
    # Generate Sdb
    genomes = list(Cdb['genome'].unique())
    Sdb = score_genomes(genomes, Gdb, **kwargs)
    logging.debug("Sdb finished")

    # Generate Wdb
    Wdb = pick_winners(Sdb,Cdb)
    logging.debug("Wdb finished")

    return Sdb, Wdb

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

def score_genomes(genomes, Gdb, **kwargs):
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

    Table = {'genome':[],'score':[]}
    for genome in genomes:
        row = Gdb[Gdb['genome'] == genome]
        score = score_row(row, **kwargs)
        Table['genome'].append(genome)
        Table['score'].append(score)

    Sdb = pd.DataFrame(Table)
    return Sdb

def score_row(row, **kwargs):
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

    Returns:
        float: score
    '''
    comW =  kwargs.get('completeness_weight',1)
    conW =  kwargs.get('contamination_weight',1)
    n50W =  kwargs.get('N50_weight',1)
    sizeW = kwargs.get('size_weight',1)
    strW =  kwargs.get('strain_heterogeneity_weight',1)

    if kwargs.get('ignoreGenomeQuality', False):
        n50 = float(row['N50'].tolist()[0])
        size = float(row['length'].tolist()[0])

        score = (np.log10(n50) * n50W) + (np.log10(size) * sizeW)
        return score

    com = float(row['completeness'].tolist()[0])
    con = float(row['contamination'].tolist()[0])
    n50 = float(row['N50'].tolist()[0])
    size = float(row['length'].tolist()[0])

    if 'strain_heterogeneity' in row:
        strh = float(row['strain_heterogeneity'].tolist()[0])
    else:
        strh = 0

    score = (com * comW) - (con * conW) + (strW * (con * (strh/100))) \
        + (np.log10(n50) * n50W) + (np.log10(size) * sizeW)
    return score

def _validate_choose_arguments(wd, **kwargs):
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

    return wd.get_db('Cdb')
