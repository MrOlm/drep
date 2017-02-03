#!/usr/bin/env python3

import logging
import glob
import pandas as pd
import os
import sys
import shutil
import numpy as np

import drep.WorkDirectory
import drep as dm
import drep
import drep.d_cluster
import drep.d_filter


'''
##################################################################
                A NOTE ABOUT PYTHON VERSIONS


For this to work you need to be able to call both python3 and python2.
To set this up on pyenv, I ran:

pyenv local anaconda2-4.1.0 anaconda3-4.1.0
##################################################################
'''

'''

'''
def d_choose_wrapper(wd,**kwargs):
    # Load the WorkDirectory.
    logging.info("Loading work directory")
    workDirectory = drep.WorkDirectory.WorkDirectory(wd)
    logging.debug(str(workDirectory))

    # Make sure you have Cdb
    Cdb = validate_arguments(workDirectory, **kwargs)

    # If Chdb already exists, validate it
    if workDirectory.hasDb('Chdb'):
        Chdb = workDirectory.get_db('Chdb')
        validate_Chdb(Chdb, Cdb)

        logging.info("CheckM db found- will use")

    # If Chdb doesn't exist, make it
    else:
        logging.info("Chdb (CheckM information) not found- running CheckM now")

        bdb = workDirectory.get_db('Bdb')
        Chdb = drep.d_filter.run_checkM_wrapper(bdb, workDirectory, **kwargs)

    # Call a method with Cdb and Chdb, returning Sdb (scored db) and Wdb (winner db)
    Sdb, Wdb = choose_winners(Cdb,Chdb,**kwargs)

    # Save Sdb and Wdb
    workDirectory.store_db(Sdb,'Sdb',overwrite=kwargs.get('overwrite',False))
    workDirectory.store_db(Wdb,'Wdb',overwrite=kwargs.get('overwrite',False))

    # Make a "winning genomes" folder and populate it
    if workDirectory.hasDb('Bdb'):
        logging.debug('saving dereplicated genomes')
        Bdb = workDirectory.get_db('Bdb')

        output_folder = workDirectory.location + '/dereplicated_genomes/'
        if os.path.exists(output_folder):
            logging.debug("{0} already exists: removing and remaking".format(output_folder))
            shutil.rmtree(output_folder)
        dm.make_dir(output_folder,dry=kwargs.get('dry',False),\
                        overwrite=kwargs.get('overwrite',False))

        for genome in Wdb['genome'].unique():
            loc = Bdb['location'][Bdb['genome'] == genome].tolist()[0]
            shutil.copy2(loc, "{0}{1}".format(output_folder,genome))

        logging.info("Done! Dereplicated genomes saved at {0}".format(output_folder))

    else:
        logging.info("Don't have Bdb, so can't populate genomes")


def choose_winners(Cdb,Chdb,**kwargs):
    # Generate Sdb
    genomes = list(Cdb['genome'].unique())
    Sdb = score_genomes(genomes, Chdb, **kwargs)
    logging.debug("Sdb finished")

    # Generate Wdb
    Wdb = pick_winners(Sdb,Cdb)
    logging.debug("Wdb finished")

    return Sdb, Wdb

def pick_winners(Sdb, Cdb):
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


def score_genomes(genomes, Chdb, **kwargs):
    Table = {'genome':[],'score':[]}
    for genome in genomes:
        row = Chdb[Chdb['Bin Id'] == genome]
        score = score_row(row, **kwargs)
        Table['genome'].append(genome)
        Table['score'].append(score)

    Sdb = pd.DataFrame(Table)
    return Sdb

def score_row(row, **kwargs):
    comW =  kwargs.get('completeness_weight',1)
    conW =  kwargs.get('contamination_weight',1)
    n50W =  kwargs.get('N50_weight',1)
    sizeW = kwargs.get('size_weight',1)
    strW =  kwargs.get('strain_heterogeneity_weight',1)

    com = float(row['Completeness'].tolist()[0])
    con = float(row['Contamination'].tolist()[0])
    n50 = float(row['N50 (scaffolds)'].tolist()[0])
    size = float(row['Genome size (bp)'].tolist()[0])
    str = float(row['Strain heterogeneity'].tolist()[0])

    score = (com * comW) - (con * conW) + (np.log10(n50) * n50W) + (np.log10(size) * sizeW) - (strW * str)
    return score


def validate_arguments(wd, **kwargs):
    if not wd.hasDb('Cdb'):
        logging.error("Can't find Cdb- quitting")
        logging.error("Cdb is not found in the work directory- you must run cluster before you ",
                + "can choose")
        sys.exit()

    return wd.get_db('Cdb')

def validate_Chdb(Chdb, Cdb):
    for genome in Cdb['genome'].unique():
        if genome not in Chdb['Bin Id'].tolist():
            logging.error("{0} is missing from Chdb- I'm going to crash now".format(genome))
            sys.exit()

def test_choose():
    print("Write this you lazy bum")

if __name__ == '__main__':
	test_choose()
