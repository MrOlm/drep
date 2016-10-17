#!/usr/bin/env python3

import logging
import glob
import pandas as pd
import os
import sys
import shutil

# !!! This is just for testing purposes, obviously
import sys
sys.path.append('/home/mattolm/Programs/drep/')
import drep_modules.WorkDirectory
import drep_modules as dm
import drep_modules.d_cluster


'''
##################################################################
                A NOTE ABOUT PYTHON VERSIONS


For this to work you need to be able to call both python3 and python2.
To set this up on pyenv, I ran:

pyenv local anaconda2-4.1.0 anaconda3-4.1.0
##################################################################
'''

'''
genomes         = list of genomes to filter

length          = minimum length requirement

completeness     = minimum completeness measurement
contamination   = maximum contamination of checkM

Chdb = already done Chdb
'''

def d_filter_wrapper(wd,**kwargs):
    
    # Load the WorkDirectory.
    logging.info("Loading work directory")
    workDirectory = drep_modules.WorkDirectory.WorkDirectory(wd)
    logging.info(str(workDirectory))
    
    # Validate arguments; figure out what you're going to filter
    bdb, saveAs = validate_arguments(workDirectory, **kwargs)
    
    # Filter by size
    if kwargs.get('length', 0) > 1:
        bdb = filter_bdb_length(bdb, kwargs['length'], verbose=True)
    
    # If a contaminant or completeness threshold exist...
    if (kwargs.get('completeness',0) > 0) or (kwargs.get('contamination',1000) < 1000):
        
        # Run checkM
        if kwargs.get('Chdb',None) != None:
            print("Loading provided CheckM data...")
            Chdb = kwargs.get('Chdb')
            validate_chdb(Chdb, bdb)
        
        elif workDirectory.hasDb('Chdb'):
            print("Loading CheckM data from work directory...")
            Chdb = workDirectory.get_db('Chdb')
            validate_chdb(Chdb, bdb)
        
        else:        
            print("Running CheckM...")
            # Make the folder to house the checkM info
            checkM_loc = workDirectory.location + '/data/checkM/'
            dm.make_dir(checkM_loc,dry=kwargs.get('dry',False),\
                        overwrite=kwargs.get('overwrite',False))
        
            # Make a 'genomes' folder for input to checkM
            checkM_genomes = checkM_loc + 'genomes/'
            dm.make_dir(checkM_genomes,dry=kwargs.get('dry',False),\
                        overwrite=kwargs.get('overwrite',False))
            copy_bdb_loc(bdb, checkM_genomes, extension='.fna')
        
            # Run checkM
            checkM_outfolder = checkM_loc + 'checkM_outdir/'
            Chdb = run_checkM(checkM_genomes, checkM_outfolder, **kwargs)
            validate_chdb(Chdb, bdb)
            
            # Save checkM run
            workDirectory.store_db(Chdb,'Chdb')
            
        # Filter bdb
        bdb = filter_bdb(bdb, Chdb, **kwargs)
    
    # Save Bdb or the new Wdb, depending on arguments above
    loc = workDirectory.location + '/data_tables/'
    if saveAs == 'Bdb':
        bdb.to_csv(loc + 'Bdb.csv')
        workDirectory.data_tables['Bdb'] = bdb
    elif saveAs == 'Wdb':
        bdb.to_csv(loc + 'Wdb.csv')
        workDirectory.data_tables['Wdb'] = bdb
    
    
def filter_bdb(bdb, chdb, **kwargs):
    min_comp = kwargs.get('completeness',0)
    max_con = kwargs.get('contamination',1000)
    start_genomes = list(bdb['genome'].unique())
    assert len(start_genomes) > 0
    
    db = chdb[(chdb['Completeness'] >= min_comp) & (chdb['Contamination'] <= max_con)]
    keep_genomes = list(db['Bin Id'].unique())
    bdb = bdb[bdb['genome'].isin(keep_genomes)]
    
    print("{0}% of genomes passed checkM filtering".format((len(keep_genomes)/len(start_genomes))*100))
    
    return bdb

def filter_bdb_length(bdb, min_length, verbose=False):
    start = len(bdb['location'].unique())
    bdb['length'] =  bdb['location'].map(dm.fasta_length)
    x = bdb[bdb['length'] >= min_length]
    end = len(x['location'].unique())
    
    print("{0}% of genomes passed length filtering".format((end/start)*100))
    
    return x
    
def validate_chdb(Chdb, bdb):
    quit = False
    b_genomes = bdb['genome'].tolist()
    for genome in b_genomes:
        if genome not in Chdb['Bin Id'].tolist():
            print("{0} is not in checkM db".format(genome))
            quit = True
    if quit:
        print("New checkM db needs to be made")
        sys.exit()
    return
    
def validate_arguments(wd,**kwargs):
    '''
    Make sure you either have a genome list, Bdb, or Wdb to filter
    
    If filtering Bdb or a genome list, crash if Cdb exists
    
    Make/get bdb- this is a dataframe LIKE Bdb, except it can be made by a genomelist,
    by loading Bdb, or made from Wdb. This is what we'll filter, and then at the
    end use it to either make a new Bdb or filter the Wdb
    '''

    # Figure out what you're going to filter, and return that db
    if wd.hasDb('Wdb'):
        print("Going to filter Wdb")
        print("Matt- write this part")
        sys.exit()
        return bdb_from_wdb(wd.get_db('Wdb')), 'Wdb'
    elif wd.hasDb('Bdb'):
        if kwargs.get('genomes',None) != None:
            print("Both Bdb and a genome list are found- either don't include\
                    a genome list or start a new work directory!")
            sys.exit()
        if wd.hasDb('Cdb'):
            print("You can't filter this work directory- it's already clustered.\
                    Either choose a winner and filter the winners, or make a new\
                    work directory")
            sys.exit()
        print("Will filter Bdb")
        return wd.get_db('Bdb'), 'Bdb'
    else:
        if kwargs.get('genomes',None) == None:
            print("I don't have anything to filter! Give me a genome list")
            sys.exit()
        print("Will filter the genome list")
        bdb = drep_modules.d_cluster.load_genomes(kwargs['genomes'])
        return bdb, 'Bdb'

def run_checkM(genome_folder,checkm_outf,**kwargs):
    t = str(kwargs.get('processors','6'))
    check_exe = '/home/mattolm/.pyenv/versions/anaconda2-4.1.0/bin/checkm'
    
    # Run checkM initial
    cmd = [check_exe,'lineage_wf',genome_folder,checkm_outf,'-f',\
            checkm_outf + '/results.tsv','--tab_table','-t',str(t),'--pplacer_threads',t]
    logging.info("Running CheckM with command: {0}".format(cmd))
    dm.run_cmd(cmd,shell=False,quiet=False)
    
    # Run checkM again for the better table
    lineage = checkm_outf + 'lineage.ms'
    desired_file = checkm_outf + 'Chdb.tsv'
    cmd = [check_exe,'qa', lineage, checkm_outf, '-f', desired_file, '-t',\
            str(t), '--tab_table','-o', '2']
    logging.info("Running CheckM with command: {0}".format(cmd))
    dm.run_cmd(cmd,shell=False,quiet=False)
    
    # Load table and return it
    chdb = pd.read_table(desired_file,sep='\t')
    return chdb
    
def copy_bdb_loc(bdb,loc,extension=''):
    for genome in bdb['location'].unique():
        shutil.copy2(genome, "{0}{1}{2}".format(loc,os.path.basename(genome),extension))
    
def test_filtering():
    # Get test genomes
	test_genomes = glob.glob(str(os.getcwd()) + '/../test/genomes/*')
	names = sorted([os.path.basename(g) for g in test_genomes])
	assert names == ['Enterococcus_casseliflavus_EC20.fasta', 
	                'Enterococcus_faecalis_T2.fna', 
	                'Enterococcus_faecalis_TX0104.fa', 
	                'Enterococcus_faecalis_YI6-1.fna', 
	                'Escherichia_coli_Sakai.fna']
	
	# Set test directory
	test_directory = str(os.getcwd()) + '/../test/test_backend/'
	dm.make_dir(test_directory,overwrite=True)
	
	# Perform functional test
	d_filter_wrapper(test_directory,genomes=test_genomes,length= 100000,\
	                completeness = 0.75, overwrite=True)
    
if __name__ == '__main__':
	test_filtering()