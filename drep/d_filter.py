#!/usr/bin/env python3

import logging
import glob
import pandas as pd
import os
import sys
import shutil
import multiprocessing


import drep.WorkDirectory
import drep as dm
import drep.d_cluster
import drep.d_choose


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
    logging.debug("Loading work directory in filter")
    workDirectory = drep.WorkDirectory.WorkDirectory(wd)
    logging.debug(str(workDirectory))

    # Validate arguments; figure out what you're going to filter
    logging.debug("Validating filter arguments")
    bdb, saveAs = validate_arguments(workDirectory, **kwargs)

    # Filter by size
    logging.debug("Filtering genomes by size")
    if kwargs.get('length', 0) > 1:
        bdb = filter_bdb_length(bdb, kwargs['length'], verbose=True)

    # If not skipping CheckM...
    if (not kwargs.get('skipCheckM',False)):

        # Run checkM
        if kwargs.get('Chdb',None) != None:
            logging.debug("Loading provided CheckM data...")
            Chdb = kwargs.get('Chdb')
            validate_chdb(Chdb, bdb)

        elif workDirectory.hasDb('Chdb'):
            logging.debug("Loading CheckM data from work directory...")
            Chdb = workDirectory.get_db('Chdb')
            validate_chdb(Chdb, bdb)

        else:
            logging.debug("Running CheckM")
            Chdb = run_checkM_wrapper(bdb, workDirectory, **kwargs)

        # Filter bdb
        bdb = filter_bdb(bdb, Chdb, **kwargs)

    # Save Bdb or the new Wdb, depending on arguments above
    loc = workDirectory.location + '/data_tables/'
    if saveAs == 'Bdb':
        bdb.to_csv(loc + 'Bdb.csv',index=False)
        workDirectory.data_tables['Bdb'] = bdb

    '''
    elif saveAs == 'Wdb':
        bdb.to_csv(loc + 'Wdb.csv')
        workDirectory.data_tables['Wdb'] = bdb
    '''

def filter_bdb(bdb, chdb, **kwargs):
    min_comp = kwargs.get('completeness',False)
    max_con = kwargs.get('contamination',False)
    min_strain_htr = kwargs.get('strain_htr',False)
    start_genomes = list(bdb['genome'].unique())
    assert len(start_genomes) > 0

    db = chdb.copy()

    if min_comp != False:
        db = db[(db['Completeness'] >= min_comp)]
    if max_con != False:
        db = db[db['Contamination'] <= max_con]
    if min_strain_htr != False:
        db = db[db['Strain heterogeneity'] <= min_strain_htr]
    keep_genomes = list(db['Bin Id'].unique())
    bdb = bdb[bdb['genome'].isin(keep_genomes)]

    logging.info("{0:.2f}% of genomes passed checkM filtering".format((len(keep_genomes)/len(start_genomes))*100))

    return bdb

def filter_bdb_length(bdb, min_length, verbose=False):
    start = len(bdb['location'].unique())
    bdb['length'] =  bdb['location'].map(dm.fasta_length)
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
        logging.warning("NOTE: Wdb already exists! This will not be filtered! Be sure you know what you're doing")
        #sys.exit()
        #return bdb_from_wdb(wd.get_db('Wdb')), 'Wdb'

    if wd.hasDb('Bdb'):
        if kwargs.get('genomes',None) != None:
            logging.error("Both Bdb and a genome list are found- either don't include "\
                    + "a genome list or start a new work directory!")
            sys.exit()
        if wd.hasDb('Cdb'):
            logging.warning("NOTE: Clustering already exists! This will not be filtered! Be sure you know what you're doing")
            #sys.exit()
        logging.info("Will filter Bdb")
        return wd.get_db('Bdb'), 'Bdb'

    else:
        if kwargs.get('genomes',None) == None:
            logging.error("I don't have anything to filter! Give me a genome list")
            sys.exit()
        logging.info("Will filter the genome list")
        bdb = drep.d_cluster.load_genomes(kwargs['genomes'])
        return bdb, 'Bdb'

def run_checkM_wrapper(bdb, workDirectory, **kwargs):
    # Make the folder to house the checkM info
    checkM_loc = workDirectory.location + '/data/checkM/'
    dm.make_dir(checkM_loc,dry=kwargs.get('dry',False),\
                overwrite=kwargs.get('overwrite',False))

    # Run prodigal
    prod_folder = workDirectory.location + '/data/prodigal/'
    if not os.path.exists(prod_folder):
        os.makedirs(prod_folder)
    logging.info("Running prodigal")
    run_prodigal(bdb, prod_folder, **kwargs)

    # Run checkM
    checkM_outfolder = checkM_loc + 'checkM_outdir/'
    Chdb = run_checkM(prod_folder, checkM_outfolder, **kwargs)
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
        g2s[genome] = dm.fasta_length(loc)
        g2n[genome] = calc_n50(loc)
    Chdb['Genome size (bp)'] = Chdb['Bin Id'].map(g2s)
    Chdb['N50 (scaffolds)'] = Chdb['Bin Id'].map(g2n)

    return Chdb

def calc_n50(loc):
    lengths = []
    sequence = []
    with open(loc) as handle:
        for line in handle:
            if line.startswith('>'):
                lengths.append(len(''.join(sequence)))
                sequence = []
            else:
                sequence += line.strip()
    lengths.append(len(''.join(sequence)))

    n50 = sorted(lengths)[int(len(lengths)/2)]
    return n50

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
        drep.d_cluster.thread_mash_cmds_status(cmds,t=int(t))
    else:
        logging.info("Past prodigal runs found- will not re-run")

def run_checkM(genome_folder,checkm_outf,**kwargs):
    t = str(kwargs.get('processors','6'))
    #check_exe = '/home/mattolm/.pyenv/versions/anaconda2-4.1.0/bin/checkm'
    loc = shutil.which('checkm')
    if loc == None:
        logging.error('Cannot locate the program {0}- make sure its in the system path'\
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
    dm.run_cmd(cmd,shell=False,quiet=False)

    # Run checkM again for the better table
    if checkm_method == 'taxonomy_wf':
        lineage = checkm_outf + 'Bacteria.ms'
    else:
        lineage = checkm_outf + 'lineage.ms'
    desired_file = checkm_outf + 'Chdb.tsv'
    cmd = [check_exe,'qa', lineage, checkm_outf, '-f', desired_file, '-t',\
            str(t), '--tab_table','-o', '2']
    logging.debug("Running CheckM with command: {0}".format(cmd))
    dm.run_cmd(cmd,shell=False,quiet=False)

    # Load table and return it
    try:
        chdb = pd.read_table(desired_file,sep='\t')
    except:
        logging.error("!!! checkM failed !!!\nIf using pyenv, make sure both python2 and " +\
            "python3 are available (for example: pyenv pyenv global 3.5.1 2.7.9)")
        sys.exit()
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
