#!/usr/bin/env python3

import pandas as pd
import os
import glob
import shutil

# !!! This is just for testing purposes, obviously
import sys
sys.path.append('/home/mattolm/Programs/drep/')
import drep_modules as dm

"""
Bdb = pandas DataFrame with the columns genome, location, and infant
Mdb = pandas containing MASH information
"""

def cluster_genomes(Bdb,output_folder,MASH_ANI=90,ANIn=99,dry=False):
    """
    Takes a number of command line arguments and returns a couple pandas dataframes

    The main method of the program
    """
    
    # Make a folder for MASH output
    MASH_folder = output_folder + 'MASH_files/'
    if not dry:
        os.makedirs(MASH_folder)
        if not os.path.exists(MASH_folder):
            os.makedirs(MASH_folder)
    
    # Run MASH clustering
    Mdb, Bdb = perform_mash_clustering(Bdb,MASH_folder,MASH_ANI=90,dry=False)
    
    # Run ANIn clustering
    
    print(Mdb)
    

def perform_mash_clustering(Bdb, MASH_folder, MASH_ANI=90, dry=False):
    """
    Run MASH pairwise within all samples in each infant
    
    Cluster genomes within each infant using MASH
    """
    
    # Run MASH
    Mdb = all_vs_all_MASH(Bdb, MASH_folder, MASH_ANI=90, dry=False)
    
    # Cluster MASH
    
    return Mdb, Bdb
    
def all_vs_all_MASH(Bdb, MASH_folder, MASH_ANI=90, dry=False):
    """
    Run MASH pairwise within all samples in Bdb
    """
    
    # Set up folders
    sketch_folder = MASH_folder + 'sketches/'
    if not dry:
        if os.path.exists(sketch_folder):
	        shutil.rmtree(sketch_folder)
        os.makedirs(sketch_folder)
    
    # Make the MASH sketches
    for fasta in Bdb['location'].unique():
        genome = Bdb['genome'][Bdb['location'] == fasta].tolist()[0]
        cmd = ['/opt/bin/bio/mash','sketch',fasta,'-o',sketch_folder + genome]
        cmd = ' '.join(cmd)
        dm.run_cmd(cmd,dry,True)
        
    # Combine MASH sketches
    cmd = ['/opt/bin/bio/mash','paste',MASH_folder + 'ALL.msh',sketch_folder + '*']
    cmd = ' '.join(cmd)
    dm.run_cmd(cmd,dry,True)
    
    # Calculate distances
    all_file = MASH_folder + 'ALL.msh'
    cmd = ['/opt/bin/bio/mash','dist',all_file,all_file,'>',MASH_folder + 'MASH_table.tsv']
    cmd = ' '.join(cmd)
    dm.run_cmd(cmd,dry,True)
    
    # Make Mdb
    Mdb = pd.DataFrame()
    file = MASH_folder + 'MASH_table.tsv'
    
    table = pd.read_csv(file,sep='\t',header = None)
    table.columns = ['genome1','genome2','dist','p','kmers']
    table['genome1'] = table['genome1'].apply(os.path.basename).str.replace('.fasta', '')
    table['genome1'] = table['genome1'].apply(os.path.basename).str.replace('.fa', '')
    table['genome2'] = table['genome2'].apply(os.path.basename).str.replace('.fasta', '')
    table['genome2'] = table['genome2'].apply(os.path.basename).str.replace('.fa', '')
    Mdb = pd.concat([Mdb,table],ignore_index=True)
    
    return Mdb
    
    

def load_genomes(genome_list):
    Table = {'genome':[],'location':[]}
    
    for genome in genome_list: 
        Table['genome'].append('.'.join(os.path.basename(genome).split('.')[:-1]))
        Table['location'].append(genome)
        
    return pd.DataFrame(Table)

def run_mash():
    pass
    
def test_clustering():
    # Get test genomes
	test_genomes = glob.glob(str(os.getcwd()) + '/../test/genomes/*.fna')
	names = [os.path.basename(g) for g in test_genomes]
	assert names == ['Enterococcus_faecalis_T2.fna', 'Enterococcus_faecalis_TX0104.fna', 'Escherichia_coli_Sakai.fna']
	Bdb = load_genomes(test_genomes)
	
	# Set test directory
	test_directory = str(os.getcwd()) + '/../test/test_backend/'
	if os.path.exists(test_directory):
	    shutil.rmtree(test_directory)
	os.makedirs(test_directory)
	
	# Perform functional test
	cluster_genomes(Bdb,test_directory)

if __name__ == '__main__':
	
	test_clustering()
	
	