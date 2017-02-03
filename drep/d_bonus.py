#!/usr/bin/env python3

import logging
import os
import pandas as pd
import sys
import shutil

import drep.WorkDirectory
import drep as dm
import drep.d_filter as d_filter

def d_bonus_wrapper(wd,**kwargs):
    logging.info("Loading work directory")
    wd = drep.WorkDirectory.WorkDirectory(wd)
    logging.debug(str(wd))

    if kwargs.get('run_tax'):
        logging.info('Running tax')
        run_taxonomy(wd,**kwargs)

def run_taxonomy(wd, **kwargs):
    # Validate arguments- make sure you have everything you need
    Bdb = wd.get_db('Bdb')
    prod_dir = wd.get_dir('prodigal')
    cent_dir = wd.get_dir('centrifuge')
    if wd.hasDb('Tdb') and (kwargs.get('overwrite',False) == False):
        logging.error('Tdb already exists- run with overwrite to overwrite')
        sys.exit()

    # Run prodigal
    d_filter.run_prodigal(Bdb, prod_dir, **kwargs)

    # Run centrifuge
    run_centrifuge(Bdb, prod_dir, cent_dir, **kwargs)

    # Call a centrifuge parser that returns Tdb
    Tdb = parse_centrifuge(Bdb, cent_dir, **kwargs)

    # Add taxonomic info to Bdb
    Bdb = add_taxonomy(Bdb,Tdb)

    # Save Tdb and Bdb
    wd.store_db(Tdb,'Tdb',overwrite=kwargs.get('overwrite',False))
    wd.store_db(Bdb,'Bdb',overwrite=True)

def run_centrifuge(Bdb, prod_dir, cent_dir, **kwargs):
    t = kwargs.get('processors','6')

    cmds = []
    files = []
    for genome in Bdb['genome'].unique():
        genes = "{0}{1}.fna".format(prod_dir, genome)
        cent = "{0}{1}".format(cent_dir, genome)
        if not (os.path.exists("{0}_hits.tsv".format(cent)) and \
                    os.path.exists("{0}_report.tsv".format(cent))):
            cmds.append(gen_centrifuge_cmd(genes,cent,**kwargs))

    if len(cmds) > 1:
        logging.info('Running Centrifuge')
        for cmd in cmds:
            logging.debug(' '.join(cmd))
        drep.d_cluster.thread_mash_cmds_status(cmds,t=int(t))

    else:
        logging.info('Past centrifuge runs found- will not re-run')

def gen_read2bin(gene_files):
    r2b = {}
    for f in gene_files:
        genome = os.path.basename(f)[:-4]
        with open(f) as handle:
          for line in handle:
            if line.startswith('>'):
              r2b.setdefault(genome,[]).append(line.strip()[1:].split(' ')[0])
    return r2b

def parse_centrifuge(Bdb, cent_dir, **kwargs):

    Tdb = pd.DataFrame()
    for genome in Bdb['genome'].unique():

        hits = parse_raw_centrifuge("{0}{1}_hits.tsv".format(cent_dir,genome), \
                    "{0}{1}_report.tsv".format(cent_dir,genome))

        if hits.empty:
            logging.debug("No centrifuge hits found for {0}- skipping".format(genome))
            continue

        x = gen_phylo_db(hits)
        x['genome'] = genome
        Tdb = pd.concat([x,Tdb], ignore_index=True)

    return Tdb

def parse_raw_centrifuge(hits, report):
    tax = pd.read_table(report)
    t2l = tax.set_index('taxID')['taxRank'].to_dict()
    t2n = tax.set_index('taxID')['name'].to_dict()

    hits = pd.read_table(hits)
    hits['gene'] = hits['readID']
    hits['scaffold'] = hits['readID'].map(get_scaff)
    hits['level'] = hits['taxID'].map(t2l)
    hits['name'] = hits['taxID'].map(t2n)
    del hits['readID']

    return hits

def add_taxonomy(Bdb,Tdb):
    g2t = {}
    for genome in Tdb['genome'].unique():
        d = Tdb[Tdb['genome'] == genome]
        tax = d['taxonomy'][d['tax_confidence'] == d['tax_confidence'].max()].tolist()[0]
        g2t[genome] = tax
    Bdb['taxonomy'] = Bdb['genome'].map(g2t)
    return Bdb

def gen_phylo_db(hits):
    Table = {'tax_confidence':[],'taxonomy':[],'tax_level':[],'tax_ID':[]}

    skip = ['uncultured bacterium','Mus musculus', 'Vitis vinifera', 'Homo sapiens']

    levels = ['leaf','subspecies','species','genus','family']
    total_ORFS = len(hits['gene'].unique())
    for level in levels:

        # Restrict to the level in question
        d = hits[hits['level'] == level]

        # Get rid of hits to shitty things
        d = d[~d['name'].isin(skip)]

        # Only get the top hit for each gene
        d = d.sort_values(by='score')
        d = d[~d.duplicated('gene')]

        # Find percentage of top tax
        try:
            top_tax = d['name'].mode()[0]
        except:
            Table['tax_level'].append(level)
            Table['tax_ID'].append(None)
            Table['taxonomy'].append(None)
            Table['tax_confidence'].append(None)
            continue

        x = d[d['name'] == top_tax]
        top_perc = (len(x['gene'].unique()) / total_ORFS) * 100
        tax_ID = x['taxID'].unique()[0]

        Table['tax_level'].append(level)
        Table['tax_ID'].append(tax_ID)
        Table['taxonomy'].append(top_tax)
        Table['tax_confidence'].append(top_perc)

    return pd.DataFrame(Table)

def get_scaff(read):
    return "_".join(read.split('_')[:-1])

def gen_centrifuge_cmd(genes,cent,**kwargs):
    cent_exe = shutil.which('centrifuge')
    if cent_exe == None:
        logging.error("Can't find centrifuge- make sure it's in your system path")
        sys.exit()

    cent_indicies = kwargs.get('cent_index', False)
    if cent_indicies == False:
        logging.error("Can't find centrifuge index- must provide for taxonomy")
        sys.exit()

    cmd = [cent_exe, '-f', '-x', cent_indicies, genes, '-S', "{0}_hits.tsv".format(cent),\
            '-p','1','--report-file',"{0}_report.tsv".format(cent)]

    return cmd

'''
def gen_giant_centrifuge_cmd(genes,base,**kwargs):
    loc = shutil.which('centrifuge')
    if loc == None:
        print("Can't find centrifuge- make sure it's in your system path")
        sys.exit()
    ind = kwargs.get('tax_db','ncbi')
    if ind == 'ncbi':
        cent_indicies = kwargs.get('cent_indicies', '/data3/Human/NIH_4/CentrifugeIndex/nt')
    elif ind == 'bac_only':
        cent_indicies = kwargs.get('cent_exe', '/home/mattolm/download/centrifuge/indices/b+h+v')
    cent_exe = kwargs.get('cent_exe', '/home/mattolm/download/centrifuge/centrifuge')
    p = kwargs.get('processors')

    cmd = [cent_exe, '-f', '-x', cent_indicies, '-U', ','.join(genes), '-S',\
        "{0}_hits.tsv".format(base),'-p',str(p),'--report-file',"{0}_report.tsv".format(base)]

    return cmd
'''

def test_bonus():
    print("Write this you lazy bum")

if __name__ == '__main__':
	test_bonus()
