#!/usr/bin/env python3

import logging
import os
import pandas as pd
import sys
import shutil
import subprocess
import time

import drep
import drep.WorkDirectory
import drep.d_filter

def d_bonus_wrapper(wd,**kwargs):
    logging.info("Loading work directory")
    wd = drep.WorkDirectory.WorkDirectory(wd)
    logging.debug(str(wd))

    if kwargs.get('check_dependencies'):
        logging.info('Checking dependencies')
        check_dependencies()

    if kwargs.get('run_tax'):
        logging.info('Running tax')
        run_taxonomy(wd,**kwargs)

def run_taxonomy(wd, **kwargs):
    # Validate arguments- make sure you have everything you need
    Bdb, prod_dir, cent_dir = validate_arguments(wd, **kwargs)

    # Run prodigal
    drep.d_filter.run_prodigal(Bdb['location'].tolist(), prod_dir, **kwargs)

    # Run centrifuge
    run_centrifuge(Bdb, prod_dir, cent_dir, wd=wd, **kwargs)

    # # Call a centrifuge parser that returns Tdb
    # Tdb = parse_centrifuge(Bdb, cent_dir, **kwargs)
    #
    # # Add taxonomic info to Bdb
    # Bdb = add_taxonomy(Bdb,Tdb)

    # Parse taxonomy
    Tdb, Bdb = parse_taxonomy(Bdb, cent_dir, **kwargs)

    # Save Tdb and Bdb
    wd.store_db(Tdb,'Tdb',overwrite=True)
    wd.store_db(Bdb,'Bdb',overwrite=True)

def parse_taxonomy(Bdb, cent_dir, **kwargs):
    '''
    take the centrifuge directory and Bdb, return Tdb and Bdb with an added 'taxonomy' column
    '''

    method = kwargs.get('tax_method')

    if method == 'max':
        Tdb = parse_centrifuge(Bdb, cent_dir, **kwargs)

    elif method == 'percent':
        Tdb = parse_centrifuge_percent(Bdb, cent_dir, **kwargs)

    else:
        logging.error("dont recognize method {0}, quitting".format(method))
        sys.exit()

    Bdb = add_taxonomy(Bdb,Tdb)
    return Tdb, Bdb

def validate_arguments(wd, **kwargs):
    '''
    make sure you have everything you need
    '''
    if wd.hasDb('Bdb'):
        if kwargs.get('genomes',None) != None:
            logging.error("Both Bdb and a genome list are found- either don't include "\
                    + "a genome list or start a new work directory!")
            sys.exit()
        Bdb = wd.get_db('Bdb')

    else:
        if kwargs.get('genomes',None) == None:
            logging.error("I don't have anything to determine the taxonomy of! Give me a genome list")
            sys.exit()
        Bdb = drep.d_cluster.utils.load_genomes(kwargs['genomes'])

    prod_dir = wd.get_dir('prodigal')
    cent_dir = wd.get_dir('centrifuge')
    # if wd.hasDb('Tdb') and (kwargs.get('overwrite',False) == False):
    #     logging.error('Tdb already exists- run with overwrite to overwrite')
    #     sys.exit()

    return Bdb, prod_dir, cent_dir

def check_dependencies(print_out=False):
    '''
    For all possible dependencies, see if you can find them
    '''
    for dep in ['mash', 'nucmer', 'checkm', 'ANIcalculator', 'prodigal', 'centrifuge',
                'nsimscan', 'fastANI']:
        loc, works = find_program(dep)
        works_message = {True:'all good', False:'!!! ERROR !!!'}[works]
        message = '{0:.<40} {1:15} (location = {2})'.format(dep, works_message, loc)
        if print_out:
            print(message)
        else:
            logging.info(message)

def find_program(dep):
    '''
    return location of progrgam, works = True/False (based on calling the help)
    '''
    # find the location of the program
    loc = shutil.which(dep)

    # make sure the help on the program works
    works = False
    if loc != None:
        try:
            o = subprocess.check_output([loc, '-h'],stderr= subprocess.STDOUT)
            works = True
        except:
            pass

    return loc, works

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

    if len(cmds) >= 1:
        logging.info('Running Centrifuge')
        for cmd in cmds:
            logging.debug(' '.join(cmd))

        if 'wd' in kwargs:
            logdir = kwargs.get('wd').get_dir('cmd_logs')
        else:
            logdir = False
        drep.thread_cmds(cmds, shell=False, logdir=logdir, t=int(t))
        #drep.d_cluster.thread_mash_cmds_status(cmds,t=int(t))

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

    # Find the best hits
    g2t = {}
    for genome in Tdb['genome'].unique():
        d = Tdb[Tdb['genome'] == genome]
        taxID = d['tax_ID'][d['tax_confidence'] == d['tax_confidence'].max()].tolist()[0]
        g2t[genome] = taxID
    Tdb['best_hit'] = [True if g2t[g] == t else False for g,t in zip(Tdb['genome'], Tdb['tax_ID'])]

    # Try and add full taxonomy string
    try:
        Tdb['full_tax'] = [lineage_from_taxId(t) if b else False for t, b in zip(\
                Tdb['tax_ID'], Tdb['best_hit'])]
    except:
        logging.info("problem determing full tax string with ete3 - skipping")

    return Tdb

def parse_centrifuge_percent(Bdb, cent_dir, **kwargs):
    min_perc = int(kwargs.get('percent'))
    min_score = kwargs.get('min_score', 250)

    Tdb = pd.DataFrame()
    for genome in Bdb['genome'].unique():
        hits = parse_raw_centrifuge("{0}{1}_hits.tsv".format(cent_dir,genome), \
                    "{0}{1}_report.tsv".format(cent_dir,genome))
        tdb = tdb_from_hits(hits[hits['score'] > min_score], minPerc= int(min_perc))
        tdb['genome'] = genome

        Tdb = pd.concat([Tdb, tdb])
    # Find the best hit

    # THIS IS HARDER BECAUSE YOU HAVE TO PICK THE LOWEST LEVEL OVER the percent
    # # Find the best hits
    # g2t = {}
    # for genome in Tdb['genome'].unique():
    #     d = Tdb[Tdb['genome'] == genome]
    #     taxID = d['tax_ID'][d['tax_confidence'] == d['tax_confidence'].max()].tolist()[0]
    #     g2t[genome] = taxID
    # Tdb['best_hit'] = [True if g2t[g] == t else False for g,t in zip(Tdb['genome'], Tdb['tax_ID'])]
    #
    # # Try and add full taxonomy string
    # try:
    #     Tdb['full_tax'] = [lineage_from_taxId(t) if b else False for t, b in zip(\
    #             Tdb['tax_ID'], Tdb['best_hit'])]
    # except:
    #     logging.info("problem determing full tax string with ete3 - skipping")

    return Tdb

def tdb_from_hits(hits, minPerc= 50, testing=False):
    '''
    Determines the lowest taxonomic level with at least minPerc certainty

    For every hit:
        reconstruct the lineage (kingdom, phylum, class, ect.)
        add a count to every rank in the lineage

    For every rank:
        see if the number of hits matching one taxa at that rank is above the minPerc
        the denominator for this equation is the number of hits that have a phyla rank

    * Note: this is complicated because some lower ranks don't have higher ranks
        For example, species [Eubacterium] rectale (taxID 39491) has no genus
        Also, species [artifical construct] (taxID 32630) has no anything but species

    '''

    from ete3 import NCBITaxa
    ncbi = NCBITaxa()

    Levels = ['superkingdom','phylum','class','order','family','genus','species']

    # generate nested dictionary for levels
    countDic = {}
    for level in Levels:
        countDic[level] = {}

    # fill in nested dictionary
    for t in hits['taxID'].tolist():
        if t == 0:
            continue

        # This try / except thing is trying to catch sporatic errors of:
        # sqlite3.OperationalError: disk I/O error
        try:
            lin = ncbi.get_lineage(t)
            lin2name  = ncbi.get_taxid_translator(lin)
            name2rank = ncbi.get_rank(lin)
        except:
            time.sleep(1)
            lin = ncbi.get_lineage(t)
            lin2name  = ncbi.get_taxid_translator(lin)
            name2rank = ncbi.get_rank(lin)

        for i in lin:
            rank = name2rank[i]
            name = lin2name[i]
            if rank in countDic:
                countDic[rank][i] = countDic[rank].get(i,0) + 1

    # make the table
    total = sum(countDic['phylum'].values())
    table = {'tax_ID':[], 'tax_confidence':[], 'tax_level':[], 'taxonomy':[]}
    count = None

    for level in Levels:
        dic = countDic[level]
        for name in sorted(dic, key=dic.get, reverse= True):
            count = dic[name]
            break

        if count == None:
            table['tax_ID'].append(None)
            table['tax_confidence'].append(0)
            table['tax_level'].append(level)
            table['taxonomy'].append('unk')

        else:
            lin = ncbi.get_lineage(name)
            lin2name  = ncbi.get_taxid_translator(lin)
            name2rank = ncbi.get_rank(lin)
            rank2name = {v: k for k, v in name2rank.items()}
            tax = (lin2name[rank2name[level]])

            table['tax_ID'].append(name)
            table['tax_confidence'].append(((count/total) *100))
            table['tax_level'].append(level)
            table['taxonomy'].append(tax)

        count = None
    tdb = pd.DataFrame(table)

    # find and mark the best hit
    best = tdb['tax_ID'][tdb['tax_confidence'] >= minPerc].tolist()[-1]
    tdb['best_hit'] = [True if i == best else False for i in tdb['tax_ID']]

    # get the full taxonomy for the best hit
    tdb['full_tax'] = [lineage_from_taxId(t) if b else False for t, b in zip(\
            tdb['tax_ID'], tdb['best_hit'])]

    return tdb

def lineage_from_taxId(t):
    from ete3 import NCBITaxa
    ncbi = NCBITaxa()

    Levels = ['superkingdom','phylum','class','order','family','genus','species']
    name = []

    lin = ncbi.get_lineage(t)

    lin2name  = ncbi.get_taxid_translator(lin)
    name2rank = ncbi.get_rank(lin)
    rank2name = {v: k for k, v in name2rank.items()}

    for level in Levels:
        if level in rank2name:
            name.append(lin2name[rank2name[level]])
        else:
            name.append('unk')

    return '|'.join([str(int(t))] + name)

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
        d = Tdb[(Tdb['genome'] == genome) & (Tdb['best_hit'] == True)]
        if len(d) != 1:
            raise ValueError('Bug')

        tax = d['taxonomy'].tolist()[0]
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
