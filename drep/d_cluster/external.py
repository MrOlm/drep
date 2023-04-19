import glob
import logging
import os
import shutil
import sys

import numpy as np
import pandas as pd

import drep
import drep.d_cluster.utils
# from drep.d_cluster.utils import _get_genome_name_from_fasta, process_deltafiles, _randomString, gen_gANI_cmd, \
#     process_gani_files, gen_goANI_cmd, process_goani_files

def run_pairwise_ANImf(genome_list, ANIn_folder, **kwargs):
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

    # Make folder if doesnt exist
    if not os.path.exists(ANIn_folder):
        os.makedirs(ANIn_folder)

    # Gen commands
    cmds = []
    files = []
    for g1 in genomes:
        # Make it so each reference is it's own folder, to spread out .delta files
        cur_folder = os.path.join(ANIn_folder, drep.d_cluster.utils._get_genome_name_from_fasta(g1))
        if not os.path.exists(cur_folder):
            os.makedirs(cur_folder)

        for g2 in genomes:
            file_name = "{0}/{1}_vs_{2}".format(cur_folder,
                        drep.d_cluster.utils._get_genome_name_from_fasta(g1),
                        drep.d_cluster.utils._get_genome_name_from_fasta(g2))
            files.append(file_name)

            # If the file doesn't already exist, add it to what needs to be run
            if not os.path.isfile(file_name + '.delta.filtered'):
                cmds.append(gen_animf_cmd(file_name,g1,g2))

    # Run commands
    if len(cmds) > 0:
        for c in cmds:
            logging.debug(c)

        if ('wd' in kwargs) and (kwargs.get('debug', False)):
            logdir = kwargs.get('wd').get_dir('cmd_logs')
        else:
            logdir = False
        drep.thread_cmds(cmds, shell=True, logdir=logdir, t=int(p))

    # Make dictionary of genome lengths
    org_lengths = {}
    for genome in genomes:
        org_lengths[drep.d_cluster.utils._get_genome_name_from_fasta(genome)] = drep.d_filter.calc_fasta_length(genome)

    # Parse output
    deltafiles = ["{0}.delta.filtered".format(f) for f in files]
    df = drep.d_cluster.utils.process_deltafiles(deltafiles, org_lengths, **kwargs)

    return df


def run_pairwise_fastANI(genome_list, outdir, **kwargs):
    p = kwargs.get('processors',6)
    code = drep.d_cluster.utils._randomString(stringLength=10)

    # Make folders
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    tmp_dir = os.path.join(outdir, 'tmp/')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # Make genome list
    glist = os.path.join(tmp_dir, 'genomeList')
    glist = _make_glist(genome_list, glist)

    # Gen command
    exe_loc = drep.get_exe('fastANI')
    out_base = os.path.join(outdir, 'fastANI_out_{0}'.format(code))
    cmd = [exe_loc, '--ql', glist, '--rl', glist, '-o', out_base, '--matrix', '-t', str(p), "--minFraction", str(0)]
    logging.debug(' '.join(cmd) + ' ' + code)

    # Run command
    if ('wd' in kwargs) and (kwargs.get('debug', False)):
        logdir = kwargs.get('wd').get_dir('cmd_logs')
    else:
        logdir = False
    drep.thread_cmds([cmd], shell=False, logdir=logdir, t=1)

    # Load results
    fdb = load_fastani(out_base)

    # fix missing ones
    try:
        fdb = _fix_fastani(fdb)
        return fdb

    # handle broken self
    except:
        logging.error("CRITICAL ERROR WITH SECONDARY CLUSTERING CODE {0}; SKIPPING".format(code))
        return pd.DataFrame()

def fastani_one_vs_many(one, many, genome_rep_file, outdir, **kwargs):
    p = kwargs.get('processors', 6)
    code = drep.d_cluster.utils._randomString(stringLength=10)
    tmp_dir = kwargs.get('tmp_dir')
    logdir = kwargs.get('logdir')
    exe_loc = kwargs.get('current_exe')
    redo = kwargs.get('redo', False)

    # Gen command
    out_base = os.path.join(outdir, 'fastANI_out_{0}'.format(code))
    cmd = [exe_loc, '-q', one, '--rl', genome_rep_file, '-o', out_base, '--matrix', '-t', str(p), "--minFraction", str(0)]
    logging.debug(' '.join(cmd) + ' ' + code)

    # Run command
    drep.thread_cmds([cmd], shell=False, logdir=logdir, t=1)

    # Load results
    fdb = load_fastani(out_base)

    return fdb

def load_fastani(file):
    fdb = pd.read_csv(file, names=['reference', 'querry', 'ANI', 'j1', 'j2'], delim_whitespace=True)
    for c in ['reference', 'querry']:
        fdb[c] = [drep.d_cluster.utils._get_genome_name_from_fasta(x) for x in fdb[c]]
    fdb = fdb.rename(columns={'ANI':'ani'})
    fdb['alignment_coverage'] = [(j1/j2) for j1, j2 in zip(fdb['j1'], fdb['j2'])]
    fdb = fdb[['reference', 'querry', 'ani', 'alignment_coverage']]
    fdb['ani'] = [x/100 for x in fdb['ani']]

    return fdb


def _fix_fastani(odb):
    # Add back missing genomes
    fdb = odb.pivot(index="reference", columns="querry", values="ani")
    fdb.reset_index(level=0, inplace=True)
    fdb.fillna(0, inplace=True)
    fdb = fdb.melt(id_vars=['reference']).rename(
            columns={'value':'ani'})

    # Add back alignment coverage
    fdb = pd.merge(fdb, odb[['reference', 'querry', 'alignment_coverage']], on=['reference', 'querry'], how='outer')
    fdb['alignment_coverage'] = fdb['alignment_coverage'].fillna(0)

    assert len(fdb['reference'].unique()) == len(fdb['querry'].unique())
    assert len(fdb) == (len(fdb['reference'].unique()) * len(fdb['querry'].unique()))

    return fdb


def _make_glist(genomes, floc):
    o = open(floc, 'w')
    for g in genomes:
        loc = os.path.abspath(g)
        assert os.path.isfile(g)
        o.write(loc + '\n')
    o.close()
    return floc


def run_pairwise_gANI(bdb, gANI_folder, prod_folder, **kwargs):
    '''
    Run pairwise gANI on a list of Genomes

    Args:
        bdb: DataFrame with ['genome', 'location']
        gANI_folder: folder to store gANI output
        prod_folder: folder containing prodigal output from genomes (will run if needed)

    Keyword arguments:
        debug: log all of the commands
        wd: if you want to log commands, you also need the wd
        processors: threads to use

    Returns:
        DataFrame: Ndb for gANI
    '''
    p = kwargs.get('processors',6)
    gANI_exe = drep.get_exe('ANIcalculator')
    genomes = bdb['location'].tolist()

    # Make folders
    if not os.path.exists(gANI_folder):
        os.makedirs(gANI_folder)
    if not os.path.exists(prod_folder):
        os.makedirs(prod_folder)

    # Remove crap folders- they shouldn't exist and if they do it messes things up
    crap_folders = glob.glob(gANI_folder + '*.gANITEMP')
    for crap_folder in crap_folders:
        if os.path.exists(crap_folder):
            shutil.rmtree(crap_folder)

    # Run prodigal
    logging.debug("Running prodigal...")
    drep.d_filter.run_prodigal(bdb['location'].tolist(), prod_folder, **kwargs)

    # Gen gANI commands
    logging.debug("Running gANI...")
    cmds = []
    files = []
    for i, g1 in enumerate(genomes):
        # Make it so each reference is it's own folder, to spread out .delta files
        cur_folder = os.path.join(gANI_folder, drep.d_cluster.utils._get_genome_name_from_fasta(g1))
        if not os.path.exists(cur_folder):
            os.makedirs(cur_folder)

        for j, g2 in enumerate(genomes):
            if i > j:
                name1= drep.d_cluster.utils._get_genome_name_from_fasta(g1)
                name2= drep.d_cluster.utils._get_genome_name_from_fasta(g2)
                file_name = "{0}/{1}_vs_{2}.gANI".format(cur_folder,
                            name1, name2)
                files.append(file_name)

                # If the file doesn't already exist, add it to what needs to be run
                if not os.path.isfile(file_name):
                    fna1 = "{0}.fna".format(os.path.join(prod_folder,name1))
                    fna2 = "{0}.fna".format(os.path.join(prod_folder,name2))
                    cmds.append(drep.d_cluster.utils.gen_gANI_cmd(file_name,fna1,fna2,gANI_exe))

    # Run commands
    if len(cmds) > 0:
        logging.debug('Running gANI commands: {0}'.format('\n'.join([' '.join(x) for x in cmds])))
        if ('wd' in kwargs) and (kwargs.get('debug', False) == True):
            logdir = kwargs.get('wd').get_dir('cmd_logs')
        else:
            logdir = False
            #logdir = "/home/mattolm/Programs/drep/tests/test_backend/logs/"
        drep.thread_cmds(cmds, logdir=logdir, t=int(p))

    else:
        logging.debug("gANI already run- will not re-run")

    # Pdrep.d_cluster.utils.arse output
    df = drep.d_cluster.utils.process_gani_files(files)

    # Add self-comparisons if there is only one genome
    if len(genomes) == 1:
        Table = {'querry':[],'reference':[],'ani':[],'alignment_coverage':[]}
        for g in genomes:
            Table['reference'].append(drep.d_cluster.utils._get_genome_name_from_fasta(g))
            Table['querry'].append(drep.d_cluster.utils._get_genome_name_from_fasta(g))
            Table['ani'].append(1)
            Table['alignment_coverage'].append(1)
        d = pd.DataFrame(Table)
        df = pd.concat([df,d],ignore_index=True)

    return df


def run_pairwise_goANI(bdb, goANI_folder, prod_folder, **kwargs):
    '''
    Run pairwise goANI on a list of Genomes

    Args:
        bdb: DataFrame with ['genome', 'location']
        goANI_folder: folder to store gANI output
        prod_folder: folder containing prodigal output from genomes (will run if needed)

    Keyword arguments:
        debug: log all of the commands
        wd: if you want to log commands, you also need the wd
        processors: threads to use

    Returns:
        DataFrame: Ndb for gANI
    '''
    p = kwargs.get('processors',6)
    nsimscan_exe = drep.get_exe('nsimscan')
    genomes = bdb['location'].tolist()

    # Make folders
    if not os.path.exists(goANI_folder):
        os.makedirs(goANI_folder)
    if not os.path.exists(prod_folder):
        os.makedirs(prod_folder)

    # Run prodigal
    logging.debug("Running prodigal...")
    drep.d_filter.run_prodigal(bdb['location'].tolist(), prod_folder, **kwargs)

    # Gen gANI commands
    logging.debug("Running goANI...")
    cmds = []
    files = []
    for i, g1 in enumerate(genomes):
        # Make it so each reference is it's own folder, to spread out .delta files
        cur_folder = os.path.join(goANI_folder, drep.d_cluster.utils._get_genome_name_from_fasta(g1))
        if not os.path.exists(cur_folder):
            os.makedirs(cur_folder)

        for j, g2 in enumerate(genomes):
            if i != j:
                name1= drep.d_cluster.utils._get_genome_name_from_fasta(g1)
                name2= drep.d_cluster.utils._get_genome_name_from_fasta(g2)
                file_name = "{0}/{1}_vs_{2}.sim".format(cur_folder,
                            name1, name2)
                files.append(file_name)

                # If the file doesn't already exist, add it to what needs to be run
                if not os.path.isfile(file_name):
                    fna1 = "{0}.fna".format(os.path.join(prod_folder,name1))
                    fna2 = "{0}.fna".format(os.path.join(prod_folder,name2))
                    cmds.append(drep.d_cluster.utils.gen_goANI_cmd(file_name,fna1,fna2,nsimscan_exe))

    # Run commands
    if len(cmds) > 0:
        logging.debug('Running goANI commands: {0}'.format('\n'.join([' '.join(x) for x in cmds])))
        if ('wd' in kwargs) and (kwargs.get('debug', False) == True):
            logdir = kwargs.get('wd').get_dir('cmd_logs')
        else:
            logdir = False
            #logdir = "/home/mattolm/Programs/drep/tests/test_backend/logs/"
        drep.thread_cmds(cmds, logdir=logdir, t=int(p))

    else:
        logging.debug("goANI already run- will not re-run")

    # Parse output
    df = drep.d_cluster.utils.process_goani_files(files)

    # Add self-comparisons if there is only one genome
    if len(genomes) == 1:
        Table = {'querry':[],'reference':[],'ani':[],'alignment_coverage':[]}
        for g in genomes:
            Table['reference'].append(drep.d_cluster.utils._get_genome_name_from_fasta(g))
            Table['querry'].append(drep.d_cluster.utils._get_genome_name_from_fasta(g))
            Table['ani'].append(1)
            Table['alignment_coverage'].append(1)
        d = pd.DataFrame(Table)
        df = pd.concat([df,d],ignore_index=True)

    return df


def parse_gani_file(file):
    '''
    Parse gANI file, return dictionary of results

    Args:
        file: location of gANI file

    Returns:
        dict: results in the gANI file
    '''
    try:
        x = pd.read_table(file)
    except:
        logging.error('gANI file {0} does not exist. The most likely reason is that '.format(file)\
            + 'one of the genomes has a .fasta header than gANI doesnt like.'\
            + ' Known issues include having a header over 160 characaters long, or any '\
            + 'kind of special character besides "_" (including ".", ":", and "-").'\
            + ' To fix this error either use ANIn (which is not as picky with fasta headers)'\
            + ', or fix the .fasta headers to conform to those rules.')
        sys.exit()
    x = x.rename(columns={'GENOME1':'reference','GENOME2':'querry','AF(1->2)':'rq_coverage',
                        'AF(2->1)':'qr_coverage','ANI(1->2)':'rq_ani','ANI(2->1)':'qr_ani'})
    dict = x.to_dict(orient='list')
    dict = {k:v[0] for k,v in dict.items()}
    dict['reference'] = dict['reference'][:-4]
    dict['querry'] = dict['querry'][:-4]
    return dict


def parse_nsim_file(file):
    '''
    Parse nsim file, return dictionary of results and gene datatable

    Args:
        file: location of nsimscan file

    Returns:
        dict: results in the gANI file
        db: gene-based alignment results
    '''
    # Load
    db = pd.read_csv(file, sep='\t')

    if len(db) == 0:
        logging.warning("File {0} is empty, indicating a nsimscan failure! Run with --debug and check the log folder for details".format(file))

        return _summarize_nsimsan(db)

    db = db.rename(columns={'#qry_id':'qry_id', '#Q_id':'qry_id', 'S_id':'sbj_id', 'trg_len':'sbj_len'})

    # Filter
    try:
        db = _filter_nsimscan(db)
    except:
        print("ERROR! FILE {0}".format(file))
        print(db)
        return _summarize_nsimsan(pd.DataFrame())

    # Make summary results
    x = _summarize_nsimsan(db)

    return x


def _filter_nsimscan(db1, minAF=0.7, minANI=70):
    '''
    Filter a single nsim scan file
    '''
    db1 = db1.sort_values(['al_len','qry_id', 'sbj_id'], ascending=False)

    # Filter like gANI
    db1['af'] = [a/min(o,t) for a,o,t in zip(db1['al_len'],
                                             db1['qry_len'],
                                             db1['sbj_len'])]
    db1 = db1[(db1['af'] >= minAF) & (db1['p_inden'] >= minANI)]

    return db1


def _summarize_nsimsan(db):
    '''
    Take all those aligned genes and return a summary ANI
    '''
    table = {}
    if len(db) > 0:
        table['ani'] = sum([p * l for p, l in zip(db['p_inden'], db['al_len'])]) \
                        / db['al_len'].sum()
        table['af'] = db['al_len'].sum() / db['qry_len'].sum()
    else:
        table['ani'] = 0
        table['af'] = 0

    return table


def gen_nucmer_cmd(prefix,ref,querry,c='65',noextend=False,maxgap='90',method='mum'):
    '''
    Generate command to run with nucmer

    Args:
        prefix: desired output name (will have .delta appended)
        ref: location of reference genome
        querry: location of querry genomes

    Keyword args:
        c: c value
        noextend: either True or False
        maxgap: maxgap
        method: detault is 'mum'

    Returns:
        list: command to run number
    '''
    cmd = ['nucmer','--' + method,'-p',prefix,'-c',str(c),'-g',str(maxgap)]
    if noextend: cmd.append('--noextend')
    cmd += [ref,querry]
    return cmd


def gen_animf_cmd(prefix,ref,querry, **kwargs):
    '''
    return animf command. It will be a single string, with the format:
    "nucmer cmd ; filter cmd"

    Args:
        prefix: desired file output name
        ref: location of reference genome
        querry: location of querry genome

    Keyword args:
        all: passed on to gen_nucmer_cmd

    Returns:
        string: format is "nucmer cmd ; filter cmd"
    '''

    nucmer_cmd = gen_nucmer_cmd(prefix,ref,querry, **kwargs)

    delta = prefix + '.delta'
    out = delta + '.filtered'

    filter_cmd = gen_filter_cmd(prefix + '.delta', out)
    return ' '.join(nucmer_cmd) + '; ' + ' '.join(filter_cmd)


def gen_filter_cmd(delta, out):
    '''
    return delta-filter command

    Args:
        delta: desired .delta file to filter
        out: desired output file

    Returns:
        list: cmd to run
    '''
    cmd = ["delta-filter", '-r', '-q', delta, '>', out]
    return cmd


def _gen_nomash_cdb(Bdb):
    '''
    From Bdb, just add a column of 'primary_cluster' = 0

    Args:
        Bdb: dataframe with [genome, location]

    Returns:
        DataFrame: Cdb
    '''
    Cdb = Bdb.copy()
    Cdb['primary_cluster'] = 0
    return Cdb


def add_avani(db):
    '''
    add a column titled 'av_ani' to the passed in dataframe

    dataframe must have rows reference, querey, and ani

    Args:
        db: dataframe
    '''

    logging.debug('making dictionary for average_ani')
    combo2value = {}
    for i, row in db.iterrows():
        combo2value["{0}-vs-{1}".format(row['querry'], row['reference'])] \
            = row['ani']

    logging.debug('list comprehension for average_ani')
    db['av_ani'] = [np.mean([combo2value["{0}-vs-{1}".format(q, r)],
                        combo2value["{0}-vs-{1}".format(r, q)]]) if r != q else 1\
                        for q, r in zip(db['querry'].tolist(),
                        db['reference'].tolist())]

    logging.debug('averageing done')


def _nucmer_preset(preset):
    '''
    Return the values from a nucmer 'preset'

    Args:
        preset: either "tight" or "normal"

    Returns:
        list: [c, n_maxgap, n_noextend, n_method]
    '''
    assert preset in ['tight','normal']

    if preset == 'tight':
        return 65, 1, True, 'mum'

    elif preset == 'normal':
        return 65, 90, False, 'mum'