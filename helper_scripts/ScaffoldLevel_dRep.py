#!/usr/bin/env python

# Version 0.1
# Matt Olm
# mattolm@berkeley.edu
# 02.15.19
# https://biotite.berkeley.edu/j/user/mattolm/notebooks/OtherProjects/OneOffs/AlexaPhage_1_writtingProgram.ipynb


import os
import sys
import gzip
import glob
import textwrap
import shutil
import distutils
import argparse
import pandas as pd

from shutil import copyfile
import subprocess
from subprocess import call
from collections import defaultdict

from tqdm import tqdm

from Bio import SeqIO

import drep

def version():
    versionFile = open(os.path.join(drep.__path__[0], 'VERSION'))
    return versionFile.read().strip()

__author__ = "Matt Olm"
__license__ = "MIT"
__version__ = version()


class abstractCommand():
    '''
    abstract class that makes other commands implement these methods
    '''
    def excecute(self):
        raise NotImplementedError

    def verify(self):
        raise NotImplementedError

class ANImCommand(abstractCommand):
    '''
    Hold all of the info necessary to run ANIm
    '''

    def __init__(self, querry=None, reference=None, out_dir=None, settings='BASIC'):
        '''
        Declare all settings
        '''
        # Load settings
        self.config = get_default_args(settings)

        # Store parameters
        self.config['querry'] = querry
        self.config['reference'] = reference
        self.config['out_dir'] = out_dir

        # Find program executable
        self.config['exe'] = find_executable('nucmer')

        # Set that the command has not been run
        self.run = False

    def verify(self):
        '''
        Autocomplete when possbile; return True if this command will likely run without problem
        '''

        # If any arguments are missing, throw exception
        for arg in self.config:
            if self.config[arg] is None:
                raise AttributeError("{0} must be set".format(arg))

        # Make sure required files exist
        for g in ['querry', 'reference']:
            if not os.path.isfile(self.config[g]):
                raise AttributeError("{0} doesn't exist ({1})".format(g,
                    self.config[g]))

        # Make sure output directory exists; try and make it if not
        for d in ['out_dir']:
            if not os.path.isdir(self.config[d]):
                os.makedirs(self.config[d])

        # Get the name of the file
        self.config['prefix'] = self.gen_prefix()

        return True

    def excecute(self, dry=False):
        ''' Run command and return a touple of (ani, coverage)
        '''
        assert self.verify()

        # 1) call the program
        cmd = gen_mummer_cmd(**self.config)
        #print(' '.join(cmd))
        if not dry:
            code = call(cmd)
            if code != 0:
                print("!!! nucmer failed with exit code {0}. Im going to crash now !!!".format(code))
                print(cmd)
                raise Exception("!!! nucmer failed with exit code {0}. Im going to crash now !!!".format(code))

        # 2) Mark that this command has been run
        self.run = True

        # 3) make sure command was successful
        assert os.path.isfile(get_delta_loc(self.gen_prefix()))
        return True

    def parse_results(self, method='simple', **kwargs):
        ''' parse results of the mummer run
        simple = (ani, coverage)
        '''
        # Make sure the command has been run
        if not self.run:
            raise AttributeError('Command has not been run')

        if method == 'simple':
            # Get the organism lengths
            qlen = fasta_length(self.config['querry'])
            rlen = fasta_length(self.config['reference'])

            # Parse the reults
            delta = get_delta_loc(self.config['prefix'])
            return get_ani_simple(delta, qlen, rlen, **kwargs)

        elif method == 'snps':
            return snps_from_delta(get_delta_loc(self.config['prefix']))

        elif method == 'snp_splits':
            # Get the filtered delta location
            f_delta = filter_delta(get_delta_loc(self.config['prefix']))

            # Get the snps location
            snp_loc = get_snps(f_delta)

            # Call the splits
            snp_splits = calc_snp_splits(f_delta, snp_loc)
            return snp_splits

    def gen_prefix(self):
        f = "{0}_vs_{1}".format(os.path.basename(self.config['reference']),
            os.path.basename(self.config['querry']))
        return os.path.join(self.config['out_dir'], f)

    def __str__(self):
        ''' Show the command parameters '''
        return pprint.pformat(self.config)

def gen_mummer_cmd(**kwargs):
    '''
    from a dictionary of arguments, return the ANIm command as an array of strings
    '''

    cmd = [kwargs['exe'],'--' + kwargs['method'],'-p',kwargs['prefix'], '-c', \
        kwargs['c'], '-g', kwargs['maxgap'], '-t', str(kwargs['p'])]

    if kwargs['noextend'] == 'True':
        cmd.append('--noextend')
    cmd += [kwargs['reference'],kwargs['querry']]

    return cmd

def gen_snp_cmd(**kwargs):
    '''
    from a dictionary of arguments, return the show-snp command
    '''
    snp_exe = kwargs.get('snp_exe')
    delta = kwargs.get('delta')
    snp_loc = kwargs.get('snp_loc', get_snp_loc(delta))

    cmd = [snp_exe, '-T', delta, '>', snp_loc]
    return ' '.join(cmd)

def gen_filter_cmd(**kwargs):
    '''
    from a dictionary of arguments, return the delta-filter command
    '''
    filter_exe = kwargs.get('filter_exe')
    delta = kwargs.get('delta')
    filtered_loc = kwargs.get('filtered_loc', get_filtered_loc(delta))

    cmd = [kwargs['filter_exe'], '-r', '-q', delta, '>', filtered_loc]

    return ' '.join(cmd)

def get_ani_simple(delta, qlen, rlen, **kwargs):
    ''' return (ani, coverage) from the delta file
    '''

    aln_length, sim_errors = parse_delta(delta)
    ani = 1 - (sim_errors / aln_length)
    cov = max((aln_length/qlen), aln_length/rlen)

    return (ani, cov)

def parse_delta(filename):
    '''
    parse the delta file and return the raw results
    Returns (alignment length, similarity errors) tuple from passed .delta.
    - filename - path to the input .delta file
    Extracts the aligned length and number of similarity errors for each
    aligned uniquely-matched region, and returns the cumulative total for
    each as a tuple.
    '''

    aln_length, sim_errors = 0, 0
    for line in [l.strip().split() for l in open(filename, 'rU').readlines()]:
        if line[0] == 'NUCMER' or line[0].startswith('>'):  # Skip headers
            continue
        # We only process lines with seven columns:
        if len(line) == 7:
            aln_length += abs(int(line[1]) - int(line[0]))
            sim_errors += int(line[4])
    return aln_length, sim_errors

def parse_delta_extra(filename):
    '''
    Return a pandas table of the details of a .delta file
    '''

    Table = {'r_start':[],'r_end':[], 'r_aln_length':[], \
             'q_start':[],'q_end':[], 'q_aln_length':[], \
             'total_errors':[], 'sim_errors':[],\
             'r_scaff':[],'q_scaff':[],'filename':[],\
             'reference':[],'querry':[]}

    for line in [l.strip().split() for l in open(filename, 'r').readlines()]:
        # Skip top part
        if len(line) == 2:
            reference = line[0]
            querry = line[1]

        # Read header
        elif line[0].startswith('>'):
            r_scaff = line[0][1:]
            q_scaff = line[1]

            r_len = line[2]
            q_len = line[3]

        # Get alignment of alignment lines
        elif len(line) == 7:
            r_start = int(line[0])
            r_end = int(line[1])

            q_start = int(line[2])
            q_end = int(line[3])

            total_errors = int(line[4])
            sim_errors = int(line[5])

            # should always be 0
            assert line[6] == '0'

            Table['r_start'].append(r_start)
            Table['q_start'].append(q_start)
            Table['r_end'].append(r_end)
            Table['q_end'].append(q_end)
            Table['r_aln_length'].append(abs(r_end - r_start))
            Table['q_aln_length'].append(abs(q_end - q_start))
            Table['total_errors'].append(total_errors)
            Table['sim_errors'].append(sim_errors)
            Table['r_scaff'].append(r_scaff)
            Table['q_scaff'].append(q_scaff)
            Table['filename'].append(filename)
            Table['reference'].append(reference)
            Table['querry'].append(querry)

    ddb = pd.DataFrame(Table)

    # Optimize dataframe size
    for cat in ['filename', 'reference', 'querry', 'r_scaff', 'q_scaff']:
        ddb[cat] = ddb[cat].astype('category')

    for num in ['r_start', 'q_start', 'r_end', 'q_end', 'r_aln_length',
            'q_aln_length', 'total_errors', 'sim_errors']:
        ddb[num] = ddb[num].astype(int)

    return ddb

def snps_from_delta(delta, snp_exe=False, filter_exe=False):
    '''
    make a .snp file from a delta file
    first filter, then get snps
    '''
    # filter the delta file
    filter_loc = filter_delta(delta)

    # call snps on filtered_delta
    snp_loc = get_snps(filter_loc)

    return parse_snp(snp_loc)

def filter_delta(delta, filter_exe=False):
    '''
    filter the .delta file
    '''
    if filter_exe is False:
        filter_exe = find_executable('delta-filter')

    filter_loc = get_filtered_loc(delta)
    cmd = gen_filter_cmd(delta=delta, filter_exe=filter_exe)
    #print(cmd)
    call(cmd, shell=True)

    return filter_loc

def get_snps(delta, snp_exe=False):
    '''
    call snps from delta file
    '''

    if snp_exe is False:
        snp_exe = find_executable('show-snps')

    snp_loc = get_snp_loc(delta)
    cmd = gen_snp_cmd(delta=delta, snp_exe=snp_exe)
    #print(cmd)
    call(cmd, shell=True)

    return snp_loc

def parse_snp(loc):
    '''
    return file as pandas dataframe
    '''
    h = ['rl', 'rb', 'qb', 'ql', 'buff', 'dist', 'R', 'Q', 'FRM', 'TAGS', 'rscaff', 'qscaff']
    x = pd.read_table(loc, skiprows=4, header=None)
    x.columns = h

    x['ql'] = x['ql'].astype(int)
    x['rl'] = x['rl'].astype(int)
    x['buff'] = x['buff'].astype(int)
    x['dist'] = x['dist'].astype(int)
    x['R'] = x['R'].astype(int)
    x['Q'] = x['Q'].astype(int)
    x['FRM'] = x['FRM'].astype(int)
    x['TAGS'] = x['TAGS'].astype(int)
    x['rscaff'] = x['rscaff'].astype('category')
    x['qscaff'] = x['qscaff'].astype('category')

    return x

def filter_snp(db):
    '''
    filter the input snp db such that indels only count as one line
    '''
    d = db.drop_duplicates(subset=['ql', 'qscaff'])
    d = d.drop_duplicates(subset=['rl', 'rscaff'])
    return d

def calc_snp_splits(delta_file, snp_file, **kwargs):
    # 0) get arguments
    split_len = kwargs.get('split_len', 1000)

    # 1) parse both files
    delta_db = parse_delta_extra(delta_file)
    snp_db = parse_snp(snp_file)

    # 2) calculate split identities
    split_dic = make_split_dic(delta_db, snp_db, split_len)

    return split_dic

def make_split_dic(delta_db, snp_db, split_len):
    '''
    Return a table of splits
    It is essential that each part of the reference only be aligned to one region
    of the querry. This is because I have a list of SNPs, but I don't know what
    alignment they're to. So if a SNP is in the reference sequence in only one
    of two alignments, I would have no way of knowing that, and will mark that
    reference location as a SNP in each alignment it goes to.
    '''
    Sdb = pd.DataFrame()

    snp_db = filter_snp(snp_db)
    for scaff in delta_db['r_scaff'].unique():
        table = {'scaffold':[], 'start':[], 'end':[], 'snps':[]}
        found_locs = [] # keep track of snps from this scaffold
        ddb = delta_db[delta_db['r_scaff'] == scaff].sort_values('r_start')
        locs = sorted(snp_db['rl'][(snp_db['rscaff'] == scaff)].tolist())

        # for every alignment in this scaffold...
        for start, end in zip(ddb['r_start'], ddb['r_end']):
            assert start < end
            cur_locs = [] # keep track of all snps in this alignment

            # set up the modulo dictionary
            split2snps = {}

            # for every snp in this scaffold...
            for x in locs:

                # if the snps is in the current alignment...
                if start <= x < end:
                    cur_locs.append(x)
                    found_locs.append(x)

                    # add it to the split
                    sp = get_split(start, x, split_len)
                    split2snps.setdefault(sp, []).append(x)

            # add the results from this alignment to the table
            for split, snps in split2snps.items():
                st, ed = get_range_from_split(start, split, split_len, end)
                table['scaffold'].append(scaff)
                table['start'].append(st)
                table['end'].append(ed)
                table['snps'].append(len(snps))

        # make sure every snp in the scaffold was added to an alignment and split
        sdb = pd.DataFrame(table)
        assert len(found_locs) == len(locs) == sdb['snps'].sum()

        # add to the table
        Sdb = pd.concat([Sdb, sdb], ignore_index=True)

    Sdb = _finalize_table(Sdb)
    return Sdb

def find_executable(prog):
    try:
        exe = distutils.spawn.find_executable(prog)
    except:
        exe = False
    return exe

def get_split(start, x, split_len):
    '''
    return the split to which x belongs
    '''
    return int((x-start)/split_len)

def get_range_from_split(start, sp, split_len, end):
    '''
    return the range of values for this split
    '''
    st = start + (sp * split_len)
    ed = min(st + split_len, end)
    return st, ed

def _finalize_table(sdb):
    sdb['length'] = sdb['end'] - sdb['start']
    sdb['ani'] = [(l-s)/l for s, l in zip(sdb['snps'], sdb['length'])]

    sdb['scaffold'] = sdb['scaffold'].astype('category')
    sdb['end'] = sdb['end'].astype(int)
    sdb['start'] = sdb['start'].astype(int)
    sdb['snps'] = sdb['snps'].astype(int)
    sdb['length'] = sdb['length'].astype(int)

    return sdb

def get_delta_loc(prefix):
    '''
    return the .delta filename based on the prefix in the original command
    '''
    return prefix + '.delta'

def get_snp_loc(delta):
    '''
    return the .snp filename based on the prefix in the original command
    '''
    return delta + '.snp'

def get_filtered_loc(delta):
    '''
    return the .snp filename based on the prefix in the original command
    '''
    return delta + '.filtered'

def main(args):
    '''
    Main entry point
    '''
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)
        # print("outdir has to exist. You gave me {0}, which does not".format(
        # args.outdir))
        # sys.exit()

    dereplicate_scaffolds(args.fasta_files, args.outdir, **vars(args))

def get_default_args(settings):
    '''
    Load default arguments from file and return dictionary of results
    '''
    if settings == 'BASIC':
        return {'c': '65', 'noextend': 'False', 'maxgap': '90', 'method': 'mum'}
    else:
        return {}


def dereplicate_scaffolds(fastas, output_folder, **kwargs):
    '''
    De-replicate several lists of scaffolds

    Arguments:
        fastas = List of .fasta files to be de-replicated on a scaffold level, in priority order
        output_folder = Folder to store output files
    '''
    verbose = kwargs.get('verbose', False)

    # Do some quick validating
    if output_folder[-1] != '/':
        output_folder = output_folder + '/'
    assert os.path.isdir(output_folder)

    if verbose:
        print("Step 1- Load scaffold to length for .fasta files")
    scaffDb = load_fastas(fastas)

    assert len(fastas) >= 1
    selfCompairson = False
    if len(fastas) == 1:
        kwargs['IgnoreSameScaffolds'] = True
        selfCompairson = True

        # Duplicate the scaffold
        scaffDb2 = scaffDb.copy()
        scaffDb2['priority'] = 1
        scaffDb = pd.concat([scaffDb, scaffDb2]).reset_index(drop=True)

    if verbose:
        print("Step 2- Run comparisons and dereplicate scaffolds")
    final_file, RMdb, Sdb = compare_scaffolds(scaffDb, output_folder, **kwargs)

    if verbose:
        print("Step 3- Generate report")
    gen_report(output_folder, final_file, RMdb, Sdb, scaffDb, verbose=verbose,
                selfCompairson=selfCompairson)

def load_fastas(fastas):
    '''
    Create a pandas dataframe of scaffold to length
    '''
    Sdb = pd.DataFrame()
    for i, fasta in enumerate(fastas):
        db = load_scaffold2length(fasta)
        db['loc'] = fasta
        db['priority'] = i
        db['sort_name'] = os.path.basename(fasta)
        Sdb = pd.concat([Sdb, db])
        #Sdb = Sdb.append(db)

    Sdb['loc'] = Sdb['loc'].astype('category')
    Sdb['sort_name'] = Sdb['sort_name'].astype('category')
    Sdb['priority'] = Sdb['priority'].astype(int)
    Sdb['length'] = Sdb['length'].astype(int)

    return Sdb

def compare_scaffolds(scaffDb, output_folder, **kwargs):
    '''
    De-replicate scaffolds

    This follows a couple of steps:
    '''
    verbose = kwargs.get('verbose')

    if verbose:
        print("Step 2.1- Prepare a directory")
    temp_folder = os.path.join(output_folder + 'comparisons/')
    if not os.path.exists(temp_folder):
        os.mkdir(temp_folder)

    if verbose:
        print("Step 2.2- Run comparisons")

    # Current is the .fasta that you differ to. Initially this is the first .fasta in the list
    indecies = sorted(list(scaffDb['priority'].unique()))[1:]
    current = scaffDb[scaffDb['priority'] == 0]['loc'].unique()[0]

    # Keep track of removed scaffolds
    RMdb = pd.DataFrame()

    for i in indecies:
        # Scaffolds that are in new and current will be removed from new
        new = scaffDb[scaffDb['priority'] == i]['loc'].unique()[0]

        if verbose:
            print('comparing {0} and {1}'.format(current, new))
        Ddb, Sdb = compare_wrapper(current, new, scaffDb, temp_folder, **kwargs)

        if verbose:
            print('De-replicating scaffolds')

        # Make a new file where scaffolds that were in new and current are removed from new
        newfile = os.path.join(temp_folder, "deReplication_v{0}.fasta".format(i))
        rmdb = reconsile(Sdb, current, new, newfile, temp_folder, **kwargs)

        # Keep track of what was removed
        #RMdb = RMdb.append(rmdb)
        RMdb = pd.concat([RMdb, rmdb])

        current = newfile

    return current, RMdb, Sdb

def gen_report(output_folder, final_file, RMdb, Sdb, scaffDb, verbose=False, selfCompairson=False):
    '''
    Generate report
    '''
    if verbose:
        print('generating report')

    # Save final file
    ffile = os.path.join(output_folder + 'dereplicated_scaffolds.fasta')
    copyfile(final_file, ffile)

    # Save scaffold report
    rfile = os.path.join(output_folder + 'DereplicationInfo.csv')
    RMdb.to_csv(rfile, index=False)

    # Print report
    print("""
    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    $$$                                              $$$
    $$$    All finished! Here's your final report    $$$
    $$$                                              $$$
    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    """)

    if selfCompairson:
        indecies = [0]
    else:
        indecies = sorted(list(scaffDb['priority'].unique()))

    for i in indecies:
        db = scaffDb[scaffDb['priority'] == i]
        loc = db['loc'].unique()[0]
        ddb = RMdb[RMdb['file_removed_from'] == os.path.basename(loc)]
        print("{0} - {2} of {1} scaffolds removed".format(os.path.basename(loc),
                len(db), len(ddb)))

    print('\nDereplicated scaffolds located at {0}'.format(ffile))
    print('Details on scaffolds removed located at {0}'.format(rfile))

    if selfCompairson:
        # Save the ANI table
        rfile = os.path.join(output_folder + 'Scaffold_ANI_values.csv')
        Sdb.to_csv(rfile, index=False)
        print('Details on scaffolds removed located at {0}'.format(rfile))


def load_scaffold2length(fasta):
    '''
    Return a database of scaffold to length
    '''
    table = defaultdict(list)

    if fasta[-3:] == '.gz':
        with gzip.open(fasta, "rt") as handle:
            for seq_record in SeqIO.parse(handle, "fasta"):
                table['scaffold'].append(seq_record.id)
                table['length'].append(len(seq_record))
    else:
        for seq_record in SeqIO.parse(str(fasta), "fasta"):
            table['scaffold'].append(seq_record.id)
            table['length'].append(len(seq_record))
    return pd.DataFrame(table)

def add_cov_ani(Ddb, s2l):
    '''
    From an unfiltered Delta, report the coverage and ANI of each alignment.

    Based on length of reference sequence
    '''
    # Add coverage and ANI at the alignment level
    Ddb['coverage'] = [x/s2l[y] for x,y in zip(Ddb['r_aln_length'], Ddb['r_scaff'])]
    Ddb['ani'] = [(l -e)/l for l,e in zip(Ddb['r_aln_length'], Ddb['total_errors'])]

    # Add coverage and ANI at the scaffold level
    table = defaultdict(list)
    for rscaff, db in Ddb.groupby('r_scaff'):
        for qscaff in db['q_scaff'].unique():
            d = db[db['q_scaff'] == qscaff]
            if len(d) == 0:
                continue

            #print(qscaff)
            table['r_scaff'].append(rscaff)
            table['q_scaff'].append(qscaff)
            table['r_aln_length'].append(d['r_aln_length'].sum())
            table['total_errors'].append(d['total_errors'].sum())
            table['alignments'].append(len(d))


    Sdb = pd.DataFrame(table)
    Sdb['coverage'] = [x/s2l[y] for x,y in zip(Sdb['r_aln_length'], Sdb['r_scaff'])]
    Sdb['ani'] = [(l -e)/l for l,e in zip(Sdb['r_aln_length'], Sdb['total_errors'])]

    return Ddb, Sdb

def split_into_chunks(file, max_len, temp_folder):
    '''
    Split the file "file" into numerous files of max_len size
    '''
    newfile_base = os.path.join(temp_folder, os.path.basename(file))
    chunks = []

    chunk_number = 0
    records = []
    chunk_len = 0
    for r in SeqIO.parse(file, "fasta"):
        records.append(r)
        chunk_len += len(r)

        if chunk_len > max_len:
            # Save this chunk
            new_loc = newfile_base + '_CHUNK_{0}.fasta'.format(chunk_number)
            SeqIO.write(records, new_loc, "fasta")
            chunks.append(new_loc)

            chunk_len = 0
            chunk_number += 1
            records = []

    # Save the final chunk
    if len(records) > 0:
        new_loc = newfile_base + '_CHUNK_{0}.fasta'.format(chunk_number)
        SeqIO.write(records, new_loc, "fasta")
        chunks.append(new_loc)

    return chunks

def compare_wrapper(current, new, scaffDb, temp_folder, **kwargs):
    '''
    A wrapper around compare that can split and paralellize if needed
    '''
    p = kwargs.get('processes')
    verbose = kwargs.get('verbose')
    max_len = kwargs.get('max_len')

    # Get length of current
    if current in set(scaffDb['loc']):
        clen = scaffDb[scaffDb['loc'] == current]['length'].sum()
    else:
        print("YOU HAVENT CODED THE ABILITY TO SPLIT NEW .FASTAS YET!!!")
        clen = max_len - 1

    # Get length of new
    if new in set(scaffDb['loc']):
        nlen = scaffDb[scaffDb['loc'] == new]['length'].sum()
    else:
        print("YOU HAVENT CODED THE ABILITY TO SPLIT NEW .FASTAS YET!!!")
        nlen = max_len - 1

    # See if either of them need to be split
    if clen > max_len:
        c_chunks = split_into_chunks(current, max_len, temp_folder)
    else:
        c_chunks = [current]

    if nlen > max_len:
        n_chunks = split_into_chunks(new, max_len, temp_folder)
    else:
        n_chunks = [new]

    # Pairwise comparisons
    ddbs = []
    sdbs = []
    pbar = tqdm(desc='Running comparisons: ', total=(len(n_chunks) * len(c_chunks)))
    for c, c_chunk in enumerate(c_chunks):
        for n, n_chunk in enumerate(n_chunks):
            ddb, sdb = compare(c_chunk, n_chunk, scaffDb, temp_folder, p=p, verbose=verbose)
            ddb['c_chunk'] = c
            ddb['n_chunk'] = n
            sdb['c_chunk'] = c
            sdb['n_chunk'] = n
            ddbs.append(ddb)
            sdbs.append(sdb)
            pbar.update(1)

    pbar.close()
    return pd.concat(ddbs).reset_index(drop=True), pd.concat(sdbs).reset_index(drop=True)

def compare(current, new, scaffDb, temp_folder, p=1, verbose=False):
    '''
    Run mummer comparisons between .fasta files
    '''

    #Make mummer command
    delta_name = "{0}-vs-{1}".format(os.path.basename(current), os.path.basename(new))

    destf = "{0}{1}.delta".format(temp_folder, delta_name)
    ANIm_cmd = ANImCommand()
    ANIm_cmd.config['querry'] = current
    ANIm_cmd.config['reference'] = new
    ANIm_cmd.config['out_dir'] = temp_folder
    ANIm_cmd.config['p'] = p
    ANIm_cmd.config['exe'] = shutil.which('nucmer')
    ANIm_cmd.verify()

    #Run original
    cmd = gen_mummer_cmd(**ANIm_cmd.config)
    delta = get_delta_loc(ANIm_cmd.gen_prefix())
    ANIm_cmd.excecute()

    #Run filter
    filter_exe = shutil.which('delta-filter')
    f_result = filter_delta(delta, filter_exe)

    #Load and parse result
    s2l = scaffDb.set_index('scaffold')['length'].to_dict()
    Ddb = parse_delta_extra(f_result)
    if len(Ddb) > 0:
        Ddb['q_scaff'] = Ddb['q_scaff'].astype('category')
        Ddb['r_scaff'] = Ddb['r_scaff'].astype('category')
        Ddb, Sdb = add_cov_ani(Ddb, s2l)

        #Store
        Sdb.to_csv("{0}{1}.summary.csv".format(temp_folder, delta_name), index=False)
        Ddb.to_csv("{0}{1}.allAlignments.csv".format(temp_folder, delta_name), index=False)

    else:
        Sdb = pd.DataFrame({'coverage':[], 'ani':[], 'q_scaff':[], 'r_scaff':[]})

    #Return
    return Ddb, Sdb

def reconsile(Sdb, current, new, newFile, temp_folder, **kwargs):
    '''
    Make a new .fasta file combining current and new, based on Ddb and Sdb
    '''
    IgnoreSameScaffolds = kwargs.get('IgnoreSameScaffolds')
    minCov = kwargs.get('minCov')
    minANI = kwargs.get('minANI')
    verbose = kwargs.get('verbose')

    Rdb = Sdb.copy()
    Rdb = Rdb[(Rdb['coverage'] >= minCov) & (Rdb['ani'] >= minANI)]
    if IgnoreSameScaffolds:
        Rdb = Rdb[Rdb['q_scaff'] != Rdb['r_scaff']]

    # All querry should be keepers
    keepers = set(Rdb['q_scaff'])

    # All ref should be removered
    to_rm = set(Rdb['r_scaff'])

    # Make a new ref
    new_ref = os.path.join(temp_folder, os.path.basename(new) + '.dr')
    write_new(new, new_ref, to_rm, verbose=verbose)

    # Merge with current
    new_current = newFile
    with open(new_current, 'w') as outfile:
        for fname in [current, new_ref]:
            if fname[-3:] == '.gz':
                with gzip.open(fname,'r') as infile:
                    for line in infile:
                        outfile.write(line.decode('utf-8'))
            else:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

    # Make an ourput to return
    Rdb = Rdb[['r_scaff', 'q_scaff', 'coverage', 'ani']]
    Rdb = Rdb.rename(columns={'q_scaff':'scaffold_retained', 'r_scaff':'scaffold_removed'})
    Rdb['file_removed_from'] = os.path.basename(new)

    return Rdb


def write_new(ori_loc, new_loc, to_rm, verbose=False):
    '''
    Write a new .fasta file removing certain scaffolds
    '''
    records = []
    if ori_loc[-3:] == '.gz':
        with gzip.open(ori_loc, "rt") as handle:
            for r in SeqIO.parse(handle, "fasta"):
                records.append(r)
    else:
        for r in SeqIO.parse(ori_loc, "fasta"):
            records.append(r)

    ori_len = len(records)
    records = (r for r in records if r.id not in to_rm)
    count = SeqIO.write(records, new_loc, "fasta")
    if verbose:
        print("Removed {0} scaffolds from {1}".format(ori_len - count, ori_loc))


class test_dRepScaffolds():
    '''
    Tests
    '''
    def setUp(self):
        self.fastas = glob.glob('/home/mattolm/Bio_scripts/test_suite/test_files/*_min1000.fa') \
        + ['/home/mattolm/Bio_scripts/test_suite/test_files/N5_271_010G1_scaffold_3.fasta',
          '/home/mattolm/Bio_scripts/test_suite/test_files/NZ_BAGX02000038.1.fasta']
        self.test_dir = '/home/mattolm/Bio_scripts/test_suite/test_backend/testdir/'

        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)
        os.mkdir(self.test_dir)

    def tearDown(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run(self):
        self.setUp()
        self.test1()
        self.tearDown()

        self.setUp()
        self.test2()
        self.tearDown()

        self.setUp()
        self.test3()
        self.tearDown()

        self.setUp()
        self.test4()
        self.tearDown()

        print('all clear')

    def test1(self):
        '''
        Test load .fasta files
        '''
        scaffDb = load_fastas(self.fastas)
        assert len(scaffDb['priority'].unique() == 2)

    def test2(self):
        '''
        Test Compare
        '''
        scaffDb = load_fastas(self.fastas)
        Ddb, Sdb = dereplicate_scaffolds(self.fastas, self.test_dir, verbose=True)

        print(Ddb.head())
        print(Sdb.head())

    def test3(self):
        '''
        Test write_new
        '''
        to_rm = ['N5_271_010G1_scaffold_119']
        new_loc = self.test_dir + 'new_file.fasta'

        ori_len = len([r for r in SeqIO.parse(self.fastas[0], "fasta")])
        write_new(self.fastas[0], new_loc, to_rm, verbose=False)
        new_len = len([r for r in SeqIO.parse(new_loc, "fasta")])

        assert new_len + 1 == ori_len

    def test4(self):
        '''
        Test all
        '''
        dereplicate_scaffolds(self.fastas, self.test_dir, verbose=True)


#test_dRepScaffolds().run()

if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=textwrap.dedent('''\
         \n
         dRep_scaffoldLevel combines many lists of scaffolds into a single,
         dereplicated list. \n\n
         -f Takes the list of .fasta files. The first .fasta file in the list
         will take priority, and all scaffolds in the second file will be removed.
         The third file will then be compared to a combination of the first and second
         files, ect.

         If only one .fasta file is given it will be compared with itself, and
         --IgnoreSameScaffolds will be automatically enabled

         -o Is the directory where all output will be created. Check in the "comparisons" folder
         for raw information

         -minCov Is the minimum alignment coverage between two scaffolds for them
         to be considered the same. minCov and minANI must both be satisfied to remove
         a scaffold

         -minANI Is in the minimum average nucleotide identity (ANI) for two scaffolds
         to be the same
         '''))

    parser.add_argument("-f", "--fasta_files", nargs='*', help='fasta files to dereplicate on a per-scaffold level')
    parser.add_argument("-o", "--outdir", help='output directory')
    parser.add_argument("-c", "--minCov", help='minimum alignment coverage', default=0.5, type=float)
    parser.add_argument("-a", "--minANI", help='minimum alignment coverage', default=0.95, type=float)
    parser.add_argument("-m", "--max_len", help='maximum length (in bp) of each .fasta file. Split .fasta files above this length. Making this number lower will increase run-time and decrease RAM usage', default=1000000000, type=int)
    parser.add_argument("--IgnoreSameScaffolds", help='Dont dereplicate scaffolds with same name', action='store_true', default=False)
    parser.add_argument("-p", "--processes", help='number of processes to use', default=6, type=int)
    parser.add_argument('-d', "--debug", help='print extra output', action='store_true', default=False)

    # Specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main(args)
