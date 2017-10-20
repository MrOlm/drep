#!/usr/bin/env python3

from subprocess import call
import os
from Bio import SeqIO
import shutil
import multiprocessing
import datetime


# def run_cmd(cmd,dry=False,shell=True,quiet= True):
#     devnull = open(os.devnull, 'w')
#     if shell:
#         if not dry:
#             if quiet:
#                 call(cmd,shell=True,stdout=devnull, stderr=devnull)
#             else:
#                 call(cmd,shell=True)
#         else:
#             print(cmd)
#     else:
#         if not dry:
#             if quiet:
#                 call(cmd,stdout=devnull, stderr=devnull)
#             else:
#                 call(cmd)
#         else: print(' '.join(cmd))
#     return

def run_cmd(cmd, dry=False, shell=True, logdir=False):
    # if dry, just print cmd
    if dry:
        if shell:
            print(cmd)
        else:
            print(' '.join(cmd))
        return

    # figure out what files you're going to store the output
    if logdir == False:
        devnull = open(os.devnull, 'w')
        sto = devnull
        ste = devnull
    else:
        uniq_filename = str(datetime.datetime.now().date()) + '_' + \
            str(datetime.datetime.now().time()).replace(':', '.')
        sto = open(os.path.join(logdir + uniq_filename + '.STDOUT'), 'w')
        ste = open(os.path.join(logdir + uniq_filename + '.STDERR'), 'w')

        # log the command
        now = open(os.path.join(logdir + uniq_filename + '.CMD'), 'w')
        if shell:
            now.write(str(cmd) + '\n')
        else:
            now.write(' '.join(cmd) + '\n')
        now.close()

    # run the command
    if shell:
        call(cmd,shell=True,stdout=sto, stderr=ste)
    else:
        call(cmd,stdout=sto, stderr=ste)
    return

def thread_cmd_wrapper(tup):
    run_cmd(*tup)

def thread_cmds(cmds, dry=False, shell=False, logdir=False, t=10):
    pool = multiprocessing.Pool(processes=t)
    tups = [(cmd, dry, shell, logdir) for cmd in cmds]
    pool.map(thread_cmd_wrapper, tups)
    pool.close()
    pool.join()
    return

def make_dir(outdirname,dry=False,overwrite=False):
    if dry:
        return
    if os.path.exists(outdirname):
        if overwrite:
            pass
            #shutil.rmtree(outdirname, ignore_errors=True)
            #os.makedirs(outdirname)
        else:
            assert False, "{0} already exsists! Will not overwrite with current settings".\
                            format(outdirname)
    else:
        os.makedirs(outdirname)

def clobber_dir(outdirname,dry=False,overwrite=False):
    if dry:
        return
    if os.path.exists(outdirname):
        if overwrite:
            shutil.rmtree(outdirname, ignore_errors=True)
            os.makedirs(outdirname)
        else:
            assert False, "{0} already exsists! Will not overwrite with current settings".\
                            format(outdirname)
    else:
        os.makedirs(outdirname)

def fasta_length(fasta):
    total = 0
    for seq_record in SeqIO.parse(fasta, "fasta"):
        total += len(seq_record)
    return total
