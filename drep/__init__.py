#!/usr/bin/env python3

from subprocess import call
import os
from Bio import SeqIO
import shutil
import multiprocessing

def run_cmd(cmd,dry=False,shell=True,quiet= True):
    devnull = open(os.devnull, 'w')
    if shell:
        if not dry: 
            if quiet:
                call(cmd,shell=True,stdout=devnull, stderr=devnull)
            else:
                call(cmd,shell=True)
        else: 
            print(cmd)
    else: 
        if not dry: 
            if quiet:
                call(cmd,stdout=devnull, stderr=devnull)
            else:
                call(cmd)
        else: print(' '.join(cmd))
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
    
def thread_cmds(cmds,t=10):
    pool = multiprocessing.Pool(processes=t)
    pool.map(run_cmd,cmds)
    pool.close()
    pool.join()
    return
    
def fasta_length(fasta):
    total = 0
    for seq_record in SeqIO.parse(fasta, "fasta"):
        total += len(seq_record)
    return total