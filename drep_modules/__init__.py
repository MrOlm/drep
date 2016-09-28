#!/usr/bin/env python3

from subprocess import call

def run_cmd(cmd,dry=False,shell=True):
    if shell:
        if not dry: call(cmd,shell=True)
        else: print(cmd)
    else: 
        if not dry: call(cmd)
        else: print(' '.join(cmd))
    return