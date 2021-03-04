import os
import glob
import shutil
import pandas as pd
import importlib
import logging
import subprocess
import pytest

import tests.test_utils as test_utils

class Empty():
    pass

@pytest.fixture()
def self():
    # Set up
    self = Empty()
    self.genomes = test_utils.load_test_genomes()
    self.wd_loc = test_utils.load_test_wd_loc()
    self.s_wd_loc = test_utils.load_solutions_wd()
    self.zipped_genome = test_utils.load_zipped_genome()

    importlib.reload(logging)
    if os.path.isdir(self.wd_loc):
        shutil.rmtree(self.wd_loc)
    os.mkdir(self.wd_loc)

    yield self

    importlib.reload(logging)
    if os.path.isdir(self.wd_loc):
        shutil.rmtree(self.wd_loc)
    os.mkdir(self.wd_loc)


# class test_bonus():
#     def __init__(self):
#         pass
#
#     def setUp(self):
#         self.genomes = test_utils.load_test_genomes()
#         self.wd_loc = test_utils.load_test_wd_loc()
#         self.s_wd_loc = test_utils.load_solutions_wd()
#         self.zipped_genome = test_utils.load_zipped_genome()
#
#         importlib.reload(logging)
#         if os.path.isdir(self.wd_loc):
#             shutil.rmtree(self.wd_loc)
#         os.mkdir(self.wd_loc)
#
#     def run(self):
#         self.setUp()
#         self.test_drep_scaffold_level()
#         self.tearDown()
#
#         self.setUp()
#         self.test_drep_scaffold_level2()
#         self.tearDown()
#
#         self.setUp()
#         self.test_drep_scaffold_level3()
#         self.tearDown()
#
#         self.setUp()
#         self.test_drep_scaffold_level4()
#         self.tearDown()
#
#         self.setUp()
#         self.test_drep_scaffold_level5()
#         self.tearDown()
#
#         self.setUp()
#         self.test_parse_stb()
#         self.tearDown()
#
#         self.setUp()
#         self.test_parse_stb2()
#         self.tearDown()

def test_drep_scaffold_level(self):
    '''
    test ScaffoldLevel_dRep.py basicaly
    '''
    script_loc = os.path.join(str(os.getcwd()),'../helper_scripts/ScaffoldLevel_dRep.py')

    cmd = "{0} -f {1} {1} -o {2}".format(script_loc, self.genomes[0], self.wd_loc)
    print(cmd)
    subprocess.call(cmd, shell=True)

    outs = glob.glob(self.wd_loc + '/*')
    assert(len(outs) == 3)
    f = [o for o in outs if 'DereplicationInfo.csv' in o][0]
    Rdb = pd.read_csv(f)
    assert len(Rdb) == 5

def test_drep_scaffold_level_2(self):
    '''
    test ScaffoldLevel_dRep.py with --IgnoreSameScaffolds
    '''
    script_loc = os.path.join(str(os.getcwd()),'../helper_scripts/ScaffoldLevel_dRep.py')

    cmd = "{0} -f {1} {1} -o {2} --IgnoreSameScaffolds".format(script_loc, self.genomes[0], self.wd_loc)
    print(cmd)
    subprocess.call(cmd, shell=True)

    outs = glob.glob(self.wd_loc + '/*')
    assert(len(outs) == 3)
    f = [o for o in outs if 'DereplicationInfo.csv' in o][0]
    Rdb = pd.read_csv(f)
    assert len(Rdb) == 0

def test_drep_scaffold_level_3(self):
    '''
    test ScaffoldLevel_dRep.py with a single fasta file
    '''
    script_loc = os.path.join(str(os.getcwd()),'../helper_scripts/ScaffoldLevel_dRep.py')

    cmd = "{0} -f {1} -o {2}".format(script_loc, self.genomes[0], self.wd_loc)
    print(cmd)
    subprocess.call(cmd, shell=True)

    outs = glob.glob(self.wd_loc + '/*')
    assert(len(outs) == 4)
    f = [o for o in outs if 'DereplicationInfo.csv' in o][0]
    Rdb = pd.read_csv(f)
    assert len(Rdb) == 0

def test_drep_scaffold_level_4(self):
    '''
    test for graceful crash
    '''
    script_loc = os.path.join(str(os.getcwd()),'../helper_scripts/ScaffoldLevel_dRep.py')

    cmd = "{0} -f {1} -o {2}".format(script_loc, self.zipped_genome, self.wd_loc)
    print(cmd)
    subprocess.call(cmd, shell=True)

    outs = glob.glob(self.wd_loc + '/*')
    assert(len(outs) == 1), outs

def test_drep_scaffold_level_5(self):
    '''
    test breaking things into chunks
    '''
    script_loc = os.path.join(str(os.getcwd()),'../helper_scripts/ScaffoldLevel_dRep.py')

    cmd = "{0} -f {1} -o {2} -m 1000000".format(script_loc, [g for g in self.genomes if 'Enterococcus_faecalis_YI6' in g][0], self.wd_loc)
    print(cmd)
    subprocess.call(cmd, shell=True)

    outs = glob.glob(self.wd_loc + '/*')
    assert(len(outs) == 4), outs

def test_parse_stb(self):
    '''
    test parse_stb.py basicaly
    '''
    script_loc = os.path.join(str(os.getcwd()),'../helper_scripts/parse_stb.py')
    out_loc = self.wd_loc + '/test.stb'

    cmd = "{0} --reverse -f {1} -o {2}".format(script_loc, ' '.join(self.genomes), out_loc)
    print(cmd)
    subprocess.call(cmd, shell=True)

    db = pd.read_csv(out_loc, sep='\t', names=['scaffold', 'bin'])
    assert len(db) == 124
    assert len(db['bin'].unique()) == 5

    # Do it the regular way
    new_out = out_loc + '_genomes_'
    cmd = "{0} -f {1} -s {2} -o {3}".format(script_loc, self.genomes[0], out_loc, new_out)
    print(cmd)
    subprocess.call(cmd, shell=True)

    out_genomes = glob.glob(new_out + '*')
    print(out_genomes)
    assert len(out_genomes) == 1

def test_parse_stb_2(self):
    '''
    test parse_stb.py on zipped genomes
    '''
    script_loc = os.path.join(str(os.getcwd()),'../helper_scripts/parse_stb.py')
    out_loc = self.wd_loc + '/test.stb'

    cmd = "{0} --reverse -f {1} -o {2}".format(script_loc, self.zipped_genome, out_loc)
    print(cmd)
    subprocess.call(cmd, shell=True)

    db = pd.read_csv(out_loc, sep='\t', names=['scaffold', 'bin'])
    assert len(db) == 1
    print(db)

    # Do it the regular way
    new_out = out_loc + '_genomes_'
    cmd = "{0} -f {1} -s {2} -o {3}".format(script_loc, self.zipped_genome, out_loc, new_out)
    print(cmd)
    subprocess.call(cmd, shell=True)

    out_genomes = glob.glob(new_out + '*')
    assert len(out_genomes) == 1, out_genomes
