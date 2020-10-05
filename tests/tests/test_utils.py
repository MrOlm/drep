import os
import glob
import pandas as pd

def load_solutions_wd():
    loc = os.path.join(str(os.getcwd()),'../tests/test_solutions/ecoli_wd')
    return loc

def load_test_wd_loc():
    loc = os.path.join(str(os.getcwd()),'../tests/test_backend/ecoli_wd')
    return loc

def load_test_wd_loc_2():
    loc = os.path.join(str(os.getcwd()),'../tests/test_backend/ecoli_wd2')
    return loc

def load_test_genomes():
    return glob.glob(os.path.join(str(os.getcwd()) + '/genomes/*'))

def load_zipped_genomes():
    return glob.glob(os.path.join(str(os.getcwd()),'../tests/test_backend/zipped/*.gz'))

def load_broken_genome():
    return glob.glob(os.path.join(str(os.getcwd()),'../tests/test_backend/other/broken_genome.fasta'))[0]

def load_zipped_genome():
    return glob.glob(os.path.join(str(os.getcwd()),'../tests/test_backend/zipped/Enterococcus_casseliflavus_EC20.fasta.gz'))[0]

def load_random_test_dir():
    loc = os.path.join(str(os.getcwd()),'../tests/test_backend/test_dir')
    return loc

def load_solutions_taxonomy_wd():
    loc = os.path.join(str(os.getcwd()),'../tests/test_solutions/ecoli_taxonomy')
    return loc

def load_large_genome_set():
    loc = '/Users/mattolm/Programs/testing_house/test_genomes/'
    genomes = glob.glob(loc + '*.fna')
    return genomes

def compare_dfs(db1, db2, round=3, verbose=False):
    '''
    Return True if dataframes are equal (order of dataframes doesn't matter)
    '''

    db1 = db1.fillna(0).round(round)
    db2 = db2.fillna(0).round(round)

    df = pd.concat([db1, db2], sort=True)
    df = df.reset_index(drop=True)
    df_gpby = df.groupby(list(df.columns))
    idx = [x[0] for x in df_gpby.groups.values() if len(x) == 1]

    identicle = (len(idx) == 0)
    if ((not identicle) and verbose):
        print("index: ", idx)
        print("db1: ",db1)
        print("db2: ",db2)
        print("df_gpby: ", str(df_gpby))

    return identicle

def compare_dfs2(db1, db2, verbose=False):
    """
    Return True if dataframes are equal (order of dataframes doesn't matter)
    """
    try:
        # noinspection PyProtectedMember
        pd._testing.assert_frame_equal(db1, db2)
        return True
    except AssertionError as e:
        if verbose:
            print(e)
        return False

def report_diff(x):
    return x[0] if x[0] == x[1] else '{} | {}'.format(*x)

def sanity_check(Swd):
    '''
    Make sure the work directory passed in is correct
    '''
    f = open(Swd.location + '/log/amiinsane.txt')
    l = f.readlines()[0].strip()
    f.close()
    assert l == "No, you're not", l

    return

def ensure_identicle(Swd, wd, skip = None):
    '''
    Compare two work directories
    '''
    if skip == None:
        skip = []

    # Compare datatables
    for d in Swd.data_tables:
        if d in skip:
            continue

        db1 = Swd.get_db(d, return_none=False)
        db2 =  wd.get_db(d, return_none=False)

        if d == 'Ndb':
            db1 = db1[['reference', 'querry', 'ani']].sort_values(['reference', 'querry'])
            db2 = db2[['reference', 'querry', 'ani']].sort_values(['reference', 'querry'])

        assert compare_dfs(db1, db2, verbose=True), "{0} is not the same!".format(d)

    # Compare the clustering files
    pass

    # Compare the graphs
    pass