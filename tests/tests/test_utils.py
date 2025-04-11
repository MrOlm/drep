import os
import glob
import importlib
import pandas as pd
import logging
import shutil
import uuid

def load_solutions_wd():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), '../test_solutions/ecoli_wd')

def get_unique_test_dir():
    """Generate a unique test directory path"""
    unique_id = str(uuid.uuid4())
    base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../test_backend')
    return os.path.join(base_dir, f'test_dir_{unique_id}')

def load_test_wd_loc():
    return get_unique_test_dir()

def load_test_wd_loc_2():
    return get_unique_test_dir()

def load_test_genomes():
    return glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../genomes/*'))

def load_zipped_genomes():
    return glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../test_backend/zipped/*.gz'))

def load_broken_genome():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), '../test_backend/other/broken_genome.fasta')

def load_zipped_genome():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), '../test_backend/zipped/Enterococcus_casseliflavus_EC20.fasta.gz')

def load_random_test_dir():
    return get_unique_test_dir()

def load_solutions_taxonomy_wd():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), '../test_solutions/ecoli_taxonomy')

def load_large_genome_set():
    loc = '/Users/mattolm/Programs/testing_house/test_genomes/'
    genomes = glob.glob(loc + '*.fna')
    return genomes

def load_test_backend():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), '../test_backend/')

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
            print("DataFrames are not equal:")
            print("\nShape difference:")
            print(f"DataFrame 1 shape: {db1.shape}")
            print(f"DataFrame 2 shape: {db2.shape}")
            
            print("\nColumn differences:")
            cols1 = set(db1.columns)
            cols2 = set(db2.columns)
            if cols1 != cols2:
                print(f"Only in df1: {cols1 - cols2}")
                print(f"Only in df2: {cols2 - cols1}")
            
            print("\nDetailed differences:")
            print(e)
            
            if not db1.empty and not db2.empty:
                # Show a sample of differing rows
                try:
                    diff_mask = (db1 != db2).any(axis=1)
                    print("\nSample of differing rows:")
                    print("DataFrame 1:")
                    print(db1[diff_mask].head())
                    print("\nDataFrame 2:")
                    print(db2[diff_mask].head())
                except:
                    pass # Skip detailed diff if shapes don't match
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

class TestingClass():
    def setup(self):
        """Initialize unique test directories for this test instance"""
        self.wd_loc = get_unique_test_dir()
        self.test_dir = get_unique_test_dir()
        self.working_wd_loc = get_unique_test_dir()
        
        # Create the directories
        os.makedirs(self.wd_loc, exist_ok=True)
        os.makedirs(self.test_dir, exist_ok=True)
        os.makedirs(self.working_wd_loc, exist_ok=True)

    def teardown(self):
        importlib.reload(logging)
        # Clean up our unique test directories
        for path in [self.wd_loc, self.test_dir, self.working_wd_loc]:
            if os.path.isdir(path):
                shutil.rmtree(path)

def load_common_self():
    self = TestingClass()
    self.setup()  # Create unique directories
    
    # Set up
    self.genomes = load_test_genomes()
    self.zipped_genome = load_zipped_genome()
    self.test_dir = load_random_test_dir()

    self.wd_loc = load_test_wd_loc()
    self.s_wd_loc = load_solutions_wd()
    self.working_wd_loc = load_test_wd_loc_2()

    self.extra_weights_loc = os.path.join(load_test_backend(), 'extra_weights.tsv')
    self.stinker_genome = os.path.join(load_test_backend(), 'other/Enterococcus_faecalis_TX0104.fa')

    importlib.reload(logging)
    if os.path.isdir(self.wd_loc):
        shutil.rmtree(self.wd_loc)
    os.mkdir(self.wd_loc)

    if os.path.isdir(self.test_dir):
        shutil.rmtree(self.test_dir)
    os.mkdir(self.test_dir)

    if os.path.isdir(self.working_wd_loc):
        shutil.rmtree(self.working_wd_loc)

    # copy over the data from solutions directory
    os.mkdir(self.working_wd_loc)
    shutil.copytree(os.path.join(self.s_wd_loc, 'data'), \
                    os.path.join(self.working_wd_loc, 'data'))
    shutil.copytree(os.path.join(self.s_wd_loc, 'data_tables'), \
                    os.path.join(self.working_wd_loc, 'data_tables'))
    shutil.copytree(os.path.join(self.s_wd_loc, 'log'), \
                    os.path.join(self.working_wd_loc, 'log'))
    importlib.reload(logging)

    return self

