import glob
import logging
import os
import sys

import pandas as pd

import drep
import drep.d_cluster.cluster_utils
import drep.d_cluster.external
import drep.d_cluster.utils
import drep.d_cluster.greedy_clustering

class genomeChunk():
    """
    This is an object that just holds stuff related to a chunk of genomes for Mash clustering
    """

    def __init__(self, locations, number, sketch_folder, genome_names, no_create=False):
        '''
        Initialize this genome chunk
        '''
        self.name = "chunk_{0}".format(number)
        self.genome_locations = locations
        self.chunk_folder = os.path.join(sketch_folder, self.name)
        self.genome_names = genome_names

        if not no_create:
            if not os.path.exists(self.chunk_folder):
                os.makedirs(self.chunk_folder)

    def gen_sketch_cmds(self, mash_exe, MASH_s):
        cmds = []
        for location, name in zip(self.genome_locations, self.genome_names):
            file = os.path.join(self.chunk_folder, name)
            if not os.path.isfile(file + '.msh'):
                cmd = [mash_exe, 'sketch', location, '-s', str(MASH_s), '-o',
                       file]
                cmds.append(cmd)
        return cmds

    def gen_paste_cmd(self, mash_exe):
        all_file = os.path.join(self.chunk_folder, 'chunk_all.msh')
        cmd = [mash_exe, 'paste', all_file] \
              + glob.glob(os.path.join(self.chunk_folder, '*'))
        self.all_file = all_file
        return cmd

    def gen_dist_cmd(self, mash_exe, mash_folder, p):
        dist_file = os.path.join(mash_folder, '{0}_MASH_table.tsv'.format(self.name))
        self.dist_file = dist_file
        cmd = [mash_exe, 'dist', '-p', str(p), self.all_file, self.all_file, '>', dist_file]
        cmd = ' '.join(cmd)

        return cmd

    def load_mash_table(self):
        Mdb = drep.d_cluster.utils.parse_mash_table(self.dist_file)

        # Filter out those genomes that are in the MASH folder but shouldn't be in Mdb
        gs = set(self.genome_names)
        Mdb = Mdb[Mdb['genome1'].isin(gs)]
        Mdb = Mdb[Mdb['genome2'].isin(gs)]

        # Reorder categories to be correct
        for g in ['genome1', 'genome2']:
            Mdb[g] = Mdb[g].cat.remove_unused_categories()
            Mdb[g] = Mdb[g].cat.reorder_categories(sorted((Mdb[g].unique())), ordered=True)

        self.Mdb = Mdb

    def cluster_mash_table(self, **kwargs):
        if len(self.Mdb) > 1:
            Cdb, cluster_ret = cluster_mash_database(self.Mdb, **kwargs)
            self.Cdb = Cdb
        else:
            self.Cdb = pd.DataFrame({'primary_cluster':[1], 'genome':[self.Mdb['genome1'].tolist()[0]]})

    def get_winning_genomes(self, g2l, **kwargs):
        Cdb = self.Cdb
        Cdb['length'] = Cdb['genome'].map(g2l)
        return self.Cdb.sort_values('length').drop_duplicates(subset=['primary_cluster'], keep='last')['genome'].tolist()

def all_vs_all_MASH(Bdb, data_folder, **kwargs):
    """
    Run MASH pairwise within all samples in Bdb

    Args:
        Bdb: dataframe with genome, location
        data_folder: location to store temporary output files

    Keyword Args:
        MASH_sketch: size of mash sketches
        dry: dont actually run anything
        processors: number of processors to multithread with
        mash_exe: location of mash excutible (will try and find with shutil if not provided)
        groupSize: max number of mash sketches to hold in each folder
        debug: if True, log all of the commands
        wd: if you want to log commands, you also need the wd
    """
    # Set up the mash folder structure
    logdir, MASH_folder, sketch_folder, mash_exe = prepare_mash(data_folder, **kwargs)

    # Set up chunks of genomes
    genome_chunks = prepare_genome_chunks(Bdb, sketch_folder, MASH_folder,  **kwargs)
    if len(genome_chunks) > 1:
        logging.info(f"  Will split genomes into {len(genome_chunks)} groups for primary clustering")

    # Process the chunks individually
    genome_chunks = run_mash_on_genome_chunks(genome_chunks, mash_exe, sketch_folder, MASH_folder, logdir,  **kwargs)

    # If there's only one chunk, we're done here
    if len(genome_chunks) == 1:
        Mdb = genome_chunks[0].Mdb
        Cdb, cluster_ret = cluster_mash_database(Mdb, **kwargs)
        return Mdb, Cdb, cluster_ret

    # If there's multiple chunks, run a second round
    logging.info("  Final step: comparing between all groups")
    return run_second_round_clustering(Bdb, genome_chunks, data_folder, verbose=True, **kwargs)

def prepare_mash(data_folder, **kwargs):
    """
    Make some folders and things
    """
    append = kwargs.get('v2', '')

    # set up logdir
    if ('wd' in kwargs) and (kwargs.get('debug', False) == True):
        logdir = kwargs.get('wd').get_dir('cmd_logs')
    else:
        logdir = False

    # Find mash excutable
    mash_exe = kwargs.get('exe_loc', None)
    if mash_exe == None:
        mash_exe = drep.get_exe('mash')

    # Make a folder to hold this information
    MASH_folder = os.path.join(data_folder, 'MASH_files{0}/'.format(append))
    if not os.path.exists(MASH_folder):
        os.makedirs(MASH_folder)

    # Make a folder in there to store sketches
    sketch_folder = os.path.join(MASH_folder, 'sketches{0}/'.format(append))
    if not os.path.exists(sketch_folder):
        os.makedirs(sketch_folder)

    return logdir, MASH_folder, sketch_folder, mash_exe

def prepare_genome_chunks(Bdb, sketch_folder, MASH_folder, **kwargs):
    groupSize = kwargs.get('primary_chunksize', 5000)
    l2g = Bdb.set_index('location')['genome'].to_dict()

    locations = list(Bdb['location'].unique())
    chunks = [locations[x:x + groupSize] for x in range(0, len(locations), groupSize)]

    genome_chunks = []
    for i, chunk in enumerate(chunks):
        genome_chunks.append(genomeChunk(chunk, i, sketch_folder, [l2g[l] for l in chunk]))

    return genome_chunks

def run_mash_on_genome_chunks(genome_chunks, mash_exe, sketch_folder, MASH_folder, logdir, **kwargs):
    dry = kwargs.get('dry', False)
    p = kwargs.get('processors', 6)
    MASH_s = kwargs.get('MASH_sketch', 1000)
    multi_round = kwargs.get('multiround_primary_clustering', True)

    # Step 1) Create Mash sketches
    cmds = []
    for GC in genome_chunks:
        cmds += GC.gen_sketch_cmds(mash_exe, MASH_s)
    if (not dry) & (len(cmds) > 0):
        drep.thread_cmds(cmds, logdir=logdir, t=int(p))

    # Step 2) Combine MASH sketches within chunks
    cmds = [GC.gen_paste_cmd(mash_exe) for GC in genome_chunks]
    if (not dry) & (len(cmds) > 0):
        drep.thread_cmds(cmds, logdir=logdir, t=int(p))

    # Merge the pasted chunks and make a new genomeChunk if thats what you want
    if (not multi_round) & (len(genome_chunks) > 1):
        cmd, new_gc = drep.d_cluster.utils.merge_genome_chunks(mash_exe, genome_chunks, sketch_folder, MASH_folder)
        genome_chunks = [new_gc]
        drep.run_cmd(cmd, dry, shell=False, logdir=logdir)

    # Step 3) Run Mash on each chunk
    cmds = [GC.gen_dist_cmd(mash_exe, MASH_folder, p) for GC in genome_chunks]
    for j, cmd in enumerate(cmds):
        if not dry:
            if len(cmds) > 1:
                logging.info(f"  Comparing group {j+1} of {len(cmds)}")
            drep.run_cmd(cmd, dry, shell=True, logdir=logdir)

    # Step 4) Load the Mash tables of each chunk
    for GC in genome_chunks:
        GC.load_mash_table()

    return genome_chunks

def run_second_round_clustering(Bdb, genome_chunks, data_folder, **kwargs):
    verbose = kwargs.get('verbose', False)

    kwargs_copy = kwargs.copy()
    kwargs_copy['multiround_primary_clustering'] = False
    kwargs_copy['v2'] = '_v2'

    mdbs = []

    # Step 1) Create a merged Cdb file
    dbs = []
    for gc in genome_chunks:
        gc.cluster_mash_table(**kwargs_copy)
        cdb = gc.Cdb
        cdb['subcluster'] = ["{0}_{1}".format(gc.name, x) for x in cdb['primary_cluster']]
        dbs.append(cdb)

        mdb = gc.Mdb
        mdb['genome_chunk'] = gc.name
        mdbs.append(mdb)

    Cdb = pd.concat(dbs)

    # Step 2) Pick winners
    g2l = Bdb.set_index('genome')['length'].to_dict()
    Cdb['length'] = Cdb['genome'].map(g2l)
    second_round_genomes = Cdb.sort_values('length').drop_duplicates(subset=['subcluster'], keep='last')['genome'].tolist()

    if verbose:
        logging.info(f"Comparing {len(second_round_genomes):,} genomes")

    # Step 3) Run a second round
    logdir, MASH_folder, sketch_folder, mash_exe = prepare_mash(data_folder, **kwargs_copy)
    genome_chunks = prepare_genome_chunks(Bdb[Bdb['genome'].isin(second_round_genomes)], sketch_folder, MASH_folder, **kwargs_copy)
    genome_chunks = run_mash_on_genome_chunks(genome_chunks, mash_exe, sketch_folder, MASH_folder, logdir, **kwargs_copy)

    # Step 4) Get results
    assert len(genome_chunks) == 1

    mdb = genome_chunks[0].Mdb
    mdb['genome_chunk'] = 'v2'
    mdbs.append(mdb)
    Mdb = pd.concat(mdbs).reset_index(drop=True)

    Cdb2, cluster_ret = cluster_mash_database(mdb, **kwargs)
    Cdb2['primary_representitive'] = True

    # Step 5) Merge the new Cdb back in with the old
    del Cdb['primary_cluster']
    Cdb = pd.merge(Cdb, Cdb2, on='genome', how='outer')
    o2n = Cdb[Cdb['primary_representitive'] == True].set_index('subcluster')['primary_cluster'].to_dict()
    Cdb['primary_cluster'] = Cdb['subcluster'].map(o2n).astype(int)

    return Mdb, Cdb, cluster_ret

def cluster_mash_database(db, **kwargs):
    '''
    From a Mash database, cluster and return Cdb

    Args:
        db: Mdb (all_vs_all Mash results)

    Keyword arguments:
        clusterAlg: how to cluster database (default = single)
        P_ani: threshold to cluster at (default = 0.9)

    Returns:
        list: [Cdb, [linkage, linkage_db, arguments]]
    '''
    logging.debug('Clustering MASH database')

    # Load key words
    P_Lmethod = kwargs.get('clusterAlg','single')
    P_Lcutoff = 1 - kwargs.get('P_ani',.9)

    # Do the actual clustering
    db['dist'] = 1 - db['similarity']
    linkage_db = db.pivot("genome1","genome2","dist")
    Cdb, linkage = drep.d_cluster.cluster_utils.cluster_hierarchical(linkage_db, linkage_method= P_Lmethod, \
                                                                     linkage_cutoff= P_Lcutoff)
    Cdb = Cdb.rename(columns={'cluster':'primary_cluster'})
    Cdb['primary_cluster'] = Cdb['primary_cluster'].astype(int)

    # Preparing clustering for return
    arguments = {'linkage_method':P_Lmethod,'linkage_cutoff':P_Lcutoff,\
                    'comparison_algorithm':'MASH'}
    cluster_ret = [linkage, linkage_db, arguments]

    return Cdb, cluster_ret

def secondary_clustering(Bdb, Cdb, algorithm, data_folder, **kwargs):
    if kwargs.get('greedy_secondary_clustering', False) != True:
        Ndb = pd.DataFrame()
        for bdb, name in iteratre_clusters(Bdb, Cdb, id='primary_cluster'):
            logging.debug('running cluster {0}'.format(name))
            # logging.debug('total memory - {0:.2f} Mbp'.format(int(process.memory_info().rss)/1000000))
            ndb = compare_genomes(bdb, algorithm, data_folder, **kwargs)

            if len(ndb) == 0:
                logging.error("CRITICAL ERROR WITH PRIMARY CLUSTER {0}; TRYING AGAIN".format(name))
                ndb = compare_genomes(bdb, algorithm, data_folder, **kwargs)

            if len(ndb) > 0:
                ndb['primary_cluster'] = name
                Ndb = Ndb.append(ndb)
            else:
                logging.error("DOUBLE CRITICAL ERROR AGAIN WITH PRIMARY CLUSTER {0}; SKIPPING".format(name))

        # Run clustering on Ndb
        Cdb, c2ret = drep.d_cluster.utils._cluster_Ndb(Ndb, comp_method=algorithm, **kwargs)

        return Ndb, Cdb, c2ret

    else:
        return drep.d_cluster.greedy_clustering.greedy_secondary_clustering(Bdb, Cdb, algorithm, data_folder, **kwargs)


def iteratre_clusters(Bdb, Cdb, id='MASH_cluster'):
    Bdb = pd.merge(Bdb, Cdb)
    for cluster in Bdb[id].unique():
        d = Bdb[Bdb[id] == cluster]
        yield d, cluster

def compare_genomes(bdb, algorithm, data_folder, **kwargs):
    '''
    Compare a list of genomes using the algorithm specified

    This method takes in bdb (a table with the columns location and genome), runs
    pair-wise comparisons between all genomes in the sample, and returns a table
    with at least the columns 'reference', 'querry', 'ani','coverage', depending
    on what algorithm is called

    Args:
        bdb: DataFrame with ['genome', 'location'] (drep.d_filter.load_genomes)
        algorithm: options are ANImf, ANIn, gANI
        data_folder: location to store output files

    Keyword Arguments:
        wd: either this or prod_folder needed for gANI
        prod_folder: either this or wd needed for gANI

    Return:
        DataFrame: Ndb (['reference', 'querry', 'ani','coverage'])
    '''
    # To handle other versions of this method which passed in a WorkDirectory
    # instead of data_folder string
    if isinstance(data_folder, drep.WorkDirectory.WorkDirectory):
        data_folder = data_folder.get_dir('data')

    if not kwargs.get('greedy_secondary_clustering', False):
        if algorithm == 'ANImf':
            genome_list = bdb['location'].tolist()
            working_data_folder = os.path.join(data_folder, 'ANImf_files/')
            df = drep.d_cluster.external.run_pairwise_ANImf(genome_list, working_data_folder, **kwargs)
            return df

        elif algorithm == 'ANIn':
            genome_list = bdb['location'].tolist()
            working_data_folder = os.path.join(data_folder, 'ANIn_files/')
            df = drep.d_cluster.utils.run_pairwise_ANIn(genome_list, working_data_folder, **kwargs)
            return df

        elif algorithm == 'fastANI':
            genome_list = bdb['location'].tolist()
            working_data_folder = os.path.join(data_folder, 'fastANI_files/')
            df = drep.d_cluster.external.run_pairwise_fastANI(genome_list, working_data_folder, **kwargs)
            return df

        elif algorithm == 'gANI':
            # Figure out prodigal folder
            wd = kwargs.get('wd', False)
            if not wd:
                prod_folder = kwargs.pop('prod_folder', False)
                assert prod_folder != False
            else:
                prod_folder = wd.get_dir('prodigal')

            working_data_folder = os.path.join(data_folder, 'gANI_files/')
            df = drep.d_cluster.external.run_pairwise_gANI(bdb, working_data_folder, \
                                                           prod_folder=prod_folder, **kwargs)
            return df

        elif algorithm == 'goANI':
            # Figure out prodigal folder
            wd = kwargs.get('wd', False)
            if not wd:
                prod_folder = kwargs.pop('prod_folder', False)
                assert prod_folder != False
            else:
                prod_folder = wd.get_dir('prodigal')

            working_data_folder = os.path.join(data_folder, 'goANI_files/')
            df = drep.d_cluster.external.run_pairwise_goANI(bdb, working_data_folder, \
                                                            prod_folder=prod_folder, **kwargs)
            return df

        else:
            logging.error("{0} not supported".format(algorithm))
            sys.exit()

    else:
        SUPPORTED = ['fastANI']
        if algorithm not in SUPPORTED:
            message = f"{algorithm} is not supported for greedy secondary clustering!\nChoose one of the following supported S_algorithm options: {' '.join(SUPPORTED)}"
            logging.error(message)
            print(message)
            raise NameError

        working_data_folder = os.path.join(data_folder, 'greedy_clustering/')
        return drep.d_cluster.greedy_clustering.compare_genomes_greedy(bdb, algorithm, working_data_folder, **kwargs)


