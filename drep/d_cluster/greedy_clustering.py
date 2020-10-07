import glob
import logging
import os
import shutil
import sys
import pandas as pd

import drep.d_cluster.external
import drep.d_cluster.compare_utils

def greedy_secondary_clustering(Bdb, Cdb, algorithm, data_folder, **kwargs):
    ndbs = []
    cdbs = []
    c2ret = {}
    for bdb, name in drep.d_cluster.compare_utils.iteratre_clusters(Bdb, Cdb, id='primary_cluster'):
        logging.debug('running cluster {0} with {1} genomes'.format(name, len(bdb)))
        # logging.debug('total memory - {0:.2f} Mbp'.format(int(process.memory_info().rss)/1000000))
        ndb, cdb, ret = drep.d_cluster.compare_utils.compare_genomes(bdb, algorithm, data_folder, **kwargs)

        if len(ndb) == 0:
            logging.error("CRITICAL ERROR WITH PRIMARY CLUSTER {0}; TRYING AGAIN".format(name))
            ndb, cdb, ret = drep.d_cluster.compare_utils.compare_genomes(bdb, algorithm, data_folder, **kwargs)

        if len(ndb) > 0:
            ndb['primary_cluster'] = name
            cdb['primary_cluster'] = name
            cdb['secondary_cluster'] = ["{0}{1}".format(name, y) for y in cdb['secondary_cluster']]
            ndbs.append(ndb)
            cdbs.append(cdb)
            c2ret[name] = ret
        else:
            logging.error("DOUBLE CRITICAL ERROR AGAIN WITH PRIMARY CLUSTER {0}; SKIPPING".format(name))

    Ndb = pd.concat(ndbs).reset_index(drop=True)
    Cdb = pd.concat(cdbs).reset_index(drop=True)
    return Ndb, Cdb, c2ret


def compare_genomes_greedy(bdb, algorithm, data_folder, **kwargs):
    ani_thresh = float(kwargs.get('S_ani', .99))
    cov_thresh = float(kwargs.get('cov_thresh', 0.5))

    assert (ani_thresh <= 1) & (ani_thresh > 0)

    # Set genome order
    odb = order_genomes_for_greedy(bdb, **kwargs)

    # Set up
    cluster = kwargs.get('cluster', '')
    genome_rep_file = os.path.join(data_folder + 'representative_genome_locations.txt')
    genome_reps = []
    rep2cluster = {}
    genome2cluster = {}
    ndbs = []

    if os.path.exists(genome_rep_file):
        os.remove(genome_rep_file)

    # Prepare for greedy clustering
    kwargs = prepare_for_greedy(algorithm, data_folder, **kwargs)

    # Iterate
    j = 1
    for i, row in odb.iterrows():
        if len(genome_reps) == 0:
            make_new_cluster = True

        else:
            ndb = genome_vs_reps(row['location'], genome_reps, genome_rep_file, algorithm, data_folder, **kwargs)
            ndbs.append(ndb)

            cluster_rep = get_cluster_rep(ndb, ani_thresh, cov_thresh)

            if cluster_rep is not False:
                make_new_cluster = False
                assert cluster_rep in rep2cluster, [cluster_rep, rep2cluster]
                genome2cluster[row['genome']] = rep2cluster[cluster_rep]
            else:
                make_new_cluster = True

        if make_new_cluster:
            new_cluster = "{0}_{1}".format(cluster,j)
            j += 1

            genome_reps.append(row['location'])
            rep2cluster[row['genome']] = new_cluster
            genome2cluster[row['genome']] = new_cluster
            with open(genome_rep_file, "a") as myfile:
                myfile.write(row['location'] + '\n')

    if len(ndbs) > 0:
        Ndb = pd.concat(ndbs)

    else:
        # Add self-comparisons if there is only one genome
        Table = {'querry': [], 'reference': [], 'ani': [], 'alignment_coverage': []}
        for g in odb['location'].tolist():
            Table['reference'].append(drep.d_cluster.utils._get_genome_name_from_fasta(g))
            Table['querry'].append(drep.d_cluster.utils._get_genome_name_from_fasta(g))
            Table['ani'].append(1)
            Table['alignment_coverage'].append(1)
        Ndb = pd.DataFrame(Table)

    Cdb, cluster_ret = generate_greedy_cdb(bdb, rep2cluster, genome2cluster, algorithm, ani_thresh, cov_thresh, **kwargs)
    return Ndb, Cdb, cluster_ret

def genome_vs_reps(new_genome, genome_reps, genome_rep_file, algorithm, data_folder, **kwargs):
    if algorithm == 'fastANI':
        # Return Ndb
        return drep.d_cluster.external.fastani_one_vs_many(new_genome, genome_reps, genome_rep_file, data_folder, **kwargs)
    else:
        logging.error("{0} algorithm is not yet supported for greedy clustering; sorry!")
        assert False


def prepare_for_greedy(algorithm, data_folder, **kwargs):
    if algorithm == 'fastANI':
        # Make folders
        if not os.path.exists(data_folder):
            os.makedirs(data_folder)
        tmp_dir = os.path.join(data_folder, 'tmp/')
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

        # Handle log directory
        if ('wd' in kwargs) and (kwargs.get('debug', False)):
            logdir = kwargs.get('wd').get_dir('cmd_logs')
        else:
            logdir = False

        kwargs['tmp_dir'] = tmp_dir
        kwargs['logdir'] = logdir
        kwargs['current_exe'] = drep.get_exe('fastANI')

    return kwargs

def order_genomes_for_greedy(bdb, **kwargs):
    return bdb.sort_values('length', ascending=False)

def get_cluster_rep(ndb, ani_thresh, cov_thresh):
    fdb = ndb[(ndb['ani'] >= ani_thresh) & (ndb['alignment_coverage'] >=  cov_thresh)]
    if len(fdb) > 0:
        return fdb['querry'].iloc[0]
    else:
        return False

def generate_greedy_cdb(bdb, rep2cluster, genome2cluster, algorithm, ani_thresh, c_thresh, **kwargs):
    reps = set(list(rep2cluster.keys()))

    cdb = bdb[['genome']]
    cdb['secondary_cluster'] = cdb['genome'].map(genome2cluster)
    cdb['threshold'] = 1 - ani_thresh
    cdb['cluster_method'] = 'greedy'
    cdb['comparison_algorithm'] = algorithm
    cdb['greedy_representative'] = [g in reps for g in cdb['genome']]

    arguments = {'linkage_method': 'greedy', 'linkage_cutoff': ani_thresh,
                 'comparison_algorithm': algorithm, 'minimum_coverage': c_thresh}
    cluster_ret = [None, pd.DataFrame(), arguments]

    return cdb, cluster_ret

