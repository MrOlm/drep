#!/usr/bin/env python3

import logging
import glob
import pandas as pd
import os
import shutil
import numpy as np
import pickle

import drep.WorkDirectory
import drep as dm
import drep.d_cluster as dClust
import drep.d_filter as d_filter
import drep.d_choose as dChoose

def d_adjust_wrapper(wd,**kwargs):
    # Load the WorkDirectory.
    logging.info("Loading work directory")
    wd = drep.WorkDirectory.WorkDirectory(wd)
    logging.debug(str(wd))

    if kwargs.get('cluster') != None:
        logging.info('adjusting cluster {0}'.format(kwargs.get('cluster')))
        adjust_cluster_wrapper(wd, **kwargs)

    elif kwargs.get('rm_cluster') != None:
        logging.info('removing cluster(s) {0}'.format(' '.join(kwargs.get('rm_cluster'))))
        remove_cluster_wrapper(wd, **kwargs)

def remove_cluster_wrapper(wd, **kwargs):
    # Get information
    Cdb = wd.get_db('Cdb')
    Wdb = wd.get_db('Wdb')

    for cluster in kwargs.get('rm_cluster'):
        type = cluster_type(cluster)

        if type == 'primary_cluster':
            # Make sure cluster is in Cdb
            if int(cluster) not in Cdb[type].tolist():
                logging.error("{0} is not in Cdb- quitting".format(cluster))
                sys.exit()

            # Make sure cluster is Wdb
            W_Pclusters = [x.split('_')[0] for x in Wdb['cluster'].tolist()]
            if cluster not in W_Pclusters:
                logging.error("{0} is not in Wdb- quitting".format(cluster))
                sys.exit()

            # Remove the cluster
            remove_primary_cluster(cluster, wd, **kwargs)

        elif type == 'secondary_cluster':

            # Make sure cluster is in Cdb
            if cluster not in Cdb[type].tolist():
                logging.error("{0} is not in Cdb- quitting".format(cluster))
                sys.exit()

            # Make sure cluster is in Wdb
            if cluster not in Wdb['cluster'].tolist():
                logging.error("{0} is not in Wdb- quitting".format(cluster))
                sys.exit()

            # Remove the cluster
            remove_secondary_cluster(cluster, wd, **kwargs)


def remove_primary_cluster(Rcluster, wd, **kwargs):
    # Get things to alter
    Cdb = wd.get_db('Cdb')
    Wdb = wd.get_db('Wdb')

    # If this is a singleton, just remove the associated secondary cluster
    if len(Cdb['genome'][Cdb['primary_cluster'] == int(Rcluster)].unique()) == 1:
        Rsecondary = "{0}_0".format(Rcluster)
        logging.error("{0} is a primary singleton- will remove {1} instead".format(Rcluster, Rsecondary))
        remove_secondary_cluster(Rsecondary, wd, **kwargs)
        return

    # Find the genome(s) to remove
    Wdb['pclust'] = [x.split('_')[0] for x in Wdb['cluster'].tolist()]
    Rgenomes = Wdb['genome'][Wdb['pclust'] == Rcluster].tolist()
    for Rgenome in Rgenomes:
        assert os.path.isfile('{0}/dereplicated_genomes/{1}'.format(wd.location, Rgenome))

    # Find the pickle to remove
    Rpickle = "{0}/data/Clustering_files/secondary_linkage_cluster_{1}.pickle".format(\
            wd.location, Rcluster)
    assert os.path.isfile(Rpickle)

    logging.info("will remove cluster {0}, genomes {1}, and pickle {2}".format(Rcluster, \
            ' '.join(Rgenomes), os.path.basename(Rpickle)))

    # Remove cluster from Wdb
    newWdb = Wdb[Wdb['pclust'] != Rcluster]
    del newWdb['pclust']

    # Remove cluster from Cdb
    newCdb = Cdb[Cdb['primary_cluster'] != int(Rcluster)]

    # Remove genomes
    for Rgenome in Rgenomes:
        os.remove('{0}/dereplicated_genomes/{1}'.format(wd.location, Rgenome))

    # Remove pickle
    os.remove(Rpickle)

    # Save removed stuff
    wd.store_db(newCdb,'Cdb',overwrite=True)
    wd.store_db(newWdb,'Wdb',overwrite=True)

    logging.info("done")


def remove_secondary_cluster(Rcluster, wd, **kwargs):
    # Get things to alter
    Cdb = wd.get_db('Cdb')
    Wdb = wd.get_db('Wdb')

    # Find the genome to remove
    Rgenome = Wdb['genome'][Wdb['cluster'] == Rcluster].tolist()[0]
    assert os.path.isfile('{0}/dereplicated_genomes/{1}'.format(wd.location, Rgenome))

    logging.info("will remove {0} and genome {1}".format(Rcluster, Rgenome))

    # Remove cluster from Wdb
    newWdb = Wdb[Wdb['cluster'] != Rcluster]

    # Remove cluster from Cdb
    newCdb = Cdb[Cdb['secondary_cluster'] != Rcluster]

    # Remove genome
    os.remove('{0}/dereplicated_genomes/{1}'.format(wd.location, Rgenome))

    # Save removed stuff
    wd.store_db(newCdb,'Cdb',overwrite=True)
    wd.store_db(newWdb,'Wdb',overwrite=True)

    logging.info("done")

def cluster_type(cluster):
    if '_' in cluster:
        return 'secondary_cluster'
    else:
        return 'primary_cluster'

def adjust_cluster_wrapper(wd, **kwargs):
    # Validate arguments
    cluster = kwargs.get('cluster')
    if cluster == None:
        logging.error("Must specify a cluster")
        sys.exit()
    comp_method = kwargs.get('clustering_method')
    clust_method = kwargs.get('clusterAlg')
    threshold = kwargs.pop('threshold',None)
    cov_thresh = float(kwargs.get('minimum_coverage'))
    if threshold != None: threshold = 1- float(threshold)

    # Make a bdb listing the genomes to cluster
    Cdb = wd.get_db('Cdb')
    Bdb = wd.get_db('Bdb')
    genomes = Cdb['genome'][Cdb['primary_cluster'] == int(cluster)].tolist()
    bdb = Bdb[Bdb['genome'].isin(genomes)]

    # Make the comparison database
    Xdb = dClust.compare_genomes(bdb,comp_method,wd,**kwargs)

    # Remove values without enough coverage
    Xdb.loc[Xdb['alignment_coverage'] <= cov_thresh, 'ani'] = 0

    # Make it symmetrical
    Xdb['av_ani'] = Xdb.apply(lambda row: dClust.average_ani (row,Xdb),axis=1)
    Xdb['dist'] = 1 - Xdb['av_ani']
    db = Xdb.pivot("reference","querry","dist")

    # Cluster it
    if threshold == None:
        threshold = float(0)
    cdb, linkage = dClust.cluster_hierarchical(db, linkage_method = clust_method, \
                            linkage_cutoff = threshold)

    # Save the pickle
    data_folder = wd.location + '/data/Clustering_files/'
    arguments = {'linkage_method':clust_method,'linkage_cutoff':threshold,\
                        'comparison_algorithm':comp_method,'minimum_coverage':cov_thresh}
    pickle_name = "secondary_linkage_cluster_{0}.pickle".format(cluster)
    logging.debug('Saving secondary_linkage pickle {1} to {0}'.format(data_folder,\
                                                        pickle_name))
    with open(data_folder + pickle_name, 'wb') as handle:
        pickle.dump(linkage, handle)
        pickle.dump(db,handle)
        pickle.dump(arguments,handle)

    # Change the cdb
    cdb['cluster_method'] = clust_method
    cdb['comparison_algorithm'] = comp_method
    cdb['threshold'] = str(threshold)
    cdb['primary_cluster'] = int(cluster)
    cdb['secondary_cluster'] = ["{0}_{1}".format(cluster, x) for x in cdb['cluster']]
    del cdb['cluster']

    # Incorporate it into the existing Cdb
    Cdb = Cdb[Cdb['primary_cluster'] != int(cluster)]
    Cdb = pd.concat([Cdb,cdb], ignore_index=True)

    # Save it
    wd.store_db(Cdb,'Cdb',overwrite=True)

    if (not kwargs.get('skip_winner',False)) and wd.hasDb('Wdb'):

        # Remake Wdb based on Sdb
        oriWdb = wd.get_db('Wdb')
        newWdb = dChoose.pick_winners(wd.get_db('Sdb'),Cdb)
        change = accounce_changes(newWdb, oriWdb)
        wd.store_db(newWdb,'Wdb',overwrite=True)

        if change:
            # Change the ./dereplicated genomes thing
            logging.info("Remaking ./dereplicated genomes...")
            output_folder = wd.location + '/dereplicated_genomes/'
            dm.clobber_dir(output_folder,dry=kwargs.get('dry',False),\
                            overwrite=True)

            for genome in newWdb['genome'].unique():
                loc = Bdb['location'][Bdb['genome'] == genome].tolist()[0]
                shutil.copy2(loc, "{0}{1}".format(output_folder,genome))

def accounce_changes(newWdb, oriWdb):
    changes = False

    for cluster in newWdb['cluster'].unique():
        w1 = newWdb['genome'][newWdb['cluster'] == cluster].unique()[0]
        if cluster not in oriWdb['cluster'].tolist():
            change = "{0} is the winner of the new cluster {1}".format(w1,cluster)
            logging.info(change)
            changes = True
            continue
        w2 = oriWdb['genome'][oriWdb['cluster'] == cluster].unique()[0]
        if w1 != w2:
            change = "the winner of cluster {0} is now {1} (was {2})".format(\
                        cluster, w1, w2)
            logging.info(change)
            changes = True
            continue

    for cluster in oriWdb['cluster'].unique():
        if cluster not in newWdb['cluster'].tolist():
            change = "cluster {0} no longer exists".format(cluster)
            logging.info(change)
            changes = True

    if not changes:
        change = "No secondary clusters were made or destroyed, and none have new winners"
        logging.info(change)

    return changes


def test_adjust():
    print("Write this you lazy bum")

if __name__ == '__main__':
	test_adjust()
