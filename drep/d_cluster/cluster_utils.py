import logging
import os
import sys

import numpy as np
import pandas as pd
import scipy.cluster
from scipy.spatial import distance as ssd
import networkx as nx

import drep.d_cluster.utils

def genome_hierarchical_clustering(Ndb, **kwargs):
    '''
    Cluster ANI database

    Args:
        Ndb: result of secondary clustering

    Keyword arguments:
        clusterAlg: how to cluster the database (default = single)
        S_ani: thershold to cluster at (default = .99)
        cov_thresh: minumum coverage to be included in clustering (default = .5)
        cluster: name of the cluster
        comp_method: comparison algorithm used

    Returns:
        list: [Cdb, {cluster:[linkage, linkage_db, arguments]}]
    '''
    logging.debug('Clustering ANIn database')

    S_Lmethod = kwargs.get('clusterAlg', 'single')
    S_Lcutoff = 1 - kwargs.get('S_ani', .99)
    cov_thresh = float(kwargs.get('cov_thresh',0.5))
    cluster = kwargs.get('cluster','')
    comp_method = kwargs.get('comp_method', 'unk')

    Table = {'genome':[],'secondary_cluster':[]}

    # Handle the case where there's only one genome
    if len(Ndb['reference'].unique()) == 1:
        Table['genome'].append(os.path.basename(Ndb['reference'].unique().tolist()[0]))
        Table['secondary_cluster'].append("{0}_0".format(cluster))
        cluster_ret = []

    else:
        # Make linkage Ndb
        Ldb = drep.d_cluster.utils.make_linkage_Ndb(Ndb, **kwargs)

        # 3) Cluster the linkagedb
        Gdb, linkage = cluster_hierarchical(Ldb, linkage_method= S_Lmethod, \
                                    linkage_cutoff= S_Lcutoff)

        # 4) Extract secondary clusters
        for clust, d in Gdb.groupby('cluster'):
            for genome in d['genome'].tolist():
                Table['genome'].append(genome)
                Table['secondary_cluster'].append("{0}_{1}".format(cluster,clust))

        # 5) Save the linkage
        arguments = {'linkage_method':S_Lmethod,'linkage_cutoff':S_Lcutoff,\
                    'comparison_algorithm':comp_method,'minimum_coverage':cov_thresh}
        cluster_ret = [linkage, Ldb, arguments]

    # Return the database
    Gdb = pd.DataFrame(Table)
    Gdb['threshold'] = S_Lcutoff
    Gdb['cluster_method'] = S_Lmethod
    Gdb['comparison_algorithm'] = comp_method

    return Gdb, cluster_ret


def iteratre_clusters(Bdb, Cdb, id='primary_cluster'):
    '''
    An iterator: Given Bdb and Cdb, yeild smaller Bdb's in the same cluster

    Args:
        Bdb: [genome, location]
        Cdb: [genome, id]
        id: what to iterate on (default = 'primary_cluster')

    Returns:
        list: [d(subset of b), cluster(name of cluster)]
    '''
    Bdb = pd.merge(Bdb,Cdb)
    for cluster, d in Bdb.groupby(id):
        yield d, cluster


def cluster_threshold_graph_optimized(db, linkage_cutoff=0.10, linkage_method='single'):
    # Filter distances below threshold first
    filtered_edges = db[db['dist'] <= linkage_cutoff]
    
    # Create graph directly from filtered edges
    G = nx.Graph()
    G.add_edges_from(filtered_edges[['genome1', 'genome2']].values)
    
    # Log that we're using the optimized method
    logging.debug("Using low-RAM optimized clustering method with {0} edges".format(len(filtered_edges)))
    
    # Find connected components
    clusters = {}
    for cluster_id, component in enumerate(nx.connected_components(G)):
        for genome in component:
            clusters[genome] = cluster_id + 1
            
    # Add isolated nodes (if needed)
    all_genomes = set(db['genome1']).union(set(db['genome2']))
    for genome in all_genomes:
        if genome not in clusters:
            clusters[genome] = len(clusters)
            
    # Return clusters and a special value indicating we used the optimized method
    return clusters, "optimized_method_used"

def cluster_hierarchical(db, linkage_method= 'single', linkage_cutoff= 0.10, low_ram=False):
    '''
    Perform hierarchical clustering on a symmetrical distiance matrix

    Args:
        db: result of db.pivot usually
        linkage_method: passed to scipy.cluster.hierarchy.fcluster
        linkage_cutoff: distance to draw the clustering line (default = .1)
        low_ram: whether to use the memory-efficient algorithm

    Returns:
        list: [Cdb, linkage]
    '''
    # Save names
    names = list(db.columns)

    if low_ram:
        # Convert to long format for optimized clustering
        db_long = db.stack().reset_index()
        db_long.columns = ['genome1', 'genome2', 'dist']
        clusters, _ = cluster_threshold_graph_optimized(db_long, linkage_cutoff, linkage_method)
        
        # Convert clusters to Cdb format
        Cdb = pd.DataFrame({'genome': list(clusters.keys()), 
                          'cluster': list(clusters.values())})
        return Cdb, _
    else:
        # Generate linkage dataframe
        arr =  np.asarray(db)
        try:
            arr = ssd.squareform(arr)
        except:
            logging.error("The database passed in is not symmetrical!")
            logging.error(arr)
            logging.error(names)
            sys.exit()
        linkage = scipy.cluster.hierarchy.linkage(arr, method= linkage_method)

        # Form clusters
        fclust = scipy.cluster.hierarchy.fcluster(linkage,linkage_cutoff, \
                        criterion='distance')
        # Make Cdb
        Cdb = drep.d_cluster.utils._gen_cdb_from_fclust(fclust,names)

        return Cdb, linkage