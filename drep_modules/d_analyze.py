#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import os
import seaborn as sns
from matplotlib import pyplot as plt

"""#############################################################################
                            MODULE ARCHITECTURE

****************************    Clustering   *******************************
Graphs - Heatmaps

*   plot_Mdb_heatmap(Mdb)
    -   Make a single heat-map of Mdb
    -   Ability to color clusters?
*   plot_MASH_clusters_ANIn(Ndb)
    -   Make a heat-map of ANIn for every MASH cluster
*   plot_MASH_clusters_cov(Ndb)
    -   Make a heat-map of ANI_alignment_coverage for every MASH cluster
*   plot_ANIn_clusters(Ndb, Cdb)
    -   Make a heat-map for every Ndb cluster
    
Graphs - Scatterplots

*   plot_MASH_vs_ANIn_ani(Mdb, Ndb)
    - Plot MASH_ani vs. ANIn_ani (including correlation)
*   plot_MASH_vs_ANIn_cov(Mdb, Ndb)
    - Plot MASH_ani vs. ANIn_cov (including correlation)
*   plot_ANIn_vs_ANIn_cov(Mdb, Ndb)
    - Plot ANIn vs. ANIn_cov (including correlation)
*   plot_MASH_vs_len(Mdb, Ndb)
    - Plot MASH_ani vs. length_difference (including correlation)
*   plot_ANIn_vs_len(Ndb)
    - Plot ANIn vs. length_difference (including correlation)
    
Graphs - Custom

*   plot_cluster_tightness(Ndb)
    - Come up with some way of visualizing the variation within ANIn clusters,
      versus the variation between clusters
      - Show the average and max tightness within and between all clusters
################################################################################
"""


"""
HEAT MAPS
"""

def plot_Mdb_heatmap(Mdb):
    db = Mdb.pivot("genome1","genome2","similarity")
    g = sns.clustermap(db)
    g.fig.suptitle("MASH ANI")
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    return g
    
def plot_MASH_clusters_ANIn(Ndb):
    gs = []
    for Mcluster in Ndb['MASH_cluster'].unique():
        db = Ndb[Ndb['MASH_cluster'] == Mcluster]
        d = db.pivot("reference","querry","ani")
        g = sns.clustermap(d)
        g.fig.suptitle("MASH cluster {0} - ANIn".format(Mcluster))
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        gs.append(g)
    return gs
    
def plot_MASH_clusters_cov(Ndb):
    gs = []
    for Mcluster in Ndb['MASH_cluster'].unique():
        db = Ndb[Ndb['MASH_cluster'] == Mcluster].copy()
        d = db.pivot("reference","querry","alignment_coverage")
        g = sns.clustermap(d)
        g.fig.suptitle("MASH cluster {0} - Alignment Coverage".format(Mcluster))
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        gs.append(g)
    return gs
    
"""
SCATTER PLOTS
"""

def plot_MASH_vs_ANIn_ani(Mdb,Ndb,exclude_zero_MASH=True):
    mdb = Mdb.copy()
    mdb.rename(columns={'genome1':'querry','genome2':'reference',
                        'similarity':'MASH_ANI'},inplace=True)
    if exclude_zero_MASH:
        mdb= mdb[mdb['MASH_ANI'] > 0]
        
    db = pd.merge(mdb,Ndb)
    db.rename(columns={'ani':'ANIn'},inplace=True)
    g = sns.jointplot(x='ANIn',y='MASH_ANI',data=db)
    return g

def plot_MASH_vs_ANIn_cov(Mdb,Ndb,exclude_zero_MASH=True):
    mdb = Mdb.copy()
    mdb.rename(columns={'genome1':'querry','genome2':'reference',
                        'similarity':'MASH_ANI'},inplace=True)
    if exclude_zero_MASH:
        mdb= mdb[mdb['MASH_ANI'] > 0]
        
    db = pd.merge(mdb,Ndb)
    db.rename(columns={'alignment_coverage':'ANIn_cov'},inplace=True)
    g = sns.jointplot(x='ANIn_cov',y='MASH_ANI',data=db)
    return g
    
def plot_ANIn_vs_ANIn_cov(Ndb):
    db = Ndb.copy()
    db.rename(columns={'alignment_coverage':'ANIn_cov','ani':'ANIn'},inplace=True)
    g = sns.jointplot(x='ANIn_cov',y='ANIn',data=db)
    return g
    
"""
OTHER
"""

def test_clustering():
    print("You should make some test cases here!")

if __name__ == '__main__':
	test_clustering()