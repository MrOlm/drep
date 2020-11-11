#!/usr/bin/env python3

import logging
import os
import pandas as pd
import sys

import drep
import drep.WorkDirectory
import drep.d_cluster.compare_utils
import drep.d_cluster.controller
import drep.d_filter
import drep.d_cluster
import drep.d_analyze

def d_evaluate_wrapper(wd,**kwargs):

    logging.debug("Loading work directory")
    wd = drep.WorkDirectory.WorkDirectory(wd)
    logging.debug(str(wd))

    # Run tertiary clustering if that's what you're here to do
    if kwargs.get('run_tertiary_clustering', False):
        logging.info("Running tertiary clustering on genome representatives")
        run_tertiary_clustering(wd, **kwargs)

    # Determine what to evaluate
    evs = kwargs.get('evaluate')
    options = ['1','2','3']
    to_eval = drep.d_analyze._parse_plot_options(options, evs)
    logging.debug("evaluating {0}".format(to_eval))

    # 1) Evaluate de-replicated genome similarity
    if '1' in to_eval:
        logging.info('will compare winners')
        Wmdb, Wndb = compare_winners(wd,**kwargs)

        # Save databases
        wd.store_db(Wmdb,'Wmdb',overwrite=True)
        wd.store_db(Wndb,'Wndb',overwrite=True)

    # 2) Throw warnings for clusters that were almost different
    if '2' in to_eval:
        logging.info('will provide warnings about clusters')
        warnings = evaluate_warnings(wd, **kwargs)

        # Save it
        warn_log = os.path.join(wd.location, 'log/warnings.txt')
        with open(warn_log, 'w') as file:
            file.write('\n'.join(warnings))
            file.write('\n')
        logging.info("{0} warnings generated: saved to {1}".format(len(warnings),warn_log))

    # 3) Generate a database of information on winning genomes
    if '3' in to_eval:
        logging.info('will produce Widb (winner information db)')
        Widb = evaluate_winners(wd, **kwargs)

        # Save it
        wd.store_db(Widb,'Widb',overwrite=True)
        loc = wd.location + 'data_tables/Widb.csv'
        logging.info("Winner database saved to {0}".format(loc))

    return

def compare_winners(wd, **kwargs):
    Bdb = wd.get_db('Bdb')
    Wdb = wd.get_db('Wdb')
    Bdb = Bdb[Bdb['genome'].isin(Wdb['genome'].tolist())]
    data_folder = os.path.join(wd.location, 'data/')
    comp_method = kwargs.get('comp_method','ANIn')

    # Generate MASH db
    Wmdb = drep.d_cluster.compare_utils.all_vs_all_MASH(Bdb, data_folder)

    # Generate ANIn db
    Wndb = drep.d_cluster.compare_utils.compare_genomes(Bdb, comp_method, wd, **kwargs)

    return Wmdb, Wndb

def evaluate_warnings(wd, **kwargs):
    Cdb = wd.get_db('Cdb')
    Ndb = wd.get_db('Ndb')
    warn_dist = float(kwargs.get('warn_dist',.25))

    if 'Blank' in Ndb:
        logging.error("Ndb is blank (was run with --SkipSecondary) - "\
                + "will skip cluster evaluation")

    if wd.hasDb('Wmdb'):
        Wmdb = wd.get_db('Wmdb')
        Wndb = wd.get_db('Wndb')
        warn_sim = float(kwargs.get('warn_sim',.98))
        warn_aln = float(kwargs.get('warn_aln',.25))

    warnings = []

    # Check for cases where there is a close clustering
    # That is, a case where a cluster was almost split or not-split
    for cluster in sorted(Cdb['primary_cluster'].unique()):
        d = Cdb[Cdb['primary_cluster'] == cluster]

        if 'Blank' in Ndb:
            continue

        # Skip if it's a singleton
        if len(d['genome'].unique()) == 1:
            continue

        # Load clustering information
        linkI = wd.get_cluster("secondary_linkage_cluster_{0}".format(cluster))
        db = linkI['db']
        linkage = linkI['linkage']
        args = linkI['arguments']
        threshold = args['linkage_cutoff']
        threshold = (1-threshold) * 100
        names = list(db.columns)

        if linkage is None:
            continue

        # Find cases where the cluster is close to the treshold used
        for node in linkage:
            t = node[2]
            x = (1-t)*100
            if abs(x - threshold) < warn_dist:

                # This means that it wasn't split, but was almost
                if x < threshold:
                    #c1 = Cdb['secondary_cluster'][Cdb['genome'] == names[node[0]]].tolist()[0]
                    #c2 = Cdb['secondary_cluster'][Cdb['genome'] == names[node[1]]].tolist()[0]

                    warning = "CLUSTERING WARNING: Primary cluster {0} was almost not split".\
                            format(cluster)
                    warnings.append(warning)

                # This means that a cluster was split, but almost not
                if x > threshold:
                    #print(node)
                    #c1 = Cdb['secondary_cluster'][Cdb['genome'] == names[int(node[0])]].tolist()[0]
                    #c2 = Cdb['secondary_cluster'][Cdb['genome'] == names[int(node[1])]].tolist()[0]
                    #assert c1 == c2

                    warning = "CLUSTERING WARNING: Primary cluster {0} was almost split".\
                            format(cluster)
                    warnings.append(warning)

        # Find cases where you have a high self comparison
        self_thresh = (1-drep.d_analyze.get_highest_self(Ndb, names))*100
        if self_thresh <= threshold + warn_dist:
            warning = "CLUSTERING WARNING: Primary cluster {0} has a high self-comparison value".\
                            format(cluster)
            warnings.append(warning)

    # Check for cases where winners are very similar
    # Either based on MASH or ANIn
    if wd.hasDb('Wmdb'):
        # See if any MASH comparisons are too similar
        Wmdb = Wmdb[(Wmdb['genome1'] != Wmdb['genome2']) & (Wmdb['similarity'] > warn_sim)\
                   & (Wmdb['genome1'] > Wmdb['genome2'])]
        for i,row in Wmdb.iterrows():
            c1 = Cdb['secondary_cluster'][Cdb['genome'] == row['genome1']].tolist()[0]
            c2 = Cdb['secondary_cluster'][Cdb['genome'] == row['genome2']].tolist()[0]

            warning = "WINNER WARNING: Genomes {0} ({3}) and {1} ({4}) have a high MASH score ({2:.2f}%)".format(\
                        row['genome1'], row['genome2'], row['similarity']*100, c1, c2)
            warnings.append(warning)

        # See if any secondary comparisons are too similar
        Wndb = Wndb[(Wndb['reference'] != Wndb['querry']) & (Wndb['ani'] > warn_sim)\
                   & (Wndb['reference'] > Wndb['querry']) & (Wndb['alignment_coverage'] > warn_aln)]
        for i,row in Wndb.iterrows():
            c1 = Cdb['secondary_cluster'][Cdb['genome'] == row['reference']].tolist()[0]
            c2 = Cdb['secondary_cluster'][Cdb['genome'] == row['querry']].tolist()[0]

            warning = "WINNER WARNING: Genomes {0} ({3}) and {1} ({4}) have a high ANIn score ({2:.2f}% ANI ".format(\
                        row['reference'], row['querry'], row['ani']*100, c1, c2) + "; {0:.2f}% aligned)".format(row['alignment_coverage'])
            warnings.append(warning)

    return warnings

def comp_str(val):
    if val == 100:
        return('perfect')
    if val > 90:
        return("near")
    if val > 70:
        return("substantial")
    if val > 50:
        return("moderate")
    if val <= 50:
        return("partial")

def con_str(val):
    if val == 0:
        return ("none")
    if val <= 5:
        return("low")
    if val < 10:
        return ("medium")
    if val < 15:
        return("high")
    if val >= 15:
        return("very high")


def evaluate_winners(wd, **kwrags):
    Wdb = wd.get_db('Wdb')
    Cdb = wd.get_db('Cdb')
    Ndb = wd.get_db('Ndb')

    # For every winning genome, give some key stats based on the information available
    Table = {'genome':[],'score':[],'completeness':[],'contamination':[],'strain_heterogeneity':[],\
            'size':[],'N50':[],'cluster':[], 'taxonomy':[], 'tax_confidence':[],\
            'cluster_members':[],'closest_cluster_member':[],'furthest_cluster_member':[],\
            'completeness_metric':[], 'contamination_metric':[]}

    for i, row in Wdb.iterrows():
        # Add info in Wdb
        Table['genome'].append(row['genome'])
        Table['score'].append(row['score'])
        Table['cluster'].append(row['cluster'])

        # Add clustering info
        if 'Blank' not in Ndb:
            d = Cdb[Cdb['secondary_cluster'] == row['cluster']]
            ndb = Ndb[(Ndb['reference'] == row['genome']) & (Ndb['querry'].isin(d['genome'].tolist()))\
                     & (Ndb['reference'] != Ndb['querry'])]
            members = len(d['genome'].unique())

            if members > 1:
                Table['closest_cluster_member'].append("{0:.2f}".format(ndb['ani'].max()*100))
                Table['furthest_cluster_member'].append("{0:.2f}".format(ndb['ani'].min()*100))
                Table['cluster_members'].append(len(d['genome'].unique()))
            else:
                Table['closest_cluster_member'].append("NA")
                Table['furthest_cluster_member'].append("NA")
                Table['cluster_members'].append(len(d['genome'].unique()))
        else:
            d = Cdb[Cdb['secondary_cluster'] == row['cluster']]
            Table['closest_cluster_member'].append("NA")
            Table['furthest_cluster_member'].append("NA")
            Table['cluster_members'].append(len(d['genome'].unique()))

        # Add checkM info
        if wd.hasDb('Chdb'):
            Chdb = wd.get_db('Chdb')
            d = Chdb[Chdb['Bin Id'] == row['genome']]
            Table['completeness'].append(d['Completeness'].tolist()[0])
            Table['contamination'].append(d['Contamination'].tolist()[0])
            Table['strain_heterogeneity'].append(d['Strain heterogeneity'].tolist()[0])
            Table['size'].append(d['Genome size (bp)'].tolist()[0])
            Table['N50'].append(d['N50 (scaffolds)'].tolist()[0])
            Table['completeness_metric'].append(comp_str(d['Completeness'].tolist()[0]))
            Table['contamination_metric'].append(con_str(d['Contamination'].tolist()[0]))
        else:
            Table['completeness'].append("NA")
            Table['contamination'].append("NA")
            Table['strain_heterogeneity'].append("NA")
            Table['size'].append("NA")
            Table['N50'].append("NA")
            Table['completeness_metric'].append("NA")
            Table['contamination_metric'].append("NA")

        # Add taxonomy info
        if wd.hasDb('Tdb'):
            Tdb = wd.get_db('Tdb')
            d = Tdb[Tdb['genome'] == row['genome']]
            Table['taxonomy'].append(d['taxonomy'][d['tax_confidence'] == d['tax_confidence'].max()].tolist()[0])
            Table['tax_confidence'].append(d['tax_confidence'].max())
        else:
            Table['taxonomy'].append('NA')
            Table['tax_confidence'].append('NA')

    Widb = pd.DataFrame(Table)
    return Widb

def run_tertiary_clustering(wd, **kwargs):
    # Create a new workdir inside data
    new_wd_loc = wd.get_dir('data') + 'tertiary_clustering'
    nWd = drep.WorkDirectory.WorkDirectory(new_wd_loc)

    # Make a copy of the kwargs
    kwargs_copy = kwargs.copy()
    if 'genomes' in kwargs_copy:
        del kwargs_copy['genomes']
    kwargs_copy['P_ani'] = kwargs_copy['P_ani'] - 0.05
    kwargs_copy['multiround_primary_clustering'] = False
    kwargs_copy['greedy_secondary_clustering'] = False

    # Figure out what genomes you're going to compare and make a new Bdb
    Bdb = wd.get_db('Bdb')
    Cdb = wd.get_db('Cdb')
    Wdb = wd.get_db('Wdb')
    BBdb = Bdb[Bdb['genome'].isin(Wdb['genome'].tolist())]
    nWd.store_db(BBdb, 'Bdb')

    # Copy over genome info for choose
    for name in ['genomeInfo', 'Chdb']:
        if wd.hasDb(name):
            nWd.store_db(wd.get_db(name), name)

    # Run the comparison
    drep.d_cluster.controller.d_cluster_wrapper(new_wd_loc, **kwargs_copy)

    # Reconsile results for Cdb
    nWd = drep.WorkDirectory.WorkDirectory(new_wd_loc)
    nCdb = nWd.get_db('Cdb')

    Cdb['original_secondary_cluster'] = Cdb['secondary_cluster']
    del Cdb['secondary_cluster']
    Cdb = pd.merge(Cdb, nCdb[['genome', 'secondary_cluster']], on='genome', how='left')

    rep_cdb = Cdb.dropna()
    old2new = rep_cdb.set_index('original_secondary_cluster')['secondary_cluster'].to_dict()
    Cdb['secondary_cluster'] = Cdb['original_secondary_cluster'].map(old2new)

    # Rename clusters
    old2new_names = {}
    for clust, db in Cdb.groupby('secondary_cluster'):
        if len(db['original_secondary_cluster'].unique()) == 1:
            old2new_names[clust] = db['original_secondary_cluster'].iloc[0]
        else:
            new_name =  f"{'.'.join(sorted(list(set([x.split('_')[0] for x in list(db['original_secondary_cluster'].unique())]))))}" + \
                        f"_{'.'.join(sorted(list(set([x.split('_')[1] for x in list(db['original_secondary_cluster'].unique())]))))}"
            old2new_names[clust] = new_name
            logging.debug(f"Clusters {list(db['original_secondary_cluster'].unique())} were merged into {new_name}")
    Cdb['secondary_cluster'] = Cdb['secondary_cluster'].map(old2new_names)

    # Store new Cdb
    wd.store_db(Cdb, 'Cdb')

    # Re-run choose
    drep.d_choose.d_choose_wrapper(wd.location, **kwargs_copy)
