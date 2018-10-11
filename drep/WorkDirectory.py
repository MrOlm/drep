#!/usr/bin/env python3
"""
This module provides access to the workDirectory

The directory layout::

    workDirectory
    ./data
    ...../MASH_files/
    ...../ANIn_files/
    ...../gANI_files/
    ...../Clustering_files/
    ...../checkM/
    ........./genomes/
    ........./checkM_outdir/
    ...../prodigal/
    ./figures
    ./data_tables
    ...../Bdb.csv  # Sequence locations and filenames
    ...../Mdb.csv  # Raw results of MASH comparisons
    ...../Ndb.csv  # Raw results of ANIn comparisons
    ...../Cdb.csv  # Genomes and cluster designations
    ...../Chdb.csv # CheckM results for Bdb
    ...../Sdb.csv  # Scoring information
    ...../Wdb.csv  # Winning genomes
    ./dereplicated_genomes
    ./log
    ...../logger.log
    ...../cluster_arguments.json

"""

import os
import logging
import pandas as pd
import pickle
import json
import sys
import shutil
import glob
import numpy as np

import drep

class WorkDirectory(object):
    '''
    Object to interact with the workDirectory

    Args:
        location (str): location to make the workDirectory


    '''
    firstLevels = ['data','figures','data_tables','dereplicated_genomes','log']

    def __init__(self, location):
        self.location = os.path.abspath(location)
        self.data_tables = {}
        self.clusters = {}
        self.arguments = {}
        self.overwrite = True
        self.name = None

        self.make_fileStructure()
        self.load_cached()

    def __str__(self):
        string = "Located: {0}\nDatatables: {1}\nCluster files: {2}\nArguments: {3}".format(\
                    self.location,list(self.data_tables.keys()),list(self.clusters.keys()),\
                    list(self.arguments.keys()))

        return string

    def make_fileStructure(self):
        '''
        Make the top level file structure
        '''
        location = self.location

        if not os.path.exists(location):
            os.makedirs(location)

        for l in WorkDirectory.firstLevels:
            loc = location + '/' + l
            if not os.path.exists(loc):
                os.makedirs(loc)

    def load_cached(self):
        '''
        The wrapper to load everything it has into attributes
        '''
        # Import data_tables
        loc = self.get_dir('data_tables')
        self.import_data_tables(loc)

        # Import pickles
        loc = self.location + '/data/Clustering_files/'
        if not os.path.exists(loc):
            os.makedirs(loc)
        self.import_clusters(loc)

        # Import arguments
        loc = self.location + '/log/'
        if not os.path.exists(loc):
            os.makedirs(loc)
        self.import_arguments(loc)

    def import_data_tables(self, loc):
        '''
        Given the location of the datatables, load them
        '''
        tables = [os.path.join(loc, t) for t in os.listdir(loc) if \
                os.path.isfile(os.path.join(loc, t))]
        for t in tables:
            self.data_tables[os.path.basename(t).replace('.csv','').\
                replace('.pickle', '')] = t

    def import_clusters(self,loc):
        '''
        Given the location of the cluster files, load them
        '''
        pickles = [os.path.join(loc, t) for t in os.listdir(loc) if os.path.isfile(os.path.join(loc, t))]
        for p in pickles:
            assert p.endswith('.pickle'), "{0} is incorrectly in the data/Clustering_files folder".format(t)

            f = open(p,'rb')
            try:
                linkage = pickle.load(f)
                db = pickle.load(f)
                args = pickle.load(f)
                name = os.path.basename(p).replace('.pickle','')
                self.clusters[name] = {'linkage':linkage,'db':db,'arguments':args}
            except:
                print("{0} is an imporperly made pickle- skipping import".format(p))

    def import_arguments(self,loc):
        '''
        Given the location of the log directory, load it
        '''
        args = [os.path.join(loc, t) for t in os.listdir(loc) if os.path.isfile(os.path.join(loc, t))]
        for j in args:
            if j.endswith('_arguments.json'):
                with open(j) as data_file:
                    data = json.load(data_file)
                self.arguments[os.path.basename(j).replace('_arguments.json','')] = data

    def hasDb(self,db):
        '''
        If db is in the data_tables, return True
        '''
        return (db in self.data_tables)

    def get_cluster(self, name):
        '''
        Get the cluster passed in

        Args:
            name: name of the cluster

        Returns:
            cluster
        '''
        # If the whole cluster name is passed in
        if name in self.clusters:
            return self.clusters[name]

        else:
            n = 'secondary_linkage_cluster_{0}'.format(name)
            if n in self.clusters:
                return self.clusters[n]

        print("Could not find cluster {0} - could it be a singleton?".format(name))
        print("The clusters I have are: {0}".format(list(self.clusters.keys())))
        print("Quitting")
        sys.exit()

    def get_primary_linkage(self):
        '''
        Get the primary linkage cluster
        '''
        return self.clusters['primary_linkage']

    def store_db(self,db,name,overwrite=None):
        '''
        Store a dataframe in the workDirectory

        Will make a physical copy in the datatables folder

        Args:
            db: pandas dataframe to store
            name: name to store it under (will add .csv automatically)
            overwrite: if True, overwrite if DataFrame with same name already exists
        '''
        loc = self.get_dir('data_tables')

        if overwrite == None:
            overwrite = self.overwrite

        if os.path.isfile(loc + name + '.csv'):
            assert overwrite == True, "data_table {0} already exists".format(name)

        if name == 'Mdb':
            floc = loc + name + '.csv'
            db.to_csv(floc, index=False)
            self.data_tables[name] = floc

        else:
            floc = loc + name + '.csv'
            db.to_csv(floc, index=False)
            self.data_tables[name] = floc

    def get_db(self, name, return_none=True, forPlotting=False):
        '''
        Get database from self.data_tables

        Args:
            name: name of dataframe
            return_none: if True will return None if database not found; otherwise assert False
            forPlotting: if True don't do fancy dType loading; it messes with order of names for dendrograms
        '''
        if name in self.data_tables:
            if name == 'Mdb':
                if forPlotting:
                    return pd.read_csv(self.data_tables[name])
                else:
                    dTypes={'genome1':'category', 'genome2':'category', 'dist':np.float32,\
                        'similarity':np.float32}
                    return pd.read_csv(self.data_tables[name], dtype=dTypes)

                print('hioh')

            else:
                return pd.read_csv(self.data_tables[name])
        else:
            if return_none:
                return None
            else:
                assert False, "Datatable {0} is not in the work directory {1}".format(\
                                name, self.location)

    def get_dir(self, dir):
        '''
        Get the location of one of the named directory types

        Args:
            dir: Name of directory to find

        Returns:
            string: Location of requested directory
        '''
        d = None
        if dir == 'data_tables':
            d = self.location + '/data_tables/'
        if dir == 'prodigal':
            d = self.location + '/data/prodigal/'
        elif dir == 'centrifuge':
            d = self.location + '/data/centrifuge/'
        elif dir == 'ESOM':
            d = self.location + '/data/ESOM/'
        elif dir == 'log':
            d = self.location + '/log/'
        elif dir == 'cmd_logs':
            d = self.location + '/log/cmd_logs/'
        elif dir == 'MASH':
            d = self.location + '/data/MASH_files/'
        elif dir == 'checkM':
            d = self.location + '/data/checkM/'
        elif dir == 'data':
            d = self.location + '/data/'
        elif dir == 'clustering':
            d = self.location + '/data/Clustering_files/'
        elif dir == 'dereplicated_genomes':
            d = self.location + '/dereplicated_genomes/'
        elif dir == 'figures':
            d = self.location + '/figures/'

        if d == None:
            assert False, "{0} is not a directory I know about".format(dir)

        if not os.path.exists(d):
            os.makedirs(d)

        return d

    def get_loc(self, what):
        '''
        Get the location of Things

        Args:
            what: string of what to get the location of

        Returns:
            string: location of what
        '''
        if what == 'log':
            return self.location + '/log/logger.log'

    def store_special(self, name, thing):
        '''
        Store special items in the work directory

        Args:
            name: what to store
            thing: actual thing to store
        '''
        if name == 'primary_linkage':
            assert len(thing) == 3
            store_loc = self.get_dir('clustering')
            logging.debug('Saving primary_linkage pickle to {0}'.format(store_loc))
            with open(store_loc + 'primary_linkage.pickle', 'wb') as handle:
                pickle.dump(thing[0], handle, protocol=4)
                pickle.dump(thing[1], handle, protocol=4)
                pickle.dump(thing[2], handle, protocol=4)

        elif name == 'secondary_linkages':
            assert type(thing) == dict
            store_loc = self.get_dir('clustering')
            for name, cluster_ret in thing.items():
                #print(name)
                if len(cluster_ret) == 0:
                    continue
                assert len(cluster_ret) == 3

                pickle_name = "secondary_linkage_cluster_{0}.pickle".format(name)
                logging.debug('Saving secondary_linkage pickle {1} to {0}'.format(pickle_name,\
                                                                    store_loc))
                with open(store_loc + pickle_name, 'wb') as handle:
                    pickle.dump(cluster_ret[0],handle, protocol=4)
                    pickle.dump(cluster_ret[1],handle, protocol=4)
                    pickle.dump(cluster_ret[2],handle, protocol=4)

        elif name == 'dereplicated_genomes':
            output_folder = self.get_dir('dereplicated_genomes')
            if os.path.exists(output_folder):
                #logging.debug("{0} already exists: removing and remaking".format(output_folder))
                shutil.rmtree(output_folder)
            os.makedirs(output_folder)

            for loc in thing:
                genome = os.path.basename(loc)
                shutil.copy2(loc, "{0}{1}".format(output_folder,genome))
            #logging.info("Done! Dereplicated genomes saved at {0}".format(output_folder))

        elif name == 'cluster_log':
            cluster_log = os.path.join(self.get_dir('log') + 'cluster_arguments.json')
            with open(cluster_log, 'w') as fp:
                json.dump(thing, fp)
            fp.close()

    def _wipe_secondary_clusters(self):
        '''
        Wipe any and all secondary clusters present in the workDirectory
        '''
        # Clear out clustering folder
        c_folder = self.get_dir('clustering')
        for fn in glob.glob(c_folder + 'secondary_linkage_cluster*'):
            os.remove(fn)
