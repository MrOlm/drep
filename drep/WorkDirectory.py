#!/usr/bin/env python3

import os
import logging
import pandas as pd
import pickle
import json
import sys

"""
This is the file that maintains the integrity of the work directory

The directory layout
########################################
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


The data tables
########################################
Bdb = sequence names and sequence locations

"""

class WorkDirectory(object):

    firstLevels = ['data','figures','data_tables','dereplicated_genomes','log']

    '''
    # Make this is Singleton object
    def __new__(cls, loc):
        if not hasattr(cls,'instance'):
            cls.instance = super(WorkDirectory, cls).__new__(cls)
            cls.instance.__initialized = False
        return cls.instance
    '''

    def __init__(self, location):
        #if(self.__initialized): return
        #self.__initialized = True

        self.location = os.path.abspath(location)
        self.data_tables = {}
        self.clusters = {}
        self.arguments = {}
        self.overwrite = False
        self.name = None

        self.make_fileStructure()
        self.load_cached()

    def __str__(self):
        string = "Located: {0}\nDatatables: {1}\nCluster files: {2}\nArguments: {3}".format(\
                    self.location,list(self.data_tables.keys()),list(self.clusters.keys()),\
                    list(self.arguments.keys()))

        return string

    def make_fileStructure(self):
        location = self.location

        if not os.path.exists(location):
            os.makedirs(location)

        for l in WorkDirectory.firstLevels:
            loc = location + '/' + l
            if not os.path.exists(loc):
                os.makedirs(loc)

    def load_cached(self):
        # Import data_tables
        loc = self.location + '/data_tables/'
        if not os.path.exists(loc):
            os.makedirs(loc)
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

    def import_data_tables(self,loc):
        tables = [os.path.join(loc, t) for t in os.listdir(loc) if os.path.isfile(os.path.join(loc, t))]
        for t in tables:
            assert t.endswith('.csv'), "{0} is incorrectly in the data_tables folder".format(t)
            self.data_tables[os.path.basename(t).replace('.csv','')] = pd.read_csv(t)

    '''
    def import_pickles(self,loc):
        pickles = [os.path.join(loc, t) for t in os.listdir(loc) if os.path.isfile(os.path.join(loc, t))]
        for p in pickles:
            assert p.endswith('.pickle'), "{0} is incorrectly in the data/Clustering_files folder".format(t)
            self.pickles[os.path.basename(p).replace('.pickle','')] = pickle.load(open(p,'rb'))
    '''
    def import_clusters(self,loc):
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
        args = [os.path.join(loc, t) for t in os.listdir(loc) if os.path.isfile(os.path.join(loc, t))]
        for j in args:
            if j.endswith('_arguments.json'):
                with open(j) as data_file:
                    data = json.load(data_file)
                self.arguments[os.path.basename(j).replace('_arguments.json','')] = data

    def hasDb(self,db):
        return (db in self.data_tables)

    def get_cluster(self, name):

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
    '''
    def get_ANIn_linkages(self):
        to_return = {}
        for pickle in self.pickles:
            if pickle.startswith('ANIn_linkage_cluster_'):
                to_return[pickle.replace('ANIn_linkage_cluster_','')] = self.pickles[pickle]
        return to_return

    def get_gANI_linkages(self):
        to_return = {}
        for pickle in self.pickles:
            if pickle.startswith('gANI_linkage_cluster_'):
                to_return[pickle.replace('gANI_linkage_cluster_','')] = self.pickles[pickle]
        return to_return
    '''
    def get_primary_linkage(self):
        return self.clusters['primary_linkage']

    def store_db(self,db,name,overwrite=False):
        loc = self.location + '/data_tables/'

        if os.path.isfile(loc + name + '.csv'):
            assert overwrite == True, "data_table {0} already exists".format(name)

        db.to_csv(loc + name + '.csv', index=False)
        self.data_tables[name] = db

    def get_db(self, name, return_none=True):
        if name in self.data_tables:
            return self.data_tables[name]
        else:
            if return_none:
                return None
            else:
                assert False, "Datatable {0} is not in the work directory {1}".format(\
                                name, self.location)

    def get_dir(self, dir):
        d = None
        if dir == 'prodigal':
            d = self.location + '/data/prodigal/'
        if dir == 'centrifuge':
            d = self.location + '/data/centrifuge/'
        if dir == 'ESOM':
            d = self.location + '/data/ESOM/'

        if not os.path.exists(d):
            os.makedirs(d)

        return d

    def get_loc(self, what):
        if what == 'log':
            return self.location + '/log/logger.log'
