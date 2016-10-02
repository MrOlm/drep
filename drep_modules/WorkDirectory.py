#!/usr/bin/env python3

import os
import logging
import pandas as pd

"""
This is the file that maintains the integrity of the work directory

The directory layout
########################################
workDirectory
./data
...../MASH_files/
...../ANIn_files/
./figures
./data_tables
...../Bdb.csv # Sequence locations and filenames
...../Mdb.csv # Raw results of MASH comparisons
...../Ndb.csv # Raw results of ANIn comparisons
...../Cdb.csv # Genomes and cluster designations
./log
...../logger.log
...../cluster_arguments.txt


The data tables
########################################
Bdb = sequence names and sequence locations

"""


class WorkDirectory:

    def __init__(self, location):
        self.location = os.path.abspath(location)
        self.data_tables = {}
        self.overwrite = False
        self.name = None
        
        # Import data_tables
        loc = self.location + '/data_tables/'
        if not os.path.exists(loc):
            os.makedirs(loc)
        self.import_data_tables(loc)
        
    def __str__(self):
        string = "Located: {0}\nDatatables: {1}".format(self.location,list(
                                                        self.data_tables.keys()))
        
        return string
        
    def import_data_tables(self,loc):
        tables = [os.path.join(loc, t) for t in os.listdir(loc) if os.path.isfile(os.path.join(loc, t))]
        for t in tables:
            assert t.endswith('.csv'), "{0} is incorrectly in the data_tables folder".format(t)
            self.data_tables[os.path.basename(t).replace('.csv','')] = pd.read_csv(t)
            
    def hasDb(self,db):
        return (db in self.data_tables) 
        
    def store_db(self,db,name):
        loc = self.location + '/data_tables/'
        
        if os.path.isfile(loc + name + '.csv'):
            assert self.overwrite == True, "data_table {0} already exists".format(name)
            
        db.to_csv(loc + name + '.csv')
        self.data_tables[name] = db