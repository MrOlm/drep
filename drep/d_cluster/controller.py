import logging
import shutil

import pandas as pd

import drep
import drep.d_cluster.external
import drep.d_cluster.utils
import drep.d_cluster.compare_utils

class GenomeClusterController(object):
    """
    Handle the logic of comparing and clustering genomes
    """
    def __init__(self, wd, **kwargs):
        self.wd = drep.WorkDirectory.WorkDirectory(wd)
        self.kwargs = kwargs

        # Handle special kwargs
        self.debug = kwargs.get('debug', False)

    def main(self, store_output=True):
        """
        Main entrypoint for the dRep cluster operation

        Formerly called d_cluster_wrapper
        """
        # Get the arguments
        self.parse_cluster_arguments()

        # Run primary clustering
        self.run_primary_clustering()

        # Run secondary clustering
        self.run_secondary_clustering()

        # Save the output
        if store_output:
            self.store_output()
        else:
            self.return_output()

    def parse_cluster_arguments(self):
        """
        Load the genomes and store Bdb in the wd
        """
        # Make sure you have the required program installed
        loc = shutil.which('mash')
        if loc is None:
            logging.error('Cannot locate the program {0}- make sure its in the system path' \
                          .format('mash'))

        # If genomes are provided, load them
        if self.kwargs.get('genomes', None) is not None:
            assert self.wd.hasDb("Bdb") == False, \
                "Don't provide new genomes- you already have them in the work directory"
            Bdb = drep.d_cluster.utils.load_genomes(self.kwargs['genomes'])

        # If genomes are not provided, don't load them
        if self.kwargs.get('genomes', None) is None:
            assert self.wd.hasDb("Bdb") != False, \
                "Must either provide a genome list, or run the 'filter' operation with the same work directory"
            Bdb = self.wd.get_db('Bdb')

        # Make sure people weren't dumb with their cutoffs
        for v in ['P_ani', 'S_ani']:
            if self.kwargs.get(v) > 1:
                logging.warning("{0} is set to {1}- this should be \
                            between 0-1, not 1-100".format(v, self.kwargs.get(v)))

        # Load length and N50 if you need it
        if self.kwargs.get('multiround_primary_clustering', False) | (self.kwargs.get('greedy_secondary_clustering', False)):
            Bdb = drep.d_filter._add_lengthN50(Bdb, Bdb)

        # Store the genomes
        self.Bdb = Bdb

    def run_primary_clustering(self):
        """
        Run primary clustering and end with an Mdb and an MCdb in this object
        """
        logging.info("Running primary clustering")
        cached = (self.debug and self.wd.hasDb('Mdb') and self.wd.hasDb('CdbF'))

        if self.kwargs.get('SkipMash', False):
            logging.info("Nevermind! Skipping Mash")
            # Make a "Cdb" where all genomes are in the same cluster
            Cdb = drep.d_cluster.external._gen_nomash_cdb(self.Bdb)
            # Make a blank "Mdb" for storage anyways
            Mdb = pd.DataFrame({'Blank': []})

        elif cached:
            logging.info('Nevermind! Loading cached primary clustering')
            Mdb = self.wd.get_db('Mdb')
            Cdb = self.wd.get_db('CdbF')
            logging.info('2. Primary clustering cache loaded')

        else:
            logging.info("Running pair-wise MASH clustering")
            Mdb, Cdb, cluster_ret = drep.d_cluster.compare_utils.all_vs_all_MASH(self.Bdb, self.wd.get_dir('MASH'), **self.kwargs)

            if self.debug:
                logging.debug("Debug mode on - saving Mdb ASAP")
                self.wd.store_db(Mdb, 'Mdb')

                logging.debug("Debug mode on - saving CdbF ASAP")
                self.wd.store_db(Cdb, 'CdbF')

            # Store the primary clustering results
            self.wd.store_special('primary_linkage', cluster_ret)

        logging.info("{0} primary clusters made".format(len(Cdb['primary_cluster'].unique())))
        self.Mdb = Mdb
        self.MCdb = Cdb

    def run_secondary_clustering(self):
        logging.info("Running secondary clustering")

        # Get the arguments
        algorithm = self.kwargs.get('S_algorithm', 'ANImf')
        cached = (self.debug and self.wd.hasDb('Ndb') and self.wd.hasDb('Cdb'))
        p = self.kwargs.get('processors', 6)
        self.deal_with_nucmer_presets()

        # Wipe any old secondary clusters
        self.wd._wipe_secondary_clusters()

        if not self.kwargs.get('SkipSecondary', False):
            if cached:
                logging.info('3. Loading cached secondary clustering')
                Ndb = self.wd.get_db('Ndb')
                Cdb = self.wd.get_db('Cdb')

                # Get rid of broken ones
                base = len(Ndb)
                Ndb = Ndb.dropna(subset=['reference'])
                logging.info(f'!!! {len(Ndb) - base} lines from Ndb failed!')

                logging.info('3. Secondary clustering cache loaded')

            # Run comparisons, make Ndb
            else:
                drep.d_cluster.utils._print_time_estimate(self.Bdb, self.MCdb, algorithm, p)
                Ndb, Cdb, c2ret = drep.d_cluster.compare_utils.secondary_clustering(self.Bdb, self.MCdb, algorithm, self.wd.get_dir('data'), wd=self.wd, **self.kwargs)
                if self.debug:
                    logging.debug("Debug mode on - saving Ndb ASAP")
                    self.wd.store_db(Ndb, 'Ndb')
                    self.wd.store_db(Cdb, 'Cdb')

                # Store the secondary clustering results
                self.wd.store_special('secondary_linkages', c2ret)

        else:
            logging.info("3. Nevermind! Skipping secondary clustering")
            Cdb = drep.d_cluster.utils._gen_nomani_cdb(self.MCdb, data_folder=self.wd.get_dir('data'), **self.kwargs)
            Ndb = pd.DataFrame({'Blank': []})

        logging.info(
            "Step 4. Return output")

        self.Cdb = Cdb
        self.Ndb = Ndb

    def store_output(self):
        logging.debug("Main program run complete- saving output")
        self.wd.store_db(self.Cdb, 'Cdb')
        self.wd.store_db(self.Mdb, 'Mdb')
        self.wd.store_db(self.Ndb, 'Ndb')
        if not self.wd.hasDb('Bdb'):
            self.wd.store_db(self.Bdb, 'Bdb')

        # Log arguments
        self.wd.store_special('cluster_log', self.kwargs)

    def return_output(self):
        return self.Cdb, self.Mdb, self.Ndb

    def deal_with_nucmer_presets(self):
        if self.kwargs.get('n_preset', None) != None:
            self.kwargs['n_c'], self.kwargs['n_maxgap'], self.kwargs['n_noextend'], self.kwargs['n_method'] \
                = drep.d_cluster.external._nucmer_preset(self.kwargs['n_PRESET'])

def d_cluster_wrapper(workDirectory, **kwargs):
    GenomeClusterController(workDirectory, **kwargs).main()


