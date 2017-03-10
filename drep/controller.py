#!/usr/bin/env python3

###############################################################################
#
# dRep - main program entry point
#
###############################################################################

'''
Controller- takes input from argparse and calls correct modules
'''


__author__ = "Matt Olm"
__license__ = "MIT"
__email__ = "mattolm@gmail.com"
__status__ = "Development"

import argparse
import logging
import os
import sys

import drep
from drep.WorkDirectory import WorkDirectory
import drep.d_cluster
import drep.d_analyze
import drep.d_filter
import drep.d_choose
import drep.d_adjust
import drep.d_bonus
import drep.d_evaluate
import drep.d_workflows

class Controller():
    def __init__(self):
        self.logger = logging.getLogger()

    def filter_operation(self, **kwargs):
        logging.debug("Starting the filter operation")
        drep.d_filter.d_filter_wrapper(kwargs['work_directory'],**kwargs)
        logging.debug("!!! Finished the filter operation !!!")

    def cluster_operation(self, **kwargs):
        if (kwargs['P_ani'] > 1) or (kwargs['S_ani'] > 1):
            logging.error("Can't assign a MASH or ANIn value over 1")
            sys.exit()

        logging.debug("Starting the clustering operation")
        drep.d_cluster.d_cluster_wrapper(kwargs['work_directory'],**kwargs)
        logging.debug("!!! Finished the clustering operation !!!")

    def analyze_operation(self, **kwargs):
        logging.debug("Starting the analyze operation")
        drep.d_analyze.d_analyze_wrapper(kwargs['work_directory'],**kwargs)
        logging.debug("!!! Finished the analyze operation !!!")

    def choose_operation(self, **kwargs):
        logging.debug("Starting the choose operation")
        drep.d_choose.d_choose_wrapper(kwargs['work_directory'],**kwargs)
        logging.debug("!!! Finished the choose operation !!!")

    def adjust_operation(self, **kwargs):
        logging.debug("Starting the adjust operation")
        drep.d_adjust.d_adjust_wrapper(kwargs['work_directory'],**kwargs)
        logging.debug("!!! Finished the adjust operation !!!")

    def bonus_operation(self, **kwargs):
        logging.debug("Starting the bonus operation")
        drep.d_bonus.d_bonus_wrapper(kwargs['work_directory'],**kwargs)
        logging.debug("!!! Finished the bonus operation !!!")

    def evaluate_operation(self, **kwargs):
        logging.debug("Starting the evaluate operation")
        drep.d_evaluate.d_evaluate_wrapper(kwargs['work_directory'],**kwargs)
        logging.debug("!!! Finished the evaluate operation !!!")

    def dereplicate_wf_operation(self, **kwargs):
        logging.debug("Starting the dereplicate_wf operation")
        drep.d_workflows.dereplicate_wrapper(kwargs['work_directory'],**kwargs)
        logging.debug("Finished the dereplicate_wf operation!")

    def compare_wf_operation(self, **kwargs):
        logging.debug("Starting the compare_wf operation")
        drep.d_workflows.compare_wrapper(kwargs['work_directory'],**kwargs)
        logging.debug("!!! Finished the compare_wf operation !!!")

        '''
    def makeload_logger(wd):
        wd = str(os.path.abspath(wd))
        if not os.path.exists(wd):
            os.makedirs(wd)

        log_dir = wd + '/log/'
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)

        logging.basicConfig(filename=log_dir + 'logger.log',level=logging.DEBUG,\
            format='%(asctime)s %(message)s')
        logging.info("***Logger started up at {0}***".format(log_dir + 'logger.log'))
        '''

    def setup_logger(self,loc):
        ''' set up logger such that DEBUG goes only to file, rest go to file and console '''

        # set up logging everything to file
        logging.basicConfig(level=logging.DEBUG,
                           format='%(asctime)s %(levelname)-8s %(message)s',
                           datefmt='%m-%d %H:%M',
                           filename=loc)

        # set up logging of INFO or higher to sys.stderr
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        formatter = logging.Formatter('%(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)

        logging.debug("!"*80)
        logging.debug("***Logger started up at {0}***".format(loc))
        logging.debug("Command to run dRep was: {0}\n".format(sys.argv))
        logging.debug("!"*80 + '\n')

    def parseArguments(self, args):
        ''' Parse user options and call the correct pipeline'''

        # Load the workDirectory
        wd_loc = str(os.path.abspath(args.work_directory))
        wd = WorkDirectory(wd_loc)

        # Set up the logger
        self.setup_logger(wd.get_loc('log'))
        logging.debug(str(args))

        # Call the appropriate workflow
        if args.operation == "dereplicate_wf":
            self.dereplicate_wf_operation(**vars(args))
        if args.operation == "compare_wf":
            self.compare_wf_operation(**vars(args))

        if args.operation == "filter":
            self.filter_operation(**vars(args))
        if args.operation == "cluster":
            self.cluster_operation(**vars(args))
        if args.operation == "analyze":
            self.analyze_operation(**vars(args))
        if args.operation == "choose":
            self.choose_operation(**vars(args))
        if args.operation == "adjust":
            self.adjust_operation(**vars(args))
        if args.operation == "bonus":
            self.bonus_operation(**vars(args))
        if args.operation == "evaluate":
            self.evaluate_operation(**vars(args))

    def loadDefaultArgs(self):
        pass
