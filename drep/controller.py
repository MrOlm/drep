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
import drep.d_cluster.controller
from drep.WorkDirectory import WorkDirectory
import drep.d_cluster
import drep.d_analyze
import drep.d_filter
import drep.d_choose
import drep.d_adjust
import drep.d_bonus
import drep.d_evaluate
import drep.d_workflows

def version():
    versionFile = open(os.path.join(drep.__path__[0], 'VERSION'))
    return versionFile.read().strip()

VERSION = version()

class Controller():
    def __init__(self):
        self.logger = logging.getLogger()

    def dereplicate_operation(self, **kwargs):
        logging.debug("Starting the dereplicate operation")
        drep.d_workflows.dereplicate_wrapper(kwargs['work_directory'],**kwargs)
        logging.debug("Finished the dereplicate operation!")

    def compare_operation(self, **kwargs):
        logging.debug("Starting the compare operation")
        drep.d_workflows.compare_wrapper(kwargs['work_directory'],**kwargs)
        logging.debug("!!! Finished the compare operation !!!")

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
        logging.debug("Command to run dRep was: {0}\n".format(' '.join(sys.argv)))
        logging.debug("dRep version {0} was run \n".format(VERSION))
        logging.debug("!"*80 + '\n')

    def parseArguments(self, args):
        ''' Parse user options and call the correct pipeline'''
        if args.operation == 'check_dependencies':
            drep.d_bonus.check_dependencies(print_out=True)
            return

        # Load the workDirectory
        wd_loc = str(os.path.abspath(args.work_directory))
        wd = WorkDirectory(wd_loc)

        # Set up the logger
        self.setup_logger(wd.get_loc('log'))
        logging.debug(str(args))

        # Do some testing
        if args.run_tertiary_clustering:
            if args.operation != "dereplicate":
                raise ValueError("Can only run tertiary clustering with dereplicate")


        # Call the appropriate workflow
        if args.operation == "dereplicate":
            self.dereplicate_operation(**vars(args))
        if args.operation == "compare":
            self.compare_operation(**vars(args))

    def loadDefaultArgs(self):
        pass
