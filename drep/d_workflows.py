#!/usr/bin/env python3

import logging
import os
import pandas as pd
import sys

import drep.WorkDirectory
import drep.d_filter
import drep.d_cluster
import drep.d_choose
import drep.d_bonus
import drep.d_evaluate

def dereplicate_wrapper(wd,**kwargs):
    validate_dereplicate(wd, **kwargs)

    message = """\
***************************************************
    ..:: dRep dereplicate Step 1. Filter ::..
***************************************************
    """
    logging.info(message)

    # Pop these arguments they're not there for future operations
    genomes = kwargs.pop('genomes',None)
    Chdb = kwargs.pop('Chdb',None)
    drep.d_filter.d_filter_wrapper(wd, genomes = genomes, Chdb = Chdb, **kwargs)

    message = """\
***************************************************
    ..:: dRep dereplicate Step 2. Cluster ::..
***************************************************
    """
    logging.info(message)
    drep.d_cluster.d_cluster_wrapper(wd, **kwargs)

    message = """\
***************************************************
    ..:: dRep dereplicate Step 3. Choose ::..
***************************************************
    """
    logging.info(message)
    drep.d_choose.d_choose_wrapper(wd, **kwargs)

    message = """\
***************************************************
    ..:: dRep dereplicate Step 4. Bonus ::..
***************************************************
    """
    logging.info(message)
    drep.d_bonus.d_bonus_wrapper(wd, **kwargs)

    message = """\
***************************************************
    ..:: dRep dereplicate Step 5. Evaluate ::..
***************************************************
    """
    logging.info(message)
    drep.d_evaluate.d_evaluate_wrapper(wd, evaluate = '23', **kwargs)

    message = """\
***************************************************
    ..:: dRep dereplicate Step 6. Analyze ::..
***************************************************
    """
    logging.info(message)
    drep.d_analyze.d_analyze_wrapper(wd, plots = 'a', **kwargs)

    loc = drep.WorkDirectory.WorkDirectory(wd).location
    message = """\

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    ..:: dRep dereplicate finished ::..

Dereplicated genomes................. {0}
Dereplicated genomes information..... {2}
Figures.............................. {1}
Warnings............................. {3}

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    """.format(loc + '/dereplicated_genomes/', loc + '/figures/', \
                loc + '/data_tables/Widb.csv', loc + '/log/warnings.txt')
    logging.info(message)

def compare_wrapper(wd,**kwargs):
    validate_compare(wd, **kwargs)

    message = """\
***************************************************
    ..:: dRep compare Step 1. Cluster ::..
***************************************************
    """
    logging.info(message)
    drep.d_cluster.d_cluster_wrapper(wd, **kwargs)

    message = """\
***************************************************
    ..:: dRep compare Step 2. Bonus ::..
***************************************************
    """
    logging.info(message)
    drep.d_bonus.d_bonus_wrapper(wd, **kwargs)

    message = """\
***************************************************
    ..:: dRep compare Step 3. Evaluate ::..
***************************************************
    """
    logging.info(message)
    drep.d_evaluate.d_evaluate_wrapper(wd, evaluate = '2', **kwargs)

    message = """\
***************************************************
    ..:: dRep compare Step 4. Analyze ::..
***************************************************
    """
    logging.info(message)
    drep.d_analyze.d_analyze_wrapper(wd, plots = '1234', **kwargs)

    loc = drep.WorkDirectory.WorkDirectory(wd).location
    message = """\

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    ..:: dRep compare finished ::..

Genome comparison data............... {0}
Figures.............................. {1}
Warnings............................. {3}

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    """.format(loc + '/data_tables/', loc + '/figures/', \
                loc + '/data_tables/Widb.csv', loc + '/log/warnings.txt')
    logging.info(message)

def validate_dereplicate(wd, **kwargs):
    pass

def validate_compare(wd, **kwargs):
    pass
