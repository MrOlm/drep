#!/usr/bin/env python

###############################################################################
#
# test_suite.py - process several E. coli genomes to verify operation of dRep
#
###############################################################################

import tests.test_rerun
import tests.test_unit
import tests.test_dereplicate
import tests.test_filter
import tests.test_cluster
import tests.test_choose
import tests.test_analyze
import tests.test_bonus
import tests.test_taxonomy

if __name__ == '__main__':
    # tests.test_rerun.test_rerun().run()
    # tests.test_unit.test_unit().run()
    tests.test_dereplicate.test_dereplicate().run()
    # tests.test_filter.test_filter().run()
    # tests.test_cluster.test_cluster().run()
    # tests.test_choose.test_choose().run()
    # tests.test_analyze.test_analyze().run()
    # tests.test_bonus.test_bonus().run()
    # tests.test_taxonomy.test_taxonomy().run()

    print("Everything seems to be working swimmingly!")
