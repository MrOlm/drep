#!/usr/bin/env python3

import argparse
from argparse import RawTextHelpFormatter


class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()  
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(formatter_class=SmartFormatter)
	subparsers = parser.add_subparsers(help='Desired operation',dest='operation')
	
	# Arguments for all operations
	parser.add_argument("working_directory",help="R|Directory where data and output is stored\
	\n*** Use the same working directory for all commands ***")
	
	# Arguments for clustering operations
	cluster_parser = subparsers.add_parser("cluster")
	cluster_parser.add_argument("-ma","--MASH_ani",help="Threshold to form MASH clusters",\
	    default=0.9)
	
	# Arguments for filtering operations
	filter_parser = subparsers.add_parser("filter")
	filter_parser.add_argument("-l",'--length',help="Minimum length")
	
	args = parser.parse_args()
	
	print("You chose {0}".format(args.operation))