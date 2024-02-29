#!/usr/bin/env python

import os
import sys
import gzip
import textwrap
import argparse

def extract_bins(fasta, stb_file, out_base):
    if (out_base != '') & (out_base[-1] != '/'):
        out_base += '_'

    # Make scaffold to bin dictionary
    stb = {}
    with open(stb_file) as stb_reader:
        for line in stb_reader:
            line = line.strip()

            if line.startswith('#') or line.startswith('scaffold_name'):
                continue

            scaffold, bin = line.split('\t')[:2]
            stb[scaffold.strip()] = bin.strip()

    # Parse the FASTA file and write to bins
    if fasta.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'

    with open_func(fasta, mode) as handle:
        record_id = ''
        record_seq = ''
        for line in handle:
            if line.startswith('>'):
                if record_id and record_id in stb:
                    bin = stb[record_id]
                    with open(f"{out_base}{bin}.fa", 'a') as nw:
                        nw.write(f">{record_id}\n{record_seq}\n")
                record_id = line[1:].strip().split(' ')[0]
                record_seq = ''
            else:
                record_seq += line.strip()

        # Write the last record
        if record_id and record_id in stb:
            bin = stb[record_id]
            with open(f"{out_base}{bin}.fa", 'a') as nw:
                nw.write(f">{record_id}\n{record_seq}\n")

import os
import gzip

def gen_stb(fastas):
    if ((len(fastas) == 1) & (not fastas[0].endswith('.gz'))):
        # See if this is a text file, not a fasta file
        text_list = True
        genomes = []
        with open(fastas[0], 'r') as o:
            for line in o.readlines():
                if line.startswith('>'):
                    text_list = False
                    break
                else:
                    genomes.append(line.strip())
        if text_list:
            print("Treating .fasta input as list")
            fastas = genomes

    stb = {}
    for fasta in fastas:
        bin = os.path.basename(fasta)
        if fasta.endswith('.gz'):
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'

        with open_func(fasta, mode) as handle:
            for line in handle:
                if line.startswith('>'):
                    id = line[1:].strip().split(' ')[0]
                    stb[id] = bin

    return stb


def print_stb(stb, output):
    file = open(output,'w')
    for key in sorted(stb, key=stb.get):
        file.write("{0}\t{1}\n".format(key,stb[key]))
    file.close()


class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=textwrap.dedent('''\
         \n
         The program has two uses related to scaffold to bin (.stb) files.
         .stb files should be tab-separated, with no header, and two columns: scaffold and bin

         Use 1) Pass a list of genomes to generate a .stb file.

         Example:
         parse_stb.py --reverse -f dereplicate_genomes/* -o representitve_genomes.stb

         Use 2) Pass a single .fasta file and a scaffold to bin file (.stb) to generate a number of
         fasta files based on the .stb file.

         Example:
         parse_stb.py -f concat_genomes.fasta -s scaffold_to_bin.tsv -o genomeList_1
         '''))

    parser.add_argument('-s','--stb',help='scaffold to bin file')
    parser.add_argument('-f','--fasta',help='fasta file to extract scaffolds from. Will treat as compressed if ends in .gz. This can also be a single text file with a genome on each line',nargs='*')
    parser.add_argument('-o','--output',help='output base name', default = '')

    parser.add_argument('--reverse',help='generate a stb from a list of genomes',\
                        action = "store_true")
    args = parser.parse_args()

    if args.reverse == False:
        if len(args.fasta) != 1:
            print("must give one and only one fasta file")
            sys.exit()
        extract_bins(args.fasta[0], args.stb, args.output)

    if args.reverse == True:
        stb = gen_stb(args.fasta)
        print_stb(stb, args.output)
