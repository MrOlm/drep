#!/usr/bin/env python

import os
import sys
import gzip
import logging
import textwrap
import argparse

from Bio import SeqIO

def extract_bins(fasta, stb_file, out_base):
    if (out_base != '') & (out_base[-1] != '/'):
        out_base += '_'

    # Make scaffold to bin dictionary
    stb = {}
    stb_reader = open(stb_file)
    for line in stb_reader:
        line = line.strip()

        if line.startswith('#') or line.startswith('scaffold_name'):
            continue

        scaffold = line.split('\t')[0].strip()
        bin = line.split('\t')[1].strip()
        stb[scaffold] = bin

    # Do it with append
    if fasta[-3:] == '.gz':
        handle = gzip.open(fasta, "rt")
        seqs = SeqIO.parse(handle, "fasta")

    else:
        seqs = SeqIO.parse(fasta, "fasta")

    for seq_record in seqs:
        id = str(seq_record.id).strip()
        if id not in stb:
            # print("{0} not in stb".format(id))
            continue
        fasta = stb[id]

        nw = open("{0}{1}.fa".format(out_base, fasta), 'a')
        nw.write('\n'.join([">{0}".format(id), str(seq_record.seq), '']))
        nw.close()

    if fasta[-3:] == '.gz':
        handle.close()

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
        if bin[-3:] == '.gz':
            with gzip.open(fasta, "rt") as handle:
                for seq_record in SeqIO.parse(handle, "fasta"):
                    id = str(seq_record.id).strip()
                    stb[id] = bin
        else:
            iter = SeqIO.parse(fasta, "fasta")
            for seq_record in iter:
                id = str(seq_record.id).strip()
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
