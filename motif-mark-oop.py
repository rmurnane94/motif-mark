#!/usr/bin/env python

import argparse


#argparse
def get_args():
    parser = argparse.ArgumentParser(description="to find motifs")
    parser.add_argument("-f", "--fasta_file", help="fasta file to find motifs", type=str, required=True)
    parser.add_argument("-m", "--motifs_file", help="file of motif sequences to identify in fasta", type=str, required=True)
    #parser.add_argument("-o", "--output_file", help="designates file to deduplicated sam file", type=str, required=True)
    #parser.add_argument("-u", "--umi_file", help="designates file containing the list of UMIs", type=str, required=True)
    return parser.parse_args()

args = get_args()

#getting motifs from the motifs file. reading it in and making a set with all the motifs. one from each line of the file.
#setting variable to get the file
motifs_file = args.motifs_file

#making set to hold the motifs
motifs = set()

#reading through the file and saving each line to the set
with open(motifs_file, "r") as open_motifs:
    for line in open_motifs:
        line = line.strip()
        motifs.add(line)



