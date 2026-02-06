#!/usr/bin/env python

#import dependencies
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


#CLASSES

#Creating class for genes. this will help keep track of entries in the fasta file and then will be able to relate exons and motifs as well.
#the length defines the total line length for the final drawing

class Gene:
    """For each entry into the fasta file"""
    # A class variable, shared by all instances
    #species = "canine"

    # The constructor method to initialize new objects
    def __init__(self, header, sequence):
        """Takes sequence and header to initialize. creates length from sequence"""
        self.header = header  # Instance variable unique to each instance
        self.sequence = sequence # Instance variable unique to each instance
        self.length = len(sequence)


class Exon:
    """For each Exon in a Gene""" #need starting location in gene sequence and length to draw the rectangle on the line at the right spot.

    # The constructor method to initialize new objects
    def __init__(self, header, sequence):
        """Takes sequence and header to initialize. creates length from sequence"""
        #self.header = header  # Instance variable unique to each instance
        #self.sequence = sequence # Instance variable unique to each instance
        self.length = len(sequence)


class Motif:
    """For each Motif in a Gene"""

    # The constructor method to initialize new objects
    def __init__(self, header, sequence):
        """Takes sequence and header to initialize. creates length from sequence"""
        self.header = header  # Instance variable unique to each instance
        self.sequence = sequence # Instance variable unique to each instance
        self.length = len(sequence)




#FUNCTIONS (WE OUT HERE TRYNA)

#find genes in FASTA file
def grab_genes(gene):
    """Goes through FASTA file and grabs all gene entries with headers and sequences"""

#find exons in each gene. takes the input gene from the saved gene classes in the full gene list for the fasta file.
def find_exons(gene):
    """Goes through the gene sequence and pulls out exon locations and lengths so they can be drawn"""

#find motifs in each gene. takes the input gene from the saved gene classes in the full gene list for the fasta file.
def find_motifs(gene):
    """Goes through gene sequence and identifies motifs with locations so they can be drawn"""

#drawing function
def draw_annotated_gene(gene, exons, motifs):
    """This takes genes and all necessary elements and draws the pictures for each"""
    print('fuckme')



#getting motifs from the motifs file. reading it in and making a set with all the motifs. one from each line of the file.
#making set to hold the motifs
motifs = set()
#reading through the file and saving each line to the set
with open(args.motifs_file, "r") as open_motifs:
    for line in open_motifs:
        line = line.strip()
        motifs.add(line)


#reading through fasta file. pulling out each gene as its own string with corresponding header and make a class for each gene with header and sequence
#genes are saved to list which keeps track of all genes for the file and will help draw the final output

#initializng list of genes to track
gene_list = []


#initializing these variables which will keep track of the current header and gene for each entry in the fasta file
current_header = 'no'
current_seq = ''

#going through fasta file to pull out genes with the headers to make classes for each gene.
with open(args.fasta_file, "r") as open_fasta:
    for line in open_fasta:
        line = line.strip()
        
        #for each time a header line is read in the file, resets the current header and sequence while adding them to the list.
        if line[0] == '>':

            #the if statement skips over the initial part of reading through the file and avoid adding an empty entry to the list.
            if current_header != 'no':
                #adds the current header with its sequence for each gene to the list as a gene class.
                gene_list.append(Gene(current_header, current_seq))

            #resets the current header when a new one is reached
            current_header = line
            
            #resets the current sequence as well
            current_seq = ''

            #skips back so the sequence lines can be read.
            continue

        #for every line that is not a header, the sequence from the line is added as a string for each gene.
        current_seq += line


print(gene_list)       
print(gene_list[0].length)


#draw all genes
for x in gene_list:
    draw_annotated_gene(x, x.exons, x.motifs)