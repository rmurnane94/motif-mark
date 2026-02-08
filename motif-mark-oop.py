#!/usr/bin/env python

#import dependencies
import argparse
import re


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

    # The constructor method to initialize new Genes
    def __init__(self, header, sequence):
        """Takes sequence and header to initialize."""
        self.header = header  # stores gene header
        self.sequence = sequence # stores gene sequence
        self.length = len(sequence) # gets length of gene for drawing
        self.exons = [] #creates empty list to store exons
        self.motifs = [] #creates empty list to store motifs

    def find_exons(self): #adds exons to the gene
        """Goes through the gene sequence and pulls out exon locations and lengths so they can be drawn"""
        #setting pattern to find all instances of one or more capital letters. sequential capital letters are exons in this case.
        pattern = r'[A-Z]+'
        # Use re.finditer() to locate capitalized exon strings with locations
        matches = re.finditer(pattern, self.sequence)
        #go through the located exons to get information and make classes
        for match in matches:
            # Get the exon string
            exon_seq = match.group()
            # Get the start location in the gene
            exon_start = match.start()
            #store exon information for the gene as exon classes
            self.exons.append(Exon(exon_seq, exon_start))
       
        # print(self.exons) 

        print('fuckme')


    #find motifs in each gene. takes the input gene from the saved gene classes in the full gene list for the fasta file.
    def find_motifs(self):
        """Goes through gene sequence and identifies motifs with locations so they can be drawn"""
        for motif in all_motifs:
            matches = re.finditer(motif, self.sequence)


            for match in matches:
                # Get the exon string
                motif_seq = match.group()
                # Get the start location in the gene
                motif_start = match.start()
                #store exon information for the gene as exon classes
                self.motifs.append(Motif(motif_seq, motif_start))


        print('suck it')
        #can maybe find something similar to the regex above for all instances of the various sequences??
        #def regex tf out of this
        #how 



class Exon:
    """For each Exon in a Gene""" #need starting location in gene sequence and length to draw the rectangle on the line at the right spot.

    # The constructor method to initialize new objects
    def __init__(self, sequence, start):
        """Takes sequence and header to initialize. creates length from sequence"""
    
        self.sequence = sequence 
        self.start = start
        self.length = len(sequence)


class Motif:
    """For each Motif in a Gene"""

    # The constructor method to initialize new objects
    def __init__(self, sequence, start):
        """Takes sequence and header to initialize. creates length from sequence"""
        
        self.sequence = sequence 
        self.length = len(sequence)
        self.start = start




#FUNCTIONS (WE OUT HERE TRYNA)

#getting motifs from the motifs file
def grab_motifs(motif_file):
    """Goes through motifs file and grabs all motifs"""
    #creates empty set for storing motifs
    motifs = set()
    #reading through the file and saving each line to the set
    with open(motif_file, "r") as open_motifs:
        for line in open_motifs:
            line = line.strip()
            motifs.add(line)
    #returns set of motifs
    return motifs

def convert_motifs(motifs):
    #making dictionary to convert motif sequences to regex patterns
    translator = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'N': '[ATGC]',  # Any nucleotide
        'R': '[AG]',    # puRine
        'Y': '[CT]',    # pYrimidine
        'S': '[GC]',    # Strong
        'W': '[AT]',    # Weak
        'K': '[GT]',    # Keto
        'M': '[AC]',    # aMino
        'B': '[CGT]',   # Not A
        'D': '[AGT]',   # Not C
        'H': '[ACT]',   # Not G
        'V': '[ACG]'    # Not T
    }

    #list to hold translated motifs
    converted_motifs = []
    #initialize variable to hold current converted motif
    converted_motif = ''

    #loop through each motif and convert the motifs to regex patterns
    for motif in motifs:
        for letter in motif: 
            if letter.isupper(): #for each uppercase letter, can use the translator as is
                converted_motif += translator.get(letter, letter) 

            else:
                converted_motif += translator.get(letter.upper(), letter).lower() #for lowercase letters, have to make uppercase to use translator

        converted_motifs.append(converted_motif) #add translated motif
        converted_motif = '' #reset the current motif
        
    return converted_motifs #returns list of translated motifs
     

        

#find genes in FASTA file
# def grab_genes(fasta_file, current_header = 'no'):
def grab_genes(fasta_file: 'str')->'list':
    """Goes through FASTA file and grabs all gene entries with headers and sequences"""

    #create final gene list that the function will output
    gene_list = []
    
    #initializing header and sequence variables that will keep track of current header and sequence from the file
    current_header = False
    current_seq=''


    #going through fasta file to pull out genes with sequences and headers to make classes for each gene.
    with open(fasta_file, "r") as open_fasta:
        for line in open_fasta:
            line = line.strip()
            
            #for each time a header line is read in the file, resets the current header and sequence while adding them to the list.
            if line[0] == '>':

                #the if statement skips over the initial part of reading through the file and avoid adding an empty entry to the list.
                if current_header:
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
    
    #adds last gene (the loop does not add the last one so this is needed)
    gene_list.append(Gene(current_header, current_seq))

    #function returns final gene list for the file.
    return gene_list




#drawing function
def draw_annotated_gene(gene, exons, motifs):
    """This takes genes and all necessary elements and draws the pictures for each"""
    print('fuckme')
    #once you get the genes, with the exons and the motifs you can draw them. how tf do we lay out the genes per fasta file??






#creates list of all genes from the fasta file with each header and its sequence saved in Gene classes
gene_list = grab_genes(args.fasta_file)

#creates set of all motifs from motifs file
all_motifs = grab_motifs(args.motifs_file)

all_converted_motifs = convert_motifs(all_motifs)

print(all_motifs)
print(all_converted_motifs)


for x in gene_list:
    x.find_exons()
    x.find_motifs()
# find_exons(gene_list[0])

# for x in gene_list:
#     find_exons(x) #is this correct????????
#     find_motifs(x)


# print(gene_list)       
# print(gene_list[0].length)


#draw all genes
# for x in gene_list:
#     draw_annotated_gene(x, x.exons, x.motifs)