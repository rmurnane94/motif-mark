#!/usr/bin/env python

from __future__ import annotations

#import dependencies
import argparse
import re
import cairo

#argparse
def get_args():
    parser = argparse.ArgumentParser(description="to find motifs")
    parser.add_argument("-f", "--fasta_file", help="fasta file to find motifs", type=str, required=True)
    parser.add_argument("-m", "--motifs_file", help="file of motif sequences to identify in fasta", type=str, required=True)
    return parser.parse_args()

args = get_args()

#CLASSES
#Creating class for genes. this will help keep track of entries in the fasta file
# the gene class will also hold its associated motifs and exons
class Gene:
    """For each entry into the fasta file. Also creates lists of its exons and motifs when instantiated""" 
    def __init__(self, 
                 header: str, #the header from the fasta file for the gene
                 sequence: str, #the sequence of the gene from the fasta file
                 motif_list: list[MotifType] #takes the list of MotifTypes which are the motifs to look for taken from the motif file.
                 ):
        """Takes sequence and header to initialize."""
        self.header = header  # stores gene header to be written in the output figure
        self.sequence = sequence # stores gene sequence so it can be searched for exons and motifs
        self.length = len(sequence) # gets length of gene for drawing it 
        self.exons = self.find_exons(sequence)#creates list of exons found in gene sequence
        self.motifs = self.find_motifs(sequence, motif_list) #creates list of motif occurences found in gene sequence

    def find_exons(self, sequence): #adds exons objects to list in the gene
        """Goes through the gene sequence and pulls out exon locations and lengths so they can be drawn"""
        #making a list to store results
        exons = []
        #setting pattern to find all instances of one or more capital letters. sequential capital letters are exons in this case.
        pattern = r'[A-Z]+'
        # Use re.finditer() to locate capitalized exon strings with locations
        matches = re.finditer(pattern, sequence)
        #go through the located exons to get information and make classes
        for match in matches:
            # Get the exon string
            exon_seq = match.group()
            # Get the start location in the gene
            exon_start = match.start()
            #store exon information for the gene as exon classes
            exons.append(Exon(exon_seq, exon_start))
       #returns list of exons with sequence and position for each
        return exons

    #find motifs in each gene. takes the input gene from the saved gene classes in the full gene list for the fasta file.
    def find_motifs(self, sequence, motif_list):
        """Goes through gene sequence and identifies motifs with locations so they can be drawn"""
        found_motifs = [] #make list to hold motifs 
        for motif in motif_list: #going through all motifs from the motifs file to search for instances in the gene
            converted_motif = convert_motif(motif.sequence) #converts motif to a regex pattern
            matches = re.finditer(f'(?={converted_motif})', sequence, re.IGNORECASE) #searches for motif matches in the gene sequence. f'(?={converted_motif})' allows for overlaps
            for match in matches: #for each matching sequence for the current motif
                # Get the motif string
                #motif_seq = match.group()
                # Get the start location in the gene
                motif_start = match.start()
                #store motif information for the gene as motif classes
                found_motifs.append(MotifFound(motif_start, motif.length, motif.color, motif.id))

        return found_motifs #returns list of found motif classes for the gene
    
    def draw(self, ctx, draw_index):
        #set the location for the gene drawing
        gene_location = 75 + 125*draw_index

        #Write the header for the gene
        ctx.set_source_rgb(0, 0, 0) # Black text
        ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ctx.set_font_size(14)
        ctx.move_to(10,gene_location-60)
        ctx.show_text(self.header)
        
        #Draw the LINE for the whole Gene
        ctx.set_line_width(10) #line width
        ctx.set_source_rgba(0, 0, 0, .8) #line color to black
        ctx.move_to(0,gene_location)   #sets line start based on gene
        ctx.line_to(self.length,gene_location) #line finish
        ctx.stroke()

        for exon in self.exons:
            exon.draw(ctx, gene_location)

        for motif in self.motifs:
            motif.draw(ctx, gene_location)
    

class MotifType:
    """For each type of Motif pulled from the file"""
    def __init__(self, sequence, id):
        """Taking sequence from file. adds color for drawing"""
        self.sequence = sequence
        self.id = id
        self.color = self.get_color(id)
        self.length = len(sequence)

    def get_color(self, color_code):
        """Sets the color tuple for each MotifType for drawing the picture"""
        #dictionary helps convert the color code to the color type
        color_dict = {
            0: (1, 0.7, 0.5, .5),
            1: (0,.5,1,.5),
            2: (.5,0,1,.5),
            3: (0,1,0,.5),
            4: (1,0,.7,.7)
            }
        #convert the color code and return the tuple to save to the MotifType object
        color = color_dict[color_code]
        return color
        

class Exon:
    """For each Exon in a Gene""" #need starting location in gene sequence and length to draw the rectangle on the gene line at the right spot.
    def __init__(self, sequence, start):
        """Takes sequence and header to initialize. creates length from sequence"""
        self.sequence = sequence #sequence of the exon
        self.start = start #start location of exon
        self.length = len(sequence) #length of exon

    def draw(self, ctx, gene_location):
        #Parameters: x, y, width, height (top-left corner coordinates and size)
        x = self.start
        y = gene_location-15
        width = self.length
        height = 30
        ctx.rectangle(x, y, width, height) #(x0,y0,x1,y1)
        ctx.set_source_rgba(0, 0, 0, .8) #set exon color to black
        ctx.fill() #draw rectangle


class MotifFound:
    """For each Motif in a Gene""" #need location and length for drawing. colorkey helps generate color and motif spacing for the final drawing
    def __init__(self, start, length, color, id):
        """Takes sequence and header to initialize. creates length from sequence"""
        #self.sequence = sequence  #sequence of motif
        self.length = length #length of motif
        self.start = start #location in gene of motif
        self.color = color #color for that motif.
        self.id = id
    
    def draw(self, ctx, gene_location):
        """Draws the output"""
        
        ctx.set_line_width(self.length)
        ctx.set_source_rgba(*self.color) #unpacks the color tuple with the *. sets color for each motif
        motif_stagger = self.id*10#creates a variable length for each motif to stagger the location around the gene
        ctx.move_to(self.start,gene_location-45+motif_stagger)        #(x,y)
        ctx.line_to(self.start,gene_location+5+motif_stagger)
        ctx.stroke()
    

#FUNCTIONS (WE OUT HERE TRYNA)

#getting motifs from the motifs file
def grab_motifs(motif_file):
    """Goes through motifs file and grabs all motifs"""
    #creates empty list for storing motifs
    motifs = []
    #reading through the file and saving motif from each line to the set. one motif per line in this case.
    with open(motif_file, "r") as open_motifs:
        for i,line in enumerate(open_motifs): #using enumerate allows keeping track of the index which is used to set the color for each MotifType
            line = line.strip()
            motifs.append(MotifType(line, i))
    #returns list of motifs
    return motifs


def convert_motif(motif: str):
    """Convert Motifs to regex friendly patterns for searching gene sequences"""
    #making dictionary to convert motif sequences to regex patterns. accounts for ambiguity
    translator = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': '[TU]', 'U': '[TU]',
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
    #initialize variable to hold current converted motif
    converted_motif = ''
    #loop through each motif and convert each letter to its matching regex pattern
    for letter in motif: 
        if letter.isupper(): #for each uppercase letter, can use the translator as is
            converted_motif += translator[letter] 
        else:
            converted_motif += translator[letter.upper()].lower() #for lowercase letters, have to make uppercase to use translator. then lowers the result

    return converted_motif #returns converted motif
     

#find genes in FASTA file
def grab_genes(fasta_file: 'str', motif_list)->'list':
    """Goes through FASTA file and grabs all gene entries with headers and sequences. Saves a list of Gene class objects with associated exons and motifs"""
    gene_list = [] #create final gene list that the function will output
    
    #initializing header and sequence variables that will keep track of current header and sequence from the file
    current_header = ''
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
                    gene_list.append(Gene(current_header, current_seq, all_motifs))
                #resets the current header when a new one is reached
                current_header = line
                #resets the current sequence as well
                current_seq = ''
                #skips back so the sequence lines can be read.
                continue

            #for every line that is not a header, the sequence from the line is added as a string for each gene.
            current_seq += line
    
    #adds last gene (the loop does not add the last one so this is needed)
    gene_list.append(Gene(current_header, current_seq, motif_list))

    #function returns final gene list for the file.
    return gene_list


#create drawing surface function
def create_context(gene_list):
    """Create pycairo surface and context to be used for creating output"""
    WIDTH = 1200
    HEIGHT = 150*len(gene_list)
    surface = cairo.ImageSurface(cairo.FORMAT_RGB24, WIDTH, HEIGHT)
    ctx = cairo.Context(surface)

    #make background white
    ctx.set_source_rgb(1.0, 1.0, 1.0)
    ctx.mask_surface(surface, 0, 0)

    return ctx, surface

#drawing function
def draw_annotated_genes(ctx, gene_list):
    """This takes genes and all necessary elements and draws the pictures for each"""
    #draw the annotated genes
    for i, gene in enumerate(gene_list):
        gene.draw(ctx, i)
       

def draw_figure_key(ctx, all_motifs):  
    """Draws Figure Key including MotifType sequences and their colors for visual identification"""
    x, y = 1010, 25
    box_size = 15
    for motif in all_motifs:
        # Color boxes
        ctx.set_source_rgba(*motif.color)
        ctx.rectangle(x, y, box_size, box_size)
        ctx.fill()
        
        # Text labels
        ctx.set_source_rgb(0, 0, 0) # Black text
        ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ctx.set_font_size(14)
        ctx.move_to(x + box_size + 10, y + box_size - 3)
        ctx.show_text(motif.sequence)
        
        y += box_size + 10 # Move down for next item

def output_figure(surface):
    fasta_name = args.fasta_file
    name_to_use = fasta_name.split('.')[0]
    surface.write_to_png(f"{name_to_use}.png")


#CALLING FUNCTIONS

#creates list of all motifs from motifs file
all_motifs = grab_motifs(args.motifs_file)

#creates list of all genes as classes from the gene file. Gene classes store exon and motif classes for each gene in addition to gene info.
gene_list = grab_genes(args.fasta_file, all_motifs)

ctx, surface = create_context(gene_list)

#draws the final drawing. goes through each gene and draws it along with its associated exons and motifs. includes color key for motifs.
draw_annotated_genes(ctx, gene_list)

#draws the key for the figure
draw_figure_key(ctx, all_motifs)

#save the figure
output_figure(surface)

print('finished')