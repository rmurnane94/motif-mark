#!/usr/bin/env python

#import dependencies
import argparse
import re
import cairo


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
        self.exons = self.find_exons()#creates list of exons for the gene
        self.motifs = self.find_motifs() #creates list to store motifs

    def find_exons(self): #adds exons to the gene
        """Goes through the gene sequence and pulls out exon locations and lengths so they can be drawn"""
        #making a list to store results
        exons = []
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
            exons.append(Exon(exon_seq, exon_start))
       #returns list of exons with sequence and position for each
        return exons


    #find motifs in each gene. takes the input gene from the saved gene classes in the full gene list for the fasta file.
    def find_motifs(self):
        """Goes through gene sequence and identifies motifs with locations so they can be drawn"""

        found_motifs = [] #make list to hold motifs 

        for i, motif in enumerate(all_motifs): #going through all motifs from the motifs file to search for instances in the gene

            converted_motif = convert_motif(motif) #converts motif to a regex pattern
            
            matches = re.finditer(converted_motif, self.sequence) #searches for motif matches in the gene sequence
            for match in matches:
                # Get the motif string
                motif_seq = match.group()
                # Get the start location in the gene
                motif_start = match.start()
                #store motif information for the gene as motif classes
                found_motifs.append(Motif(motif_seq, motif_start, i))


        return found_motifs #returns list of motif classes for the gene
    


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
    def __init__(self, sequence, start, color_key):
        """Takes sequence and header to initialize. creates length from sequence"""
        
        self.sequence = sequence 
        self.length = len(sequence)
        self.start = start
        self.color = color_dict[color_key]

    # def select_color(self, color_key):
        
    #     return color_dict[color_key]


#FUNCTIONS (WE OUT HERE TRYNA)

#getting motifs from the motifs file
def grab_motifs(motif_file):
    """Goes through motifs file and grabs all motifs"""
    #creates empty set for storing motifs
    motifs = []
    #reading through the file and saving each line to the set
    with open(motif_file, "r") as open_motifs:
        for line in open_motifs:
            line = line.strip()
            motifs.append(line)
    #returns set of motifs
    return motifs

def convert_motif(motif):
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
    # converted_motifs = []
    #initialize variable to hold current converted motif
    converted_motif = ''

    #loop through each motif and convert the motifs to regex patterns
    # for motif in motifs:
    for letter in motif: 
        if letter.isupper(): #for each uppercase letter, can use the translator as is
            converted_motif += translator.get(letter, letter) 

        else:
            converted_motif += translator.get(letter.upper(), letter).lower() #for lowercase letters, have to make uppercase to use translator

    # converted_motifs.append(converted_motif) #add translated motif
    # converted_motif = '' #reset the current motif
    
    return converted_motif #returns converted motif
     

        

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
def draw_annotated_genes(gene_list):
    """This takes genes and all necessary elements and draws the pictures for each"""
    print('hi')

    #Create pycairo surface
    WIDTH, HEIGHT = 1000, 1000
    surface = cairo.ImageSurface(cairo.FORMAT_RGB24, WIDTH, HEIGHT)
    ctx = cairo.Context(surface)

    #make background white
    ctx.set_source_rgb(1.0, 1.0, 1.0)
    ctx.mask_surface(surface, 0, 0)

    for i, gene in enumerate(gene_list):
        print(i, gene)

        #set the location for the gene drawing
        gene_location = 50 + 100*i
        
        #Draw the LINE for the whole Gene
        ctx.set_line_width(10) #line width
        ctx.set_source_rgb(0, 0, 0) #line color to black
        ctx.move_to(0,gene_location)   #sets line start based on gene
        ctx.line_to(gene.length,gene_location) #line finish
        ctx.stroke()

        #Draw RECTANGLES for Exons in Gene
        for exon in gene.exons:
            # print('exon')
        
        
            #Parameters: x, y, width, height (top-left corner coordinates and size)
            x = exon.start
            y = gene_location-25
            width = exon.length
            height = 50
            ctx.rectangle(x, y, width, height) #(x0,y0,x1,y1)
            ctx.set_source_rgb(0, 0, 0) #set exon color to black
            ctx.fill() #draw rectangle

    

        #Draw lines for motifs in Gene
        for motif in gene.motifs:
            ctx.set_line_width(motif.length)

            
            ctx.set_source_rgb(*motif.color) #unpacks the color tuple with the *. sets color for each motif
          

            ctx.move_to(motif.start,gene_location-25)        #(x,y)
            ctx.line_to(motif.start,gene_location+25)
            ctx.stroke()

    
    # Draw Key
    x, y = 800, 800
    box_size = 15
    for i, motif in enumerate(all_motifs):
        # Color box
        ctx.set_source_rgb(*color_dict[i])
        ctx.rectangle(x, y, box_size, box_size)
        ctx.fill()
        
        # Text label
        ctx.set_source_rgb(0, 0, 0) # Black text
        ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ctx.set_font_size(14)
        ctx.move_to(x + box_size + 10, y + box_size - 3)
        ctx.show_text(motif)
        
        y += box_size + 10 # Move down for next item

    #save the thing
    surface.write_to_png("my_figure.png")
    print('fuckme')
    #once you get the genes, with the exons and the motifs you can draw them. how tf do we lay out the genes per fasta file??


#creates set of all motifs from motifs file
all_motifs = grab_motifs(args.motifs_file)

color_dict = {0: (1,0,0),
            1: (0,1,0),
            2: (0,0,1),
            3: (1,1,0),
            4: (1,0,1)
            }

#creates list of all genes from the fasta file with each header and its sequence saved in Gene classes
gene_list = grab_genes(args.fasta_file)

draw_annotated_genes(gene_list)

# draw_annotated_genes(gene_list)
# print(gene_list[3].motifs)

# draw_annotated_genes(gene_list)

# for x in gene_list:
#     for y in x.motifs:
#         print(type(y))
#         print(y.sequence)
#     print('------')


# print(all_motifs)


# print([convert_motif(x) for x in all_motifs])


# print(gene_list[0].exons)


# find_exons(gene_list[0])

# for x in gene_list:
#     find_exons(x) #is this correct????????
#     find_motifs(x)


# print(gene_list)       
# print(gene_list[0].length)


#draw all genes
# for x in gene_list:
#     draw_annotated_gene(x, x.exons, x.motifs)