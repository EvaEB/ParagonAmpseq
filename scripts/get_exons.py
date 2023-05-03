"""
script that extracts exon positions of all genes covered by amplicons

usage: python get_exons.py reference.gff amplicon_positions.csv output_file.csv
- reference.gff is a gff file containing features of the reference genome. It 
  should contain all genes and exons of interest and is typically downloaded
  from a database (e.g. PlasmoDB)
- amplicon_positions.csv must contain the following columns: Amplicon,chromosome,left_start,right_end
- output will be writted to output_file.csv, with columns: ID,proteinName,chromosome,exonStart,exonEnd. When protein 
  name is not given in the .gff file, the name of the first amplicon within this gene will be used.
"""
import sys

import pandas as pd 

#----------------------
# Define files to use
#----------------------
gff_file = sys.argv[1]
amplicon_pos_file = sys.argv[2]
exon_file = sys.argv[3]

#----------------------
# Read in gff file
#----------------------
# Find the number of lines to skip in the gff file (comments at start of file)
skip = 0
with open(gff_file) as f:
    while f.readline()[0] == '#':
        skip += 1

# (rownames as defined in GFF file specification: https://www.ensembl.org/info/website/upload/gff3.html )
rownames = ['seqid', # name of the chromosome or scaffold
            'source', # name of the program that generated this feature, or the data source (database or project name)
            'type', # type of feature. Must be a term or accession from the SOFA sequence ontology
            'start', # Start position of the feature, with sequence numbering starting at 1.
            'end', # End position of the feature, with sequence numbering starting at 1.
            'score', # A floating point value.
            'strand', # defined as + (forward) or - (reverse).
            'phase', # One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
            'attributes' # A semicolon-separated list of tag-value pairs, providing additional information about each feature
           ]  
gff = pd.read_csv(gff_file,sep='\t',skiprows=skip,names=rownames)

#---------------------
# read in amplicon positions file
#---------------------
amplicon_positions = pd.read_csv(amplicon_pos_file,sep='\t')

#----------------------
# Extract genes & exons
#----------------------
gene = gff[gff.type == 'protein_coding_gene'] # gff entries that are genes
exons = gff[gff.type == 'CDS'] # gff entries that are coding seqs
already_done = [] #will store genes that are already processed to prevent double entries

with open(exon_file,'w') as f:
    f.write('\t'.join(['ID','proteinName','chromosome','exonStart','exonEnd'])+'\n') #write the header to the csv
    for _,row in amplicon_positions.iterrows():
        # find genes that overlap with amplicon
        right_chrom = gene[gene.seqid == row.chrom] 
        overlap = right_chrom[(right_chrom.start < row.right_end) & (right_chrom.end > row.left_start)]
        
        # extract attributes
        attributes_list = overlap.attributes.iloc[0].split(';')
        attributes = {}
        for attr in attributes_list:
            desc,value = attr.split('=')
            attributes[desc] = value
        
        # find exons and store if new gene
        if attributes['ID'] not in already_done:
            # define name if not given in attributes
            if 'Name' not in attributes:
                attributes['Name'] = row.Amplicon
            exons_this_gene = exons[exons['attributes'].str.contains('Parent='+attributes['ID'])]
            for _,exon in exons_this_gene.iterrows():
                f.write('\t'.join([attributes['ID'], attributes['Name'], row.chrom, str(exon.start),str(exon.end)])+'\n')
            already_done.append(attributes['ID'])
    