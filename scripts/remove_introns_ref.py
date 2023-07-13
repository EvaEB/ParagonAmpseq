"""
removes introns from the reference sequence in a primerfile

usage: python remove_introns_ref.py amplicon_file.csv exon_file.csv primer_file.csv
"""

import sys

import pandas as pd
import numpy as np

amplicon_file = sys.argv[1]
exon_file = sys.argv[2]
primer_file = sys.argv[3]

amplicons = pd.read_csv(amplicon_file,sep='\t')
amplicons = amplicons.set_index('Amplicon')
exons = pd.read_csv(exon_file,sep='\t')

first = True
with open(primer_file) as f:
    for line in f.readlines():
        if first:
            print(line)
            first=False
        else:
            marker,fw,rv,ref = line.strip().split(';')
            this_amplicon = amplicons.loc[marker]
            exons_right_chromosome = exons[exons.chromosome == this_amplicon.chrom]
            exons_right_positions = exons_right_chromosome[(exons_right_chromosome.exonStart < this_amplicon.right_start) & (exons_right_chromosome.exonEnd > this_amplicon.left_end)]
            exon_positions = (np.hstack([np.arange(i.exonStart,i.exonEnd-1) for _,i in exons_right_positions.iterrows()]))
            amplicon_positions = np.arange(this_amplicon.left_end,this_amplicon.right_start-1)
            positions_to_keep = [i for i in amplicon_positions if i in exon_positions]
            new_ref = ''.join([ref[i-this_amplicon.left_end] for i in positions_to_keep])
            print(';'.join([marker,fw,rv,new_ref]))
            
