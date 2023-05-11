"""
removes primers & intronic regions from reads

usage: python extract_amplicons.py reads.bam exon_positions.csv amplicon_positions.csv marker out.fastq
"""
import sys 

import pysam
import pandas as pd
import numpy as np


# define files
bamfile = sys.argv[1]
exon_positions = sys.argv[2]
amplicon_positions = sys.argv[3]
marker = sys.argv[4]
outfile = sys.argv[5]

# read in amplicon file and find amplicon positions for this marker
amplicons = pd.read_csv(amplicon_positions,sep='\t')
this_amplicon = amplicons[amplicons.Amplicon == marker].iloc[0]

# read in exon file and derive relevant introns
exons = pd.read_csv(exon_positions,sep='\t')
introns = []
for protein, exons_this_protein in exons.groupby('proteinName'): 
    if len(exons_this_protein) > 1:
        # introns are the spaces between exons --> so positions between end of one exon and start of the next
        intron_starts = exons_this_protein.exonEnd.iloc[:-1]
        intron_ends = exons_this_protein.exonStart.iloc[1:]
        introns.append(pd.DataFrame({'protein': protein,
                                     'chromosome': exons_this_protein.chromosome.iloc[0],
                                     'intronStart': list(intron_starts),
                                     'intronEnd': list(intron_ends)}))
introns = pd.concat(introns).reset_index(drop=True)

#iterate over the reads, extracting only the parts relevant for the amplicon (no primers, no introns)
right_chromosome = introns[introns.chromosome == this_amplicon.chrom]
introns_in_amplicon = right_chromosome[(right_chromosome.intronStart < this_amplicon.right_end) &  (right_chromosome.intronEnd > this_amplicon.left_start)]

with open(outfile, 'w') as f:
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam:
            ref_positions = read.get_reference_positions()
            sequence = read.seq
            quals = read.query_qualities
            
            current_position=0
            # deal with insertions, clipping etc.
            for cigar in read.cigartuples:
                if cigar[0] == 1: #insertion:
                    insert = [ref_positions[current_position]+1]*cigar[1]
                    ref_positions = ref_positions[:current_position] + insert + ref_positions[current_position:]
                elif cigar[0] == 4: #soft clipping
                    sequence = sequence[:current_position]+sequence[current_position+cigar[1]:]
                    quals = quals[:current_position]+quals[current_position+cigar[1]:]

                elif cigar[0] in [0,2,5]: #match,deletion,hard clipping
                    pass
                else:
                    raise 
                if cigar[0] != 2:
                    current_position += cigar[1]
            
            ref_positions = np.array(ref_positions)
            
            #find primer regions to remove
            to_remove = [ref_positions<this_amplicon.left_end,ref_positions>=this_amplicon.right_start-1]

            # find intron regions to remove
            for _,intron in introns_in_amplicon.iterrows():
                to_remove.append((intron.intronStart < ref_positions) & (intron.intronEnd > ref_positions))
            
            to_remove = np.logical_or.reduce(to_remove)
            
            #remove primers and introns from sequence and qualities
            primers_introns_removed = np.array([i for i in sequence])[~to_remove]


            new_seq = ''.join(primers_introns_removed)
            new_qualities = (np.array(quals)[~to_remove])
            
            # format as fastq and print
            if len(new_seq)>0:
                fastq_string = "@%s\n%s\n+\n%s\n" % (read.query_name, new_seq, pysam.qualities_to_qualitystring(new_qualities))
                f.write(fastq_string)
        