"""
SNP_positions.py

Author: Eva Bons

Usage: SNP_positions.py [haplotype_file] [haplotype_seqs] [primer_file] [amplicon_postion_file] [exon_positions] [marker] [sample]"
"""


import os
import numpy as np
import pandas as pd

import sys

DNA_to_protein = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

rev_complement = {'A':'T',
                  'T':'A',
                  'C':'G',
                  'G':'C'}

def determine_haplotype(seq,ref,marker,AmpliconStart_genome_0_NT):
    """finds the SNPs in a given sequence
    """
    diffs = [(i,ref[i],seq[i]) for i in range(len(seq)) if seq[i]!=ref[i]]
    exons_in_this_chromosome = exons[(exons.chromosome == chromosome)]
    strand = exons_in_this_chromosome[exons_in_this_chromosome.exonEnd >= AmpliconStart_genome_0_NT].iloc[0].strand
    gene = exons_in_this_chromosome[exons_in_this_chromosome.exonEnd >= AmpliconStart_genome_0_NT].iloc[0].proteinName

    # see if amplicon starts in an intron (amplicon_start is not found in any exon)
    if len(exons_in_this_chromosome[(exons_in_this_chromosome.exonStart < AmpliconStart_genome_0_NT) & (exons_in_this_chromosome.exonEnd > AmpliconEnd_genome_0_NT)]) == 0:
        if strand == '+':
            start_after = (exons_in_this_chromosome[exons_in_this_chromosome.exonStart > AmpliconStart_genome_0_NT])
            AmpliconStart_genome_0_NT = start_after.iloc[start_after.exonStart.argmin()].exonStart - 1
            #ampliconStart needs to be updated to match the sequence (where introns have been removed)
        elif strand == '-':
            raise NotImplementedError('amplicon starting in intron on negative strand not yet implemented')

    mutIDs = []
    mut_details = []
    for dif in diffs:
        SNPpos_amplicon_0_NT = dif[0]
        SNPpos_genome_0_NT = AmpliconStart_genome_0_NT + SNPpos_amplicon_0_NT
        if strand == '+':
            #find the exon where this amplicon start
            before_and_including = exons_in_this_chromosome[exons_in_this_chromosome.exonStart<=SNPpos_genome_0_NT] 
            exonWithSNPStart_genome_0_NT = before_and_including.exonStart.max() - 1 #-1 because exon-file is 1-indexed
            if len(before_and_including) == 1: #no introns before the amplicons starts
                SNPpos_gene_0_NT = SNPpos_genome_0_NT - exonWithSNPStart_genome_0_NT
            else: #there might be introns before the amplicon starts
                before = before_and_including[before_and_including.exonStart < exonWithSNPStart_genome_0_NT] #extract exons before the exon where the amplicon starts
                before = before[before.proteinName == gene] #only include exons of the same gene

                #position in the gene is the length of all preceding exons + the position in the current exon
                SNPpos_gene_0_NT = ((before.exonEnd - (before.exonStart - 1)).sum()+(SNPpos_genome_0_NT-exonWithSNPStart_genome_0_NT))
        else: #negative strand gene, same logic but now using exons after (in + sense)
            AmpliconLength = len(seq)
            after_and_including = exons_in_this_chromosome[exons_in_this_chromosome.exonEnd>=SNPpos_genome_0_NT]
            exonWithSNPEnd_genome_0_NT = after_and_including.exonEnd.min() - 1
            after_and_including = after_and_including[after_and_including.proteinName == gene]

            if len(after_and_including) == 1:
                SNPpos_gene_0_NT = exonWithSNPEnd_genome_0_NT - SNPpos_genome_0_NT + 1
            else:
                print(after_and_including)
                raise NotImplementedError('reverse strand genes with introns not yet implemented')
            
            SNPpos_amplicon_0_NT = AmpliconLength - SNPpos_amplicon_0_NT #convert to - strand
            
        # find position in terms of AA
        SNPpos_gene_0_AA = SNPpos_gene_0_NT//3
        
        gene_base_position = SNPpos_gene_0_NT%3
        offset = int((SNPpos_amplicon_0_NT-gene_base_position)%3)
        SNPpos_amplicon_0_AA = (SNPpos_amplicon_0_NT-offset)//3
    
        if strand == '-': #reverse complement sequence and reference
            seq = ''.join(rev_complement[i] for i in seq[::-1])
            ref = ''.join(rev_complement[i] for i in ref[::-1])

        # find codons of sequence and reference & translate
        codons = [seq[i+offset:i+offset+3] for i in range(0,len(seq),3)]
        codons_ref = [ref[i+offset:i+offset+3] for i in range(0,len(seq),3)]
        while len(codons[-1]) != 3:
            codons = codons[:-1]
            codons_ref = codons_ref[:-1]
        translation = [DNA_to_protein[i] for i in codons]
        translation_ref = [DNA_to_protein[i] for i in codons_ref]
        
        # find the actual SNPs in terms of AA
        if SNPpos_amplicon_0_AA < len(translation_ref):
            AA_old = translation_ref[SNPpos_amplicon_0_AA]
            AA_new = translation[SNPpos_amplicon_0_AA]
        else:
            AA_old = '-'
            AA_new = '-'
        if strand == '+':
            NT_old = dif[1]
            NT_new = dif[2]
        else:
            NT_old = rev_complement[dif[1]]
            NT_new = rev_complement[dif[2]]
        mutIDs.append(f'{SNPpos_amplicon_0_NT+1}{NT_old}{NT_new}/{SNPpos_amplicon_0_AA+1}{AA_old}{AA_new}')
        mut_details.append([SNPpos_genome_0_NT, 
                            SNPpos_gene_0_NT, SNPpos_amplicon_0_NT, NT_old, NT_new,
                            SNPpos_gene_0_AA, SNPpos_amplicon_0_AA, AA_old, AA_new])
    
    if len(mutIDs)>0:
        haplotypeID = marker + '_' + '-'.join(mutIDs)
    else:
        haplotypeID = marker + '_WT'
    return haplotypeID, mutIDs, mut_details

# command line args
Haplotypes = sys.argv[1]
haplotype_seqs = sys.argv[2]
primer_file = sys.argv[3] 
amplicon_postion_file = sys.argv[4]
exon_positions = sys.argv[5]
marker = sys.argv[6] 
sample = sys.argv[7]
master_haplotypes_file = sys.argv[8] 
master_SNPs_file = sys.argv[9]

# read in relevant files
try:
  haplotypes = pd.read_csv(Haplotypes,sep='\t')
except pd.errors.EmptyDataError:
  sys.stderr.write("No data in file - terminating")
  exit()  

total_reads = haplotypes.Reads.sum() ###THIS MIGHT NEED TO CHANGE TO EXCLUDE INDELS/CHIMERA

names = []
seqs = []

with open(haplotype_seqs) as f:
    for line in f.readlines():
        if '>' in line:
            names.append(line.strip('>\n'))
            seqs.append('')
        else:
            seqs[-1] += line.strip()
      
SNP_names = ['SNP_ID','chromosome','genome_pos','gene_pos_nt','amplicon_position_nt','nt_old','nt_new','gene_position_AA','amplicon_position_AA','AA_old','AA_new']


with open(primer_file) as f:
    for line in f.readlines():
        fields = line.strip().split(';')
        if fields[0] == marker:
            ref = fields[-1]
            
with open(amplicon_postion_file) as f:
    for line in f.readlines():
        fields = line.strip().split('\t')
        if fields[0] == marker:
            AmpliconStart_genome_0_NT = int(fields[3])
            AmpliconEnd_genome_0_NT = int(fields[5])   
            chromosome = fields[1]
            break
exons = pd.read_csv(exon_positions,sep='\t')

# analyze each sequence
for seq,name in zip(seqs,names):
    try:
        n_reads = haplotypes[haplotypes.Haplotype == name].iloc[0].Reads
    except IndexError:
        continue
    freq = n_reads/total_reads
    
    gene_details = []
    haplotype_ID, SNP_IDs,SNP_details = determine_haplotype(seq,ref,marker,AmpliconStart_genome_0_NT)
    out = '\t'.join([sample,marker,haplotype_ID]) + '\t'

    if (len(SNP_IDs) > 0) and (SNP_IDs[0] != marker+'_WT'):
        for s,SNP_ID in enumerate(SNP_IDs):
            SNP_detail = SNP_details[s]
            SNP_names = ['SNP_ID','chromosome','genome_pos','gene_pos_nt','amplicon_position_nt','nt_old','nt_new','gene_position_AA','amplicon_position_AA','AA_old','AA_new']
            gene_details.append([f'{marker}_{SNP_ID}',chromosome] + SNP_detail)
    else:
        gene_details.append(['NA']*11)
        gene_details[-1][1] = chromosome
        
    for gene in gene_details: 
        for i in [2,3,4,7,8]:
            if gene[i] != 'NA':
                gene[i] += 1 #switch to one-based indexing for output

        print(out + '\t'.join([str(i) for i in gene]) + '\t' + str(freq) + '\t' + str(total_reads))
        
