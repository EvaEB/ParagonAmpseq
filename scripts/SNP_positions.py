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

def determine_haplotype(seq,ref,marker,ref_translation,translation_offset):
    diffs = [(i,ref[i],seq[i]) for i in range(len(seq)) if seq[i]!=ref[i]]
    mut_IDs,mut_details = get_SNP_IDs(diffs,seq,ref_translation,translation_offset)
    if len(mut_IDs)>0:
        haplotypeID = marker + '_' + '-'.join(mut_IDs)
    else:
        haplotypeID = marker + '_WT'
    mut_IDs = [marker+'_'+i for i in mut_IDs]
    return haplotypeID, mut_IDs,mut_details


def get_SNP_IDs(muts,seq,ref_translation,translation_offset):
    codons = [seq[i+translation_offset:i+3+translation_offset] for i in range(0,len(seq),3)]
    if len(codons[-1]) != 3:
        codons = codons[:-1]
    seq_translation = [DNA_to_protein[i] for i in codons]
    
    mut_IDs = []
    mut_details = []
    for mut in muts:
        amplicon_position_nt = mut[0]
        nt_old = mut[1]
        nt_new = mut[2]
        
        amplicon_position_AA = (mut[0]-translation_offset)//3
        AA_old = ref_translation[amplicon_position_AA]
        AA_new = seq_translation[amplicon_position_AA]
        
        mut_IDs.append(f'{amplicon_position_nt}{nt_old}{nt_new}/{amplicon_position_AA}{AA_old}{AA_new}')
        mut_details.append((amplicon_position_nt,nt_old,nt_new,amplicon_position_AA, AA_old, AA_new))
    return mut_IDs, mut_details

def get_SNP_details(mut_ID, mut_detail,this_gene,amplicon_start):
    amplicon_position_nt,nt_old,nt_new,amplicon_position_AA, AA_old, AA_new = mut_detail
    genome_start_nt = amplicon_start+amplicon_position_nt
    
    gene_position_nt = 0
    for _,i in this_gene.iterrows():
        if i.exonEnd > genome_start_nt:
            gene_position_nt += genome_start_nt-i.exonStart 
            break
        else:
            gene_position_nt += i.exonEnd - i.exonStart + 1
    gene_position_AA = (gene_position_nt//3) + 1
    
    chromosome = this_gene.chromosome.iloc[0]
    
    return mut_ID,chromosome,genome_start_nt,gene_position_nt,amplicon_position_nt,nt_old,nt_new, gene_position_AA,amplicon_position_AA,AA_old,AA_new

Haplotypes = sys.argv[1]
haplotype_seqs = sys.argv[2]
primer_file = sys.argv[3] 
amplicon_postion_file = sys.argv[4]
exon_positions = sys.argv[5]
marker = sys.argv[6] 
sample = sys.argv[7]
master_haplotypes_file = sys.argv[8] 
master_SNPs_file = sys.argv[9]


haplotypes = pd.read_csv(Haplotypes,sep='\t')
total_reads = haplotypes.Reads.sum()

names = []
seqs = []

with open(haplotype_seqs) as f:
    for line in f.readlines():
        if '>' in line:
            names.append(line.strip('>\n'))
            seqs.append('')
        else:
            seqs[-1] += line.strip()
            

with open(primer_file) as f:
    for line in f.readlines():
        fields = line.strip().split(';')
        if fields[0] == marker:
            ref = fields[-1]
            
with open(amplicon_postion_file) as f:
    for line in f.readlines():
        fields = line.strip().split('\t')
        if fields[0] == marker:
            amplicon_start = int(fields[3])
            amplicon_end = int(fields[5])
            chromosome = fields[1]
            
exons = pd.read_csv(exon_positions,sep='\t')
exons = exons[exons.chromosome == chromosome]
this_exon = exons[(exons.exonStart < amplicon_start) & (amplicon_start < exons.exonEnd)]
if len(this_exon)>1:
    raise ValueError("more than two exons overlap with start of amplicon -- overlapping genes? this is not implemented")
elif len(this_exon)<1: #amplicon starts in an intronic region
    print(amplicon_start)
    this_exon = exons[(exons.exonStart > amplicon_start)]
    this_exon = this_exon[this_exon.exonStart == this_exon.exonStart.min()]
    amplicon_start = this_exon.exonStart.iloc[0]
    

this_gene = exons[exons.proteinName == this_exon.proteinName.iloc[0]]

exon_start = this_exon.exonStart.iloc[0]


exons_before = this_gene[this_gene.exonStart <= exon_start]
intronic_nt = sum(np.array(exons_before.exonStart.iloc[1:]) - np.array(exons_before.exonEnd.iloc[:-1])) - 2
gene_start = exons_before.exonStart.min()
translation_offset = (amplicon_start - exon_start)%3 + 1

codons_ref = [ref[i+translation_offset:i+3+translation_offset] for i in range(0,len(ref),3)]
if len(codons_ref[-1]) != 3:
    codons_ref = codons_ref[:-1]
ref_translation = [DNA_to_protein[i] for i in codons_ref]


try:
    master_haplotypes = pd.read_csv(master_haplotypes_file,sep='\t')
except FileNotFoundError:
    master_haplotypes = pd.DataFrame(columns=['haplotype_ID','seq'])
    
SNP_names = ['SNP_ID','chromosome','genome_pos','gene_pos_nt','amplicon_position_nt','nt_old','nt_new','gene_position_AA','amplicon_position_AA','AA_old','AA_new']

try:
    master_SNPs = pd.read_csv(master_SNPs_file,sep='\t')
except FileNotFoundError:
    master_SNPs = pd.DataFrame(columns=SNP_names)
for seq,name in zip(seqs,names):
    try:
        n_reads = haplotypes[haplotypes.Haplotype == name].iloc[0].Reads
    except IndexError:
        continue
    freq = n_reads/total_reads
    
    gene_details = []
    if any(seq == master_haplotypes.seq):
        haplotype_ID = (master_haplotypes[seq == master_haplotypes.seq].iloc[0].haplotype_ID)
        SNP_IDs = haplotype_ID.split('_')[1].split('-')
        SNP_IDs = [marker+'_'+i for i in SNP_IDs]
    else:
        haplotype_ID, SNP_IDs,SNP_details = determine_haplotype(seq,ref,marker,ref_translation,translation_offset)
        master_haplotypes = master_haplotypes.append({'haplotype_ID':haplotype_ID,'seq':seq},ignore_index=True)
    out = '\t'.join([sample,marker,haplotype_ID]) + '\t'

    if (len(SNP_IDs) > 0) and (SNP_IDs[0] != marker+'_WT'):
        for s,SNP_ID in enumerate(SNP_IDs):
            if any(SNP_ID == master_SNPs.SNP_ID):
                gene_details.append(list(master_SNPs[SNP_ID == master_SNPs.SNP_ID].iloc[0]))
            else:
                try:
                    SNP_detail = SNP_details[s]
                except NameError:
                    raise NameError(f"{SNP_ID}\n{haplotype_ID}\n{master_haplotypes}\n{master_SNPs}")
                gene_details.append(get_SNP_details(SNP_ID,SNP_detail,this_gene,amplicon_start))
                master_SNPs = master_SNPs.append({i:j for i,j in zip(SNP_names,gene_details[-1])},ignore_index=True)

    else:
        gene_details.append(['NA']*11)
        gene_details[-1][1] = chromosome
        
    for gene in gene_details: 
        print(out + '\t'.join([str(i) for i in gene]) + '\t' + str(freq))
        
master_haplotypes.to_csv(master_haplotypes_file,sep='\t',index=False)
master_SNPs.to_csv(master_SNPs_file,sep='\t',index=False)