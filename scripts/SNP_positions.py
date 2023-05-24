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

def determine_haplotype(seq,ref,marker):
    diffs = [(i,ref[i],seq[i]) for i in range(len(seq)) if seq[i]!=ref[i]]
    exons_in_this_chromosome = exons[(exons.chromosome == chromosome)]

    mutIDs = []
    mut_details = []
    for dif in diffs:
        SNPpos_amplicon_0_NT = dif[0]
        SNPpos_genome_0_NT = AmpliconStart_genome_0_NT + SNPpos_amplicon_0_NT
        before_and_including = exons_in_this_chromosome[exons_in_this_chromosome.exonStart<=SNPpos_genome_0_NT]
        exonWithSNPStart_genome_0_NT = before_and_including.exonStart.max() - 1
        if len(before_and_including) == 1:
            SNPpos_gene_0_NT = SNPpos_genome_0_NT - exonWithSNPStart_genome_0_NT
        else:
            before = before_and_including[before_and_including.exonStart < exonWithSNPStart_genome_0_NT]
            SNPpos_gene_0_NT = ((before.exonEnd - (before.exonStart - 1)).sum()+(SNPpos_genome_0_NT-exonWithSNPStart_genome_0_NT))
        
        SNPpos_gene_0_AA = SNPpos_gene_0_NT//3
        gene_base_position = SNPpos_gene_0_NT%3
        offset = int((SNPpos_amplicon_0_NT-gene_base_position)%3)
        SNPpos_amplicon_0_AA = (SNPpos_amplicon_0_NT-offset)//3

        codons = [seq[i+offset:i+offset+3] for i in range(0,len(seq),3)]
        codons_ref = [ref[i+offset:i+offset+3] for i in range(0,len(seq),3)]
        while len(codons[-1]) != 3:
            codons = codons[:-1]
            codons_ref = codons_ref[:-1]
        translation = [DNA_to_protein[i] for i in codons]
        translation_ref = [DNA_to_protein[i] for i in codons_ref]
        
        AA_old = translation_ref[SNPpos_amplicon_0_AA]
        AA_new = translation[SNPpos_amplicon_0_AA]
        NT_old = dif[1]
        NT_new = dif[2]
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
      
#try:
#    master_haplotypes = pd.read_csv(master_haplotypes_file,sep='\t')
#except FileNotFoundError:
#    master_haplotypes = pd.DataFrame(columns=['haplotype_ID','seq'])

SNP_names = ['SNP_ID','chromosome','genome_pos','gene_pos_nt','amplicon_position_nt','nt_old','nt_new','gene_position_AA','amplicon_position_AA','AA_old','AA_new']
#try:
#    master_SNPs = pd.read_csv(master_SNPs_file,sep='\t')
#except FileNotFoundError:
#    master_SNPs = pd.DataFrame(columns=SNP_names)


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

for seq,name in zip(seqs,names):
    try:
        n_reads = haplotypes[haplotypes.Haplotype == name].iloc[0].Reads
    except IndexError:
        continue
    freq = n_reads/total_reads
    
    gene_details = []
    #if any(seq == master_haplotypes.seq):
    #    haplotype_ID = (master_haplotypes[seq == master_haplotypes.seq].iloc[0].haplotype_ID)
    #    SNP_IDs = haplotype_ID.split('_')[1].split('-')
    #    SNP_IDs = [marker+'_'+i for i in SNP_IDs]
    #else:
    haplotype_ID, SNP_IDs,SNP_details = determine_haplotype(seq,ref,marker)
    #new_line = pd.DataFrame({'haplotype_ID':[haplotype_ID],'seq':[seq]})
    #master_haplotypes = pd.concat((master_haplotypes,new_line),ignore_index=True)
    out = '\t'.join([sample,marker,haplotype_ID]) + '\t'

    if (len(SNP_IDs) > 0) and (SNP_IDs[0] != marker+'_WT'):
        for s,SNP_ID in enumerate(SNP_IDs):
            #if any(SNP_ID == master_SNPs.SNP_ID):
            #    gene_details.append(list(master_SNPs[SNP_ID == master_SNPs.SNP_ID].iloc[0]))
            #else:
            #    try:
            #        SNP_detail = SNP_details[s]
            #    except NameError:
            #        raise NameError(f"{SNP_ID}\n{haplotype_ID}")
            SNP_detail = SNP_details[s]
            SNP_names = ['SNP_ID','chromosome','genome_pos','gene_pos_nt','amplicon_position_nt','nt_old','nt_new','gene_position_AA','amplicon_position_AA','AA_old','AA_new']
            gene_details.append([f'{marker}_{SNP_ID}',chromosome] + SNP_detail)
            #new_line = pd.DataFrame({i:[j] for i,j in zip(SNP_names,gene_details[-1])})
            #master_SNPs = pd.concat((master_SNPs,new_line),ignore_index=True)

    else:
        gene_details.append(['NA']*11)
        gene_details[-1][1] = chromosome
        
    for gene in gene_details: 
        for i in [2,3,4,7,8]:
            if gene[i] != 'NA':
                gene[i] += 1 #switch to one-based indexing for output

        print(out + '\t'.join([str(i) for i in gene]) + '\t' + str(freq))
        
#master_haplotypes.to_csv(master_haplotypes_file,sep='\t',index=False)
#master_SNPs.to_csv(master_SNPs_file,sep='\t',index=False)
