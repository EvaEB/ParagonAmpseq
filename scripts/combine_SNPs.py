###############################################################################
### Combine_SNPs.py 
###
### Combines all [snps].csv files in a given directory into a single file and
### combines triplicates (only keeping variants present in >2/3) when applicable
### 
### Author: Eva Bons
### 
### Usage: 
### combine_SNPs.py [directory] [triplicate(None|left|right)] [outfile]
### 
###############################################################################

import os
import pandas as pd
import numpy as np
import re
from collections import Counter
import sys


directory = sys.argv[1]
triplicate = sys.argv[2] #one of 'None' (dataset is not in triplicate), 'left' (replicate number is on left ), 'right' (replicate number is on right)
outfile = sys.argv[3] 

#read in all SNP files
first=True
for file in os.listdir(directory):
    new_SNPs = pd.read_csv(directory+file,sep='\t',names=['sample','marker','haplotype_ID', 'SNP_ID','chromosome',
                                                        'genome_pos','gene_position_nt','amplicon_position_nt',
                                                        'old_nt','new_nt','gene_position_aa','amplicon_position_aa',
                                                        'old_aa','new_aa','fraction','reads'])
    if first:
        SNPs = new_SNPs #initialize the dataframe
        first=False
    else:
        SNPs = pd.concat((SNPs, new_SNPs),ignore_index=True) #add to the dataframe


if triplicate != 'None': # combine triplicates if needed
    SNPs['sample']= SNPs['sample'].astype('string')
    # extract the replicate ID & sample name
    if triplicate == 'left':
        SNPs['replicate'] = SNPs['sample'].str[0]
        SNPs['sample'] = SNPs['sample'].str[1:]
    elif triplicate == 'right':
        SNPs['replicate'] = SNPs['sample'].str[-1]
        SNPs['sample'] = SNPs['sample'].str[:-1]
    else:
        raise ValueError(f"value for 'triplicate' should be one of 'None','left' or 'right', is {triplicate}")
    
    # reduce to replicate-sample-haplotype-fraction
    haplotypes = SNPs[['replicate','sample','haplotype_ID','fraction','reads']].drop_duplicates()
    # group by haplotype ID and sample
    grouped = haplotypes.groupby(['sample','haplotype_ID'],as_index=False)

    # combine groups (weighted average for fraction, number of replicates)
    combined = grouped.agg(combined_fraction = ("fraction",lambda x: np.average(x,weights=haplotypes.loc[x.index,"reads"])),
                        count = ("replicate","count"),reads = ("reads","sum"))

    # filter out those with count < 2, keep those with count >=2
    keep = combined[combined['count'] >= 2]
    
    haplotypes_3plicates = SNPs.merge(keep,'right',on='haplotype_ID') #combine SNPs and the filtered/grouped dataset, keeping only those relevant
    print(haplotypes_3plicates) #print for logging purposes

    haplotypes_3plicates = haplotypes_3plicates.drop(columns=['replicate','fraction','count','sample_y','reads_x']) #drop columns we don't need anymore
    haplotypes_3plicates = haplotypes_3plicates.rename(columns={'combined_fraction':'fraction',
                                                                'sample_x':'sample','reads_y': 'reads'}) #rename the combined fraction to normal fraction
    haplotypes_3plicates = haplotypes_3plicates.drop_duplicates()
    SNPs = haplotypes_3plicates 

all_haplotypeIDs = SNPs.haplotype_ID.unique() #extract the unique haplotypeIDs to rename to something more straightforward
marker_counts = Counter() #keeps track of the number of haplotypes per marker we've named
new_names = {}
for i in all_haplotypeIDs:
    if (i[-2:] != 'WT'): #if it's not the wildtype haplotype
        marker = re.match('(.*)_[0-9]+[ATCG][ATCG]/',i).groups()[0] #find the marker name
        index = marker_counts[marker] 
        marker_counts[marker] += 1  
        new_names[i] = f'{marker}_{index}'
    else:
        new_names[i] = i #keep the name if WT
        
SNPs['haplotype_ID'] = SNPs['haplotype_ID'].map(new_names) #rename the haplotypes with the just generated simpler names

SNPs.to_csv(outfile,index=False)
