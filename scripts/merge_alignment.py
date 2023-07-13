import pysam
import numpy as np
import sys
from collections import defaultdict

def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    (see https://www.biostars.org/p/306041/)
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        qname = read.query_name
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            yield read, None
        elif qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

unmerged_aligned_file = sys.argv[1]
outfile = sys.argv[2]

samfile = pysam.AlignmentFile(unmerged_aligned_file)
reads_paired = []

new_seqs = []
unpaired = 0
merged_total = 0
not_merged_total = 0

for read1, read2 in read_pair_generator(samfile):
    if read2 is None:
        new_seqs.append((read1.seq, read1.query_qualities,read1.query_name+'_unmerged'))
        unpaired += 1
    else:
        read_1_aligned = read1.get_reference_positions()
        read_2_aligned = read2.get_reference_positions()
        overlap = set(read_1_aligned).intersection(read_2_aligned)
        if len(overlap) != 0:
            seq_pos1 = [i for i,a in enumerate(read_1_aligned) if a in overlap]
            seq_1_only = read1.query_alignment_sequence[:min(seq_pos1)]
            overlap_1 = read1.query_alignment_sequence[min(seq_pos1):max(seq_pos1)]
            quals_1 = read1.query_qualities[read1.query_alignment_start:read1.query_alignment_end]
            quals_1_only = quals_1[:min(seq_pos1)]
            quals_1_overlap = quals_1[min(seq_pos1):max(seq_pos1)]
            
            seq_pos2 = [i for i,a in enumerate(read_2_aligned) if a in overlap]        
            overlap_2 = read2.query_alignment_sequence[min(seq_pos2):max(seq_pos2)]
            seq_2_only = read2.query_alignment_sequence[max(seq_pos2):]
            quals_2 = read2.query_qualities[read2.query_alignment_start:read2.query_alignment_end]
            quals_2_overlap = quals_2[min(seq_pos2):max(seq_pos2)]
            quals_2_only = quals_2[max(seq_pos2):]
            
            if overlap_1 == overlap_2:
                new_overlap = overlap_1
                new_quals_overlap = [max(a,b) for a,b in zip(quals_1_overlap,quals_2_overlap)]
            else:
                new_overlap = ''
                new_quals_overlap = []
                for b1,b2,q1,q2 in zip(overlap_1, overlap_2,quals_1_overlap, quals_2_overlap):
                    if (q1>q2):
                        new_overlap += b1 
                        new_quals_overlap.append(int(q1))
                    else:
                        new_overlap += b2 
                        new_quals_overlap.append(int(q2))
    
                
            new_seq = seq_1_only + new_overlap + seq_2_only
            new_quals = np.hstack((quals_1_only, new_quals_overlap, quals_2_only))
            new_seqs.append((new_seq, new_quals,read1.query_name+'-combined'))
            merged_total += 2
        else:
            new_seqs.append((read1.seq, read1.query_qualities,read1.query_name+'-unmerged_fw'))
            new_seqs.append((read2.seq, read2.query_qualities,read2.query_name+'-unmerged_rv'))
            not_merged_total += 2
    
print(f"{unpaired} reads with no mate or not mapped")    
print(f"{merged_total} reads could be merged ({merged_total//2} merged reads)")
print(f"{not_merged_total} reads left unmerged due to no overlap")

with open(outfile,'w') as f:
    for i in new_seqs:
        try:
            f.write(f'@{i[2]}\n{i[0]}\n+\n{"".join([chr(round(j+33)) for j in i[1]])}\n')
        except TypeError:
            print(i)
            raise
        
