import matplotlib.pyplot as plt
import matplotlib
import os 
import gzip
import numpy as np
import pandas as pd
import sys

experiment = sys.argv[1]
sample = sys.argv[2]
directory = experiment
samples = [i[:-4] for i in os.listdir(directory+'/logs/fusedReads/')]
amp_positions = pd.read_csv(f'{directory}/input/AmpliconPositions/amplicon_positions.csv',sep='\t')
longest_amplicon = max(amp_positions.right_start - amp_positions.left_end)

colors = [(0.6091898794104678, 0.0, 1.0, 1.0), (1.0, 0.18529430294136176, 0.0, 1.0), (0.4944837886014357, 1.0, 0.0, 1.0), (0.023161131617013827, 0.016177736765972346, 1.0, 1.0), (0.0, 1.0, 0.20036897885323546, 1.0), (0.0, 0.08566375661963932, 1.0, 1.0), (0.0, 1.0, 0.8952182373286364, 1.0), (0.10073339485104227, 1.0, 0.0, 1.0), (0.0, 0.201472695957991, 1.0, 1.0), (1.0, 0.6022064845594256, 0.0, 1.0), (1.0, 0.0, 0.7886036360301064, 1.0), (0.30808663713075457, 0.0, 1.0, 1.0), (0.0, 0.6878702411790647, 1.0, 1.0), (1.0, 0.0, 0.18639715147068092, 1.0), (0.7029398794104678, 1.0, 0.0, 1.0), (0.0, 1.0, 0.7099251017351961, 1.0), (1.0, 0.0, 0.8812507875007872, 1.0), (0.9113959702194997, 1.0, 0.0, 1.0), (1.0, 0.48639754522107453, 0.0, 1.0), (0.09963054632172239, 0.0, 1.0, 1.0), (0.8176459702194999, 0.0, 1.0, 1.0), (1.0, 0.9033097268391386, 0.0, 1.0), (0.5165427279397867, 0.0, 1.0, 1.0), (1.0, 0.0, 0.27904430294136173, 1.0), (0.19338054632172286, 1.0, 0.0, 1.0), (0.0, 0.38676699889935195, 1.0, 1.0), (0.0, 0.9889734834587774, 1.0, 1.0), (0.0, 1.0, 0.10772241105651532, 1.0), (0.4018366371307548, 1.0, 0.0, 1.0), (0.0, 1.0, 0.594116891989296, 1.0), (0.0, 0.8036791805174157, 1.0, 1.0), (0.0235287477934537, 1.0, 0.015442504413092598, 1.0), (0.7955870308811486, 1.0, 0.0, 1.0), (1.0, 0.09264715147068088, 0.0, 1.0), (0.9801467577202873, 0.0, 0.9772064845594255, 1.0), (0.0, 0.5025759382377031, 1.0, 1.0), (1.0, 0.7875007875007874, 0.0, 1.0), (1.0, 0.6948536360301065, 0.0, 1.0), (1.0, 0.0, 0.5801475452210745, 1.0), (0.0, 0.5952230897083839, 1.0, 1.0), (0.0, 1.0, 0.29301554664995555, 1.0), (0.7018370308811488, 0.0, 1.0, 1.0), (0.40073378860143555, 0.0, 1.0, 1.0), (0.2154394856600736, 0.0, 1.0, 1.0), (0.0, 0.29411984742867114, 1.0, 1.0), (0.309189485660074, 1.0, 0.0, 1.0), (0.9102931216901803, 0.0, 1.0, 1.0), (0.0, 1.0, 0.8025716695319163, 1.0), (0.0, 1.0, 0.4088237563958557, 1.0), (1.0, 0.3011032422797128, 0.0, 1.0), (1.0, 0.0, 0.4875003937503936, 1.0), (0.9805143738967269, 0.9764712522065463, 0.0, 1.0), (0.6102927279397868, 1.0, 0.0, 1.0), (1.0, 0.0, 0.39485324227971275, 1.0), (0.0, 0.8963263319880966, 1.0, 1.0), (1.0, 0.0, 0.0, 1.0), (1.0, 0.0, 0.09375, 1.0), (1.0, 0.3937503937503937, 0.0, 1.0), (1.0, 0.0, 0.6959564845594255, 1.0), (0.0, 1.0, 0.5014703241925758, 1.0)]

NT = 'ACGT'

cmap = matplotlib.colormaps['hsv']
colors_nt = {NT[i]: cmap(i/4) for i in range(4)}


with open(f'{directory}/logs/cutadapt/{sample}_R1.log') as f:
    for line in f.readlines():
        if 'Total reads processed:' in line:
            seq_in_r1 = int(line.split()[-1].replace(',',''))
        elif "Reads with adapters" in line:
            adapters_removed_r1 = int(line.split()[-2].replace(',',''))
        elif "Reads written (passing filters):" in line:
            total_out_r1 = int(line.split()[-2].replace(',',''))

with open(f'{directory}/logs/cutadapt/{sample}_R2.log') as f:
    for line in f.readlines():
        if 'Total reads processed:' in line:
            seq_in_r2 = int(line.split()[-1].replace(',',''))
        elif "Reads with adapters" in line:
            adapters_removed_r2 = int(line.split()[-2].replace(',',''))
        elif "Reads written (passing filters):" in line:
            total_out_r2 = int(line.split()[-2].replace(',',''))
total_in = seq_in_r1 + seq_in_r2

with open(f'{directory}/logs/fusedReads/{sample}.log') as f:
    for line in f.readlines():
        if 'Pairs' in line:
            pairs = int(line.split()[0])
        elif 'Merged' in line:
            merged = int(line.split()[0])
        elif 'Not merged' in line:
            not_merged = int(line.split()[0])
            break
          
with open(f'{directory}/logs/fusedReads2_alignmentBased/{sample}.log') as f:
    vals = []
    for line in f.readlines():
        vals.append(int(line.split(' ')[0]))
    notMapped, merged2, notMerged = vals

with open(f'{directory}/logs/aligned/{sample}.log') as f:
    for line in f.readlines():
        if 'reads mapped:' in line:
            mapped = int(line.split()[-1])
        elif 'reads unmapped' in line:
            unmapped = int(line.split()[-1])

counts_per_marker = {}
for i in os.listdir(directory+'/processed/splitByMarker/'):
    if '.fastq' in i:
        if i[:i.find('_marker')] == sample:
            length = 0
            with open(directory+'/processed/splitByMarker/'+i) as f:
                for line in f.readlines():
                    if '@' == line[0]:
                        length += 1
            marker = i.split('.fastq')[0]
            pos = marker.find('marker_')
            marker = marker[pos+7:]
            counts_per_marker[marker] = length


amplicon_counts = {}
for marker in counts_per_marker:
    count = 0
    with gzip.open(f'{directory}/processed/extractedAmplicons/{sample}_marker_{marker}.fastq.gz') as f: 
        for i,l in enumerate(f):
            if (l[:1] == b'@'):
                count += 1
        amplicon_counts[marker] = count

haplotypes_per_amplicon = {}
for marker in counts_per_marker:
    try:
        with open(f'{directory}/processed/Haplotypes/{sample}/finalHaplotypeList_Hcov3_Scov25_occ2_sens0.0100_{marker}.txt') as f:
            first=True
            haplotypes = {}
            for line in f.readlines():
                if first:
                    first = False
                else:
                    fields = line.split()
                    try:
                      haplotypes[fields[3].strip('"')] = int(fields[4])
                    except ValueError:
                      haplotypes[fields[3].strip('"')] = 0
            haplotypes_per_amplicon[marker] = haplotypes
    except FileNotFoundError:
        pass


plt.figure(figsize=[20,4])

#raw reads
pos = 0

plt.text(0,pos,f'Raw reads    ',ha= 'right',size=15)

plt.plot([0,seq_in_r1], [pos,pos],marker='|',color='black')
plt.text(0,pos+0.1,f'fw: {seq_in_r1}')
plt.plot([seq_in_r1,seq_in_r1+seq_in_r2], [pos,pos],marker='|',color='grey')
plt.text(seq_in_r1,pos+0.1,f'rv: {seq_in_r2}',color='grey')
#cutadapt
pos -= 1
plt.text(0,pos,f'Remove adapters    ',ha= 'right',size=15)
plt.plot([0,adapters_removed_r1], [pos,pos],marker='|',color='black')
plt.text(0,pos+0.1,f'reads with adapters: {adapters_removed_r1}')

plt.plot([adapters_removed_r1,total_out_r1], [pos,pos],marker='|',color='black')
plt.text(adapters_removed_r1,pos-0.1,f'reads without adapters: {total_out_r1 - adapters_removed_r1}',va='top')

plt.plot([seq_in_r1,seq_in_r1+adapters_removed_r2], [pos,pos],marker='|',color='grey')
plt.text(seq_in_r1,pos+0.1,f'reads with adapters: {adapters_removed_r2}',color='grey')

plt.plot([seq_in_r1+adapters_removed_r2,seq_in_r1+total_out_r2], [pos,pos],marker='|',color='grey')
plt.text(seq_in_r1+adapters_removed_r2,pos-0.1,f'reads without adapters: {total_out_r2 - adapters_removed_r2}',va='top',color='grey')

#read merging

pos -= 1
plt.text(0,pos,f'Pair merging    ',ha= 'right',size=15)

plt.plot([0,merged],[pos,pos],marker='|',color='black')
plt.text(0,pos+0.1,f'merged by vsearch: {merged}')

plt.plot([merged,merged+merged2],[pos,pos],marker='|',color='black')
plt.text(merged,pos-0.1,f'merged by alignment: {merged2}',va='top')

plt.plot([merged+merged2,merged+merged2+notMapped+notMerged],[pos,pos],marker='|',color='darkred')
plt.text(merged+merged2,pos+0.1,f'not merged: {notMapped+notMerged}',color='darkred')

# aligning
pos -= 1
plt.text(0,pos,f'Aligning    ',ha= 'right',size=15)

plt.plot([0, mapped],[pos,pos],marker='|',color='black')
plt.text(0,pos+0.1,f'mapped: {mapped}')

plt.plot([mapped,mapped+unmapped],[pos,pos],marker='|',color='darkred')
plt.text(mapped,pos+0.1,f'unmapped: {unmapped}',color='darkred')

# split by marker
pos -= 1
start = 0
text_pos = 'up'
plt.text(0,pos,f'Split by marker    ',ha= 'right',size=15)

for a,i in enumerate(counts_per_marker):
    plt.plot([start, start+counts_per_marker[i]],[pos,pos],marker='|',color=colors[a])
    if text_pos == 'up':
        plt.text(start,pos+0.1,f'{i}: {int(counts_per_marker[i])}',color=colors[a])
        text_pos = 'down'
    else:
        plt.text(start,pos-0.1,f'{i}: {int(counts_per_marker[i])}',va='top',color=colors[a])
        text_pos = 'up'

    start += counts_per_marker[i]

# get amplicons
pos -= 1
start = 0
text_pos = 'up'
plt.text(0,pos,f'Complete amplicons    ',ha= 'right',size=15)

for a,i in enumerate(amplicon_counts):
    plt.plot([start, start+amplicon_counts[i]],[pos,pos],marker='|',color=colors[a])
    if text_pos == 'up':
        plt.text(start,pos+0.1,f'{i}: {amplicon_counts[i]}',color=colors[a])
        text_pos = 'down'
    else:
        plt.text(start,pos-0.1,f'{i}: {amplicon_counts[i]}',va='top',color=colors[a])
        text_pos = 'up'
    start += amplicon_counts[i]

#haplotypes
pos -= 1
new_start = 0
plt.text(0,pos,f'Valid Haplotypes    ',ha= 'right',size=15)

for a,haplotypes in enumerate(haplotypes_per_amplicon):
    start = new_start
    this_haplotype = 0
    for haplotype,count in haplotypes_per_amplicon[haplotypes].items():
        if haplotype in ['Noice','Chimera','Indels','Singleton']:
            c = 'grey'
        else:
            c = colors[a]
        plt.plot([start, start+count], [pos,pos],marker='|',color=c)

        start += count
        this_haplotype += count

    if text_pos == 'up':
        plt.text(new_start,pos+0.1,f'{haplotypes}: {this_haplotype}',color=colors[a])
        text_pos = 'down'
    else:
        plt.text(new_start,pos-0.1,f'{haplotypes}: {this_haplotype}',va='top',color=colors[a])
        text_pos = 'up'
    new_start += counts_per_marker[haplotypes]
plt.axis('off')
plt.title(f'Sample: {sample}',size=20)
plt.savefig(f'{experiment}/plots/{sample}_pipeline.png')


n_markers = len(amp_positions)
plt.figure(figsize=[20,n_markers])
fignumber = 1
for file in os.listdir(directory+'/processed/SNPs/'):
    if file[:file.find('_marker')] == sample:
        marker = file[file.find('marker_')+7:-4]
        amplicon_length = (amp_positions[amp_positions.Amplicon == marker].right_start - amp_positions[amp_positions.Amplicon == marker].left_end).iloc[0]
        SNPs = pd.read_csv(directory +'/processed/SNPs/'+ file,sep='\t',names= ['sample','marker','haplotype','SNP_ID','chromosome','genome_pos','gene_pos_nt','amplicon_position_nt','nt_old','nt_new','gene_position_AA','amplicon_position_AA','AA_old','AA_new','fraction'])
        ypos = 0
        plt.subplot(n_markers,1,fignumber)
        fignumber += 1
        for haplotype, SNPs_this_haplotype in SNPs.groupby('haplotype'):
            plt.plot([0,amplicon_length],[ypos,ypos],color='grey',lw=1)
            positions = (SNPs_this_haplotype.amplicon_position_nt)
            if str(SNPs_this_haplotype.nt_new.iloc[0]) in ['A','T','G','C']:
                plt.scatter(positions, [ypos for i in positions],marker='|',c=[colors_nt[i] for i in SNPs_this_haplotype.nt_new])
            ypos -=1
        plt.axis([0,longest_amplicon,1,ypos])
        plt.axis('off')
        plt.title(marker)
plt.suptitle(f'Sample: {sample}',size=20)
plt.savefig(f'{experiment}/plots/{sample}_highlighter.png')

plt.figure(figsize=[10,len(counts_per_marker)/3])

ypos = 0
first=True
green = True
grey=True
for marker in sorted(counts_per_marker.keys()):
    plt.barh(ypos,counts_per_marker[marker],color='black',label=['__nolabel__','assigned to marker',][first])
    plt.barh(ypos,amplicon_counts[marker],color='grey',label=['__nolabel__','full amplicon recovered',][first])
    plt.text(0,ypos,marker,ha='right')
    pos_so_far = 0
    if marker not in haplotypes_per_amplicon:
      continue
    for hap,count in haplotypes_per_amplicon[marker].items():
        if count != 0:
            if marker in hap:
                color = 'green'
                edgecolor = 'darkgreen'
                label = label=['__nolabel__','haplotype called'][green]
                green=False
            else:
                color = 'lightgrey'
                edgecolor = 'lightgrey'
                label = label=['__nolabel__','haplotype discarded (indels, chimera,...)'][grey]
                grey = False
            plt.barh(ypos,count,left=pos_so_far,color=color,edgecolor=edgecolor,label=label)
            pos_so_far += count
    ypos += 1
    first=False
plt.axis('off')
plt.legend()
plt.savefig(f'{experiment}/plots/{sample}_overview_markers.png')
