try:
  experiment = config['experiment']
except KeyError:
  experiment = 'DATA'


try:
  genome = config['genome'] 
except KeyError:
  genome = 'PlasmoDB-62_Pfalciparum3D7'


with open(experiment+'/samples') as f:
 sample = f.readlines()
 sample = [i.strip() for i in sample if i[0] != '#']
 
with open(experiment+'/markers') as f:
 marker = f.readlines()
 marker = [i.strip() for i in marker if i[0] != '#']


 
 
localrules: all, index, splitByMarker, getExons, get_primer_file

rule all:
    input:
       expand("{experiment}/processed/SNPs/{sample}_marker_{marker}.csv",sample=sample,marker=marker,experiment=experiment)

rule cutadapt:
    input: 
        R1 = "{experiment}/input/raw/{sample}_R1.fastq.gz",
        R2 = "{experiment}/input/raw/{sample}_R2.fastq.gz"
    output:
        R1 = "{experiment}/processed/cutadapt/{sample}_R1.fastq.gz",
        R2 = "{experiment}/processed/cutadapt/{sample}_R2.fastq.gz"
    conda:
        "envs/cutadapt.yaml"
    shell:
        "mkdir -p $(dirname {wildcards.experiment}/logs/cutadapt/{wildcards.sample}_R1.log);"
        "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o {output.R1} {input.R1} > {wildcards.experiment}/logs/cutadapt/{wildcards.sample}_R1.log;"
        "cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o {output.R2} {input.R2} > {wildcards.experiment}/logs/cutadapt/{wildcards.sample}_R2.log"


rule fuseReads:
    input:
        R1 = "{experiment}/processed/cutadapt/{sample}_R1.fastq.gz",
        R2 = "{experiment}/processed/cutadapt/{sample}_R2.fastq.gz"
    output: 
        "{experiment}/processed/fusedReads/{sample}.fastq.gz"
    conda:
        "envs/vsearch.yaml"
    resources:
        mem_mb=10000
    shell:
        "mkdir -p $(dirname {wildcards.experiment}/logs/fusedReads/{wildcards.sample}.log);"
        "mkdir -p $(dirname {wildcards.experiment}/processed/fusedReads/{wildcards.sample}.fastq);"
        "vsearch --fastq_mergepairs {input.R1} --reverse {input.R2} --fastqout {wildcards.experiment}/processed/fusedReads/{wildcards.sample}_temp_merged.fastq "
        "--fastqout_notmerged_fwd  {wildcards.experiment}/processed/fusedReads/{wildcards.sample}_temp_not_merged_fw.fastq "
        "--fastqout_notmerged_rev  {wildcards.experiment}/processed/fusedReads/{wildcards.sample}_temp_not_merged_rv.fastq "
        "--fastq_truncqual 1 --fastq_maxns 0 --fastq_allowmergestagger  &>  {wildcards.experiment}/logs/fusedReads/{wildcards.sample}.log;"
        "cat  {wildcards.experiment}/processed/fusedReads/{wildcards.sample}_temp_*.fastq >  {wildcards.experiment}/processed/fusedReads/{wildcards.sample}.fastq;"
        "rm  {wildcards.experiment}/processed/fusedReads/{wildcards.sample}_temp_*.fastq;"
        "gzip  {wildcards.experiment}/processed/fusedReads/{wildcards.sample}.fastq"

rule index:
    input: 
        "{genome_full}.fasta"
    output:
        "{genome_full}.fasta.sa"
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa index {wildcards.genome_full}.fasta"

        
rule align:
    input: 
        reads = "{experiment}/processed/fusedReads/{sample}.fastq.gz",
        ref = expand("reference/{genome}_Genome.fasta",genome=genome),
        index = expand("reference/{genome}_Genome.fasta.sa",genome=genome)
    output:
        al = "{experiment}/processed/aligned/{sample}.sam",
        al_sorted = "{experiment}/processed/aligned/{sample}_sorted.sam",
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa mem {input.ref} {input.reads} > {output.al};"
        "samtools sort {output.al} > {output.al_sorted};"
        "samtools index {output.al_sorted};"

rule getAmpliconPosition_list:
    input:
         "{experiment}/input/AmpliconPositions/amplicon_positions.csv"
    output:
         "{experiment}/input/AmpliconPositions/amplicon_position_list.csv"
    shell:
        """
        cat {input} | awk '{{print $1\"\\t\"$2\":\"$4\"-\"$5}}' > {output}
        """

rule splitByMarker:
    input:
        al = "{experiment}/processed/aligned/{sample}_sorted.sam",
        ampliconPos = "{experiment}/input/AmpliconPositions/amplicon_position_list.csv"
    output:
        bamfile = "{experiment}/processed/splitByMarker/{sample}_marker_{marker}.bam",
        fastqfile = "{experiment}/processed/splitByMarker/{sample}_marker_{marker}.fastq"
    conda:
        "envs/bwa.yaml"
    shell:
        "position=$(cat {input.ampliconPos} | grep {wildcards.marker} | awk '{{print $2}}');"
        "samtools view -b {input.al} $position > {output.bamfile};"
        "samtools fastq {output.bamfile} > {output.fastqfile}"

rule getExons:
    input:
        gff = expand("reference/{genome}.gff",genome=genome),
        ampPos = "{experiment}/input/AmpliconPositions/amplicon_positions.csv"
    output:
        "{experiment}/input/AmpliconPositions/exon_positions.csv"
    conda: 
        "envs/AmpSeqPython.yaml"
    shell:
        "python scripts/get_exons.py {input.gff} {input.ampPos} {output}"

rule extractAmplicons:
    input:
       bam = "{experiment}/processed/splitByMarker/{sample}_marker_{marker}.bam",
       primerfile = "{experiment}/input/AmpliconPositions/amplicon_positions.csv",
       exonfile = "{experiment}/input/AmpliconPositions/exon_positions.csv"
    params:
       fastq = "{experiment}/processed/extractedAmplicons/{sample}_marker_{marker}.fastq",
    output:
       gz = "{experiment}/processed/extractedAmplicons/{sample}_marker_{marker}.fastq.gz",
    conda:
        "envs/AmpSeqPython.yaml"
    shell:
       "mkdir -p $(dirname {wildcards.experiment}/processed/extractedAmplicons/{wildcards.sample}.log);"
       "mkdir -p $(dirname {wildcards.experiment}/logs/extractedAmplicons/{wildcards.sample}.log);"
       "python scripts/extract_amplicons.py {input.bam} {input.exonfile} {input.primerfile} {wildcards.marker} {params.fastq} > {wildcards.experiment}/logs/extractedAmplicons/{wildcards.sample}_marker_{wildcards.marker}.log;"
       "gzip {params.fastq}"
 
rule get_primer_file:
    input:
        amplicon_positions = "{experiment}/input/AmpliconPositions/amplicon_positions.csv",
        genome = expand("reference/{genome}_Genome.fasta",genome=genome)
    output:
        "{experiment}/input/primer_files/primers_generated.csv"
    conda:
        "envs/bwa.yaml"
    shell:
        "bash scripts/amplicon_positions_to_primer.sh {input} > {output} "
        
rule remove_exons_from_primer:
    input:
        amplicon_file = '{experiment}/input/AmpliconPositions/amplicon_positions.csv',
        exon_file = '{experiment}/input/AmpliconPositions/exon_positions.csv',
        primer_file = '{experiment}/input/primer_files/primers_generated.csv'
    output:
        "{experiment}/input/primer_files/primers_generated_no_introns.csv"
    conda:
        "envs/AmpSeqPython.yaml"
    shell:
        "python scripts/remove_introns_ref.py {input.amplicon_file} {input.exon_file} {input.primer_file} > {output}"


rule callHaplotypes:
    input: 
        fq="{experiment}/processed/extractedAmplicons/{sample}_marker_{marker}.fastq.gz",
        primers="{experiment}/input/primer_files/primers_generated_no_introns.csv"
    output:
        haplotypeFile = "{experiment}/processed/Haplotypes/{sample}/finalHaplotypeList_Hcov3_Scov25_occ2_sens0.0100_{marker}.txt",
        haplotype_seqs = "{experiment}/processed/Haplotypes/{sample}/{marker}_HaplotypeSeq.fasta"
    envmodules: 
        "R"
    shell: 
        "mkdir -p $(dirname {wildcards.experiment}/logs/callHaplotypes/{wildcards.sample}_{wildcards.marker}.log); "
        "Rscript --vanilla scripts/CallHaplotypes.R {wildcards.marker} {input.fq} {input.primers} {wildcards.experiment}/processed/Haplotypes/{wildcards.sample} "
        "&> {wildcards.experiment}/logs/callHaplotypes/{wildcards.sample}_{wildcards.marker}.log;"
        "touch {output}"

rule call_SNPs:
    input:
        haplotype_file = "{experiment}/processed/Haplotypes/{sample}/finalHaplotypeList_Hcov3_Scov25_occ2_sens0.0100_{marker}.txt",
        haplotype_seqs = "{experiment}/processed/Haplotypes/{sample}/{marker}_HaplotypeSeq.fasta",
        primer_file = '{experiment}/input/primer_files/primers_generated_no_introns.csv',
        amplicon_postion_file = '{experiment}/input/AmpliconPositions/amplicon_positions.csv',
        exon_positions = '{experiment}/input/AmpliconPositions/exon_positions.csv'
    output:
        "{experiment}/processed/SNPs/{sample}_marker_{marker}.csv"
    conda:
        "envs/AmpSeqPython.yaml"
    shell:
        "python scripts/SNP_positions.py {input.haplotype_file} {input.haplotype_seqs} "
        "{input.primer_file} {input.amplicon_postion_file} {input.exon_positions} "
        "{wildcards.marker} {wildcards.sample} "
        "processed/master_files/{wildcards.marker}_haplotypes.csv "
        "processed/master_files/{wildcards.marker}_SNPs.csv > {output}"
