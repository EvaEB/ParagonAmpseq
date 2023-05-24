genome = 'PlasmoDB-62_Pfalciparum3D7'
with open('samples') as f:
 sample = f.readlines()
 sample = [i.strip() for i in sample]
 
with open('markers') as f:
 marker = f.readlines()
 marker = [i.strip() for i in marker]	

 
localrules: all, index, splitByMarker, getExons, get_primer_file

rule all:
    input:
       expand("processed/SNPs/{sample}_marker_{marker}.csv",sample=sample,marker=marker)

rule cutadapt:
    input: 
        R1 = "input/raw/{sample}_R1.fastq.gz",
        R2 = "input/raw/{sample}_R2.fastq.gz"
    output:
        R1 = "processed/cutadapt/{sample}_R1.fastq.gz",
        R2 = "processed/cutadapt/{sample}_R2.fastq.gz"
    conda:
        "envs/cutadapt.yaml"
    shell:
        "mkdir -p logs/cutadapt;"
        "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o {output.R1} {input.R1} > logs/{wildcards.sample}_R1.log;"
        "cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o {output.R2} {input.R2} > logs/{wildcards.sample}_R2.log"


rule fuseReads:
    input:
        R1 = "processed/cutadapt/{sample}_R1.fastq.gz",
        R2 = "processed/cutadapt/{sample}_R2.fastq.gz"
    output:
        "processed/fusedReads/{sample}.fastq.gz"
    envmodules: 
        "R"
    resources:
        mem_mb=10000
    shell:
        "mkdir -p logs/fusedReads;"
        "Rscript --vanilla scripts/FuseReads.R processed/fusedReads {input} &> logs/fusedReads/{wildcards.sample}.log"

rule index:
    input: 
        "{genome_full}.fasta"
    output:
        "{genome_full}.fasta.fai"
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa index {wildcards.genome_full}.fasta"

        
rule align:
    input: 
        reads = "processed/fusedReads/{sample}.fastq.gz",
        ref = expand("input/reference/{genome}_Genome.fasta",genome=genome),
        index = expand("input/reference/{genome}_Genome.fasta.fai",genome=genome)
    output:
        al = "processed/aligned/{sample}.sam",
        al_sorted = "processed/aligned/{sample}_sorted.sam",
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa mem {input.ref} {input.reads} > {output.al};"
        "samtools sort {output.al} > {output.al_sorted};"
        "samtools index {output.al_sorted};"

rule splitByMarker:
    input:
        al = "processed/aligned/{sample}_sorted.sam",
        ampliconPos = "input/AmpliconPositions/amplicon_position_list.csv"
    output:
        bamfile = "processed/splitByMarker/{sample}_marker_{marker}.bam",
        fastqfile = "processed/splitByMarker/{sample}_marker_{marker}.fastq"
    conda:
        "envs/bwa.yaml"
    shell:
        "position=$(cat {input.ampliconPos} | grep {wildcards.marker} | awk '{{print $2}}');"
        "samtools view -b {input.al} $position > {output.bamfile};"
        "samtools fastq {output.bamfile} > {output.fastqfile}"

rule getExons:
    input:
        gff = expand("input/reference/{genome}.gff",genome=genome),
        ampPos = "input/AmpliconPositions/amplicon_positions.csv"
    output:
        "input/AmpliconPositions/exon_positions.csv"
    conda: 
        "envs/AmpSeqPython.yaml"
    shell:
        "python scripts/get_exons.py {input.gff} {input.ampPos} {output}"

rule extractAmplicons:
    input:
       bam = "processed/splitByMarker/{sample}_marker_{marker}.bam",
       primerfile = "input/AmpliconPositions/amplicon_positions.csv",
       exonfile = "input/AmpliconPositions/exon_positions.csv"
    params:
       fastq = "processed/extractedAmplicons/{sample}_marker_{marker}.fastq",
    output:
       gz = "processed/extractedAmplicons/{sample}_marker_{marker}.fastq.gz",
    conda:
        "envs/AmpSeqPython.yaml"
    shell:
       "mkdir -p processed/extractedAmplicons;"
       "mkdir -p logs/extractedAmplicons;"
       "python scripts/extract_amplicons.py {input.bam} {input.exonfile} {input.primerfile} {wildcards.marker} {params.fastq} > logs/extractedAmplicons/{wildcards.sample}_marker_{wildcards.marker}.log;"
       "gzip {params.fastq}"
 
rule get_primer_file:
    input:
        amplicon_positions = "input/AmpliconPositions/amplicon_positions.csv",
        genome = expand("input/reference/{genome}_Genome.fasta",genome=genome)
    output:
        "input/primer_files/primers_generated.csv"
    conda:
        "envs/bwa.yaml"
    shell:
        "bash scripts/amplicon_positions_to_primer.sh {input} > {output} "
        
rule remove_exons_from_primer:
    input:
        amplicon_file = 'input/AmpliconPositions/amplicon_positions.csv',
        exon_file = 'input/AmpliconPositions/exon_positions.csv',
        primer_file = 'input/primer_files/primers_generated.csv'
    output:
        "input/primer_files/primers_generated_no_introns.csv"
    conda:
        "envs/AmpSeqPython.yaml"
    shell:
        "python scripts/remove_introns_ref.py {input.amplicon_file} {input.exon_file} {input.primer_file} > {output}"


rule callHaplotypes:
    input: 
        fq="processed/extractedAmplicons/{sample}_marker_{marker}.fastq.gz",
        primers="input/primer_files/primers_generated_no_introns.csv"
    output:
        haplotypeFile = "Haplotypes/{sample}/finalHaplotypeList_Hcov3_Scov25_occ2_sens0.0100_{marker}.txt",
        haplotype_seqs = "Haplotypes/{sample}/{marker}_HaplotypeSeq.fasta"
    envmodules: 
        "R"
    shell: 
        "Rscript --vanilla scripts/CallHaplotypes.R {wildcards.marker} {input.fq} {input.primers} Haplotypes/{wildcards.sample};"
        "touch {output}"

rule call_SNPs:
    input:
        haplotype_file = "Haplotypes/{sample}/finalHaplotypeList_Hcov3_Scov25_occ2_sens0.0100_{marker}.txt",
        haplotype_seqs = "Haplotypes/{sample}/{marker}_HaplotypeSeq.fasta",
        primer_file = 'input/primer_files/primers_generated_no_introns.csv',
        amplicon_postion_file = 'input/AmpliconPositions/amplicon_positions.csv',
        exon_positions = 'input/AmpliconPositions/exon_positions.csv'
    output:
        "processed/SNPs/{sample}_marker_{marker}.csv"
    conda:
        "envs/AmpSeqPython.yaml"
    shell:
        "mkdir -p processed/master_files;"
        "python scripts/SNP_positions.py {input.haplotype_file} {input.haplotype_seqs} "
        "{input.primer_file} {input.amplicon_postion_file} {input.exon_positions} "
        "{wildcards.marker} {wildcards.sample} "
        "processed/master_files/{wildcards.marker}_haplotypes.csv "
        "processed/master_files/{wildcards.marker}_SNPs.csv > {output}"
