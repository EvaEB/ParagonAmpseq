rule all:
    input:
        "Haplotypes/1M_1/finalHaplotypeList_Hcov3_Scov25_occ2_sens0.0100_2.cpmp2.txt"


rule cutadapt:
    input: 
        R1 = "input/raw/{sample}_R1.fastq.gz",
        R2 = "input/raw/{sample}_R2.fastq.gz"
    output:
        R1 = "processed/cutadapt/{sample}_R1.fastq.gz",
        R2 = "processed/{sample}_R2.fastq.gz"
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
    shell:
        "mkdir -p logs/fusedReads;"
        "Rscript --vanilla scripts/FuseReads.R fusedReads {input} &> logs/fusedReads/{wildcards.sample}.log"

rule index:
    input: 
        "{genome}.fasta"
    output:
        "{genome}.fasta.sa"
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa index {wildcards.genome}.fasta"

        
rule align:
    input: 
        reads = "processed/fusedReads/{sample}.fastq.gz",
        ref = "input/reference/PlasmoDB-59_Pfalciparum3D7_Genome.fasta",
        index = "input/reference/PlasmoDB-59_Pfalciparum3D7_Genome.fasta.sa"
    output:
        al = "processed/aligned/{sample}.sam",
        al_sorted = "processed/aligned/{sample}_sorted.sam",
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa mem {input.ref} {input.reads} > {output.al};"
        "samtools sort {output.al} > {output.al_sorted};"
        "samtools index {output.al_sorted}"

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

rule extractAmplicons:
    input:
       fastq = "processed/splitByMarker/{sample}_marker_{marker}.fastq",
       primerfile = "input/primer_files/primers_generated.csv"
    output:
       "processed/extractedAmplicons/{sample}_marker_{marker}.fastq.gz"
    resources:
        mem_mb = 2000
    shell:
       "mkdir -p logs/extractedAmplicons;"
       "Rscript scripts/ExtractAmplicons.R {input.fastq} {output} {wildcards.marker} {input.primerfile} 0.2 &> logs/extractedAmplicons/{wildcards.sample}_{wildcards.marker}.log"
 
rule callHaplotypes:
    input: 
        fq="processed/extractedAmplicons/{sample}_marker_{marker}.fastq.gz",
        primers="input/primer_files/primers_generated.csv"
    output:
        "Haplotypes/{sample}/finalHaplotypeList_Hcov3_Scov25_occ2_sens0.0100_{marker}.txt"
    shell: 
        "Rscript --vanilla scripts/CallHaplotypes.R {wildcards.marker} {input.fq} {input.primers} Haplotypes/{wildcards.sample};"
        "touch {output}"