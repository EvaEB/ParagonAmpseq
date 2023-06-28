**WARNING** the pipeline has recently been refactored but the documentation does not yet reflect this.
The documentation will be updated soon but for now: the steps below are not up to date and running the pipeline like this WILL FAIL


# Paragon AmpSeq Pipeline
This is a snakemake analysis pipeline designed to extract haplotypes and SNPs from multiplexed amplicon sequencing experiments using Paragon's amplicon sequencing kits, but can be used for any amplicon sequencing experiment

## Installation and setup
This pipeline works in a linux environment only. 
The pipeline needs a conda environemnt running snakemake and HaplotypR installed. All other dependencies will be installed automatically when needed.

### Install conda/mamba
If you do not yet have conda or mamba installed, install mamba by following the instructions [here](https://github.com/conda-forge/miniforge)

note: mamba and conda provide the same functionality, mamba is just much faster. I recommend using mamba, but everything will work when using conda as well.

### Set up Snakemake Environment
```
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

### Install HaplotypR	
follow instructions [here](https://github.com/lerch-a/HaplotypR)

### clone this repo
```
git clone git@github.com:EvaEB/ParagonAmpseq.git
```

## Running the pipeline
### Setup
Snakemake relies on the presence and absence of certain files to decide what components to execute next. It is therefore important to not rename files and maintain the directory structure described here

#### sequencing files
the raw sequencing sequencing files should be placed in the folder `input/raw`. If desired, input files can be sorted into seperate folders within this folder - this structure will then be maintained in all folders with processed data.
The sequencing files are assumed to already by demultiplexed by sample, and files should be named `{sample}_R1.fastq.gz` for the forward reads and `{sample}_R2.fastq.gz` for the reverse reads, where {sample} is replaced with the sample name of your choice.

Minimal example:
```
📁 raw
├── 🧬 sample1_R1.fastq.gz
├── 🧬 sample1_R2.fastq.gz
```

More complex example:
```
📁 raw
├──📁 experiment1
│   ├── 🧬 sample1_R1.fastq.gz
│   ├── 🧬 sample1_R2.fastq.gz
│   ├── 🧬 sample2_R1.fastq.gz
│   ├── 🧬 sample2_R2.fastq.gz
├──📁 experiment 2
│   ├── 📁 experiment2.1
│   │   ├── 🧬 sample1_R1.fastq.gz
│   │   ├── 🧬 sample1_R2.fastq.gz
│   │   ├── 🧬 sample2_R1.fastq.gz
│   │   ├── 🧬 sample2_R2.fastq.gz
│   ├── 📁 experiment2.2
│   │   ├── 🧬 sample1_R1.fastq.gz
│   │   ├── 🧬 sample1_R2.fastq.gz
```
#### Amplicon positions
In the folder `input/AmpliconPositions` there has to be a file `amplicon_positions.csv` that contains the marker names and the forward/reverse primer positions. This file is the same as that provided by Paragon. If needed, the provided template can be used to create this file manually.

#### Reference 
In the folder `input/reference`, the reference genome (named `{genome_name}_Genome.fasta` and the reference annotation (named `{genome_name}.gff`) should be present. For *Plasmodium*, these files can be downloaded from [PlasmoDB](https://plasmodb.org/plasmo/app/downloads)

#### markers & samples
These files contain the names of the markers and samples to be analyzed. This makes it possible to only analyze a subset of the sequenced information at one time. Snakemake will only run the pipeline for the samples and markers in these files that are not yet analyzed.

**markers** needs to contain a list of markers, one per line, names matching those specified in `amplicon_positions.csv`

**samples** needs to contain a list of samples, one per line, names matching the files in the raw input, including any directory structure (for example `experiment1/sample1`)

If any samples are missing, the pipeline will not run. Missing markers will result in empty downstream files

#### Summary
before running the pipeline, the directory structure should now look like this
(running the pipeline will create several new folders in addition to the ones listed here)

```
.
├── 📁 envs
│   ├── ...
├── 📁 input
│   ├── 📁 AmpliconPositions
│   │   ├── 🖹 amplicon_positions.csv
│   ├── 📁 raw
│   │   ├── 🧬 {sample}_R1.fastq.gz
│   │   ├── 🧬 {sample}_R2.fastq.gz
│   │   ├── 🧬 ...
│   └── 📁 reference
│       ├── 🧬 {genome_name}_Genome.fasta
│       └── 🧬 {genome_name}.gff
├── 🖹 markers
├── 🖹 samples
├── 📁 scripts
│   ├── ...
├── 🖹 Snakefile
```

### Running the pipeline
#### Local execution
Once the directory structure is set up, make sure you are in the main directory (`cd ParagonAmpseq`), then enter the following (if more than one core is available, replace `--cores 1` by however many cores available

```bash
mamba activate snakemake #or whatever you named the environment containing snakemake
snakemake --cores 1 --use-conda
``` 

#### Cluster execution




