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

### Experiment Folder
Create a folder within this directory for each independent experiment (e.g. seperate sequencing run with differing markers). This will contain the raw & processed data and any files specifically related to this experiment (e.g. what markers and samples are present)


#### sequencing files
the raw sequencing sequencing files should be placed in the folder `{experiment_name}/input/raw`. 
The sequencing files are assumed to already by demultiplexed by sample, and files should be named `{sample}_R1.fastq.gz` for the forward reads and `{sample}_R2.fastq.gz` for the reverse reads, where {sample} is replaced with the sample name of your choice.

Minimal example:
```
📁 raw
├── 🧬 sample1_R1.fastq.gz
├── 🧬 sample1_R2.fastq.gz
```

#### Amplicon positions
In the folder `input/AmpliconPositions` there has to be a file `amplicon_positions.csv` that contains the marker names and the forward/reverse primer positions. This file is the same as that provided by Paragon. If needed, the provided template can be used to create this file manually.

#### Reference 
In the folder `reference`, the reference genome (named `{genome_name}_Genome.fasta` and the reference annotation (named `{genome_name}.gff`) should be present. For *Plasmodium*, these files can be downloaded from [PlasmoDB](https://plasmodb.org/plasmo/app/downloads)

#### markers & samples
These files contain the names of the markers and samples to be analyzed. This makes it possible to only analyze a subset of the sequenced information at one time. Snakemake will only run the pipeline for the samples and markers in these files that are not yet analyzed.

**markers** needs to contain a list of markers, one per line, names matching those specified in `amplicon_positions.csv`

**samples** needs to contain a list of samples, one per line, names matching the files in the raw input

If any samples are missing, the pipeline will not run. Missing markers will result in empty downstream files.

Any lines starting with `#` will be ignored - this can be useful to comment out samples/markers you are currently not interested in.

#### Summary
before running the pipeline, the directory structure should now look like this
(running the pipeline will create several new folders in addition to the ones listed here)

```
├── 📁 envs
│   ├── ...
├── 📁 {experiment_name}
│   ├── 🖹 markers
│   ├── 🖹 samples
│   ├── 📁 input
│   │   ├── 📁 AmpliconPositions
│   │   │   ├── 🖹 amplicon_positions.csv
│   │   ├── 📁 raw
│   │   │   ├── 🧬 {sample}_R1.fastq.gz
│   │   │   ├── 🧬 {sample}_R2.fastq.gz
│   │   │   └── 🧬 ...
├── 📁 reference
│   ├── 🧬 {genome_name}_Genome.fasta
│   └── 🧬 {genome_name}.gff
├── 📁 scripts
│   └── ...
├── 🖹 Snakefile
└── 🖹 Readme.md
```

### Running the pipeline
Once the directory structure is set up, make sure you are in the main directory (`cd ParagonAmpseq`), then enter the following (if more than one core is available, replace `--cores 1` by however many cores available, 

When running locally:
```bash
mamba activate snakemake #or whatever you named the environment containing snakemake
snakemake --cores 1 --use-conda --config experiment='{experiment_name} genome='{genome_name}'
``` 

When running in a cluster environment (e.g. SciCore):
```bash
mamba activate snakemake #or whatever you named the environment containing snakemake
snakemake --cores 1 --use-conda --use-envmodules --config experiment='{experiment_name} genome='{genome_name}'
```

replacing `{experiment_name}` and `{genome_name}` by the right values. If these are not provided, they default to `DATA` and `PlasmoDB-62_Pfalciparum3D7` )


## Output & Further analysis
Once the pipeline has completed, there will be three new folders in the `{experiment_name}` folder:
```
├── 📁 {experiment_name}
│   ├── ...
│   ├── 📁 logs
│   ├── 📁 plots
│   └── 📁 processed
```
`processed` is the main output folder. Each pipeline step has its own folder, with the main output (SNP and haplotype frequencies) stored in `SNPs`: for each marker & sample, there is a .csv file containing complete information which can be used for any downstream analyses.

`logs` contains logging files from the pipeline, which are used to create some of the plots and can be used to debug if any issues arise. There will be a folder per pipeline step with one file per time this step was run.

`plots` contains some plots, one per sample
- {sample}_overview_markers.png - diagnostics for haplotype extraction: the length of the bar represents the amount of reads aligning to this marker in total, green part indicates how many of these reads could successfully be used for haplotype calling. Large black proportions indicate many incomplete reads (not fully covering the amplicon, e.g. due to failing to merge), large grey areas indicate likely presence of indels.
- {sample}_highligher.png - an overview of the SNPs present per haplotype for every marker in the sample. 
- {sample}_pipeline.png - an overview of how many reads get processed in each step of the pipeline.  



