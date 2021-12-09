# CHROME ChIPseq
This repository contains the pipeline/code that was used to analyze the ChIPseq data in GEO: GSE190413. This pipeline was written for execution on the NYU big purple HPC. This readme describes how to execute the snake make workflow for paired-end ChIP-seq data pre-processing (fastq -> peak calling), Utilizing bowtie2 for alignment and MACS2 for peak calling. Post peak calling R was utilized for differential peak enrichment analysis and peak annotation. The conda environment used for the snakemake pipeline is in ChIPseq.yml. The conda environment used for the DiffBind analysis in R is in DiffBind.yml.

To run this analysis, clone this github repo and then download the fastq.gz files listed in samples_info.tab to a directory titled 'fastq' within the cloned 'ChIPseq' directory. Clone the ChIPseq.yml. Run the snakemake work flow via 'bash snakemake_init.sh'. 'snakemake_init.sh' load the conda environment and executes the workflow. This workflow was designed to run on the NYU BigPurple HPC which uses the Slurm workload manager. This workflow will likely have to be augmented for running on other systems. To run the differential peak enrichment analysis, clone the DiffBind.yml conda environment and execute DiffBind.R.

# Description of files:
## Snakefile
This file contains the processing work flow (fastq -> peak calling). Each command and associated parameters can be found here.
## samples_info.tab
This file contains a tab deliminated table with file information and associated meta data.
## config.yaml
This file contains general configuaration info.

		1. Where to locate the samples_info.tab file
		2. Path to bowtie2 indexed genome
## cluster_config.yml
Sbatch parameters for each rule in the Snakefile workflow
## cat_rename.py
This script was utilized to concatenate samples split over multiple lanes and rename them as specified in samples_info.tab
This script:

		1. Concatenates fastq files for samples that were split over multiple sequencing lanes
		2. Renames the fastq files from the generally verbose ids given by the sequencing center to those supplied in the Samples_info.tab file.
		3. The sample name, condition, and replicate columns are concatenated and form the new sample_id_Rx.fastq.gz files
		4. This script is executed snakemake_init.sh prior to snakemake execution
## snakemake_init.sh
This bash script:

		1. loads the miniconda3/4.6.14 module
		2. Loads the conda environment. You can clone the conda environment using the ChIPseq.yml file and modify this bash script to load the env.
		3. Executes snakemake
## ChIPseq.yml
This file contains the environment used by this pipeline.
 
## FRP.py
This file computes the fraction of reads in peaks (FRP) and outputs a table with FRP, total fragments, and fragments within peaks.

## DiffBind.R
This file contains the code used for differential enrichment analysis and peak annotation

## DiffBind.yml
This file contains the R environment used by DiffBind.R