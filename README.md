# CORAL: Compact-genome Oriented RNA-based Annotation using Long reads

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.24.1-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

The CORAL protocol is a snakemake workflow design to <b>annotate compact genomes</b> using <b>long-read RNAseq data</b>.

It uses as **input** clean (primer-trimmed) pre-processed fastq files and maps them to the given genome using <i>Minimap2</i>.
Then creates non-assembled annotations for each fastq given using <i>StringTie v3.0.2</i> and looks for potential operons within those annotations (implementing <i>GAMBA v1.3.2</i>).

After identifying operon transcripts, operon-contained transcripts, and non-operon-related transcripts, it generates consensus annotations for the three sets of transcripts. Then, these sets are merged (using <i>StringTie v3.0.1</i>) to generate two final consensus annotations: 'Merge clean_andOPRNs GTF', that contains all three sets, and 'Merge clean_noOPRNs GTF', that only includes operon-contained transcripts and non-operon-related transcripts.

The quality of the annotation is assayed with <i>BUSCO (v5.8)</i>, and also with _Gffcompare_ when a reference annotation is provided. Finally, CORAL generates an expression matrix of the consensus annotation with all the transcripts (Merge clean_andOPRNs GTF) when specified in the configuration file.

Schematic pipeline:

<img width="695" height="1551" alt="Figure1_new" src="https://github.com/user-attachments/assets/1fe688f6-d507-4a5a-a065-45dd8b265965" />

## Installation
This pipeline is build on _Snakemake_, therefore, you will need to have _Snakemake_ installed (tested on v5.24.1).

Source files for executing this pipeline can be directly downloaded from the Releases page on this repository.
However, due to the presence of a submodule we recommend downloading it using the following:

    git clone --recursive https://github.com/EvoDevoGenomics-UB/CORAL.git

## How to run

To run **CORAL** you just need to modify the 'CORAL-config.yaml' with your parameters and execute it as any other _Snakemake_ file. We recomend to use _conda_ so it will create an environment where install all the dependecies specified in the 'CORAL-env.yml' and 'CORAL-env.merge.yml' files. The command will be like:

    snakemake --use-conda --snakefile CORAL.smk --configfile CORAL-config.yaml --cores 4

## Output files

The CORAL pipeline creates several folders including:
* **alignments**: contains all the reads alignments for each sample individually.
* **index**: contains the minimap2 index of the genome.
* **logs**: contais several log files of diffrent processes (i.e.: minimap2, input_files_stats, etc.).
* **sample_annotations**: contains the GTF created for each samples,
* **annotations**: contains the consensus annotations (merged annotations).
* **GAMBA_results**: contians the output of the GAMBA tool for each samples (i.e. the operons found on each sample).
* **busco_downloads**: contians the BUSCO database used for the BUSCO analysis
* **busco_analysis**: contains the BUSCO results for the main consensus annotaitons.
* **Expression_matix**: contains the ouputs generated for creat the expression matrix of the 'andOPRNs' consensus annotation.
