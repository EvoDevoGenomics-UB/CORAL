# CORAL: Compact-genome Oriented RNA-based Annotation using Long reads

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.24.1-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

The CORAL protocol is a snakemake workflow design to <b>annotate compact genomes</b> using <b>long-read RNAseq data</b>.

It uses as **input** clean (primer-trimmed) pre-processed FASTQ files, mapping them to the provided genome using _Minimap2_.
Then, it creates non-assembled annotations for each FASTQ file using <i>StringTie v3.0.2</i> and identifies potential operons within those annotations (implementing <i>GAMBA v2.0</i>).

After identifying operon transcripts, operon-contained transcripts, and non-operon-related transcripts, CORAL generates consensus annotations for each of the three sets. These sets are then merged (using <i>StringTie v3.0.1</i>) to generate two final consensus annotations:
* **Merge clean_andOPRNs GTF**: contains all three transcripts sets
* **Merge clean_noOPRNs GTF**: includes only operon-contained transcripts and non-operon-related transcripts.

The quality of the annotation is assayed with <i>BUSCO (v5.8)</i>, and optionally with _Gffcompare_ when a reference annotation is provided. Finally, CORAL can generate an expression matrix for the consensus annotation including all the transcript sets (Merge clean_andOPRNs GTF), when specified in the configuration file.

Schematic pipeline:

<img width="695" height="1551" alt="Figure1_new" src="https://github.com/user-attachments/assets/1fe688f6-d507-4a5a-a065-45dd8b265965" />

## Installation
This pipeline is build on _Snakemake_; therefore, you need to have _Snakemake_ installed (tested on v5.24.1).

Source files for runing this pipeline can be directly downloaded from the **Releases** page on this repository.
However, due to the presence of a submodule we recommend downloading it using:

    git clone --recursive https://github.com/EvoDevoGenomics-UB/CORAL.git

## How to run

To run **CORAL**, simply modify the <code>CORAL-config.yaml</code> file with your desired parameters and execute it as any other _Snakemake_ workflow. We recommend running it with <code>--use-conda</code>, which will automatically create an environment to install all the dependecies specified in the <code>CORAL-env.yml</code> and <code>CORAL-env.merge.yml</code> files. Example command:

    snakemake --use-conda --snakefile CORAL.smk --configfile CORAL-config.yaml --cores 4

## Indicating the FASTQ files to use in the Config file
There are two ways to indicate to CORAL where to find your long-read FASTQ files:
1. By using a samplesheet file

   Edit the <code>samplesheet</code> parameter to point into a **TSV** file:
   
       samplesheet: "/absolute/path/to/TSV_file.tsv"

   The TSV file should contain the '_sample names_' and their '_absolute paths_' separated by _tab_ (<code>\t</code>). It should look like this:

       Sample1    /absolute/path/to/your/sample1.fq
       Sample2    /absolute/path/to/your/sample2.part1.fastq
       Sample2    /absolute/path/to/yout/sample2.part2.fastq

   **NOTE**: This format supports multiple FASTQ files for a single '_sample name_', and accepts different FASTQ suffixes (<code>.fq</code> or <code>.fastq</code>)
  
2. Using directory and naming  parameters
  
   Edit the parameters <code>data_dir</code>, <code>samples</code>, and <code>data_suffix</code>. Example:

       data_dir: "/absolute/path/to/the/data/files/"
       samples: ["Sample1","Sample2","Sample3"]
       data_suffix: "_chip-runXXXXX.fastq"

   In this example, CORAL will use as '_sample name_' the ones provided in samples ("Sample1", "Sample2", and "Sample3"), and will interpret that the FASTQ files to use are:
   
       /absolute/path/to/the/data/files/Sample1_chip-runXXXXX.fastq
       /absolute/path/to/the/data/files/Sample2_chip-runXXXXX.fastq
       /absolute/path/to/the/data/files/Sample3_chip-runXXXXX.fastq

   This method is useful when the sample files have **highly similar names**, with just different IDs, and a **single FASTQ file** to process.
   In this case all sample files must share the **same FASTQ suffix** (either <code>.fq</code> or <code>.fastq</code>).

## Output files

The CORAL pipeline creates several folders, including:
* **alignments**: contains all the reads alignments for each sample individually.
* **index**: contains the minimap2 index of the genome.
* **logs**: contains the log files for the different processes.
* **sample_annotations**: contains the GTF annotation files created for each sample.
* **annotations**: contains the consensus annotations (merged annotations).
* **GAMBA_results**: contians the output of the GAMBA tool for each sample (i.e. the operons found on each sample).
* **busco_downloads**: contians the BUSCO database used for the BUSCO analysis.
* **busco_analysis**: contains the BUSCO results for the main consensus annotaitons.
* **Expression_matix**: contains the outputs generated for create the expression matrix of the '_andOPRNs_' consensus annotation.
