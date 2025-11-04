# CORAL: Compact-genome Oriented RNA-based Annotation using Long reads

The CORAL protocol is a set of scripts to <b>annotate compact genomes</b> using <b>long-read RNAseq data</b>.

It uses as **input** clean (primer-trimmed) pre-processed fastq files and maps them to the given genome using <i>Minimap2</i>.
Then creates non-assembled annotation for each samples using <i>Stringtie (v3.0)</i> and look for the potential operons within those annotations (<i>GAMBA v1.3.2</i>). After identifing operon transcripts, operon-conteined transcripts, and non-operon-related transcriptsit, it generates a consensus annotations for the three sets of transcirpts.

Those sets are merge together to generate the final consensus annotation. From the annotation no containg operons there is selected the longest annotation for each gene (<i>Longest_transcript_filter.py</i>) to generate a fasta file, which is use to asses the quality of the annotation with <i>BUSCO (v5.8)</i>. The transcriptome obtained from the consensus annotations are also assed with <i>BUSCO</i>.

Schematic pipeline:

<img width="695" height="1551" alt="Figure1_new" src="https://github.com/user-attachments/assets/1fe688f6-d507-4a5a-a065-45dd8b265965" />

## Installation
This pipeline is build on _Snakemake_, therefore, you will need to have _Snakemake_ installed (tested on v5.24.1).

Source and binary packages for execute this pipeline can be directly downloaded from the Releases page on this repository.
However, due to the presence of a submodule we recomend to download it using the following:

    git clone --recursive https://github.com/EvoDevoGenomics-UB/CORAL.git

## How to run

To run **CORAL** you just need to modify the 'CORAL-config.yaml' with your parameters and execute it as any other _Snakemake_ file. We recomend to use _conda_ so it will create an environment where install all the dependecies specified in the 'CORAL-env.yml' and 'CORAL-env.merg.yml' files. The command will be like:

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
