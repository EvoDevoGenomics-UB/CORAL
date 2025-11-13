import os
from os import path
import pandas as pd
from collections import defaultdict
import sys

from snakemake.utils import min_version
min_version("5.24")

configfile: path.join(path.dirname(workflow.snakefile),"CORAL-config.yaml")
#workdir: path.join(config["workdir_top"], config["pipeline"])

WORKDIR = path.join(path.dirname(workflow.snakefile),config["workdir_top"], config["pipeline"])
SNAKEDIR = path.dirname(workflow.snakefile)
env_file = path.join(path.dirname(workflow.snakefile),"envs/CORAL-env.yml")
env_file2 = path.join(path.dirname(workflow.snakefile),"envs/CORAL-env.merge.yml")

in_genome = config["genome_fasta"]
REF = config["reference_annot"]

include: "rules/common.smk"
include: "rules/alignment.smk"
include: "rules/sample-annot.smk"
include: "rules/gamba.smk"
include: "rules/consensus.smk"
include: "rules/busco.smk"
include: "rules/gffcmp.smk"
include: "rules/exp-matrix.smk"

rule all:
    input:
        get_final_output()

rule dump_versions:
    log: "logs/versions.txt"
    conda: env_file
    shell: "command -v conda > /dev/null && conda list > {log}"

## Check input files
rule input_files_stats:
    input:
        fastq = lambda wc: SAMPLE_TO_FASTQ[wc.sample]
    log:
        file="logs/{sample}_stats_input_reads.txt"
    conda: env_file
    shell:"""
        ( seqkit stats {input.fastq} ) 2>&1 | tee {log.file}
    """

#### Specific rules
## Alignments of the reads by minimap2 & samtools
rule do_alignment:
    input:
        expand("index/{specie}_genome_index.mmi", specie=config["specie"]),
        expand("alignments/{specie}/{specie}_{sample}_reads_aln_v{intron}.sorted.bam",
            specie=config["specie"], sample=SAMPLES, intron=config["minimap2_max_intron"]),
        expand("alignments/{specie}/{specie}_{sample}_reads_aln_sorted_v{intron}.stats.txt", 
            specie=config["specie"], sample=SAMPLES, intron=config["minimap2_max_intron"])


## Sample annotations
rule do_stringtie_sample_annotations:
    input:
        expand("sample_annotations/{specie}/{specie}_{sample}_guide{ref}_v{intron}.gtf", 
            specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"])


## Find and annotate operons and contained genes
rule do_operon_annotations:
    input:
        expand("GAMBA_results/{specie}/{specie}_{sample}_guide{ref}_v{intron}_operons_found_t{threshold}.tsv", specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("GAMBA_results/{specie}/{specie}_{sample}_guide{ref}_v{intron}_Operons_t{threshold}.clean.gtf", specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("GAMBA_results/{specie}/{specie}_{sample}_guide{ref}_v{intron}_OperonGenes_t{threshold}.clean.gtf", specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("GAMBA_results/{specie}/{specie}_{sample}_guide{ref}_v{intron}_opCLEAN_t{threshold}.clean.gtf", specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/{specie}/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"])

## Find and annotate operons and contained genes
rule do_consensus_annotations:
    input:
        expand("annotations/{specie}/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-andOPRNs.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.counts.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"])

## BUSCO-related rules
rule do_busco_analyses:
    input:
        expand("annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.gtf",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.fasta",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.fasta",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-andOPRNs.fasta",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_andOPRNs",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs_longest_trans_only",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_andOPRNs",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}/BUSCO_results_all_summaries_{specie}_guide{ref}_v{intron}_gambat{threshold}",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"])
        
#### OPTIONAL steps
## Comparing new annotations againts reference one
rule do_gffcompare:
    input:
        expand("Gffcompare_results/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}",
            specie=config["specie"], ref=config["stringtie_guide_opts"], intron=config["minimap2_max_intron"], threshold=config["operon_threshold"])

## Expression matrix creation
rule do_expression_matrix:
    input:
        expand("Expression_matrix/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-andOPRNs/gene_count_matrix.csv", 
            specie=config["specie"], ref=config["stringtie_guide_opts"], intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("Expression_matrix/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-andOPRNs/transcript_count_matrix.csv",
            specie=config["specie"], ref=config["stringtie_guide_opts"], intron=config["minimap2_max_intron"], threshold=config["operon_threshold"])
