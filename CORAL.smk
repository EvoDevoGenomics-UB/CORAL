import os
from os import path
import pandas as pd
from collections import defaultdict
import sys

from snakemake.utils import min_version
min_version("5.24")

configfile: path.join(path.dirname(workflow.snakefile),"CORAL-config.yaml")
workdir: path.join(config["workdir_top"], config["pipeline"])

WORKDIR = path.join(path.dirname(workflow.snakefile),config["workdir_top"], config["pipeline"])
SNAKEDIR = path.dirname(workflow.snakefile)
env_file = path.join(path.dirname(workflow.snakefile),"envs/CORAL-env.yml")
env_file2 = path.join(path.dirname(workflow.snakefile),"envs/CORAL-env.merge.yml")

in_genome = config["genome_fasta"]
REF = config["reference_annot"]

include: "{SNAKEDIR}/rules/common.smk"
include: "{SNAKEDIR}/rules/alignment.smk"
include: "{SNAKEDIR}/rules/sample-annot.smk"
include: "{SNAKEDIR}/rules/gamba.smk"
include: "{SNAKEDIR}/rules/busco.smk"

rule all:
    input:
        get_final_output()

rule dump_versions:
    log: "logs/versions.txt"
    conda: env_file
    shell: "command -v conda > /dev/null && conda list > {log}"

rule build_GAMBA:
    output: "gamba"
    params:
        snakedir = SNAKEDIR,
        workdir = WORKDIR
    conda: env_file
    log: "logs/log_build_GAMBA.log"
    shell: """
    (cd {params.snakedir}/scripts/gamba-tool
    cargo build --release
    cd {params.workdir}
    cp -r {params.snakedir}/scripts/gamba-tool/target/release/{output} {params.workdir}
    {params.workdir}/{output} --help ) 2>&1 | tee {log}
    """

## Check input files
rule input_files_stats:
    input:
        fastq = lambda wc: SAMPLE_TO_FASTQ[wc.sample]
    log:
        file="logs/{specie}_{sample}_stats_input_reads.txt"
    conda: env_file
    shell:"""
        ( seqkit stats {input.fastq} ) 2>&1 | tee {log.file}
    """

## Alignments of the reads by minimap2 & samtools
rule do_alignment:
    input:
        expand("index/{specie}_genome_index.mmi", specie=config["specie"]),
        expand("alignments/{specie}_{sample}_reads_aln_v{intron}.sorted.bam",
            specie=config["specie"], sample=SAMPLES, intron=config["minimap2_max_intron"]),
        expand("alignments/{specie}_{sample}_reads_aln_sorted_v{intron}.stats.txt", 
            specie=config["specie"], sample=SAMPLES, intron=config["minimap2_max_intron"])


## Sample annotations
rule do_stringtie_sample_annotations:
    input:
        expand("sample_annotations/{specie}_{sample}_guide{ref}_v{intron}.gtf", 
            specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"])


## Find and annotate operons and contained genes
rule do_operon_annotations:
    input:
        expand("GAMBA_results/{specie}_guide{ref}_{sample}_v{intron}_operons_found_t{threshold}.tsv", specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("GAMBA_results/{specie}_guide{ref}_{sample}_v{intron}_Operons_t{threshold}.clean.gtf", specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("GAMBA_results/{specie}_guide{ref}_{sample}_v{intron}_OperonGenes_t{threshold}.clean.gtf", specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("GAMBA_results/{specie}_guide{ref}_{sample}_v{intron}_opCLEAN_t{threshold}.clean.gtf", specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"])


## BUSCO-related rules
rule do_busco_analyses:
    input:
        expand("annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.gtf",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.fasta",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.fasta",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-and-OPRNs.fasta",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_andOPRNs",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs_longest_trans_only",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_andOPRNs",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/BUSCO_results_all_summaries_{specie}_guide{ref}_v{intron}_gambat{threshold}",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"])
        

# Obtaining coverage of final annotation
rule run_recover_coverage:
    input:
        gtf = rules.run_final_annotation.output.andOPRNs ,
        gtf2 = rules.run_final_annotation.output.noOPRNs ,
        bams = expand("alignments/{specie}_{sample}_reads_aln_v{intron}.sorted.bam", 
            specie=config["specie"], sample=SAMPLES, intron=config["minimap2_max_intron"])
    output:
        gtfFinal = "annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-and-OPRNs.counts.gtf",
        gtfFinal2 = "annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.counts.gtf"
    conda: env_file
    log: "logs/log_recover_coverage_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    shell: """ (
    stringtie -G {input.gtf2} -e -o {output.gtfFinal2} {input.bams}
    stringtie -G {input.gtf} -e -o {output.gtfFinal} {input.bams} ) 2> {log}
    """

#### OPTIONAL steps
## Comparing new annotations againts reference one
rule do_gffcompare:
    input:
        expand("Gffcompare_results/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"])

rule run_gffcompare:
    input:
        ref = REF ,
        gtf_longest = rules.run_longest_trans_filter.output.filtergtf ,
        gtf_noOPRNs = rules.run_final_annotation.output.noOPRNs ,
        gtf_andOPRNs = rules.run_final_annotation.output.andOPRNs
    output:
        gffcmp_dir = directory("Gffcompare_results/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}")
    params:
        prefix = "{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}"
    conda: env_file
    log: "logs/log_gffcomapre_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    shell:""" (
    if [[ ! -f {input.ref} ]] ; then
        echo \"Error: No reference annotation provided.\" >&2
        exit 1
    else
        mkdir -p {output.gffcmp_dir}
        gffcompare -r {input.ref} {input.gtf_longest} -o ./{output.gffcmp_dir}/{params.prefix}_longest_trans
        gffcompare -r {input.ref} {input.gtf_noOPRNs} -o ./{output.gffcmp_dir}/{params.prefix}_noOPRNs
        gffcompare -r {input.ref} {input.gtf_andOPRNs} -o ./{output.gffcmp_dir}/{params.prefix}_andOPRNs

        echo 'Performing counting of transcript types...'
        for i in ./{output.gffcmp_dir}/{params.prefix}_longest_trans ./{output.gffcmp_dir}/{params.prefix}_noOPRNs ./{output.gffcmp_dir}/{params.prefix}_andOPRNs ; do
            (for x in "=" c k m n j e o s x i y p r u ; do
                count=$(awk -v a="$x" '{{if($4==a) print $5,$4}}' $i.tracking | wc -l)
                echo "$x $count"
            done) > $i.gffcmp_trans_types.txt
            echo "File '$i.gffcmp_trans_types.txt' created."
        done
    fi ) 2>&1 | tee {log}
    """

## Expression matrix creation
rule run_expression_matrix:
    input:
        gtf = rules.run_final_annotation.output.andOPRNs ,
        bams = expand("alignments/{{specie}}_{sample}_reads_aln_v{{intron}}.sorted.bam", sample=SAMPLES)
    output:
        out_file_g = "Expression_matrix/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-and-OPRNs/gene_count_matrix.csv",
        out_file_t = "Expression_matrix/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-and-OPRNs/transcript_count_matrix.csv"
    params:
        result_dir = directory("Expression_matrix/{specie}"),
        samples = config["samples"],
        length = config["length"],
        snakedir = SNAKEDIR
    conda: env_file
    log:
        log1 = "logs/run_expression_matrix_part1.{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log",
        log2 = "logs/run_expression_matrix_part2.{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    threads: config["threads"]
    shell:"""
    mkdir -p "{params.result_dir}"
    
    samplelist=$(python {params.snakedir}/scripts/StringTie_counts.py \
     -f {input.gtf} -b {input.bams} --outdir {params.result_dir} -t {threads}  --log {log.log1})

    ( echo "Create final matrix with all counts"
    python {params.snakedir}/scripts/prepDE.py3 -l "{params.length}" -i "$samplelist" -g {output.out_file_g} -t {output.out_file_t} 
    [[ -f {output.out_file_t} ]] && echo "Expression matrix created succsesfully!" ) 2>&1 | tee {log.log2}
    """
