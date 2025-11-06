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

def validate_samplesheet(samples_df):
    """Validate sample sheet and return dict of sample -> list of FASTQ paths."""
    grouped = defaultdict(list)
    seen_paths = defaultdict(list)
    for _, row in samples_df.iterrows():
        sample = str(row.iloc[0])
        fq = str(row.iloc[1])
        if not os.path.exists(fq):
            print(f"ERROR: Missing FASTQ file for sample '{sample}': {fq}", file=sys.stderr)
            sys.exit(1)
        grouped[sample].append(fq)
        seen_paths[fq].append(sample)
    # Check for duplicates within a sample
    for sample, fqs in grouped.items():
        if len(fqs) != len(set(fqs)):
            print(f"ERROR: Sample '{sample}' lists the same FASTQ file multiple times.", file=sys.stderr)
            sys.exit(1)
    # Check for FASTQs reused across samples
    reused = {fq: s for fq, s in seen_paths.items() if len(set(s)) > 1}
    if reused:
        print("ERROR: Some FASTQ files are used by multiple samples:", file=sys.stderr)
        for fq, s in reused.items():
            print(f"  {fq}  ->  samples: {', '.join(set(s))}", file=sys.stderr)
        sys.exit(1)
    print(f"[INFO] Loaded {len(grouped)} samples, all FASTQ paths validated.")
    return grouped

# Load and validate samples
if "samplesheet" in config and config["samplesheet"]:
    samples_df = pd.read_csv(config["samplesheet"], sep="\t", header=None)
    SAMPLE_TO_FASTQ = validate_samplesheet(samples_df)
    SAMPLES = list(SAMPLE_TO_FASTQ.keys())
else:
    SAMPLES = config["samples"]
    SAMPLE_TO_FASTQ = {
        s: [os.path.join(config["data_dir"], f"{s}{config['data_sufix']}")]
        for s in SAMPLES
    }

rule_all_input_list=["versions.txt","operon-finder",
        expand("logs/{specie}_{sample}_stats_input_reads.txt", specie=config["specie"], sample=SAMPLES),
        expand("alignments/{specie}_{sample}_reads_aln_v{intron}.sorted.bam", specie=config["specie"], sample=SAMPLES, intron=config["minimap2_max_intron"]),
        expand("alignments/{specie}_{sample}_reads_aln_sorted_v{intron}.stats.txt", specie=config["specie"], sample=SAMPLES, intron=config["minimap2_max_intron"]),
        expand("sample_annotations/{specie}_{sample}_guide{ref}_v{intron}.gtf", specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"]),
        expand("GAMBA_results/{specie}_guide{ref}_{sample}_v{intron}_opCLEAN_t{threshold}.clean.gtf", specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-and-OPRNs.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.counts.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.fasta", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.fasta", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-and-OPRNs.fasta", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_andOPRNs", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/BUSCO_results_all_summaries_{specie}_guide{ref}_v{intron}_gambat{threshold}", specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"])]

if config["reference_annot"] not in (None, [], ""):
    rule_all_input_list.append(expand("busco_analysis/BUSCO_trans_{specie}_LRannot_REF", specie=config["specie"]))
    if config["run_gffcomapre"] == True :
        rule_all_input_list.append(expand("Gffcompare_results/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}",specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]))        

if config["run_expression_matrix"] == True :
    rule_all_input_list.append(expand("Expression_matrix/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-and-OPRNs/transcript_count_matrix.csv",specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]))        
                
rule all:
    input:
        rule_all_input_list

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

rule build_minimap_index:
    input:
        genome = in_genome
    output:
        index = "index/{specie}_genome_index.mmi"
    params:
        opts = config["minimap_index_opts"]
    conda: env_file
    log: "logs/log_minimap2_{specie}_genome_index.log"
    threads: config["threads"]
    shell:"""
        (minimap2 -t {threads} {params.opts} -I 1000G -d {output.index} {input.genome}
        samtools faidx {input.genome}) 2>&1 | tee {log}
    """

rule prepare_fastqs:
    input:
        fastq = lambda wc: SAMPLE_TO_FASTQ[wc.sample]
    output:
        fq = temp("tmp/{sample}_concatenated.fastq")
    conda: env_file
    log: "logs/log_prepare_fastqs_{sample}.log"
    shell: """
    (mkdir -p tmp
    if [ $(echo {input} | wc -w) -ge 2 ]; then
        echo "Concatenating FASTQ files for {wildcards.sample}..."
        seqkit seq {input} > {output.fq}
    else
        echo "Single FASTQ for {wildcards.sample}, copying..."
        ln -sf {input} {output.fq}
    fi )2> {log}
    """

rule run_minimap2:
    input:
        index = rules.build_minimap_index.output,
        fastq = rules.prepare_fastqs.output.fq
    output:
        sam = temp("alignments/{specie}_{sample}_reads_aln_v{intron}.sam")
    params:
        opts = config["minimap2_opts"],
        max_intron = config["minimap2_max_intron"]
    threads: config["threads"]
    conda: env_file
    log: "logs/{specie}_{sample}_v{intron}_minimap2_run.log"
    shell: """
    (minimap2 -t {threads} {params.opts} -G {params.max_intron} {input.index} {input.fastq} > {output.sam} ; \
    head -2 {output.sam} ) 2>&1 | tee {log}
    echo "Minimap2 alignment done: {output.sam} created"
    """

rule run_samtools:
    input:
        sam = rules.run_minimap2.output.sam,
        index = rules.build_minimap_index.output
    output:
        bam = protected("alignments/{specie}_{sample}_reads_aln_v{intron}.sorted.bam"),
        bai = "alignments/{specie}_{sample}_reads_aln_v{intron}.sorted.bam.bai"
    params:
        qual = config["minimum_mapping_quality"]
    threads: config["threads"]
    conda: env_file
    log: "logs/log_{specie}_{sample}_reads_aln_v{intron}_samtools.log"
    shell:"""
    (samtools view -h -bt {input.index} {input.sam} | seqkit bam -j {threads} -q {params.qual} -x -\
    | samtools sort -@ {threads} -O BAM -o {output.bam} -;
    echo \"BAM created: {output.bam}\"
    samtools index {output.bam}
    echo \"Index created: {output.bai}\" ) 2> {log}
    """

rule run_aling_stats:
    input:
        bam = "alignments/{specie}_{sample}_reads_aln_v{intron}.sorted.bam"
    output:
        stats="alignments/{specie}_{sample}_reads_aln_sorted_v{intron}.stats.txt"
    conda: env_file
    threads: config["threads"]
    log: "logs/log_{specie}_{sample}_reads_aln_sorted_v{intron}.stats.log"
    shell:"""
    (samtools stats -@ {threads} {input.bam} > {output.stats} ) 2> {log}
    """

# Stringite v3.0
rule do_stringtie_sample_annotations:
    input:
        expand("sample_annotations/{specie}_{sample}_guide{ref}_v{intron}.gtf", 
            specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"])

rule run_stringtie_sample_annotations:
    input: 
        bam = "alignments/{specie}_{sample}_reads_aln_v{intron}.sorted.bam"
    output:
        gtf = "sample_annotations/{specie}_{sample}_guide{ref}_v{intron}.gtf"
    params:
        opts = config["stringtie_opts"],
        strand = config["stringtie_strand"],
        ref_annot = REF
    conda: env_file
    log: "logs/log_stringtie_annotation_{specie}_{sample}_guide{ref}_v{intron}.log"
    threads: config["threads"]
    shell:"""
        (mkdir -p sample_annotations
        input_guide=\"{wildcards.ref}\"
        stringtie --version
        if [ $input_guide == "REF" ] ; then
            echo \"Comand: stringtie --fr -L -R -p {threads} {params.opts} -G {params.ref_annot} -o {output.gtf} {input.bam}\"
            stringtie {params.strand} -L -R -p {threads} {params.opts} -G {params.ref_annot} -o {output.gtf} {input.bam}
            echo \"Stringtie {wildcards.ref} guided gtf created: {output.gtf}\"
        else
            echo \"Comand: stringtie --fr -L -R -p {threads} {params.opts} -o {output.gtf} {input.bam}\"
            stringtie {params.strand} -L -R -p {threads} {params.opts} -o {output.gtf} {input.bam} ; \
            echo \"Stringtie no-guide no-assembly gtf created: {output.gtf}\"
        fi ) 2> {log}
        """

## Find and annotate operons and contained genes
rule do_operon_annotations:
    input:
        expand("GAMBA_results/{specie}_guide{ref}_{sample}_v{intron}_operons_found_t{threshold}.tsv", specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("GAMBA_results/{specie}_guide{ref}_{sample}_v{intron}_Operons_t{threshold}.clean.gtf", specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("GAMBA_results/{specie}_guide{ref}_{sample}_v{intron}_OperonGenes_t{threshold}.clean.gtf", specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("GAMBA_results/{specie}_guide{ref}_{sample}_v{intron}_opCLEAN_t{threshold}.clean.gtf", specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"])

#Rust tool GAMBA
rule run_GAMBA_and_sanatizing:
    input:
        GAMBA = rules.build_GAMBA.output,
        gtf = rules.run_stringtie_sample_annotations.output.gtf
    output:
        file = "GAMBA_results/{specie}_{sample}_guide{ref}_v{intron}_operons_found_t{threshold}.tsv",
        gtfOPRNs = "GAMBA_results/{specie}_{sample}_guide{ref}_v{intron}_Operons_t{threshold}.clean.gtf",
        gtfOpGs = "GAMBA_results/{specie}_{sample}_guide{ref}_v{intron}_OperonGenes_t{threshold}.clean.gtf",
        gtfCLEAN = "GAMBA_results/{specie}_{sample}_guide{ref}_v{intron}_opCLEAN_t{threshold}.clean.gtf"
    params:
        threshold = config["operon_threshold"]
    conda: env_file
    log: "logs/{specie}_guide{ref}_{sample}_v{intron}_gambat{threshold}_GAMBA_run.log"
    threads: config["threads"]
    shell:"""
    "sample_annotations/{specie}_{sample}_guide{ref}_v{intron}.gtf"
    gtf_name=$(basename {input.gtf} ".gtf")
    ./{input.GAMBA} -f {input.gtf} --threshold {params.threshold} -o "GAMBA_results" --log {log}
    
    awk \'{{if($4>$5) print $1,$2,$3,$5,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18 ; \
    else print $0}}\' GAMBA_results/${gtf_name}_Operons_t{params.threshold}.gtf > ${gtf_name}_Operons_t{params.threshold}.tmp ; \
    gffread --sort-alpha -F -T -o {output.gtfOPRNs} ${gtf_name}_Operons_t{params.threshold}.tmp

    awk \'{{if($4>$5) print $1,$2,$3,$5,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18 ; \
    else print $0}}\' GAMBA_results/${gtf_name}_OperonGenes_t{params.threshold}.gtf > ${gtf_name}_OperonGenes_t{params.threshold}.tmp ; \
    gffread --sort-alpha -F -T -o {output.gtfOpGs} ${gtf_name}_OperonGenes_t{params.threshold}.tmp

    awk \'{{if($4>$5) print $1,$2,$3,$5,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18 ; \
    else print $0}}\' GAMBA_results/${gtf_name}_opCLEAN_t{params.threshold}.gtf > ${gtf_name}_opCLEAN_t{params.threshold}.tmp ; \
    gffread --sort-alpha -F -T -o {output.gtfCLEAN} ${gtf_name}_opCLEAN_t{params.threshold}.tmp

    rm {params.name}*{params.threshold}.tmp
    """

gtfsoperons_samples=[]
gtfsopgenes_samples=[]
gtfsclean_samples=[]
for SAMPLE in config["samples"]:
    gtfsoperons_samples.append("GAMBA_results/{{specie}}_{}_guide{{ref}}_v{{intron}}_Operons_t{{threshold}}.clean.gtf".format(SAMPLE))
    gtfsopgenes_samples.append("GAMBA_results/{{specie}}_{}_guide{{ref}}_v{{intron}}_OperonGenes_t{{threshold}}.clean.gtf".format(SAMPLE))
    gtfsclean_samples.append("GAMBA_results/{{specie}}_{}_guide{{ref}}_v{{intron}}_opCLEAN_t{{threshold}}.clean.gtf".format(SAMPLE))

#Create oepron and operon-contained genes annotations
rule run_operon_annotation:
    input:
        gtfsoperons = gtfsoperons_samples,
        gtfsopgenes = gtfsopgenes_samples
    output:
        operongtf = "annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_OPRNs.gtf",
        opgenesgtf = "annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_OpGs.gtf",
        merge = "annotations/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted.gtf",
        def_file = "annotations/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.clean.gtf",
        opgenesgtfCLEAN = "annotations/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.OpGclean.gtf",
        excluded_file = "annotations/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.excluded.gtf"
    params:
        g_param = config["stringtie_g"],
        name = "{specie}_guide{ref}_v{intron}_gambat{threshold}",
        snakedir = SNAKEDIR
    log: 
        logOPRN = "logs/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.log",
        logSTRG = "logs/log_StrignTie_merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    conda: env_file
    shell:"""
    (mkdir -p annotations
    (for i in {input.gtfsopgenes} ; do echo $i ; done) > annotations/List_merge_OpGs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ; \
    stringtie --merge -l OpG -f 0 -F 0 -T 0 -c 0 -g '-60' -o {output.opgenesgtf} annotations/List_merge_OpGs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ; \
    grep 'StringTie	transcript' {output.opgenesgtf} | wc -l ; \
    (for i in {input.gtfsoperons} ; do echo $i ; done) > annotations/List_merge_OPRNs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ; \
    stringtie --merge -l OPRN -f 0 -F 0 -T 0 -c 0 -g 0 -o {output.operongtf} annotations/List_merge_OPRNs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ; \
    cat {output.operongtf} {output.opgenesgtf} > {params.name}.tmp.gtf ; \
    gffread --sort-alpha -F -T -o {output.merge} {params.name}.tmp.gtf ; rm {params.name}.tmp.gtf

    python {params.snakedir}/scripts/operon_validation.py -f {output.merge} --log {log.logOPRN}
    grep 'StringTie	transcript' {output.opgenesgtfCLEAN} | wc -l ; ) 2>&1 | tee {log.logSTRG}
    """

# Create gene final concensus annotations
rule run_gCLEAN_annotation:
    input:
        gtfsclean = gtfsclean_samples
    output:
        cleanfinal = "annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_mergeCLEAN.gtf"
    params:
        freq = config["stringtie_freq"],
        g_param = config["stringtie_g"],
        opts = config["stringtie_merge_opts"]
    conda: env_file
    log: "logs/log_StrignTie_merge_opCLEAN_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    shell:"""
    ((for i in {input.gtfsclean} ; do echo $i ; done ) > annotations/List_merge_opCLEAN.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ; \
    stringtie --version

    stringtie --merge annotations/List_merge_opCLEAN.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt \
     -l g -f {params.freq} {params.opts} -g {params.g_param} \
     -o {output.cleanfinal} ; \
    echo "  Final merge CLEAN done" ; \
    grep 'StringTie	transcript' {output.cleanfinal} | wc -l ) 2>&1 | tee {log}
    """

# Create Merge final concensus annotations
rule run_final_annotation:
    input:
        cleanfinal = rules.run_gCLEAN_annotation.output.cleanfinal ,
        mergegtf = rules.run_operon_annotation.output.def_file,
        excluded_file = rules.run_operon_annotation.output.excluded_file,
        opgenesgtf = rules.run_operon_annotation.output.opgenesgtfCLEAN
    output:
        noOPRNs = "annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.gtf" ,
        andOPRNs = "annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-and-OPRNs.gtf"
    params:
        freq = config["stringtie_freq"],
        g_param = config["stringtie_g"]
    conda: env_file2
    log: "logs/log_final_annotations_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    shell:"""
    (stringtie --merge {input.cleanfinal} {input.excluded_file} \
     -G {input.opgenesgtf} \
     -l g -f {params.freq} -F 0 -T 0 -c 0 -g {params.g_param} \
     -o {output.noOPRNs} ; \
    echo "  Final CLEAN-noOPRNs done" ; \
    grep 'StringTie	transcript' {output.noOPRNs} | wc -l ; \
    
    stringtie --merge {input.cleanfinal} {input.excluded_file} \
     -G {input.mergegtf} \
     -l g -f {params.freq} -F 0 -T 0 -c 0 -g {params.g_param} \
     -o {output.andOPRNs} ; \
    echo "  Final CLEAN-and-OPRNs done"
    grep 'StringTie	transcript' {output.andOPRNs} | wc -l ) 2>&1 | tee {log}
    """

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
        
rule run_longest_trans_filter:
    input:
        gtf = rules.run_final_annotation.output.noOPRNs
    output:
        filtergtf = "annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.gtf"
    params:
        snakedir = SNAKEDIR
    conda: env_file
    log: "logs/log_long_trans_filter_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.log"
    shell:"""
    (python {params.snakedir}/scripts/Longest_transcript_filter.py {input.gtf}
    touch {output.filtergtf} ) 2> {log}
    """

rule run_obtaining_fasta:
    input:
        genome = in_genome ,
        gtf = rules.run_longest_trans_filter.output.filtergtf ,
        gtf_noOPRNs = rules.run_final_annotation.output.noOPRNs ,
        gtf_andORPNs = rules.run_final_annotation.output.andOPRNs
    output:
        fasta = "busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.fasta" ,
        fasta_noOPRNs = "busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.fasta" ,
        fasta_andOPRNs = "busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-and-OPRNs.fasta"
    conda: env_file
    log: "logs/log_obtaining_fasta_GTFs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    shell:"""
    (mkdir -p busco_analysis
    gffread -g {input.genome} -w {output.fasta} {input.gtf}
    gffread -g {input.genome} -w {output.fasta_noOPRNs} {input.gtf_noOPRNs}
    gffread -g {input.genome} -w {output.fasta_andOPRNs} {input.gtf_andORPNs} ) 2> {log}
    """
rule busco_download_lineage:
    output:
        lin_dir = directory(path.join(WORKDIR, "busco_downloads/lineages/", config["lineages"]))
    params:
        lineage = config["lineages"]
    conda: env_file
    log: "logs/log_busco_download_lineage.log"
    shell:"""
        (busco --download {params.lineage} --download_path {output.lin_dir} ) 2> {log}
    """
rule run_busco_analyses:
    input:
        lin_dir = rules.busco_download_lineage.output.lin_dir ,
        fa_longtrans = rules.run_obtaining_fasta.output.fasta ,
        fa_noOPRNs = rules.run_obtaining_fasta.output.fasta_noOPRNs ,
        fa_andOPRNs = rules.run_obtaining_fasta.output.fasta_andOPRNs
    output:
        out_longtrans = directory("busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs_longest_trans_only"),
        out_noOPRNs = directory("busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs"),
        out_andOPRNs = directory("busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_andOPRNs")
    params:
        lineage = config["lineages"]
    conda: env_file
    log: "logs/log_BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    shell:""" (
        busco -i {input.fa_longtrans} -l {input.lin_dir} -o {output.out_longtrans} -m transcriptome
        busco -i {input.fa_noOPRNs} -l {input.lin_dir} -o {output.out_noOPRNs} -m transcriptome
        busco -i {input.fa_andOPRNs} -l {input.lin_dir} -o {output.out_andOPRNs} -m transcriptome ) 2> {log}
    """

rule run_busco_reference_annot:
    input:
        lin_dir = rules.busco_download_lineage.output.lin_dir ,
        genome = in_genome ,
        ref_annot = REF
    output:
        fasta = "busco_analysis/{specie}_LRannot_REF.fasta" ,
        out_ref = directory("busco_analysis/BUSCO_trans_{specie}_LRannot_REF")
    params:
        lineage = config["lineages"]
    conda: env_file
    log: "logs/log_busco_reference_annot_{specie}.log"
    shell:""" (
        mkdir -p busco_analysis
        gffread -g {input.genome} -w {output.fasta} {input.ref_annot}
        busco -i {output.fasta} -l {input.lin_dir} -o {output.out_ref} -m transcriptome ) 2> {log}
    """

busco_ref_input=[]
if config["reference_annot"] not in (None, [], ""):
    busco_ref_input.append("busco_analysis/BUSCO_trans_{specie}_LRannot_REF")

rule busco_plot:
    input:
        out_longtrans = "busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs_longest_trans_only",
        out_noOPRNs = "busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs",
        out_andOPRNs = "busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_andOPRNs",
        out_ref = busco_ref_input
    output:
        out_dir = directory("busco_analysis/BUSCO_results_all_summaries_{specie}_guide{ref}_v{intron}_gambat{threshold}")
    params:
        workdir = WORKDIR,
        snakedir = SNAKEDIR
    conda: env_file
    log: "logs/log_busco_polt_BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    shell:""" (
        mkdir -p {output.out_dir}
        cp {params.workdir}{input.out_longtrans}/short_summary.*.txt {output.out_dir}/short_summary.specific.metazoa_odb10.noOPRNs_longtrans.txt
        cp {params.workdir}{input.out_noOPRNs}/short_summary.*.txt {output.out_dir}/short_summary.specific.metazoa_odb10.noOPRNs.txt
        cp {params.workdir}{input.out_andOPRNs}/short_summary.*.txt {output.out_dir}/short_summary.specific.metazoa_odb10.andOPRNs.txt
        if [ {input.out_ref} != "" ]; then
            cp {params.workdir}{input.out_ref}/short_summary.*.txt {output.out_dir}/short_summary.specific.metazoa_odb10.REF.txt
        fi
        python3 {params.snakedir}/scripts/generate_plot.py -wd {output.out_dir} ) 2> {log}
    """

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
