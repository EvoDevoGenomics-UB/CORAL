import os
from os import path

configfile: path.join(path.dirname(workflow.snakefile),"compact-genome-annotation-config_NEWold.yaml")
workdir: path.join(config["workdir_top"], config["pipeline"])

WORKDIR = path.join(path.dirname(workflow.snakefile),config["workdir_top"], config["pipeline"])
SNAKEDIR = path.dirname(workflow.snakefile)
exeOpF = "python {SNAKEDIR}/scripts/operon_finder_v9.7.py"

in_genome = config["genome_fasta"]
env_file = "compact-genome-annotation-env.yml"
REF = config["reference_annot"]
exeOpF = "./operon-finder"

rule_all_input_list=["versions.txt","operon-finder",
        expand("logs/{specie}_{sample}_stats_input_reads.txt", specie=config["specie"], sample=config["samples"]),
        expand("alignments/{specie}_{sample}_reads_aln_v{intron}.sorted.bam", specie=config["specie"], sample=config["samples"], intron=config["minimap2_max_intron"]),
        expand("sample_annotations/{specie}_{sample}_guide{ref}_v{intron}.gtf", specie=config["specie"], sample=config["samples"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"]),
        expand("operon_finder_results/{specie}_guide{ref}_{sample}_v{intron}_opCLEAN_v9.t{threshold}.clean.gtf", specie=config["specie"], sample=config["samples"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}.sorted.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-and-OPRNs.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-noOPRNs.counts.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.fasta", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-noOPRNs.fasta", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-and-OPRNs.fasta", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_andOPRNs", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/BUSCO_results_all_summaries_{specie}_guide{ref}_v{intron}_OFv9t{threshold}", specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("logs/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_operon_finder_run_FINAL.log", specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]) ]

if config["run_gffcomapre"] == True :
    rule_all_input_list.append(expand("Gffcompare_results/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}",specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]))        
                
rule all:
    input:
        rule_all_input_list

rule dump_versions:
    output: ver = "versions.txt"
    conda: env_file
    shell: "command -v conda > /dev/null && conda list > {output.ver}"

rule build_operon_finder:
    output: "operon-finder"
    conda: env_file
    shell: """
    cd {SNAKEDIR}/scripts/operon-finder-rust
    cargo build --release
    mv {SNAKEDIR}/scripts/operon-finder-rust/target/release/{output} {SNAKEDIR}/scripts/
    {SNAKEDIR}/scripts/{output} --help
    """

## Check input files
rule input_files_stats:
    input:
        fastq = config["data_dir"] + "{sample}" + config["data_sufix"]
    output:
        file="logs/{specie}_{sample}_stats_input_reads.txt"
    conda: env_file
    shell:"""
        seqkit stats {input.fastq} -o {output.file}
    """

## Alignments of the reads by minimap2 & samtools
rule do_alignment:
    input:
        expand("index/{prefix}_genome_index.mmi", prefix=config["specie"]),
        expand("alignments/{specie}_{sample}_reads_aln_v{intron}.sorted.bam",
            specie=config["specie"], sample=config["samples"], intron=config["minimap2_max_intron"]),
        expand("alignments/{specie}_{sample}_reads_aln_sorted_v{intron}.stats.txt", 
            specie=config["specie"], sample=config["samples"], intron=config["minimap2_max_intron"])

rule build_minimap_index:
    input:
        genome = in_genome
    output:
        index = "index/{prefix}_genome_index.mmi"
    params:
        prefix = config["specie"],
        opts = config["minimap_index_opts"]
    conda: env_file
    threads: config["threads"]
    shell:"""
        minimap2 -t {threads} {params.opts} -I 1000G -d {output.index} {input.genome}
    """
rule run_minimap2:
    input:
        index = "index/{specie}_genome_index.mmi",
        fastq = config["data_dir"] + "{sample}" + config["data_sufix"]
    output:
        bam = "alignments/{specie}_{sample}_reads_aln_v{intron}.sorted.bam"
    params:
        name = "alignments/{specie}_{sample}_reads_aln_v{intron}",
        opts = config["minimap2_opts"],
        max_intron = config["minimap2_max_intron"],
        qual = config["minimum_mapping_quality"]
    threads: config["threads"]
    conda: env_file
    log: "logs/{specie}_{sample}_v{intron}_minimap2_run.log"
    shell: """
    (minimap2 -t {threads} {params.opts} -G {params.max_intron} {params.opts} {input.index} {input.fastq} > {params.name}.sam ; \
    head -2 {params.name}.sam ) 2> {log}
    echo \"{params.name}.sam created\"
    samtools view {params.name}.sam -@ {threads} -O BAM -o {params.name}.bam
    echo \"{params.name}.bam created\"
    seqkit bam -j {threads} -q {params.qual} -x {params.name}.bam > {params.name}.clean.sam
    echo \"{params.name}.clean.sam created\"
    rm {params.name}.sam
    samtools sort -@ {threads} -O BAM -o {output.bam} {params.name}.clean.sam
    samtools index {output.bam}
    rm {params.name}.clean.sam
    """
rule run_aling_stats:
    input:
        bam = "alignments/{specie}_{sample}_reads_aln_v{intron}.sorted.bam"
    output:
        stats="alignments/{specie}_{sample}_reads_aln_sorted_v{intron}.stats.txt"
    conda: env_file
    threads: config["threads"]
    shell:"""
    samtools stats -@ {threads} {input.bam} > {output.stats}
    """

#Stringite v3.0
rule do_stringtie_sample_annotations:
    input:
        expand("sample_annotations/{specie}_{sample}_guide{ref}_v{intron}.gtf", 
            specie=config["specie"], sample=config["samples"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"])

rule run_stringtie_sample_annotations:
    input: 
        bam = "alignments/{specie}_{sample}_reads_aln_v{intron}.sorted.bam"
    output:
        gtf = "sample_annotations/{specie}_{sample}_guide{ref}_v{intron}.gtf"
    params:
        opts = config["stringtie_opts"]
    conda: env_file
    threads: config["threads"]
    shell:"""
        mkdir -p sample_annotations
        input_guide=\"{wildcards.ref}\"
        if [ $input_guide == "REF" ] ; then
            echo \"Comand: stringtie --rf -L -R -p {threads} {params.opts} -G {REF} -o {output.gtf} {input.bam}\"
            stringtie --rf -L -R -p {threads} {params.opts} -G {REF} -o {output.gtf} {input.bam}
            echo \"Stringtie {wildcards.ref} guided gtf created: {output.gtf}\"
        else
            echo \"Comand: stringtie --rf -L -R -p {threads} {params.opts} -o {output.gtf} {input.bam}\"
            stringtie --rf -L -R -p {threads} {params.opts} -o {output.gtf} {input.bam} ; \
            echo \"Stringtie no-guide no-assembly gtf created: {output.gtf}\"
        fi
        """

##Find and annotate operons and contained genes
rule do_operon_annotations:
    input:
        expand("operon_finder_results/{specie}_guide{ref}_{sample}_v{intron}_operons_found_v9.t{threshold}.tsv", specie=config["specie"], sample=config["samples"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("operon_finder_results/{specie}_guide{ref}_{sample}_v{intron}_Operons_v9.t{threshold}.clean.gtf", specie=config["specie"], sample=config["samples"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("operon_finder_results/{specie}_guide{ref}_{sample}_v{intron}_OperonGenes_v9.t{threshold}.clean.gtf", specie=config["specie"], sample=config["samples"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("operon_finder_results/{specie}_guide{ref}_{sample}_v{intron}_opCLEAN_v9.t{threshold}.clean.gtf", specie=config["specie"], sample=config["samples"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}.sorted.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"])

#Python script operon_finder
rule run_operon_finder_and_sanatizing:
    input:
        operon_finder = rules.build_operon_finder.output,
        gtf = rules.run_stringtie_sample_annotations.output.gtf
    output:
        file = "operon_finder_results/{specie}_guide{ref}_{sample}_v{intron}_operons_found_v9.t{threshold}.tsv",
        gtfOPRNs = "operon_finder_results/{specie}_guide{ref}_{sample}_v{intron}_Operons_v9.t{threshold}.clean.gtf",
        gtfOpGs = "operon_finder_results/{specie}_guide{ref}_{sample}_v{intron}_OperonGenes_v9.t{threshold}.clean.gtf",
        gtfCLEAN = "operon_finder_results/{specie}_guide{ref}_{sample}_v{intron}_opCLEAN_v9.t{threshold}.clean.gtf"
    params:
        name = "operon_finder_results/{specie}_guide{ref}_{sample}_v{intron}",
        threshold = config["operon_threshold"]
    conda: env_file
    log: "logs/{specie}_guide{ref}_{sample}_v{intron}_OFv9t{threshold}_operon_finder_run.log"
    threads: config["threads"]
    shell:"""
    mkdir -p operon_finder_results
    {exeOpF} -f {input.gtf} --threshold {params.threshold} -o {params.name} --log {log}
    
    head -3 {output.gtfOPRNs}
    head -3 {output.gtfOpGs}
    head -3 {output.gtfCLEAN}
    """

gtfsoperons_samples=[]
gtfsopgenes_samples=[]
gtfsclean_samples=[]
for SAMPLE in config["samples"]:
    gtfsoperons_samples.append("operon_finder_results/{{specie}}_guide{{ref}}_{}_v{{intron}}_Operons_v9.t{{threshold}}.clean.gtf".format(SAMPLE))
    gtfsopgenes_samples.append("operon_finder_results/{{specie}}_guide{{ref}}_{}_v{{intron}}_OperonGenes_v9.t{{threshold}}.clean.gtf".format(SAMPLE))
    gtfsclean_samples.append("operon_finder_results/{{specie}}_guide{{ref}}_{}_v{{intron}}_opCLEAN_v9.t{{threshold}}.clean.gtf".format(SAMPLE))

#Create oepron and operon-contained genes annotations
rule run_operon_annotation:
    input:
        gtfsoperons = gtfsoperons_samples,
        gtfsopgenes = gtfsopgenes_samples
    output:
        operongtf = "annotations/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_OPRNs.gtf",
        opgenesgtf = "annotations/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_OpGs.gtf",
        merge = "annotations/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}.sorted.gtf",
        def_file = "annotations/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}.sorted_OPRNstatistics.clean.gtf"
    params:
        freq = config["stringtie_freq"],
        name = "{specie}_guide{ref}_v{intron}_OFv9t{threshold}"
    log: "logs/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}.sorted_OPRNstatistics.log"
    conda: env_file
    shell:"""
    mkdir -p annotations
    (for i in {input.gtfsopgenes} ; do echo $i ; done) > annotations/List_merge_OpGs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}OFv9t{wildcards.threshold}.txt ; \
    stringtie --merge -l OpG -f {params.freq} -F 0 -T 0 -c 0 -g '-50' -o {output.opgenesgtf} annotations/List_merge_OpGs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}OFv9t{wildcards.threshold}.txt ; \
    (for i in {input.gtfsoperons} ; do echo $i ; done) > annotations/List_merge_OPRNs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}OFv9t{wildcards.threshold}.txt ; \
    stringtie --merge -l OPRN -f {params.freq} -F 0 -T 0 -c 0 -g 0 -o {output.operongtf} annotations/List_merge_OPRNs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}OFv9t{wildcards.threshold}.txt ; \
    cat {output.operongtf} {output.opgenesgtf} > {params.name}.tmp.gtf ; \
    gffread --sort-alpha -F -T -o {output.merge} {params.name}.tmp.gtf ; rm {params.name}.tmp.gtf

    python3 {SNAKEDIR}/scripts/operon_statistics.py -f {output.merge} --log {log}
    """

#Create final concensus annotations
rule run_final_annotation:
    input:
        gtfsclean = gtfsclean_samples ,
        mergegtf = expand("annotations/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}.sorted.gtf",
            specie=config["specie"], ref=config["stringtie_guide_opts"], intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        opgenesgtf = expand("annotations/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_OpGs.gtf",
            specie=config["specie"], ref=config["stringtie_guide_opts"], intron=config["minimap2_max_intron"], threshold=config["operon_threshold"])
    output:
        cleanfinal = "annotations/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_mergeCLEAN.gtf" ,
        noOPRNs = "annotations/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-noOPRNs.gtf" ,
        andOPRNs = "annotations/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-and-OPRNs.gtf"
    params:
        freq = config["stringtie_freq"]
    conda: env_file
    shell:"""
    (for i in {input.gtfsclean} ; do echo $i ; done) > annotations/List_merge_opCLEAN.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}OFv9t{wildcards.threshold}.txt ; \
    stringtie --merge annotations/List_merge_opCLEAN.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}OFv9t{wildcards.threshold}.txt \
     -l g -f {params.freq} -F 0.5 -T 0 -c 10 -g '-75' \
     -o {output.cleanfinal} ; \
    echo "  Final merge CLEAN done" ; \
    stringtie --merge {output.cleanfinal} \
     -G {input.opgenesgtf} \
     -l g -f {params.freq} -F 0 -T 0 -c 0 -g '-75' \
     -o {output.noOPRNs} ; \
    echo "  Final CLEAN-noOPRNs done" ; \
    stringtie --merge {output.cleanfinal} \
     -G {input.mergegtf} \
     -l g -f {params.freq} -F 0 -T 0 -c 0 -g '-75' \
     -o {output.andOPRNs} ; \
    echo "  Final CLEAN-and-OPRNs done"
    """

##Extra things
rule do_busco_analyses:
    input:
        expand("annotations/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.gtf",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.fasta",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-noOPRNs.fasta",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-and-OPRNs.fasta",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_andOPRNs",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_noOPRNs_longest_trans_only",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_noOPRNs",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_andOPRNs",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/BUSCO_results_all_summaries_{specie}_guide{ref}_v{intron}_OFv9t{threshold}",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"])
        
rule run_longest_trans_filter:
    input:
        gtf = rules.run_final_annotation.output.noOPRNs
    output:
        filtergtf = "annotations/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.gtf"
    conda: env_file
    shell:"""
    python {SNAKEDIR}/scripts/Longest_transcript_filter.py {input.gtf}
    touch {output.filtergtf}
    """

rule run_obtaining_fasta:
    input:
        genome = in_genome ,
        gtf = rules.run_longest_trans_filter.output.filtergtf ,
        gtf_noOPRNs = rules.run_final_annotation.output.noOPRNs ,
        gtf_andORPNs = rules.run_final_annotation.output.andOPRNs
    output:
        fasta = "busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.fasta" ,
        fasta_noOPRNs = "busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-noOPRNs.fasta" ,
        fasta_andOPRNs = "busco_analysis/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-and-OPRNs.fasta"
    conda: env_file
    shell:"""
    mkdir -p busco_analysis
    gffread -g {input.genome} -w {output.fasta} {input.gtf}
    gffread -g {input.genome} -w {output.fasta_noOPRNs} {input.gtf_noOPRNs}
    gffread -g {input.genome} -w {output.fasta_andOPRNs} {input.gtf_andORPNs}
    """
rule busco_download_lineage:
    output:
        lin_dir = directory(path.join(WORKDIR, "busco_downloads/lineages/", config["lineages"]))
    params:
        lineage = config["lineages"]
    conda: env_file
    shell:"""
        busco --download {params.lineage} --download_path {output.lin_dir}
    """
rule run_busco_analyses:
    input:
        lin_dir = rules.busco_download_lineage.output.lin_dir ,
        fa_longtrans = rules.run_obtaining_fasta.output.fasta ,
        fa_noOPRNs = rules.run_obtaining_fasta.output.fasta_noOPRNs ,
        fa_andOPRNs = rules.run_obtaining_fasta.output.fasta_andOPRNs
    output:
        out_longtrans = directory("busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_noOPRNs_longest_trans_only"),
        out_noOPRNs = directory("busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_noOPRNs"),
        out_andOPRNs = directory("busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_andOPRNs")
    params:
        lineage = config["lineages"]
    conda: env_file
    shell:"""
        busco -i {input.fa_longtrans} -l {input.lin_dir} -o {output.out_longtrans} -m transcriptome
        busco -i {input.fa_noOPRNs} -l {input.lin_dir} -o {output.out_noOPRNs} -m transcriptome
        busco -i {input.fa_andOPRNs} -l {input.lin_dir} -o {output.out_andOPRNs} -m transcriptome
    """

rule busco_plot:
    input:
        out_longtrans = "busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_noOPRNs_longest_trans_only",
        out_noOPRNs = "busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_noOPRNs",
        out_andOPRNs = "busco_analysis/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_andOPRNs"
    output:
        out_dir = directory("busco_analysis/BUSCO_results_all_summaries_{specie}_guide{ref}_v{intron}_OFv9t{threshold}")
    conda: env_file
    shell:"""
        mkdir -p {output.out_dir}
        cp {WORKDIR}{input.out_longtrans}/short_summary.*.txt {output.out_dir}/short_summary.specific.metazoa_odb10.noOPRNs_longtrans.txt
        cp {WORKDIR}{input.out_noOPRNs}/short_summary.*.txt {output.out_dir}/short_summary.specific.metazoa_odb10.noOPRNs.txt
        cp {WORKDIR}{input.out_andOPRNs}/short_summary.*.txt {output.out_dir}/short_summary.specific.metazoa_odb10.andOPRNs.txt
        python3 {SNAKEDIR}/scripts/generate_plot.py -wd {output.out_dir}
        """

####FINAL steps
#Obtaining coverage of final annotation
rule run_recover_coverage:
    input:
        gtf = rules.run_final_annotation.output.andOPRNs ,
        gtf2 = rules.run_final_annotation.output.noOPRNs ,
        bams = expand("alignments/{specie}_{sample}_reads_aln_v{intron}.sorted.bam", 
            specie=config["specie"], sample=config["samples"], intron=config["minimap2_max_intron"])
    output:
        gtfFinal = "annotations/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-and-OPRNs.counts.gtf",
        gtfFinal2 = "annotations/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-noOPRNs.counts.gtf"
    conda: env_file
    shell: """
    stringtie -G {input.gtf2} -e -o {output.gtfFinal2} {input.bams}
    stringtie -G {input.gtf} -e -o {output.gtfFinal} {input.bams}
    """
rule run_final_operon_search:
    input:
        gtf = rules.run_recover_coverage.output.gtfFinal
    output:
        dir_name = directory("operon_finder_results/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_DEF")
    params:
        name = "{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_StringtieMerge.clean-and-OPRNs.counts",
        threshold = config["operon_threshold"]
    log: "logs/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}_operon_finder_run_FINAL.log"
    conda: env_file
    shell:"""
    mkdir -p {output.dir_name}
    {exeOpF} -f {input.gtf} --threshold {params.threshold} -o {output.dir_name}/{params.name} --log {log}
    """

##Comparing new annoatation againts reference one
rule do_gffcompare:
    input:
        expand("Gffcompare_results/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}",
            specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"])

rule run_gffcompare:
    input:
        ref = config["reference_annot"] ,
        gtf_longest = rules.run_longest_trans_filter.output.filtergtf ,
        gtf_noOPRNs = rules.run_final_annotation.output.noOPRNs ,
        gtf_andOPRNs = rules.run_final_annotation.output.andOPRNs
    output:
        gffcmp_dir = directory("Gffcompare_results/{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}")
    params:
        prefix = "{specie}_LRannot_guide{ref}_v{intron}_OFv9t{threshold}"
    conda: env_file
    shell:"""
    if [[ {input.ref} == "" ]] ; then
        echo \"Error: No reference annoation provided.\" >&2
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
    fi
    """
