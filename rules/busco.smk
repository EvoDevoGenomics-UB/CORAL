##Busco-related rules
rule run_longest_trans_filter:
    input:
        gtf=ancient(rules.run_recover_coverage.output.gtfFinal)
    output:
        filtergtf="annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.counts_longest_trans_only.gtf"
    log:
        "logs/{specie}/log_long_trans_filter_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.counts.log"
    conda:
        env_file
    params:
        snakedir=SNAKEDIR
    shell:
        """
    (python {params.snakedir}/scripts/Longest_transcript_filter.py {input.gtf}
    touch -c {output.filtergtf} ) 2> {log}
    """


rule run_obtaining_fasta:
    input:
        genome=ancient(rules.check_genome_format.output.genome),
        gtf=rules.run_longest_trans_filter.output.filtergtf,
        gtf_noOPRNs=rules.run_final_annotation.output.noOPRNs,
        gtf_andORPNs=rules.run_final_annotation.output.andOPRNs
    output:
        fasta="busco_analysis/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs_longest_trans_only.fasta",
        fasta_noOPRNs="busco_analysis/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs.fasta",
        fasta_andOPRNs="busco_analysis/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_andOPRNs.fasta"
    log:
        "logs/{specie}/log_obtaining_fasta_GTFs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    conda:
        env_file
    shell:
        """
    (mkdir -p busco_analysis
    gffread -g {input.genome} -w {output.fasta} {input.gtf}
    gffread -g {input.genome} -w {output.fasta_noOPRNs} {input.gtf_noOPRNs}
    gffread -g {input.genome} -w {output.fasta_andOPRNs} {input.gtf_andORPNs} ) 2> {log}
    """

rule run_fasta_reference_annot:
    input:
        genome=ancient(rules.check_genome_format.output.genome),
        ref_annot=REF
    output:
        fasta="busco_analysis/{specie}/{specie}_LRannot_REF.fasta"
    log:
        "logs/{specie}/log_fasta_reference_annot_{specie}.log"
    conda:
        env_file
    shell:
        """ (
        mkdir -p busco_analysis
        gffread -g {input.genome} -w {output.fasta} {input.ref_annot}
        ) 2>&1 | tee {log}
    """

rule busco_download_lineage:
    output:
        lin_dir=directory(path.join("busco_downloads/lineages/", config["lineages"])),
        lin_dir_root=directory("busco_downloads/")
    log:
        "logs/log_busco_download_lineage.log"
    conda:
        env_file
    params:
        lineage=config["lineages"]
    shell:
        """
        (
        busco --download {params.lineage} --download_path {output.lin_dir_root}
        ls -l {output.lin_dir} )  2>&1 | tee {log}
    """


rule run_busco_analyses:
    input:
        lin_dir=ancient(rules.busco_download_lineage.output.lin_dir),
        fasta="busco_analysis/{specie}/{specie}_{filename}.fasta"
    output:
        outdir=directory(
            "busco_analysis/{specie}/BUSCO_trans_{specie}_{filename}"
        ),
        summary=expand("busco_analysis/{{specie}}/BUSCO_trans_{{specie}}_{{filename}}/short_summary.specific.{lin}.BUSCO_trans_{{specie}}_{{filename}}.txt", lin=config["lineages"])
    log:
        "logs/{specie}/log_BUSCO_trans_{specie}_{filename}.log"
    conda:
        env_file
    params:
        lineage=config["lineages"]
    shell:
        """ (
        busco -i {input.fasta} -l {input.lin_dir} -o {output.outdir} -m transcriptome -f
         ) 2> {log}
    """


busco_ref_input = []
if config["reference_annot"] not in (None, [], ""):
    busco_ref_input.append("busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_REF")


rule busco_plot:
    input:
        out_longtrans="busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs_longest_trans_only",
        out_noOPRNs="busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs",
        out_andOPRNs="busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_andOPRNs",
        out_ref=busco_ref_input
    output:
        out_dir=directory(
            "busco_analysis/{specie}/BUSCO_results_all_summaries_{specie}_guide{ref}_v{intron}_gambat{threshold}"
        )
    log:
        "logs/{specie}/log_busco_polt_BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    conda:
        env_file
    params:
        workdir=WORKDIR,
        snakedir=SNAKEDIR
    shell:
        """ (
        mkdir -p {output.out_dir}
        cp {params.workdir}{input.out_longtrans}/short_summary.*.txt {output.out_dir}/short_summary.specific.metazoa_odb10.noOPRNs_longtrans.txt
        cp {params.workdir}{input.out_noOPRNs}/short_summary.*.txt {output.out_dir}/short_summary.specific.metazoa_odb10.noOPRNs.txt
        if [ -s {params.workdir}{input.out_andOPRNs}/short_summary.*.txt ] ; then
            cp {params.workdir}{input.out_andOPRNs}/short_summary.*.txt {output.out_dir}/short_summary.specific.metazoa_odb10.andOPRNs.txt ; \
        fi
        if [ {input.out_ref} != "" ]; then
            cp {params.workdir}{input.out_ref}/short_summary.*.txt {output.out_dir}/short_summary.specific.metazoa_odb10.REF.txt
        fi
        python3 {params.snakedir}/scripts/generate_plot.py -wd {output.out_dir} ) 2> {log}
    """
