##Busco-related rules
rule check_genome_format:
    input:
        genome = in_genome
    output:
        genome = temp("{specie}_genome.fasta")
    conda: env_file
    log: "logs/log_{specie}_genome_format.log"
    shell:"""
    (if ls {input.genome} | grep -q '.gz' ; then
        seqkit seq {input.genome} -o {output.genome}
    else
        ln -sf {input.genome} {output.genome} 
    fi ) 2> {log}
    """

rule run_longest_trans_filter:
    input:
        gtf = rules.run_final_annotation_part1.output.noOPRNs
    output:
        filtergtf = "annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.gtf"
    params:
        snakedir = SNAKEDIR
    conda: env_file
    log: "logs/{specie}/log_long_trans_filter_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.log"
    shell:"""
    (python {params.snakedir}/scripts/Longest_transcript_filter.py {input.gtf}
    touch {output.filtergtf} ) 2> {log}
    """

rule run_obtaining_fasta:
    input:
        genome = rules.check_genome_format.output.genome ,
        gtf = rules.run_longest_trans_filter.output.filtergtf ,
        gtf_noOPRNs = rules.run_final_annotation_part1.output.noOPRNs ,
        gtf_andORPNs = rules.run_final_annotation_part2.output.andOPRNs
    output:
        fasta = "busco_analysis/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.fasta" ,
        fasta_noOPRNs = "busco_analysis/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.fasta" ,
        fasta_andOPRNs = "busco_analysis/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-andOPRNs.fasta"
    conda: env_file
    log: "logs/{specie}/log_obtaining_fasta_GTFs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    shell:"""
    (mkdir -p busco_analysis
    gffread -g {input.genome} -w {output.fasta} {input.gtf}
    gffread -g {input.genome} -w {output.fasta_noOPRNs} {input.gtf_noOPRNs}
    gffread -g {input.genome} -w {output.fasta_andOPRNs} {input.gtf_andORPNs} ) 2> {log}
    """
rule busco_download_lineage:
    output:
        lin_dir = directory(path.join("busco_downloads/lineages/", config["lineages"]))
    params:
        lineage = config["lineages"]
    conda: env_file
    log: "logs/log_busco_download_lineage.log"
    shell:"""
        (busco --download {params.lineage} --download_path {output.lin_dir} 
        ls -l {output.lin_dir} ) 2> {log}
    """
rule run_busco_analyses:
    input:
        lin_dir = rules.busco_download_lineage.output.lin_dir ,
        fa_longtrans = rules.run_obtaining_fasta.output.fasta ,
        fa_noOPRNs = rules.run_obtaining_fasta.output.fasta_noOPRNs ,
        fa_andOPRNs = rules.run_obtaining_fasta.output.fasta_andOPRNs
    output:
        out_longtrans = directory("busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs_longest_trans_only"),
        out_noOPRNs = directory("busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs"),
        out_andOPRNs = directory("busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_andOPRNs")
    params:
        lineage = config["lineages"]
    conda: env_file
    log: "logs/{specie}/log_BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    shell:""" (
        busco -i {input.fa_longtrans} -l {input.lin_dir} -o {output.out_longtrans} -m transcriptome
        busco -i {input.fa_noOPRNs} -l {input.lin_dir} -o {output.out_noOPRNs} -m transcriptome
        busco -i {input.fa_andOPRNs} -l {input.lin_dir} -o {output.out_andOPRNs} -m transcriptome ) 2> {log}
    """

rule run_busco_reference_annot:
    input:
        lin_dir = rules.busco_download_lineage.output.lin_dir ,
        genome = rules.check_genome_format.output.genome ,
        ref_annot = REF
    output:
        fasta = "busco_analysis/{specie}/{specie}_LRannot_REF.fasta" ,
        out_ref = directory("busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_REF")
    params:
        lineage = config["lineages"]
    conda: env_file
    log: "logs/{specie}/log_busco_reference_annot_{specie}.log"
    shell:""" (
        mkdir -p busco_analysis
        gffread -g {input.genome} -w {output.fasta} {input.ref_annot}
        busco -i {output.fasta} -l {input.lin_dir} -o {output.out_ref} -m transcriptome ) 2>&1 | tee {log}
    """

busco_ref_input=[]
if config["reference_annot"] not in (None, [], ""):
    busco_ref_input.append("busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_REF")

rule busco_plot:
    input:
        out_longtrans = "busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs_longest_trans_only",
        out_noOPRNs = "busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs",
        out_andOPRNs = "busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_andOPRNs",
        out_ref = busco_ref_input
    output:
        out_dir = directory("busco_analysis/{specie}/BUSCO_results_all_summaries_{specie}_guide{ref}_v{intron}_gambat{threshold}")
    params:
        workdir = WORKDIR,
        snakedir = SNAKEDIR
    conda: env_file
    log: "logs/{specie}/log_busco_polt_BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
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
