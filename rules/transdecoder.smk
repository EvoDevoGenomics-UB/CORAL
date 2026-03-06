rule mmseqs2_databases:
    output:
        db_dir = directory("TD2_results/mmseq2_DBs/"),
        db_alias = expand("TD2_results/mmseq2_DBs/{db_alias}", db_alias=config["db_alias"])
    params:
        db_name = config["db_name"],
        db_alias = config["db_alias"]
    conda: env_file
    log: "logs/log_mmseq2_DBs.log"
    shell: """
    mkdir -p TD2_results
    mkdir -p {output.db_dir}

    (if [[ -f {output.db_alias} ]] ;  then
        echo "Database {params.db_alias} already exists, skipping download..."
        echo "Database {params.db_alias} already exists, skipping download..." >&2
    else
        mmseqs databases {params.db_name} {output.db_alias} tmp
    fi ) 2> {log}
    """

rule TransDecoder:
    input:
        genome = rules.check_genome_format.output.genome,
        gtf = rules.run_final_annotation.output.noOPRNs,
        swissprot = rules.mmseqs2_databases.output.swissprot
    output:
        outdir = directory("TD2_results/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs_files"),
        gff3 = "TD2_results/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.fasta.TD2.genome.gff3"
    params:
        snakedir = SNAKEDIR
    conda: env_file
    log: "logs/log_TD2_{specie}_guide{ref}_v{intron}_gambat{threshold}.log"
    shell: """
    mkdir -p TD2_results/
    mkdir -p TD2_results/{wildcards.specie}/
    ( {params.snakedir}/scripts/TransDecoder2_script.sh {input.genome} {input.gtf} \
    "TD2_results/{wildcards.specie}/" \
    {input.swissprot} "{params.snakedir}/scripts" \
    ) 2>&1 | tee {log}
    """

rule run_obtaining_pep_TD2:
    input:
        genome = rules.check_genome_format.output.genome,
        gff_TD2 = rules.TransDecoder.output.gff3 
    output:
        pep_TD2 = "busco_analysis/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.TD2.pep.fasta"
    conda: env_file
    log: "logs/{specie}/log_obtaining_pep_GTFs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    shell:"""
    (mkdir -p busco_analysis
    gffread -g {input.genome} -y {output.pep_TD2} {input.gff_TD2} ) 2> {log}
    """

rule run_busco_prot:
    input:
        lin_dir = rules.busco_download_lineage.output.lin_dir ,
        pep_TD2 = rules.run_obtaining_pep_TD2.output.pep_TD2
    output:
        out_pep_TD2 = directory("busco_analysis/{specie}/BUSCO_prot_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs")
    conda: env_file
    log: "logs/{specie}/log_BUSCO_prot_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    shell:""" (
        /Users/nuriat/anaconda3/envs/busco_env/bin/busco -i {input.pep_TD2} -l {input.lin_dir} -o {output.out_pep_TD2} -m proteins) 2> {log}
    """

rule run_busco_prot_ref:
    input:
        lin_dir = rules.busco_download_lineage.output.lin_dir,
        genome = rules.check_genome_format.output.genome,
        gtf_ref = REF
    output:
        pep_ref = "busco_analysis/{specie}/{specie}_LRannot_REF.pep.fasta",
        out_ref = directory("busco_analysis/{specie}/BUSCO_prot_{specie}_LRannot_REF")
    params:
        lineage = config["lineages"]
    conda: env_file
    log: "logs/{specie}/log_busco_prot_REF_{specie}.log"
    shell:""" (
        mkdir -p busco_analysis
        gffread -g {input.genome} -y {output.pep_ref} {input.gtf_ref}
        /Users/nuriat/anaconda3/envs/busco_env/bin/busco -i {output.pep_ref} -l {input.lin_dir} -o {output.out_ref} -m proteins ) 2>&1 | tee {log}
    """