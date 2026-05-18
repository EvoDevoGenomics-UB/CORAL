rule mmseqs2_databases:
    output:
        db_dir=directory("TD2_results/mmseq2_DBs/"),
        db_alias=expand(
            "TD2_results/mmseq2_DBs/{db_alias}", db_alias=config["db_alias"]
        )
    log:
        "logs/log_mmseq2_DBs.log"
    conda:
        env_file
    params:
        db_name=config["db_name"],
        db_alias=config["db_alias"]
    shell:
        """
    mkdir -p TD2_results
    mkdir -p {output.db_dir}

    (if [[ -f {output.db_alias} ]] ;  then
        echo "Database {params.db_alias} already exists, skipping download..."
        echo "Database {params.db_alias} already exists, skipping download..." >&2
    else
        mmseqs databases {params.db_name} {output.db_alias} tmp
    fi ) 2> {log}
    """


rule run_TransDecoder:
    input:
        genome=rules.check_genome_format.output.genome,
        gtf="GTFs_inputs/{gtf_name}.gtf",
        db_alias=rules.mmseqs2_databases.output.db_alias
    output:
        gff3="TD2_results/{specie}/{gtf_name}.fasta.TD2.genome.gff3"
    log:
        log1="logs/log_TD2_{specie}_{gtf_name}.log",
        log2="logs/{specie}/log_gffcomapre_{gtf_name}_TD2.log"
    conda:
        env_file
    params:
        scriptsDIR=path.join(SNAKEDIR, "scripts"),
        prefix="{gtf_name}_TD2",
        TD2options=config["TD2options"],
        gtf_ref=REF
    shell:
        """
    mkdir -p TD2_results/
    mkdir -p TD2_results/{wildcards.specie}/
    ({params.scriptsDIR}/TransDecoder2_script.sh \
       {input.genome} \
       {input.gtf} \
    "TD2_results/{wildcards.specie}/" \
       {input.db_alias} \
       "{params.scriptsDIR}" \
       "{params.TD2options}" \
    ) 2>&1 | tee {log.log1}

    (mkdir -p Gffcompare_results
    mkdir -p ./Gffcompare_results/TD2_annotations
        gffcompare -r {input.gtf} {output.gff3} -o ./Gffcompare_results/TD2_annotations/{params.prefix}
        gffcompare -r {output.gff3} {input.gtf} -o ./Gffcompare_results/TD2_annotations/{params.prefix}revers
    if [[ -f "{params.gtf_ref}" ]] ; then
        gffcompare -r {params.gtf_ref} {output.gff3} -o ./Gffcompare_results/TD2_annotations/{params.prefix}vsREF
    fi
    echo 'Performing counting of transcript types...'
    for i in Gffcompare_results/TD2_annotations/{params.prefix} ./Gffcompare_results/TD2_annotations/{params.prefix}revers ./Gffcompare_results/TD2_annotations/{params.prefix}vsREF ; do
        if [[ -f "$i.tracking" ]] ; then
        (for x in "=" c k m n j e o s x i y p r u ; do
            count=$(awk -v a="$x" '{{if($4==a) print $5,$4}}' $i.tracking | wc -l)
            echo "$x $count"
        done) > $i.gffcmp_trans_types.txt
        echo "File '$i.gffcmp_trans_types.txt' created."
        fi
    done
    touch -c ./Gffcompare_results/TD2_annotations/{params.prefix}.annotated.gtf
    ) 2>&1 | tee {log.log2}
    """


rule run_obtaining_pep_TD2:
    input:
        genome=rules.check_genome_format.output.genome,
        gff_TD2=rules.run_TransDecoder.output.gff3
    output:
        pep_TD2="busco_analysis/{specie}/{gtf_name}.TD2.pep.fasta"
    log:
        "logs/{specie}/log_obtaining_pep_GTFs_{gtf_name}.log"
    conda:
        env_file
    shell:
        """
    (mkdir -p busco_analysis
    gffread -g {input.genome} -y {output.pep_TD2} {input.gff_TD2} ) 2> {log}
    """


rule run_busco_prot:
    input:
        lin_dir=ancient(rules.busco_download_lineage.output.lin_dir),
        pep_TD2=rules.run_obtaining_pep_TD2.output.pep_TD2
    output:
        out_pep_TD2=directory("busco_analysis/{specie}/BUSCO_prot_{gtf_name}")
    log:
        "logs/{specie}/log_BUSCO_prot_{gtf_name}.log"
    conda:
        env_file
    shell:
        """ (
        busco -i {input.pep_TD2} -l {input.lin_dir} -o {output.out_pep_TD2} -m proteins) 2> {log}
    """


rule run_busco_prot_ref:
    input:
        lin_dir=ancient(rules.busco_download_lineage.output.lin_dir),
        genome=ancient(rules.check_genome_format.output.genome),
        gtf_ref=REF
    output:
        pep_ref="busco_analysis/{specie}/{specie}_LRannot_REF.pep.fasta",
        out_ref=directory("busco_analysis/{specie}/BUSCO_prot_{specie}_LRannot_REF")
    log:
        "logs/{specie}/log_busco_prot_REF_{specie}.log"
    conda:
        env_file
    params:
        lineage=config["lineages"]
    shell:
        """ (
        mkdir -p busco_analysis
        gffread -g {input.genome} -y {output.pep_ref} {input.gtf_ref}
        busco -i {output.pep_ref} -l {input.lin_dir} -o {output.out_ref} -m proteins ) 2>&1 | tee {log}
    """
