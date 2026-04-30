# Rust tool GAMBA
rule build_GAMBA:
    output:
        "gamba"
    log:
        "logs/log_build_GAMBA.log"
    conda:
        env_file
    params:
        snakedir=SNAKEDIR,
        workdir=WORKDIR
    shell:
        """
    (pwd
    cd {params.snakedir}/scripts/gamba-tool
    cargo build --release
    cd {params.workdir}
    cp -r {params.snakedir}/scripts/gamba-tool/target/release/{output} {params.workdir}
    pwd
    ./{output} --help ) 2>&1 | tee {log}
    """


rule run_GAMBA_and_sanatizing:
    input:
        GAMBA=rules.build_GAMBA.output,
        gtf=ancient(rules.run_stringtie_sample_annotations.output.gtf)
    output:
        file="GAMBA_results/{specie}/{specie}_{sample}_guide{ref}_v{intron}_operons_found_t{threshold}.tsv",
        gtfOPRNs="GAMBA_results/{specie}/{specie}_{sample}_guide{ref}_v{intron}_Operons_t{threshold}.clean.gtf",
        gtfOpGs="GAMBA_results/{specie}/{specie}_{sample}_guide{ref}_v{intron}_OperonGenes_t{threshold}.clean.gtf",
        gtfCLEAN="GAMBA_results/{specie}/{specie}_{sample}_guide{ref}_v{intron}_opCLEAN_t{threshold}.clean.gtf"
    log:
        "logs/{specie}/{specie}_{sample}_guide{ref}_v{intron}_gambat{threshold}_GAMBA_run.log"
    conda:
        env_file
    threads: config["threads"]
    params:
        threshold=config["operon_threshold"],
        threshold2=float(config["operon_threshold"])*3
    shell:
        """
    gtf_name=$(basename {input.gtf} ".gtf")
    echo "GTF file: $gtf_name"
    echo "Command: ./{input.GAMBA} -f {input.gtf} --threshold {params.threshold} --monoexonic-t {params.threshold2} -o "GAMBA_results/{wildcards.specie}" --log {log}"

    ./{input.GAMBA} -f {input.gtf} --threshold {params.threshold} --monoexonic-t {params.threshold2} -o "GAMBA_results/{wildcards.specie}" --log {log}
    
    awk \'{{if($4>$5) print $1,$2,$3,$5,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18 ; \
    else print $0}}\' GAMBA_results/{wildcards.specie}/${{gtf_name}}_Operons_t{params.threshold}.gtf > ${{gtf_name}}_Operons_t{params.threshold}.tmp ; \
    gffread --sort-alpha -F -T -o {output.gtfOPRNs} ${{gtf_name}}_Operons_t{params.threshold}.tmp

    awk \'{{if($4>$5) print $1,$2,$3,$5,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18 ; \
    else print $0}}\' GAMBA_results/{wildcards.specie}/${{gtf_name}}_OperonGenes_t{params.threshold}.gtf > ${{gtf_name}}_OperonGenes_t{params.threshold}.tmp ; \
    gffread --sort-alpha -F -T -o {output.gtfOpGs} ${{gtf_name}}_OperonGenes_t{params.threshold}.tmp

    awk \'{{if($4>$5) print $1,$2,$3,$5,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18 ; \
    else print $0}}\' GAMBA_results/{wildcards.specie}/${{gtf_name}}_opCLEAN_t{params.threshold}.gtf > ${{gtf_name}}_opCLEAN_t{params.threshold}.tmp ; \
    gffread --sort-alpha -F -T -o {output.gtfCLEAN} ${{gtf_name}}_opCLEAN_t{params.threshold}.tmp

    rm ${{gtf_name}}*{params.threshold}.tmp
    """


rule run_gCLEAN_filter:
    input:
        gtf = "GAMBA_results/{specie}/{specie}_{sample}_guide{ref}_v{intron}_opCLEAN_t{threshold}.clean.gtf",
        bam=ancient("alignments/{specie}/{specie}_{sample}_reads_aln_v{intron}.sorted.bam")
    output:
        gtf = "GAMBA_results/{specie}/{specie}_{sample}_guide{ref}_v{intron}_opCLEAN_t{threshold}.clean.filter.gtf"
    log:
        "GAMBA_results/{specie}/{specie}_{sample}_guide{ref}_v{intron}_opCLEAN_t{threshold}.clean.filter.log"
    conda:
        env_file
    threads: config["threads"]
    params:
        freq=config["stringtie_freq"],
        g_param='-50',
        opts=config["stringtie_merge_opts"]
    shell: """
       (stringtie --merge -p {threads} \
        -c 2 -F 1 -T 0.5 -m 200 -g {params.g_param} -f {params.freq} \
        -l {wildcards.sample}g -o {output.gtf}.tmp {input.gtf} ; \
        stringtie -G {output.gtf}.tmp -e -o {output.gtf} {input.bam} ) 2>&1 | tee {log}
        rm {output.gtf}.tmp
    """