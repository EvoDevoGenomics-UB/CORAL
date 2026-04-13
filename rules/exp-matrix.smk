SCRIPTDIR = path.join(SNAKEDIR,"scripts")

## Expression matrix creation
rule run_expression_matrix:
    input:
        gtf = ancient("GTFs_inputs/{gtf_name}.gtf"),
        bams = ancient(expand("alignments/{{specie}}/{{specie}}_{sample}_reads_aln_v{{intron}}.sorted.bam", sample=SAMPLES))
    output:
        out_file_g = "Expression_matrix/{specie}/{gtf_name}/gene_count_matrix_v{intron}.csv",
        out_file_t = "Expression_matrix/{specie}/{gtf_name}/transcript_count_matrix_v{intron}.csv"
    params:
        result_dir = lambda wildcards, output: os.path.dirname(os.path.dirname(output[0])),
        length = config["length"],
        scriptsdir = SCRIPTDIR
    conda: env_file
    log:
        log1 = "logs/{specie}/run_expression_matrix_part1.{gtf_name}_v{intron}.log",
        log2 = "logs/{specie}/run_expression_matrix_part2.{gtf_name}_v{intron}.log"
    threads: config["threads"]
    shell:"""
    mkdir -p "{params.result_dir}"
    
    samplelist=$(python {params.scriptsdir}/StringTie_counts.py \
     -f {input.gtf} -b {input.bams} --outdir {params.result_dir} -t {threads}  --log {log.log1})

    ( echo "Create final matrix with all counts"
    python {params.scriptsdir}/prepDE.py3 -l "{params.length}" -i "$samplelist" -g {output.out_file_g} -t {output.out_file_t} 
    [[ -f {output.out_file_t} ]] && echo "Expression matrix created succsesfully!" ) 2>&1 | tee {log.log2}
    """

## Expression matrix creation
rule run_expression_matrix_REF:
    input:
        gtf = REF ,
        bams = ancient(expand("alignments/{{specie}}/{{specie}}_{sample}_reads_aln_v{{intron}}.sorted.bam", sample=SAMPLES))
    output:
        out_file_g = "Expression_matrix/{specie}/ref_annotation/gene_count_matrix_v{intron}.csv",
        out_file_t = "Expression_matrix/{specie}/ref_annotation/transcript_count_matrix_v{intron}.csv"
    params:
        scriptsdir = SCRIPTDIR ,
        result_dir = lambda wildcards, output: os.path.dirname(os.path.dirname(output[0])),
        length = config["length"]
    conda: env_file
    log:
        log1 = "logs/{specie}/run_expression_matrix_part1.ref_annotation_v{intron}.log",
        log2 = "logs/{specie}/run_expression_matrix_part2.ref_annotation_v{intron}.log"
    threads: config["threads"]
    shell:"""
    mkdir -p "{params.result_dir}"
    
    samplelist=$(python {params.scriptsdir}/StringTie_counts.py \
     -f {input.gtf} -b {input.bams} --outdir {params.result_dir} -t {threads}  --log {log.log1})

    ( echo "Create final matrix with all counts"
    python {params.scriptsdir}/prepDE.py3 -l "{params.length}" -i "$samplelist" -g {output.out_file_g} -t {output.out_file_t} 
    [[ -f {output.out_file_t} ]] && echo "Expression matrix created succsesfully!" ) 2>&1 | tee {log.log2}
    """
