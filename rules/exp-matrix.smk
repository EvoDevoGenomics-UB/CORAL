## Expression matrix creation
rule run_expression_matrix:
    input:
        #gtf = rules.run_final_annotation.output.andOPRNs ,
        gtf = rules.run_gffcompare.output.gffcmp_out,
        bams = expand("alignments/{{specie}}/{{specie}}_{sample}_reads_aln_v{{intron}}.sorted.bam", sample=SAMPLES)
    output:
        out_file_g = "Expression_matrix/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs.annotated/gene_count_matrix.csv",
        out_file_t = "Expression_matrix/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs.annotated/transcript_count_matrix.csv"
    params:
        result_dir = lambda wildcards, output: os.path.dirname(os.path.dirname(output[0])),
        samples = config["samples"],
        length = config["length"],
        snakedir = SNAKEDIR
    conda: env_file
    log:
        log1 = "logs/{specie}/run_expression_matrix_part1.{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log",
        log2 = "logs/{specie}/run_expression_matrix_part2.{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    threads: config["threads"]
    shell:"""
    mkdir -p "{params.result_dir}"
    
    samplelist=$(python {params.snakedir}/scripts/StringTie_counts.py \
     -f {input.gtf} -b {input.bams} --outdir {params.result_dir} -t {threads}  --log {log.log1})

    ( echo "Create final matrix with all counts"
    python {params.snakedir}/scripts/prepDE.py3 -l "{params.length}" -i "$samplelist" -g {output.out_file_g} -t {output.out_file_t} 
    [[ -f {output.out_file_t} ]] && echo "Expression matrix created succsesfully!" ) 2>&1 | tee {log.log2}
    """

## Expression matrix creation
rule run_expression_matrix_REF:
    input:
        #gtf = rules.run_final_annotation.output.andOPRNs ,
        gtf = REF ,
        bams = expand("alignments/{{specie}}/{{specie}}_{sample}_reads_aln_v{{intron}}.sorted.bam", sample=SAMPLES)
    output:
        out_file_g = "Expression_matrix/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_REF/gene_count_matrix.csv",
        out_file_t = "Expression_matrix/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_REF/transcript_count_matrix.csv"
    params:
        result_dir = lambda wildcards, output: os.path.dirname(os.path.dirname(output[0])),
        samples = config["samples"],
        length = config["length"],
        snakedir = SNAKEDIR
    conda: env_file
    log:
        log1 = "logs/{specie}/run_expression_matrix_part1.{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.REF.log",
        log2 = "logs/{specie}/run_expression_matrix_part2.{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.REF.log"
    threads: config["threads"]
    shell:"""
    mkdir -p "{params.result_dir}"
    
    samplelist=$(python {params.snakedir}/scripts/StringTie_counts.py \
     -f {input.gtf} -b {input.bams} --outdir {params.result_dir} -t {threads}  --log {log.log1})

    ( echo "Create final matrix with all counts"
    python {params.snakedir}/scripts/prepDE.py3 -l "{params.length}" -i "$samplelist" -g {output.out_file_g} -t {output.out_file_t} 
    [[ -f {output.out_file_t} ]] && echo "Expression matrix created succsesfully!" ) 2>&1 | tee {log.log2}
    """
