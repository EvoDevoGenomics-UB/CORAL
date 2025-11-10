## Alignments of the reads by minimap2 & samtools

rule build_minimap_index:
    input:
        genome = in_genome
    output:
        index = "index/{specie}_genome_index.mmi"
    params:
        opts = config["minimap_index_opts"]
    conda: env_file
    log: "logs/{specie}/log_minimap2_{specie}_genome_index.log"
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
    fi ) 2> {log}
    """

rule run_minimap2:
    input:
        index = rules.build_minimap_index.output,
        fastq = rules.prepare_fastqs.output.fq
    output:
        sam = temp("alignments/{specie}/{specie}_{sample}_reads_aln_v{intron}.sam")
    params:
        opts = config["minimap2_opts"],
        max_intron = config["minimap2_max_intron"]
    threads: config["threads"]
    conda: env_file
    log: "logs/{specie}/{specie}_{sample}_v{intron}_minimap2_run.log"
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
        bam = protected("alignments/{specie}/{specie}_{sample}_reads_aln_v{intron}.sorted.bam"),
        bai = "alignments/{specie}/{specie}_{sample}_reads_aln_v{intron}.sorted.bam.bai"
    params:
        qual = config["minimum_mapping_quality"]
    threads: config["threads"]
    conda: env_file
    log: "logs/{specie}/log_{specie}_{sample}_reads_aln_v{intron}_samtools.log"
    shell:"""
    (samtools view -h -bt {input.index} {input.sam} | seqkit bam -j {threads} -q {params.qual} -x -\
    | samtools sort -@ {threads} -O BAM -o {output.bam} -;
    echo \"BAM created: {output.bam}\"
    samtools index {output.bam}
    echo \"Index created: {output.bai}\" ) 2> {log}
    """

rule run_aling_stats:
    input:
        bam = "alignments/{specie}/{specie}_{sample}_reads_aln_v{intron}.sorted.bam"
    output:
        stats="alignments/{specie}/{specie}_{sample}_reads_aln_sorted_v{intron}.stats.txt"
    conda: env_file
    threads: config["threads"]
    log: "logs/{specie}/log_{specie}_{sample}_reads_aln_sorted_v{intron}.stats.log"
    shell:"""
    (samtools stats -@ {threads} {input.bam} > {output.stats} ) 2> {log}
    """
