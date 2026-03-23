## Alignments of the reads by minimap2 & samtools
rule check_genome_format:
    input:
        genome = in_genome
    output:
        genome = temp("index/{specie}_genome.fasta")
    conda: env_file
    log: "logs/log_{specie}_genome_format.log"
    shell:"""
    (if ls {input.genome} | grep -q '.gz' ; then
        seqkit seq {input.genome} -o {output.genome}
    else
        ln -sf {input.genome} {output.genome} 
    fi ) 2> {log}
    """

rule build_minimap_index:
    input:
        genome = rules.check_genome_format.output.genome
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
        fq = temp("tmp/{sample}_concatenated.fastq"),
        flq = "processed_reads/{sample}_full_length_reads.clean.fq"
    params:
        pc = config["run_pychopper"],
        pc_opts = config["pychopper_opts"],
        seqkitQ = config["seqkitQ"],
        min_seqkit = config["min_seqkit"]
    conda: env_file
    log:
        log1 = "logs/log_prepare_fastqs_{sample}.log",
        log2 = "logs/{sample}_pychopper_log.txt"
    shell: """
    (
    mkdir -p tmp
    mkdir -p processed_reads
    if [ -s {output.flq} ] ; then 
        echo "{output.flq} already exists, skiped..."
    else
        if [ $(echo {input} | wc -w) -ge 2 ]; then
            echo "Concatenating FASTQ files for {wildcards.sample}..."
            seqkit seq {input} -m 50 -o {output.fq}
        else
            echo "Single FASTQ for {wildcards.sample}, copying..."
            #ln -sf {input} {output.fq}
            seqkit seq {input} -m 50 -o {output.fq}
        fi 
        if [[ "{params.pc}" == "True" ]]; then
            echo "Performing pychopper of {wildcards.sample}..."
            mkdir -p processed_reads/Reports/
            ( head -3 {output.fq} ; \
            pychopper -r processed_reads/Reports/{wildcards.sample}_report.pdf \
            -S processed_reads/Reports/{wildcards.sample}_statistics.tsv \
            -u processed_reads/{wildcards.sample}_unclassified.fq \
            -w processed_reads/{wildcards.sample}_rescued.fq \
            {params.pc_opts} {output.fq} processed_reads/{wildcards.sample}_full_length_reads.fq \
            ) 2> {log.log2}
            seqkit seq -Q {params.seqkitQ} -m {params.min_seqkit} \
            processed_reads/{wildcards.sample}_full_length_reads.fq processed_reads/{wildcards.sample}_rescued.fq > {output.flq}
            gzip processed_reads/{wildcards.sample}_full_length_reads.fq
            gzip processed_reads/{wildcards.sample}_unclassified.fq
            gzip processed_reads/{wildcards.sample}_rescued.fq
        else
            seqkit seq -Q {params.seqkitQ} -m {params.min_seqkit} {output.fq} > {output.flq}
        fi
    fi ) 2> {log.log1}
    """

rule run_minimap2:
    input:
        index = rules.build_minimap_index.output,
        fastq = rules.prepare_fastqs.output.flq
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
