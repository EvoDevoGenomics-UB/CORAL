
rule run_stringtie_sample_annotations:
    input:
        bam=ancient("alignments/{specie}/{specie}_{sample}_reads_aln_v{intron}.sorted.bam")
    output:
        gtf="sample_annotations/{specie}/{specie}_{sample}_guide{ref}_v{intron}.gtf"
    log:
        "logs/{specie}/log_stringtie_annotation_{specie}_{sample}_guide{ref}_v{intron}.log"
    conda:
        env_file
    threads: config["threads"]
    params:
        opts=config["stringtie_opts"],
        strand=config["stringtie_strand"],
        ref_annot=REF
    shell:
        """
        (mkdir -p sample_annotations ; mkdir -p sample_annotations/{wildcards.specie} ;
        input_guide=\"{wildcards.ref}\"
        stringtie --version
        if [ $input_guide == "REF" ] ; then
            echo \"Comand: stringtie {params.strand} -L -R -p {threads} {params.opts} -G {params.ref_annot} -o {output.gtf} {input.bam}\"
            stringtie {params.strand} -L -R -p {threads} {params.opts} -G {params.ref_annot} -l {wildcards.sample}g -o {output.gtf} {input.bam}
            echo \"Stringtie {wildcards.ref} guided gtf created: {output.gtf}\"
        else
            echo \"Comand: stringtie {params.strand} -L -R -p {threads} {params.opts} -l {wildcards.sample}g -o {output.gtf} {input.bam}\"
            stringtie {params.strand} -L -R -p {threads} {params.opts} -o {output.gtf} {input.bam} ; \
            echo \"Stringtie no-guide no-assembly gtf created: {output.gtf}\"
        fi ) 2> {log}
        """
