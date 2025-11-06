#Rust tool GAMBA
rule run_GAMBA_and_sanatizing:
    input:
        GAMBA = rules.build_GAMBA.output,
        gtf = rules.run_stringtie_sample_annotations.output.gtf
    output:
        file = "GAMBA_results/{specie}_{sample}_guide{ref}_v{intron}_operons_found_t{threshold}.tsv",
        gtfOPRNs = "GAMBA_results/{specie}_{sample}_guide{ref}_v{intron}_Operons_t{threshold}.clean.gtf",
        gtfOpGs = "GAMBA_results/{specie}_{sample}_guide{ref}_v{intron}_OperonGenes_t{threshold}.clean.gtf",
        gtfCLEAN = "GAMBA_results/{specie}_{sample}_guide{ref}_v{intron}_opCLEAN_t{threshold}.clean.gtf"
    params:
        threshold = config["operon_threshold"]
    conda: env_file
    log: "logs/{specie}_guide{ref}_{sample}_v{intron}_gambat{threshold}_GAMBA_run.log"
    threads: config["threads"]
    shell:"""
    "sample_annotations/{specie}_{sample}_guide{ref}_v{intron}.gtf"
    gtf_name=$(basename {input.gtf} ".gtf")
    ./{input.GAMBA} -f {input.gtf} --threshold {params.threshold} -o "GAMBA_results" --log {log}
    
    awk \'{{if($4>$5) print $1,$2,$3,$5,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18 ; \
    else print $0}}\' GAMBA_results/${gtf_name}_Operons_t{params.threshold}.gtf > ${gtf_name}_Operons_t{params.threshold}.tmp ; \
    gffread --sort-alpha -F -T -o {output.gtfOPRNs} ${gtf_name}_Operons_t{params.threshold}.tmp

    awk \'{{if($4>$5) print $1,$2,$3,$5,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18 ; \
    else print $0}}\' GAMBA_results/${gtf_name}_OperonGenes_t{params.threshold}.gtf > ${gtf_name}_OperonGenes_t{params.threshold}.tmp ; \
    gffread --sort-alpha -F -T -o {output.gtfOpGs} ${gtf_name}_OperonGenes_t{params.threshold}.tmp

    awk \'{{if($4>$5) print $1,$2,$3,$5,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18 ; \
    else print $0}}\' GAMBA_results/${gtf_name}_opCLEAN_t{params.threshold}.gtf > ${gtf_name}_opCLEAN_t{params.threshold}.tmp ; \
    gffread --sort-alpha -F -T -o {output.gtfCLEAN} ${gtf_name}_opCLEAN_t{params.threshold}.tmp

    rm {params.name}*{params.threshold}.tmp
    """

gtfsoperons_samples=[]
gtfsopgenes_samples=[]
gtfsclean_samples=[]
for SAMPLE in config["samples"]:
    gtfsoperons_samples.append("GAMBA_results/{{specie}}_{}_guide{{ref}}_v{{intron}}_Operons_t{{threshold}}.clean.gtf".format(SAMPLE))
    gtfsopgenes_samples.append("GAMBA_results/{{specie}}_{}_guide{{ref}}_v{{intron}}_OperonGenes_t{{threshold}}.clean.gtf".format(SAMPLE))
    gtfsclean_samples.append("GAMBA_results/{{specie}}_{}_guide{{ref}}_v{{intron}}_opCLEAN_t{{threshold}}.clean.gtf".format(SAMPLE))

#Create oepron and operon-contained genes annotations
rule run_operon_annotation:
    input:
        gtfsoperons = gtfsoperons_samples,
        gtfsopgenes = gtfsopgenes_samples
    output:
        operongtf = "annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_OPRNs.gtf",
        opgenesgtf = "annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_OpGs.gtf",
        merge = "annotations/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted.gtf",
        def_file = "annotations/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.clean.gtf",
        opgenesgtfCLEAN = "annotations/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.OpGclean.gtf",
        excluded_file = "annotations/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.excluded.gtf"
    params:
        g_param = config["stringtie_g"],
        name = "{specie}_guide{ref}_v{intron}_gambat{threshold}",
        snakedir = SNAKEDIR
    log: 
        logOPRN = "logs/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.log",
        logSTRG = "logs/log_StrignTie_merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    conda: env_file
    shell:"""
    (mkdir -p annotations
    (for i in {input.gtfsopgenes} ; do echo $i ; done) > annotations/List_merge_OpGs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ; \
    stringtie --merge -l OpG -f 0 -F 0 -T 0 -c 0 -g '-60' -o {output.opgenesgtf} annotations/List_merge_OpGs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ; \
    grep 'StringTie	transcript' {output.opgenesgtf} | wc -l ; \
    (for i in {input.gtfsoperons} ; do echo $i ; done) > annotations/List_merge_OPRNs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ; \
    stringtie --merge -l OPRN -f 0 -F 0 -T 0 -c 0 -g 0 -o {output.operongtf} annotations/List_merge_OPRNs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ; \
    cat {output.operongtf} {output.opgenesgtf} > {params.name}.tmp.gtf ; \
    gffread --sort-alpha -F -T -o {output.merge} {params.name}.tmp.gtf ; rm {params.name}.tmp.gtf

    python {params.snakedir}/scripts/operon_validation.py -f {output.merge} --log {log.logOPRN}
    grep 'StringTie	transcript' {output.opgenesgtfCLEAN} | wc -l ; ) 2>&1 | tee {log.logSTRG}
    """

# Create gene final concensus annotations
rule run_gCLEAN_annotation:
    input:
        gtfsclean = gtfsclean_samples
    output:
        cleanfinal = "annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_mergeCLEAN.gtf"
    params:
        freq = config["stringtie_freq"],
        g_param = config["stringtie_g"],
        opts = config["stringtie_merge_opts"]
    conda: env_file
    log: "logs/log_StrignTie_merge_opCLEAN_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    shell:"""
    ((for i in {input.gtfsclean} ; do echo $i ; done ) > annotations/List_merge_opCLEAN.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ; \
    stringtie --version

    stringtie --merge annotations/List_merge_opCLEAN.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt \
     -l g -f {params.freq} {params.opts} -g {params.g_param} \
     -o {output.cleanfinal} ; \
    echo "  Final merge CLEAN done" ; \
    grep 'StringTie	transcript' {output.cleanfinal} | wc -l ) 2>&1 | tee {log}
    """

# Create Merge final concensus annotations
rule run_final_annotation:
    input:
        cleanfinal = rules.run_gCLEAN_annotation.output.cleanfinal ,
        mergegtf = rules.run_operon_annotation.output.def_file,
        excluded_file = rules.run_operon_annotation.output.excluded_file,
        opgenesgtf = rules.run_operon_annotation.output.opgenesgtfCLEAN
    output:
        noOPRNs = "annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.gtf" ,
        andOPRNs = "annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-and-OPRNs.gtf"
    params:
        freq = config["stringtie_freq"],
        g_param = config["stringtie_g"]
    conda: env_file2
    log: "logs/log_final_annotations_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    shell:"""
    (stringtie --merge {input.cleanfinal} {input.excluded_file} \
     -G {input.opgenesgtf} \
     -l g -f {params.freq} -F 0 -T 0 -c 0 -g {params.g_param} \
     -o {output.noOPRNs} ; \
    echo "  Final CLEAN-noOPRNs done" ; \
    grep 'StringTie	transcript' {output.noOPRNs} | wc -l ; \
    
    stringtie --merge {input.cleanfinal} {input.excluded_file} \
     -G {input.mergegtf} \
     -l g -f {params.freq} -F 0 -T 0 -c 0 -g {params.g_param} \
     -o {output.andOPRNs} ; \
    echo "  Final CLEAN-and-OPRNs done"
    grep 'StringTie	transcript' {output.andOPRNs} | wc -l ) 2>&1 | tee {log}
    """

# Obtaining coverage of final annotation
rule run_recover_coverage:
    input:
        gtf = rules.run_final_annotation.output.andOPRNs ,
        gtf2 = rules.run_final_annotation.output.noOPRNs ,
        bams = expand("alignments/{specie}_{sample}_reads_aln_v{intron}.sorted.bam", 
            specie=config["specie"], sample=SAMPLES, intron=config["minimap2_max_intron"])
    output:
        gtfFinal = "annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-and-OPRNs.counts.gtf",
        gtfFinal2 = "annotations/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.counts.gtf"
    conda: env_file
    log: "logs/log_recover_coverage_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    shell: """ (
    stringtie -G {input.gtf2} -e -o {output.gtfFinal2} {input.bams}
    stringtie -G {input.gtf} -e -o {output.gtfFinal} {input.bams} ) 2> {log}
    """
