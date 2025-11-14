gtfsoperons_samples=[]
gtfsopgenes_samples=[]
gtfsclean_samples=[]
for SAMPLE in SAMPLES:
    gtfsoperons_samples.append("GAMBA_results/{{specie}}/{{specie}}_{}_guide{{ref}}_v{{intron}}_Operons_t{{threshold}}.clean.gtf".format(SAMPLE))
    gtfsopgenes_samples.append("GAMBA_results/{{specie}}/{{specie}}_{}_guide{{ref}}_v{{intron}}_OperonGenes_t{{threshold}}.clean.gtf".format(SAMPLE))
    gtfsclean_samples.append("GAMBA_results/{{specie}}/{{specie}}_{}_guide{{ref}}_v{{intron}}_opCLEAN_t{{threshold}}.clean.gtf".format(SAMPLE))

#Create oepron and operon-contained genes annotations
rule run_operon_annotation:
    input:
        gtfsoperons = gtfsoperons_samples,
        gtfsopgenes = gtfsopgenes_samples
    output:
        operongtf = temp("annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_OPRNs.gtf"),
        opgenesgtf = temp("annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_OpGs.gtf"),
        merge = "annotations/{specie}/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted.gtf",
        def_file = "annotations/{specie}/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.clean.gtf",
        opgenesgtfCLEAN = "annotations/{specie}/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.OpGclean.gtf",
        oprngtfCLEAN = "annotations/{specie}/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.OPRNclean.gtf",
        excluded_file = "annotations/{specie}/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.excluded.gtf",
        db_file = temp("annotations/{specie}/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.db")
    params:
        g_param = config["stringtie_OpGs_g"],
        name = "{specie}_guide{ref}_v{intron}_gambat{threshold}",
        snakedir = SNAKEDIR
    log: 
        logOPRN = "logs/{specie}/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.log",
        logSTRG = "logs/{specie}/log_StrignTie_merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    conda: env_file
    shell:"""
    (mkdir -p annotations ; mkdir -p annotations/{wildcards.specie} ;
    (for i in {input.gtfsopgenes} ; do echo $i ; done) > annotations/{wildcards.specie}/List_merge_OpGs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ; \
    stringtie --merge -l OpG -f 0 -F 0 -T 0 -c 0 -g {params.g_param} -o {output.opgenesgtf} annotations/{wildcards.specie}/List_merge_OpGs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ; \
    grep 'StringTie	transcript' {output.opgenesgtf} | wc -l ; \
    (for i in {input.gtfsoperons} ; do echo $i ; done) > annotations/{wildcards.specie}/List_merge_OPRNs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ; \
    stringtie --merge -l OPRN -f 0 -F 0 -T 0 -c 0 -g 0 -o {output.operongtf} annotations/{wildcards.specie}/List_merge_OPRNs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ; \
    cat {output.operongtf} {output.opgenesgtf} > {params.name}.tmp.gtf ; \
    gffread --sort-alpha -F -T -o {output.merge} {params.name}.tmp.gtf ; rm {params.name}.tmp.gtf

    python {params.snakedir}/scripts/operon_validation.py -f {output.merge} --log {log.logOPRN}
    grep 'StringTie	transcript' {output.opgenesgtfCLEAN} | wc -l ; ) 2>&1 | tee {log.logSTRG}
    """

# Create gene final consensus annotations
rule run_gCLEAN_annotation:
    input:
        gtfsclean = gtfsclean_samples,
        oprngtf = rules.run_operon_annotation.output.oprngtfCLEAN
    output:
        cleanfinal = "annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_mergeCLEAN.gtf"
    params:
        freq = config["stringtie_freq"],
        g_param = config["stringtie_g"],
        opts = config["stringtie_merge_opts"]
    conda: env_file
    log: "logs/{specie}/log_StrignTie_merge_opCLEAN_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    shell:"""
    ((for i in {input.gtfsclean} ; do echo $i ; done ) > annotations/{wildcards.specie}/List_merge_opCLEAN.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ; \
    stringtie --version

    stringtie --merge annotations/{wildcards.specie}/List_merge_opCLEAN.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt \
     -l g -f {params.freq} {params.opts} -g {params.g_param} \
     -o GTFfile.tmp ; \
    
    gffcompare -r {input.oprngtf} GTFfile.tmp -o filter
    awk '{{if($3=="=" ) print $5}}' filter.GTFfile.tmp.tmap > filter.list.tmp
    gffread --nids filter.list.tmp GTFfile.tmp -o {output.cleanfinal}
    rm filter.* ; rm GTFfile.tmp ;

    echo "  Final merge CLEAN done" ; \
    grep 'StringTie	transcript' {output.cleanfinal} | wc -l ) 2>&1 | tee {log}
    """

# Create Merge final consensus annotations
rule run_final_annotation:
    input:
        cleanfinal = rules.run_gCLEAN_annotation.output.cleanfinal,
        excluded_file = rules.run_operon_annotation.output.excluded_file,
        opgenesgtf = rules.run_operon_annotation.output.opgenesgtfCLEAN,
        oprnsgtf = rules.run_operon_annotation.output.oprngtfCLEAN
    output:
        noOPRNs = "annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.gtf",
        andOPRNs = "annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-andOPRNs.gtf"
    params:
        freq = config["stringtie_freq"],
        g_param = config["stringtie_g"],
        snakedir = SNAKEDIR
    conda: env_file
    log:
        log1 = "logs/{specie}/log_final_annotations_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_part1.log",
        log2 = "logs/{specie}/log_final_annotations_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_part2.log"
    shell:"""
    (stringtie --merge {input.cleanfinal} {input.excluded_file} \
     -G {input.opgenesgtf} \
     -l g -f {params.freq} -F 0 -T 0 -c 0 -g {params.g_param} \
     -o {output.noOPRNs} ; \
    echo "  Final CLEAN-noOPRNs done" ; \
    grep 'StringTie	transcript' {output.noOPRNs} | wc -l ; ) 2>&1 | tee {log.log1}

    (cat {input.oprnsgtf} {output.noOPRNs} > {output.andOPRNs}.1.tmp ; \
    gffread --sort-alpha -F -T -o {output.andOPRNs}.2.tmp {output.andOPRNs}.1.tmp ; 
    {params.snakedir}/scripts/add_operonID.py -f {output.andOPRNs}.2.tmp -o {output.andOPRNs} --log {output.andOPRNs}.log
    rm {output.andOPRNs}*.OPRNids.db {output.andOPRNs}*.tmp

    echo "  Final CLEAN-andOPRNs done"
    grep 'StringTie	transcript' {output.andOPRNs} | wc -l ) 2>&1 | tee {log.log2}
    """

# Obtaining coverage of final annotation
rule run_recover_coverage:
    input:
        gtf = rules.run_final_annotation.output.andOPRNs ,
        gtf2 = rules.run_final_annotation.output.noOPRNs ,
        bams = expand("alignments/{specie}/{specie}_{sample}_reads_aln_v{intron}.sorted.bam", 
            specie=config["specie"], sample=SAMPLES, intron=config["minimap2_max_intron"])
    output:
        gtfFinal = "annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-andOPRNs.counts.gtf",
        gtfFinal2 = "annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.counts.gtf"
    conda: env_file
    log: "logs/{specie}/log_recover_coverage_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    shell: """ (
    stringtie -G {input.gtf2} -e -o {output.gtfFinal2} {input.bams}
    stringtie -G {input.gtf} -e -o {output.gtfFinal} {input.bams} ) 2> {log}
    """
