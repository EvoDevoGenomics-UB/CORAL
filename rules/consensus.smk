gtfsoperons_samples = []
gtfsopgenes_samples = []
gtfsclean_samples = []
for SAMPLE in SAMPLES:
    gtfsoperons_samples.append(
        "GAMBA_results/{{specie}}/{{specie}}_{}_guide{{ref}}_v{{intron}}_Operons_t{{threshold}}.clean.gtf".format(SAMPLE)
    )
    gtfsopgenes_samples.append(
        "GAMBA_results/{{specie}}/{{specie}}_{}_guide{{ref}}_v{{intron}}_OperonGenes_t{{threshold}}.clean.gtf".format(SAMPLE)
    )
    gtfsclean_samples.append(
        "GAMBA_results/{{specie}}/{{specie}}_{}_guide{{ref}}_v{{intron}}_opCLEAN_t{{threshold}}.clean_longest_trans_only.gtf".format(SAMPLE)
    )

# fmt: off
# Create oepron and operon-contained genes annotations
rule run_operon_annotation:
    input:
        gtfsoperons=ancient(gtfsoperons_samples),
        gtfsopgenes=ancient(gtfsopgenes_samples)
    output:
        operongtf=temp("annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_OPRNs.gtf"),
        opgenesgtf=temp("annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_OpGs.gtf"),
        merge="annotations/{specie}/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted.gtf",
        def_file=touch("annotations/{specie}/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.clean.gtf"),
        opgenesgtfCLEAN=touch("annotations/{specie}/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.OpGclean.gtf"),
        oprngtfCLEAN=touch("annotations/{specie}/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.OPRNclean.gtf"),
        excluded_file=touch("annotations/{specie}/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.excluded.gtf"),
        db_file=temp("annotations/{specie}/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.db")
    log:
        logOPRN="logs/{specie}/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted_OPRNvalidation.log",
        logSTRG="logs/{specie}/log_StrignTie_merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    conda:
        env_file
    params:
        g_param=config["stringtie_OpGs_g"],
        name="{specie}_guide{ref}_v{intron}_gambat{threshold}",
        snakedir=SNAKEDIR
    shell:
        """
    (mkdir -p annotations ; mkdir -p annotations/{wildcards.specie} ;
    (for i in {input.gtfsopgenes}; do
        if [ -s "$i" ]; then echo "$i" ; fi
    done) \
    > annotations/{wildcards.specie}/List_merge_OpGs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt
    if [ -s annotations/{wildcards.specie}/List_merge_OpGs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ] ; then
        stringtie --merge -l OpG -f 0 -F 0 -T 0 -c 0 -m 200 -g {params.g_param} -o {output.opgenesgtf} annotations/{wildcards.specie}/List_merge_OpGs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ; \
        grep 'StringTie	transcript' {output.opgenesgtf} | wc -l ; \
    else
        echo "No non-empty OpG GTFs found" >&2
        touch {output.opgenesgtf}
    fi

    (for i in {input.gtfsoperons}; do
        if [ -s "$i" ]; then echo "$i" ; fi
    done) \
    > annotations/{wildcards.specie}/List_merge_OPRNs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt
    if [ -s annotations/{wildcards.specie}/List_merge_OPRNs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ] ; then
        stringtie --merge -l OPRN -f 0 -F 0 -T 0 -c 0 -m 200 -g 0 -o {output.operongtf} annotations/{wildcards.specie}/List_merge_OPRNs.{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}.txt ; \
    else
        echo "No non-empty OPRN GTFs found" >&2
        touch ./{output.operongtf}
    fi

    cat {output.operongtf} {output.opgenesgtf} > {params.name}.tmp.gtf ; \
    gffread --sort-alpha -F -T -o {output.merge} {params.name}.tmp.gtf ; rm {params.name}.tmp.gtf

    if [ -s {output.merge} ] ; then
        python {params.snakedir}/scripts/operon_validation.py -f {output.merge} --log {log.logOPRN}
        grep 'StringTie	transcript' {output.opgenesgtfCLEAN} | wc -l 
    else
        (echo "Empty OpG-OPRN GTF, validation skiped"
        touch ./{output.db_file}
        touch ./{output.opgenesgtfCLEAN}
        touch ./{output.oprngtfCLEAN}
        touch ./{output.excluded_file} ) 2>&1 | tee {log.logOPRN}
    fi
    ) 2>&1 | tee {log.logSTRG}
    """


# Create gene final consensus annotations
rule run_gCLEAN_annotation:
    input:
        gtfsclean=ancient(gtfsclean_samples),
        oprngtf=ancient(rules.run_operon_annotation.output.oprngtfCLEAN),
        opgsgtf=ancient(rules.run_operon_annotation.output.opgenesgtfCLEAN)
    output:
        cleanfinal="annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_mergeCLEAN.gtf",
    log:
        "logs/{specie}/log_StrignTie_merge_opCLEAN_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log",
    conda:
        env_file
    threads: config["threads"]
    params:
        freq=config["stringtie_freq"],
        g_param=config["stringtie_g"],
        g_param0=int(config["stringtie_OpGs_g"]) ,
        opts=config["stringtie_merge_opts"]
    shell:
        """
    var_name="{wildcards.specie}guide{wildcards.ref}v{wildcards.intron}gambat{wildcards.threshold}"
    ((for i in {input.gtfsclean} ; do echo $i ; done ) > annotations/{wildcards.specie}/List_merge_opCLEAN.$var_name.txt ; \
    stringtie --version

    stringtie --merge annotations/{wildcards.specie}/List_merge_opCLEAN.$var_name.txt \
     -l g -f {params.freq} {params.opts} -g {params.g_param} -o GTFfileLose.$var_name.tmp -p {threads}; \
    echo "GTFfileLose.$var_name.tmp done!"

    stringtie --merge annotations/{wildcards.specie}/List_merge_opCLEAN.$var_name.txt \
     -l g -f {params.freq} -c 25 -F 14.0 -T 3.0 -m 200 -g {params.g_param0} -o GTFfileStrict.$var_name.tmp -p {threads}; \
    echo "GTFfileStrict.$var_name.tmp done!"

    gffcompare -r GTFfileLose.$var_name.tmp GTFfileStrict.$var_name.tmp -o filter0_$var_name
    awk '{{if($3=="=" || $3=="c" ||$3=="j") print $2}}' filter0_$var_name.GTFfileStrict.$var_name.tmp.tmap > filter0.$var_name.list.tmp
    gffread --nids filter0.$var_name.list.tmp GTFfileLose.$var_name.tmp -T -o GTFfileLose.$var_name.clean.tmp
    echo "GTFfileLose.$var_name.clean.tmp done!"

    stringtie --merge GTFfileStrict.$var_name.tmp GTFfileLose.$var_name.clean.tmp \
      -l g -f {params.freq} -F 0 -T 0 -c 0 -m 200 -g {params.g_param0} -o GTFfile.$var_name.tmp -p {threads}
    
    if [ $(grep 'StringTie	transcript' {input.oprngtf} | wc -l ) -ge 1 ] ; then
        gffcompare -r {input.oprngtf} GTFfile.$var_name.tmp -o filter_$var_name
        awk '{{if($3=="=") print $5}}' filter_$var_name.GTFfile.$var_name.tmp.tmap > filter.$var_name.list.tmp
        gffcompare -r {input.opgsgtf} GTFfile.$var_name.tmp -o filter2_$var_name
        awk '{{if($3=="=" || $3=="c" || $3=="j") print $5}}' filter2_$var_name.GTFfile.$var_name.tmp.tmap > filter2.$var_name.list.tmp
        cat filter.$var_name.list.tmp filter2.$var_name.list.tmp > filter3.$var_name.list.tmp

        gffread --nids filter3.$var_name.list.tmp GTFfile.$var_name.tmp -T -o {output.cleanfinal}
        rm filter*${{var_name}}*
    else
        cp GTFfile.$var_name.tmp {output.cleanfinal}
    fi

    rm GTFfile.$var_name.tmp ;

    echo "  Final merge CLEAN done" ; \
    grep 'StringTie	transcript' {output.cleanfinal} | wc -l ) 2>&1 | tee {log}
    """


# Create Merge final consensus annotations
rule run_final_annotation:
    input:
        cleanfinal=rules.run_gCLEAN_annotation.output.cleanfinal,
        excluded_file=rules.run_operon_annotation.output.excluded_file,
        opgenesgtf=rules.run_operon_annotation.output.opgenesgtfCLEAN,
        oprnsgtf=rules.run_operon_annotation.output.oprngtfCLEAN
    output:
        noOPRNs="annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.gtf",
        andOPRNs="annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-andOPRNs.gtf"
    log:
        log1="logs/{specie}/log_final_annotations_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_part1.log",
        log2="logs/{specie}/log_final_annotations_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_part2.log"
    conda:
        env_file
    params:
        freq=config["stringtie_freq"],
        g_param=config["stringtie_OpGs_g"],#g_param=config["stringtie_g"],
        snakedir=SNAKEDIR
    shell:
        """
    (if [ -s {input.excluded_file} ] ; then 
        FILES="{input.cleanfinal} {input.excluded_file}"
    else
        FILES="{input.cleanfinal}"
    fi
    if [ -s {input.opgenesgtf} ] ; then 
        stringtie --merge $FILES \
        -G {input.opgenesgtf} \
        -l g -f {params.freq} -F 0 -T 0 -c 0 -g {params.g_param} \
        -o {output.noOPRNs} ; \
    else
        stringtie --merge $FILES \
        -l g -f {params.freq} -F 0 -T 0 -c 0 -g {params.g_param} \
        -o {output.noOPRNs} ; \
    fi
    echo "  Final CLEAN-noOPRNs done" ; \
    grep 'StringTie	transcript' {output.noOPRNs} | wc -l ; ) 2>&1 | tee {log.log1}

    (if [ -s {input.oprnsgtf} ] ; then 
        cat {input.oprnsgtf} {output.noOPRNs} > {output.andOPRNs}.1.tmp ; \
        gffread --sort-alpha -F -T -o {output.andOPRNs}.2.tmp {output.andOPRNs}.1.tmp ; 
        python {params.snakedir}/scripts/add_operonID.py -f {output.andOPRNs}.2.tmp -o {output.andOPRNs} --log {output.andOPRNs}.log
        rm {output.andOPRNs}*.OPRNids.db {output.andOPRNs}*.tmp
        echo "  Final CLEAN-andOPRNs done"
        grep 'StringTie	transcript' {output.andOPRNs} | wc -l 
    else
     touch {output.andOPRNs}
     echo "  Final CLEAN-andOPRNs not created because of lack of OPRNs..."
    fi ) 2>&1 | tee {log.log2}
    """
# fmt: on

# Obtaining coverage of final annotation
rule run_recover_coverage:
    input:
        gtf=rules.run_final_annotation.output.andOPRNs,
        gtf2=rules.run_final_annotation.output.noOPRNs,
        bams=ancient(
            expand(
                "alignments/{specie}/{specie}_{sample}_reads_aln_v{intron}.sorted.bam",
                specie=config["specie"],
                sample=SAMPLES,
                intron=config["minimap2_max_intron"],
            )
        )
    output:
        gtfFinal="annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-andOPRNs.counts.gtf",
        gtfFinal2="annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.counts.gtf"
    log:
        "logs/{specie}/log_recover_coverage_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    conda:
        env_file
    shell:
        """ (
    stringtie -G {input.gtf2} -e -o {output.gtfFinal2} {input.bams}
    if [ -s {input.gtf} ] ; then
        stringtie -G {input.gtf} -e -o {output.gtfFinal} {input.bams}
    else
        touch {output.gtfFinal}
    fi ) 2> {log}
    """
