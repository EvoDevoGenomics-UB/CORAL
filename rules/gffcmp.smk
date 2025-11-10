#### OPTIONAL steps
## Comparing new annotations againts reference one
rule run_gffcompare:
    input:
        ref = REF ,
        gtf_longest = rules.run_longest_trans_filter.output.filtergtf ,
        gtf_noOPRNs = rules.run_final_annotation.output.noOPRNs ,
        gtf_andOPRNs = rules.run_final_annotation.output.andOPRNs
    output:
        gffcmp_dir = directory("Gffcompare_results/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}")
    params:
        prefix = "{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}"
    conda: env_file
    log: "logs/{specie}/log_gffcomapre_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.log"
    shell:""" (
    if [[ ! -f {input.ref} ]] ; then
        echo \"Error: No reference annotation provided.\" >&2
        exit 1
    else
        mkdir -p {output.gffcmp_dir}
        gffcompare -r {input.ref} {input.gtf_longest} -o ./{output.gffcmp_dir}/{params.prefix}_longest_trans
        gffcompare -r {input.ref} {input.gtf_noOPRNs} -o ./{output.gffcmp_dir}/{params.prefix}_noOPRNs
        gffcompare -r {input.ref} {input.gtf_andOPRNs} -o ./{output.gffcmp_dir}/{params.prefix}_andOPRNs

        echo 'Performing counting of transcript types...'
        for i in ./{output.gffcmp_dir}/{params.prefix}_longest_trans ./{output.gffcmp_dir}/{params.prefix}_noOPRNs ./{output.gffcmp_dir}/{params.prefix}_andOPRNs ; do
            (for x in "=" c k m n j e o s x i y p r u ; do
                count=$(awk -v a="$x" '{{if($4==a) print $5,$4}}' $i.tracking | wc -l)
                echo "$x $count"
            done) > $i.gffcmp_trans_types.txt
            echo "File '$i.gffcmp_trans_types.txt' created."
        done
    fi ) 2>&1 | tee {log}
    """
