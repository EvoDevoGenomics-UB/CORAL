import os
from os import path
import pandas as pd
from collections import defaultdict
import sys

def validate_samplesheet(samples_df):
    """Validate sample sheet and return dict of sample -> list of FASTQ paths."""
    grouped = defaultdict(list)
    seen_paths = defaultdict(list)
    for _, row in samples_df.iterrows():
        sample = str(row.iloc[0])
        fq = str(row.iloc[1])
        if not os.path.exists(fq):
            print(f"ERROR: Missing FASTQ file for sample '{sample}': {fq}", file=sys.stderr)
            sys.exit(1)
        grouped[sample].append(fq)
        seen_paths[fq].append(sample)
    # Check for duplicates within a sample
    for sample, fqs in grouped.items():
        if len(fqs) != len(set(fqs)):
            print(f"ERROR: Sample '{sample}' lists the same FASTQ file multiple times.", file=sys.stderr)
            sys.exit(1)
    # Check for FASTQs reused across samples
    reused = {fq: s for fq, s in seen_paths.items() if len(set(s)) > 1}
    if reused:
        print("ERROR: Some FASTQ files are used by multiple samples:", file=sys.stderr)
        for fq, s in reused.items():
            print(f"  {fq}  ->  samples: {', '.join(set(s))}", file=sys.stderr)
        sys.exit(1)
    print(f"[INFO] Loaded {len(grouped)} samples, all FASTQ paths validated.")
    return grouped

# Load and validate samples
if "samplesheet" in config and config["samplesheet"]:
    samples_df = pd.read_csv(config["samplesheet"], sep="\t", header=None)
    SAMPLE_TO_FASTQ = validate_samplesheet(samples_df)
    SAMPLES = list(SAMPLE_TO_FASTQ.keys())
else:
    SAMPLES = config["samples"]
    SAMPLE_TO_FASTQ = {
        s: [f"{config['data_dir']}{s}{config['data_suffix']}"]
        for s in SAMPLES
    }

def get_gtf():
    gtf_input_list = [config["extraGTF"]]
        
    if config["run_gffcomapre"] == True:
        gtf_input_list.extend(expand(
            "Gffcompare_results/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_noOPRNs.annotated.gtf",
            specie=config["specie"],
            ref=config["stringtie_guide_opts"],
            intron=config["minimap2_max_intron"],
            threshold=config["operon_threshold"]
        ))
    else:
        gtf_input_list.extend(expand(
            "annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.gtf",
            specie=config["specie"],
            ref=config["stringtie_guide_opts"],
            intron=config["minimap2_max_intron"],
            threshold=config["operon_threshold"]
        ))

    return gtf_input_list

GTF_DICT = {
    os.path.basename(gtf).replace(".gtf", ""): gtf
    for gtf in get_gtf()
    if gtf
}

NAMES = list(GTF_DICT.keys())

def get_TD2_output():
    input_TD2_list=[expand("TD2_results/{specie}/{gtf_name}.fasta.TD2.genome.gff3", gtf_name = NAMES ,specie=config["specie"]),
                    expand("busco_analysis/{specie}/{gtf_name}.TD2.pep.fasta", gtf_name = NAMES ,specie=config["specie"]),
                    expand("busco_analysis/{specie}/BUSCO_prot_{gtf_name}", gtf_name = NAMES ,specie=config["specie"])
                    ]
    if config["reference_annot"] not in (None, [], ""):
        input_TD2_list.extend(expand("busco_analysis/{specie}/{specie}_LRannot_REF.pep.fasta", specie=config["specie"]))
    return input_TD2_list
    
def get_expMatrix_output():
    expMatrix_input_list=[expand("Expression_matrix/{specie}/{gtf_name}/gene_count_matrix_v{intron}.csv",
            specie=config["specie"], intron=config["minimap2_max_intron"], gtf_name = NAMES)]
    if config["reference_annot"] not in (None, [], ""):
        expMatrix_input_list.extend(expand("Expression_matrix/{specie}/ref_annotation/transcript_count_matrix_v{intron}.csv", specie=config["specie"],intron=config["minimap2_max_intron"]))
    
    return expMatrix_input_list

def get_final_output():
    rule_all_input_list=["logs/versions.txt","gamba",
        expand("logs/{sample}_stats_input_reads.txt", sample=SAMPLES),
        expand("alignments/{specie}/{specie}_{sample}_reads_aln_v{intron}.sorted.bam", specie=config["specie"], sample=SAMPLES, intron=config["minimap2_max_intron"]),
        expand("alignments/{specie}/{specie}_{sample}_reads_aln_sorted_v{intron}.stats.txt", specie=config["specie"], sample=SAMPLES, intron=config["minimap2_max_intron"]),
        expand("sample_annotations/{specie}/{specie}_{sample}_guide{ref}_v{intron}.gtf", specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"]),
        expand("GAMBA_results/{specie}/{specie}_{sample}_guide{ref}_v{intron}_opCLEAN_t{threshold}.clean.gtf", specie=config["specie"], sample=SAMPLES, ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/{specie}/Merge_OPRNs-OpGs_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}.sorted.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-andOPRNs.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("annotations/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.counts.gtf", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs_longest_trans_only.fasta", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-noOPRNs.fasta", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_StringtieMerge.clean-andOPRNs.fasta", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}_andOPRNs", specie=config["specie"], ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]),
        expand("busco_analysis/{specie}/BUSCO_results_all_summaries_{specie}_guide{ref}_v{intron}_gambat{threshold}", specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"])]
        
    if config["reference_annot"] not in (None, [], ""):
        rule_all_input_list.extend(expand("busco_analysis/{specie}/BUSCO_trans_{specie}_LRannot_REF", specie=config["specie"]))
        if config["run_gffcomapre"] == True :
            rule_all_input_list.extend(expand("Gffcompare_results/{specie}_LRannot_guide{ref}_v{intron}_gambat{threshold}",specie=config["specie"],ref=config["stringtie_guide_opts"],intron=config["minimap2_max_intron"], threshold=config["operon_threshold"]))
            
    if config["run_expression_matrix"] == True :
        rule_all_input_list.extend(expand("Expression_matrix/{specie}/{gtf_name}/gene_count_matrix_v{intron}.csv",
            specie=config["specie"], intron=config["minimap2_max_intron"], gtf_name = NAMES )),
        if config["reference_annot"] not in (None, [], ""):
            rule_all_input_list.extend(expand("Expression_matrix/{specie}/ref_annotation/transcript_count_matrix_v{intron}.csv", specie=config["specie"],intron=config["minimap2_max_intron"]))

    if config["run_TransDecoder"] == True :
        rule_all_input_list.extend(expand("TD2_results/{specie}/{gtf_name}.fasta.TD2.genome.gff3",
                gtf_name = NAMES , specie=config["specie"])),
        rule_all_input_list.extend(expand("busco_analysis/{specie}/{gtf_name}.TD2.pep.fasta",
                gtf_name = NAMES , specie=config["specie"])),
        rule_all_input_list.extend(expand("busco_analysis/{specie}/BUSCO_prot_{gtf_name}",
                gtf_name = NAMES , specie=config["specie"]))
        if config["reference_annot"] not in (None, [], ""):
            rule_all_input_list.extend(expand("busco_analysis/{specie}/{specie}_LRannot_REF.pep.fasta", specie=config["specie"]))
        
    return rule_all_input_list
