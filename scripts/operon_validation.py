#!/usr/bin/env python3
import gffutils # type: ignore
import sys
import os
import argparse
import logging
import re
from collections import defaultdict, Counter

# Argument parser setup
parser = argparse.ArgumentParser(
    description="Count operons and their included genes from a GTF file."
)
parser.add_argument(
    "-f", "--file", 
    required=True, 
    help="Path to the input GTF file."
)
parser.add_argument(
    "-o","--output",
    #required=True,
    help="Prefix of the output files [default: input base name]."
)
parser.add_argument(
    "--log",
    #required=True,
    help="Log file [default: input base name .log]."
)
# Parse arguments
args = parser.parse_args()
# Extract values
gtf_file = args.file

# If no output was specified, use the input base name
if args.output:
    out_prefix = args.output
else:
    out_prefix = os.path.splitext(gtf_file)[0]

#Set log output file
if args.log:
    log_file = args.log
else:
    log_file = os.path.splitext(gtf_file)[0] + "_OPRNvalidation.log"
#Define logger
logging.basicConfig(filename= log_file, 
					format='%(asctime)s %(levelname)s - %(message)s', 
					filemode='w',
                    level=logging.INFO) 

# Ensure the file exists
if not os.path.isfile(gtf_file):
    print(f"Error: File '{gtf_file}' not found.")
    logging.error(f'File {gtf_file} not found.')
    sys.exit(1)

# Create a database from the GTF file (stored in memory)
db_filename = os.path.splitext(gtf_file)[0] + "_OPRNvalidation.db"
db = gffutils.create_db(
    gtf_file,
    dbfn=db_filename,
    force=True,
    keep_order=True,
    disable_infer_transcripts=True,
    disable_infer_genes=True
)
print(f"File '{db_filename}' created.")
logging.info(f"File '{db_filename}' created.")

# Define output file
output_file = out_prefix + "_OPRNvalidation.keep.tsv"
output_file_outs = out_prefix + "_OPRNvalidation.removed.tsv"
output_file_counts = out_prefix + "_OPRNvalidation.counts.tsv"
# Dictionary to store transcripts per chromosome
chrom_transcripts = defaultdict(list)

#####################################
# Organize transcripts by chromosome
for transcript in db.features_of_type("transcript"):
    chrom_transcripts[transcript.chrom].append(transcript)

# Store detected transcript pairs
contained_pairs = []
excluded_pairs = []
# Find contained transcripts (with progress tracking)
for chrom, transcripts in chrom_transcripts.items():
    print(f"Processing {chrom} ({len(transcripts)} transcripts)...")
    logging.info(f"Processing {chrom} ({len(transcripts)} transcripts)...")
    for idx, transcript in enumerate(transcripts, 1):
        if idx % max(1, len(transcripts) // 20) == 0:  # Print progress every 5% intervals
            progress = (idx / len(transcripts)) * 100
            print(f"Chrom {chrom} Progress: {progress:.1f}%", end="\r")
        trans_gene_id = transcript.attributes['gene_id'][0]
        if re.search("^OPRN", trans_gene_id):
            #print(f"{trans_gene_id}")
            for sub_transcript in transcripts:  # Now only compares within the same chromosome
                if transcript.strand == sub_transcript.strand:
                    operon_gene_id = transcript.attributes['gene_id'][0]
                    opg_gene_id = sub_transcript.attributes['gene_id'][0]
                    if transcript.id == sub_transcript.id or operon_gene_id == opg_gene_id:
                        continue # Skip when comparing itself
                    # Check if coordinates suggest containment
                    if (transcript.start <= (sub_transcript.start + 250) < (transcript.end + 250)) and \
                    (transcript.end >= (sub_transcript.end - 250) > (transcript.start - 250)):

                        # Now check exon overlap — make sure the sub_transcript overlaps at least one operon exon
                        operon_exons = list(db.children(transcript, featuretype='exon', order_by='start'))
                        sub_exons = list(db.children(sub_transcript, featuretype='exon', order_by='start'))

                        overlap_found = False
                        for op_exon in operon_exons:
                            for sub_exon in sub_exons:
                                # Allow small tolerance (±250 bp)
                                if (op_exon.start <= (sub_exon.start + 250) < (op_exon.end + 250)) and \
                                    (op_exon.end >= (sub_exon.end - 250) > (op_exon.start - 250)):
                                #if (op_exon.start - 250 <= sub_exon.end) and (sub_exon.start <= op_exon.end + 250):
                                    overlap_found = True
                                    break
                            if overlap_found:
                                break

                        # Only keep this pair if there's true exon overlap
                        if overlap_found:
                            contained_pairs.append(
                                (transcript.chrom, transcript.strand, transcript.id, operon_gene_id, sub_transcript.id, opg_gene_id)
                            )
                        else:
                            excluded_pairs.append((operon_gene_id, opg_gene_id))
                            logging.warning(
                                f"Excluded {sub_transcript.id} (gene {opg_gene_id}) inside {transcript.id} "
                                f"because it lies entirely within an intron (no exon overlap)."
                            )

######
# Group operon/transcripts by chr and strand
chr_to_operons = defaultdict(list)
for chrom, strand, op_trans, op_gene_id, transcript_id, trans_gene_id in contained_pairs:
    chr_to_operons[chrom, strand].append((op_trans, op_gene_id, transcript_id, trans_gene_id))

prefinal_pairs = []
seen_transcripts = set()
for (chrom, strand), op_trans_list in chr_to_operons.items():
    # Sort transcripts by start position
    #op_trans_list.sort(key=lambda x: (x[1]))  # Sort by start coordinate of operons
    print(f"Checking operons on strand '{strand}' of '{chrom}'")
    for current_transcript in op_trans_list:
        op_trans, op_gene_id, transcript_id, trans_gene_id = current_transcript
        if trans_gene_id in seen_transcripts:
            continue
        prefinal_pairs.append((op_gene_id, trans_gene_id))
        seen_transcripts.add(trans_gene_id)
    
# Count how many transcripts each operon contains
operon_counts_def = Counter(operon for operon, _ in prefinal_pairs)
# Keep only operons that contain **two or more** transcripts
out_operons = {operon for operon, count in operon_counts_def.items() if count < 2}
# Filter out rows where the first column appears in the second column
final_pairs_DEF = [pair for pair in prefinal_pairs if pair[0] not in out_operons]
out_pairs_DEF = [pair for pair in prefinal_pairs if pair not in final_pairs_DEF]
for x in excluded_pairs :
    out_pairs_DEF.append(x)
##########################
# Write to output file
with open(output_file, "w") as out_file:
    out_file.write("Operon\tContained_gene\n")
    for operon_gene_id, opg_gene_id in final_pairs_DEF:
        out_file.write(f"{operon_gene_id}\t{opg_gene_id}\n")
# Write to output file
with open(output_file_outs, "w") as output_file_outs:
    output_file_outs.write("Operon\tContained_gene\n")
    for operon_gene_id, opg_gene_id in out_pairs_DEF:
        output_file_outs.write(f"{operon_gene_id}\t{opg_gene_id}\n")

print(f"Operon-Genes found saved to TSV file {output_file}")
logging.info(f"Operon-Genes found saved to TSV file {output_file}")

################
#Create the file with the right operon and operon genes.
# Store all transcripts for gene_ids
gene_to_trans_ids = []
all_trans_ids = []
removed_trans_ids = []
operon_ids = {operon for operon, _ in final_pairs_DEF}
contained_ids = {opg_gene_id for _, opg_gene_id in final_pairs_DEF}
for transcript in db.features_of_type("transcript"):
    gene_id = transcript.attributes['gene_id'][0]
    if gene_id in operon_ids or gene_id in contained_ids:
        all_trans_ids.append(transcript.id)
        if gene_id in contained_ids:
            gene_to_trans_ids.append(transcript.id)
    else:
        removed_trans_ids.append(transcript.id)

# Define output GTF filenames
operon_gtf_file = out_prefix + "_OPRNvalidation.clean.gtf"
# Build mapping: contained gene_id → operon_id
gene_to_operon = {opg_gene_id: operon_gene_id for operon_gene_id, opg_gene_id in final_pairs_DEF}
for operon_gene_id, _ in final_pairs_DEF:
    gene_to_operon[operon_gene_id] = operon_gene_id

# Write operon transcripts to GTF
with open(operon_gtf_file, "w") as operon_out:
    for trans_id in all_trans_ids:
        trans_feature = db[trans_id]
        trans_gene_id = trans_feature.attributes['gene_id'][0]
        # Write the transcript feature
        if trans_gene_id in gene_to_operon:
            op_id = gene_to_operon[trans_gene_id]
            new_attrs = f'; operon_id "{op_id}";'
        else:
            new_attrs = ""
        #operon_out.write(str(trans_feature).replace('""', '"').replace('";"', '";') + new_attrs + "\n")
        operon_out.write((str(trans_feature) + new_attrs + "\n").replace(';;', ';'))
        # Write its child features (e.g., exons)
        for feature in db.children(trans_id, featuretype='exon', order_by='start'):
            #operon_out.write(str(feature).replace('""', '') + "\n")
            operon_out.write((str(feature) + ";\n").replace('""', '"').replace('";"', '";'))

# Define output GTF filenames
opg_gtf_file = out_prefix + "_OPRNvalidation.OpGclean.gtf"
# Write operon transcripts to GTF
with open(opg_gtf_file, "w") as opg_out:
    for trans_id in gene_to_trans_ids:
        trans_feature = db[trans_id]
        trans_gene_id = trans_feature.attributes['gene_id'][0]
        # Write the transcript feature
        if trans_gene_id in gene_to_operon:
            op_id = gene_to_operon[trans_gene_id]
            new_attrs = f'; operon_id "{op_id}";'
        else:
            new_attrs = ""
        # Write the transcript feature
        #opg_out.write(str(trans_feature).replace('""', '"').replace('";"', '";') + new_attrs + "\n")
        opg_out.write((str(trans_feature) + new_attrs + "\n").replace(';;', ';'))
        # Write its child features (e.g., exons)
        for feature in db.children(trans_id, featuretype='exon', order_by='start'):
            #opg_out.write(str(feature).replace('""', '') + "\n")
            opg_out.write((str(feature) + ";\n").replace('""', '"').replace('";"', '";').replace(';;', ';'))

# Define output GTF filenames
op_gtf_file = out_prefix + "_OPRNvalidation.OPRNclean.gtf"
# Write operon transcripts to GTF
with open(op_gtf_file, "w") as op_out:
    for trans_id in all_trans_ids:
        if trans_id not in gene_to_trans_ids:
            trans_feature = db[trans_id]
            trans_gene_id = trans_feature.attributes['gene_id'][0]
            # Write the transcript feature
            if trans_gene_id in gene_to_operon:
                op_id = gene_to_operon[trans_gene_id]
                new_attrs = f'; operon_id "{op_id}";'
            else:
                new_attrs = ""
            # Write the transcript feature
            #op_out.write(str(trans_feature).replace('""', '"').replace('";"', '";') + new_attrs + "\n")
            op_out.write((str(trans_feature) + new_attrs + "\n").replace(';;', ';'))
            # Write its child features (e.g., exons)
            for feature in db.children(trans_id, featuretype='exon', order_by='start'):
                #op_out.write(str(feature).replace('""', '') + "\n")
                op_out.write((str(feature) + ";\n").replace('""', '"').replace('";"', '";').replace(';;', ';'))


# Define output GTF filenames
excuded_gtf_file = out_prefix + "_OPRNvalidation.excluded.gtf"
# Write operon transcripts to GTF
with open(excuded_gtf_file, "w") as excluded_out:
    for trans_id in removed_trans_ids:
        trans_feature = db[trans_id]
        # Write the transcript feature
        #excluded_out.write(str(trans_feature).replace('""', '"').replace('";"', '";') + ";\n")
        excluded_out.write(str(trans_feature) + ";\n")
        # Write its child features (e.g., exons)
        for feature in db.children(trans_id, featuretype='exon', order_by='start'):
            #excluded_out.write(str(feature).replace('""', '') + "\n")
            excluded_out.write((str(feature) + ";\n").replace('""', '"').replace('";"', '";').replace(';;', ';'))

################
# Counting the number of genes by operon
operon_size_counter = Counter()
true_final_pairs_DEF = []

for operon, transcript in final_pairs_DEF:
    true_final_pairs_DEF.append((operon, transcript))

for operon, transcript in true_final_pairs_DEF:
    operon_size_counter[operon] += 1

print(f"Total number of operons: {len(operon_size_counter)}")
logging.info(f"Total number of operons: {len(operon_size_counter)}")
print(f"Total number of genes inside operons: {len(true_final_pairs_DEF)}")
logging.info(f"Total number of genes inside operons: {len(true_final_pairs_DEF)}")

# Grouping of operons
summary = {
    "2 genes": 0,
    "3 genes": 0,
    "4 genes": 0,
    "5 genes": 0,
    ">5 genes": 0,
}

for count in operon_size_counter.values():
    if count == 2:
        summary["2 genes"] += 1
    elif count == 3:
        summary["3 genes"] += 1
    elif count == 4:
        summary["4 genes"] += 1
    elif count == 5:
        summary["5 genes"] += 1
    elif count > 5:
        summary[">5 genes"] += 1

# Report of the results into the terminal and log file
print("Summary of operons by gene number:")
logging.info("Summary of operons by gene number:")
for category, total in summary.items():
    print(f"{category}: {total}")
    logging.info(f"{category}: {total}")

# Find the operon with the most genes
max_operon = max(operon_size_counter, key=operon_size_counter.get)
max_count = operon_size_counter[max_operon]

print(f"Operon with most genes: {max_operon} ({max_count} genes)")
logging.info(f"Operon with most genes: {max_operon} ({max_count} genes)")

# Write to output file
with open(output_file_counts, "w") as output_file_counts:
    output_file_counts.write("Operon\tNum_contained_genes\n")
    for operon, count in operon_size_counter.items():
        output_file_counts.write(f"{operon}\t{count}\n")
    
print(f"Detail Operon-count saved to TSV file {output_file_counts}")
logging.info(f"Detail Operon-count saved to TSV file {output_file_counts}")
