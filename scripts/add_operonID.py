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
    description="Add the OPRN.ID to OpGs of a given GTF file (with OPRNs)."
)
parser.add_argument(
    "-f", "--file", 
    required=True, 
    help="Path to the input GTF file."
)
parser.add_argument(
    "-o","--output",
    required=True,
    help="Output file."
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
    output_file = args.output
    out_prefix = os.path.splitext(output_file)[0]
else:
    out_prefix = os.path.splitext(gtf_file)[0]
    output_file = out_prefix + ".OPRNids.gtf"

#Set log output file
if args.log:
    log_file = args.log
else:
    log_file = os.path.splitext(gtf_file)[0] + ".OPRNids.log"
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
db_filename = os.path.splitext(gtf_file)[0] + ".OPRNids.db"
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

# Dictionary to store transcripts per chromosome
chrom_transcripts = defaultdict(list)

#####################################
# Organize transcripts by chromosome
for transcript in db.features_of_type("transcript"):
    chrom_transcripts[transcript.chrom].append(transcript)

# Store detected transcript pairs
contained_pairs = []
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

######
# Group operon/transcripts by chr and strand
chr_to_operons = defaultdict(list)
for chrom, strand, op_trans, op_gene_id, transcript_id, trans_gene_id in contained_pairs:
    chr_to_operons[chrom, strand].append((op_trans, op_gene_id, transcript_id, trans_gene_id))

prefinal_pairs = []
seen_transcripts = set()
for (chrom, strand), op_trans_list in chr_to_operons.items():
    for current_transcript in op_trans_list:
        op_trans, op_gene_id, transcript_id, trans_gene_id = current_transcript
        if trans_gene_id in seen_transcripts:
            continue
        prefinal_pairs.append((op_gene_id, trans_gene_id))
        seen_transcripts.add(trans_gene_id)

final_pairs_DEF = [pair for pair in prefinal_pairs ]

################
#Create the file with the right operon and operon genes.
# Store all transcripts for gene_ids
gene_to_trans_ids = []
all_trans_ids = []
peron_ids = {operon for operon, _ in final_pairs_DEF}
contained_ids = {opg_gene_id for _, opg_gene_id in final_pairs_DEF}
for transcript in db.features_of_type("transcript"):
    gene_id = transcript.attributes['gene_id'][0]
    all_trans_ids.append(transcript.id)
    if gene_id in contained_ids:
        gene_to_trans_ids.append(transcript.id)

# Build mapping: contained gene_id → operon_id
gene_to_operon = {opg_gene_id: operon_gene_id for operon_gene_id, opg_gene_id in final_pairs_DEF}

# Write operon transcripts to GTF
with open(output_file, "w") as gtf_out:
    for trans_id in all_trans_ids:
        trans_feature = db[trans_id]
        trans_gene_id = trans_feature.attributes['gene_id'][0]
        # Write the transcript feature
        if trans_gene_id in gene_to_operon:
            op_id = gene_to_operon[trans_gene_id]
            new_attrs = f' operon_id "{op_id}";'
        else:
            new_attrs = ""
        gtf_out.write(str(trans_feature).replace('""', '"').replace('";"', '";') + new_attrs + "\n")
        # Write its child features (e.g., exons)
        for feature in db.children(trans_id, featuretype='exon', order_by='start'):
            gtf_out.write(str(feature).replace('""', '') + "\n")
