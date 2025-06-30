import gffutils # type: ignore
import sys
import os
import argparse
import logging
from collections import defaultdict, Counter

# Argument parser setup
parser = argparse.ArgumentParser(
    description="Detect operons from a GTF file with optional coverage threshold filtering."
)
parser.add_argument(
    "-f", "--file", 
    required=True, 
    help="Path to the input GTF file."
)
parser.add_argument(
    "-t", "--threshold", 
    type=float, 
    default=1.25, 
    help="Factor applied to modify selection threshold. " \
    "Inner transcript coverage should be equal or bigger than "
    "(operon-coverage * THRESHOLD) [default: 1.25]."
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
threshold = args.threshold
# If no output was specified, use the input base name
if args.output:
    out_prefix = args.output
else:
    out_prefix = os.path.splitext(gtf_file)[0]

#Set log output file
if args.log:
    log_file = args.log
else:
    log_file = os.path.splitext(gtf_file)[0] + "_OFv9.log"
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
db_filename = os.path.splitext(gtf_file)[0] + "_annotation_v9.t"+ str(threshold) +".db"
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
output_file = out_prefix + "_operons_found_v9.t"+ str(threshold) +".tsv"
output_file_outs = out_prefix + "_operons_found_v9.t"+ str(threshold) +".outs.tsv"
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
        gene_id = transcript.attributes['gene_id'][0]
        transcript_cov = float(transcript.attributes['cov'][0])
        
        for sub_transcript in transcripts:  # Now only compares within the same chromosome
            if transcript.id == sub_transcript.id:
                continue # Skip when comparing itself
            if transcript.strand == sub_transcript.strand:
                sub_transcript_cov = float(sub_transcript.attributes['cov'][0])
                if (transcript_cov * threshold) < sub_transcript_cov:
                    if transcript.start <= (sub_transcript.start + 250 ) and transcript.end >= (sub_transcript.end - 250):
                        exons = len(list(db.children(sub_transcript.id, featuretype='exon')))
                        if exons > 1 :
                            contained_pairs.append((transcript.chrom, transcript.strand, transcript.id, transcript.start, transcript.end, sub_transcript.id, sub_transcript.start, sub_transcript.end, float(sub_transcript.attributes['FPKM'][0])))
                            # if transcript.id=="STRG.295.2":
                            #     print(f"{sub_transcript.id}")
                        elif sub_transcript_cov > (transcript_cov * threshold * 3): #Include monoexonic if theri cov is over 3.
                            contained_pairs.append((transcript.chrom, transcript.strand, transcript.id, transcript.start, transcript.end, sub_transcript.id, sub_transcript.start, sub_transcript.end, float(sub_transcript.attributes['FPKM'][0])))
######
# Count how many transcripts each operon contains
operon_first_counts = Counter(op_trans for _, _, op_trans, _, _, _, _, _, _ in contained_pairs)
# Keep only operons that contain **two or more** transcripts
first_valid_operons = {op_trans for op_trans, count in operon_first_counts.items() if count > 1}

# Group contained transcripts by operon
operon_to_transcripts = defaultdict(list)
for chrom, strand, op_trans, op_start, op_end, transcript, trans_start, trans_end, trans_fpkm in contained_pairs:
    if op_trans in first_valid_operons:
        operon_to_transcripts[op_trans, chrom, strand, op_start, op_end].append((transcript, trans_start, trans_end, trans_fpkm))

# Remove overlapping contained transcripts **within the same operon**
prefinal_pairs = []
for (op_trans, chrom, strand, op_start, op_end), transcript_list in operon_to_transcripts.items():
    # Sort transcripts by start position
    transcript_list.sort(key=lambda x: (x[1]))  # Sort by start coordinate
    
    non_overlapping = []
    for current_transcript in transcript_list:
        transcript_id, start, end , fpkm = current_transcript
        if not non_overlapping:  # First add
            non_overlapping.append(current_transcript)
        else:
            last_id, last_start, last_end, last_fpkm = non_overlapping[-1]
            if start > (last_end - 50): #No overlap with previous (50bp tolerance)
                non_overlapping.append(current_transcript)
            else:
                # Hay solapamiento, conservar el de mayor cobertura
                if fpkm > last_fpkm:
                    non_overlapping[-1] = current_transcript
    
    # Add non-overlapping transcripts to final output
    for transcript_id, start, end , fpkm in non_overlapping:
        prefinal_pairs.append((op_trans, chrom, strand, op_start, op_end, transcript_id, start, end , fpkm))

# Count how many transcripts each operon contains
operon_second_counts = Counter(op_trans for op_trans, _, _, _, _, _, _, _, _ in prefinal_pairs)
# Keep only operons that contain **two or more** transcripts
second_valid_operons = {op_trans for op_trans, count in operon_second_counts.items() if count > 1}

chr_to_operons = defaultdict(list)
for op_trans, chrom, strand, op_start, op_end, transcript_id, start, end , fpkm in prefinal_pairs:
    if op_trans in second_valid_operons:
        chr_to_operons[chrom, strand].append((op_trans, op_start, op_end, transcript_id, start, end , fpkm))

overlapping = []
seen_transcripts = set()
for (chrom, strand), op_trans_list in chr_to_operons.items():
    # Sort transcripts by start position
    op_trans_list.sort(key=lambda x: (x[1]))  # Sort by start coordinate of operons
    print(f"Checking operons on strand '{strand}' of '{chrom}'")
    for current_transcript in op_trans_list:
        operon, op_start, op_end, transcript_id, start, end , fpkm = current_transcript
        if not overlapping:  # First add
            counter = int(1)
            op_ID = str("OPRN."+ str(counter))
            overlapping.append((operon, op_start, op_end, transcript_id, start, end , fpkm, op_ID))
            seen_transcripts.add(transcript_id)
            #print(f"'{overlapping}'")
        else:
            #print(f"'{overlapping[-1]}'")
            l_operon, l_op_start, l_op_end, l_transcript_id, l_start, l_end , l_fpkm, l_op_ID = overlapping[-1]
            if transcript_id in seen_transcripts: #skip if already added
                continue
            if op_start <= (l_op_end - 250):
                # Hay solapamiento superior a 250 bp, conservar el OP_id
                op_ID = l_op_ID
                overlapping.append((operon, op_start, op_end, transcript_id, start, end , fpkm, op_ID))
                seen_transcripts.add(transcript_id)
            else:
                counter = counter + 1
                op_ID = str("OPRN."+ str(counter))
                overlapping.append((operon, op_start, op_end, transcript_id, start, end , fpkm, op_ID))
                seen_transcripts.add(transcript_id)
    
# Group contained transcripts by operon
operon_to_transcripts_2 = defaultdict(list)
# Add overlapping transcripts to final output
for operon, _, _, transcript_id, start, end , fpkm, op_ID in overlapping:
    operon_to_transcripts_2[op_ID].append((transcript_id, start, end , fpkm, operon))

# Remove overlapping contained transcripts **within the same operon**
final_pairs = []
for op_ID, transcript_list2 in operon_to_transcripts_2.items():
    # Sort transcripts by start position
    transcript_list2.sort(key=lambda x: (x[1]))  # Sort by start coordinate
    
    non_overlapping2 = []
    for current_transcript2 in transcript_list2:
        transcript_id, start, end , fpkm, operon = current_transcript2
        if not non_overlapping2:  # First add
            non_overlapping2.append(current_transcript2)
        else:
            last_id, last_start, last_end, last_fpkm, operon = non_overlapping2[-1]
            if start > (last_end - 50): #No overlap with previous (50bp tolerance)
                non_overlapping2.append(current_transcript2)
            else:
                # Hay solapamiento, conservar el de mayor cobertura
                if fpkm > last_fpkm:
                    non_overlapping2[-1] = current_transcript2
    
    # Add non-overlapping transcripts to final output
    for transcript_id, _, _ , _, operon in non_overlapping2:
        final_pairs.append((op_ID, operon, transcript_id))

# Count how many transcripts each operon contains
operon_counts_def = Counter(operon for operon, _, _ in final_pairs)
# Keep only operons that contain **two or more** transcripts
out_operons = {operon for operon, count in operon_counts_def.items() if count < 2}
# Filter out rows where the first column appears in the second column
final_pairs_DEF = [pair for pair in final_pairs if pair[0] not in out_operons]
out_pairs_DEF = [pair for pair in final_pairs if pair[0] in out_operons]

##########################
# Write to output file
with open(output_file, "w") as out_file:
    out_file.write("Operon\tOperonTrans\tContained_transcript\n")
    for operon, op_trans, transcript in final_pairs_DEF:
        out_file.write(f"{operon}\t{op_trans}\t{transcript}\n")
# Write to output file
with open(output_file_outs, "w") as output_file_outs:
    output_file_outs.write("Operon\tOperonTrans\tContained_transcript\n")
    for operon, op_trans, transcript in out_pairs_DEF:
        output_file_outs.write(f"{operon}\t{op_trans}\t{transcript}\n")

print(f"Operon-Genes found saved to TSV file {output_file}")
logging.info(f"Operon-Genes found saved to TSV file {output_file}")

# Create sets of operon and gene transcript IDs
operon_ids = {op_trans for _, op_trans, _ in final_pairs_DEF}
contained_ids = {transcript for _, _, transcript in final_pairs_DEF}

# Store all gene_ids of operon-genes
gene_ids = []
for transcript in db.features_of_type("transcript"): #db.features_of_type("transcript")
    if transcript.id in contained_ids:
        gene_id = transcript.attributes['gene_id'][0]
        #print(gene_id)
        gene_ids.append(gene_id)

# Store all transcripts for gene_ids
gene_transcripts = []
for transcript in db.features_of_type("transcript"):
    if transcript.id not in operon_ids:
        gene_id = transcript.attributes['gene_id'][0]
        if gene_id in gene_ids:
            gene_transcripts.append(transcript.id)

# Define output GTF filenames
operon_gtf_file = out_prefix + "_Operons_v9.t" + str(threshold) + ".gtf"
contained_gtf_file = out_prefix + "_OperonGenes_v9.t" + str(threshold) + ".gtf"
containedALL_gtf_file = out_prefix + "_OperonGenesALL_v9.t" + str(threshold) + ".gtf"
clean_gtf_file = out_prefix + "_opCLEAN_v9.t" + str(threshold) + ".gtf"

# Write operon transcripts to GTF
with open(operon_gtf_file, "w") as operon_out:
    for operon_id in operon_ids:
        operon_feature = db[operon_id]
        # Write the transcript feature
        operon_out.write(str(operon_feature).replace('""', '"').replace('";"', '";') + ";\n")
        # Write its child features (e.g., exons)
        for feature in db.children(operon_id, featuretype='exon', order_by='start'):
            operon_out.write(str(feature).replace('""', '') + "\n")

# Write non-overlaped contained gene transcripts to GTF
with open(contained_gtf_file, "w") as contained_out:
    for transcript_id in contained_ids:
        transcript_feature = db[transcript_id]
        contained_out.write(str(transcript_feature).replace('""', '"').replace('";"', '";') + ";\n")
        for feature in db.children(transcript_id, featuretype='exon', order_by='start'):
            contained_out.write(str(feature).replace('""', '') + "\n")
# Write ALL contained gene transcripts to GTF
with open(containedALL_gtf_file, "w") as containedALL_out:
    for transcript_id in gene_transcripts:
        transcript_feature = db[transcript_id]
        containedALL_out.write(str(transcript_feature).replace('""', '"').replace('";"', '";') + ";\n")
        for feature in db.children(transcript_id, featuretype='exon', order_by='start'):
            containedALL_out.write(str(feature).replace('""', '') + "\n")

# Write Clean GTF no containgin OPRNs nor OpGenes
with open(clean_gtf_file, "w") as clean_out:
    for transcript in db.features_of_type("transcript"):
        transcript_cov = float(transcript.attributes['cov'][0])
        if transcript_cov <= 1:
            continue
        if transcript.id not in operon_ids and transcript.id not in gene_transcripts:
            transcript_feature = db[transcript.id]
            clean_out.write(str(transcript_feature).replace('""', '"').replace('";"', '";') + ";\n")
            for feature in db.children(transcript.id, featuretype='exon', order_by='start'):
                clean_out.write(str(feature).replace('""', '') + "\n")

print(f"GTF files saved: \n {operon_gtf_file} (operons) \n {contained_gtf_file} (non-overlaped contained genes) \n {containedALL_gtf_file} ( ALL contained genes) \n {clean_gtf_file} (clean)")
logging.info(f"GTF files saved: \n {operon_gtf_file} (operons) \n {contained_gtf_file} (non-overlaped contained genes) \n {containedALL_gtf_file} ( ALL contained genes) \n {clean_gtf_file} (clean)")

################
# Counting the number of genes by operon
operon_size_counter = Counter()
true_final_pairs_DEF = []

for operon, _, transcript in final_pairs_DEF:
    true_final_pairs_DEF.append((operon,transcript))

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
