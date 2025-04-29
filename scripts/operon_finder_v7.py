import gffutils # type: ignore
import sys
import os
import argparse
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
    "--threshold", 
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

# Ensure the file exists
if not os.path.isfile(gtf_file):
    print(f"Error: File '{gtf_file}' not found.")
    sys.exit(1)

# Create a database from the GTF file (stored in memory)
db_filename = os.path.splitext(gtf_file)[0] + "_annotation_v7.t"+ str(threshold) +".db"
db = gffutils.create_db(
    gtf_file,
    dbfn=db_filename,
    force=True,
    keep_order=True,
    disable_infer_transcripts=True,
    disable_infer_genes=True
)
print(f"File '{db_filename}' created.")

# Define output file
output_file = out_prefix + "_operons_found_v7.t"+ str(threshold) +".tsv"

# Dictionary to store transcripts per chromosome
chrom_transcripts = defaultdict(list)

# Organize transcripts by chromosome
for transcript in db.features_of_type("transcript"):
    chrom_transcripts[transcript.chrom].append(transcript)

# Store detected transcript pairs
contained_pairs = []

# Find contained transcripts (with progress tracking)
for chrom, transcripts in chrom_transcripts.items():
    print(f"Processing {chrom} ({len(transcripts)} transcripts)...")
    
    for idx, transcript in enumerate(transcripts, 1):
        if idx % max(1, len(transcripts) // 20) == 0:  # Print progress every 5% intervals
            progress = (idx / len(transcripts)) * 100
            print(f"Chrom {chrom} Progress: {progress:.1f}%", end="\r")

        for sub_transcript in transcripts:  # Now only compares within the same chromosome
            if transcript.id == sub_transcript.id:
                continue # Skip when comparing itself
            if transcript.strand == sub_transcript.strand:
                if float(transcript.attributes['cov'][0])* threshold <= float(sub_transcript.attributes['cov'][0]) :
                    if transcript.start <= (sub_transcript.start + 150 ) and transcript.end >= (sub_transcript.end + 200):
                        exons = len(list(db.children(sub_transcript.id, featuretype='exon')))
                        if exons > 1:
                            contained_pairs.append((transcript.id, sub_transcript.id, sub_transcript.start, sub_transcript.end))
                    if transcript.start <= (sub_transcript.start - 200 ) and transcript.end >= (sub_transcript.end - 150):
                        exons = len(list(db.children(sub_transcript.id, featuretype='exon')))
                        if exons > 1:
                            contained_pairs.append((transcript.id, sub_transcript.id, sub_transcript.start, sub_transcript.end))

# Create a set of contained transcript IDs
contained_transcripts = {pair[1] for pair in contained_pairs}

# Filter out rows where the first column appears in the second column
filtered_pairs = [pair for pair in contained_pairs if pair[0] not in contained_transcripts]

# Count how many transcripts each operon contains
operon_counts = Counter(operon for operon, _, _, _ in filtered_pairs)

# Keep only operons that contain **two or more** transcripts
valid_operons = {operon for operon, count in operon_counts.items() if count > 1}

# Group contained transcripts by operon
operon_to_transcripts = defaultdict(list)
for operon, transcript, start, end in filtered_pairs:
    if operon in valid_operons:
        operon_to_transcripts[operon].append((transcript, start, end))

# Remove overlapping contained transcripts **within the same operon**
final_pairs = []
for operon, transcript_list in operon_to_transcripts.items():
    # Sort transcripts by start position
    transcript_list.sort(key=lambda x: x[1])  # Sort by start coordinate
    
    non_overlapping = []
    for current_transcript in transcript_list:
        transcript_id, start, end = current_transcript
        if not non_overlapping or start > non_overlapping[-1][2]:  # No overlap with previous
            non_overlapping.append(current_transcript)

    # Add non-overlapping transcripts to final output
    for transcript_id, _, _ in non_overlapping:
        final_pairs.append((operon, transcript_id))

# Count how many transcripts each operon contains
operon_counts_def = Counter(operon for operon, _ in final_pairs)
# Keep only operons that contain **two or more** transcripts
out_operons = {operon for operon, count in operon_counts_def.items() if count < 2}
# Filter out rows where the first column appears in the second column
final_pairs_DEF = [pair for pair in final_pairs if pair[0] not in out_operons]

# Write to output file
with open(output_file, "w") as out_file:
    out_file.write("Operon\tContained_transcript\n")
    for operon, transcript in final_pairs_DEF:
        out_file.write(f"{operon}\t{transcript}\n")

print(f"Operon-Genes found saved to TSV file {output_file}")

# Create sets of operon and gene transcript IDs
operon_ids = {operon for operon, _ in final_pairs_DEF}
contained_ids = {transcript for _, transcript in final_pairs_DEF}

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
operon_gtf_file = out_prefix + "_Operons_v7.t" + str(threshold) + ".gtf"
contained_gtf_file = out_prefix + "_OperonGenes_v7.t" + str(threshold) + ".gtf"
containedALL_gtf_file = out_prefix + "_OperonGenesALL_v7.t" + str(threshold) + ".gtf"
clean_gtf_file = out_prefix + "_opCLEAN_v7.t" + str(threshold) + ".gtf"

# Write operon transcripts to GTF
with open(operon_gtf_file, "w") as operon_out:
    for operon_id in operon_ids:
        operon_feature = db[operon_id]
        # Write the transcript feature
        operon_out.write(str(operon_feature) + "\n")
        # Write its child features (e.g., exons)
        for feature in db.children(operon_id, featuretype='exon', order_by='start'):
            operon_out.write(str(feature) + "\n")

# Write non-overlaped contained gene transcripts to GTF
with open(contained_gtf_file, "w") as contained_out:
    for transcript_id in contained_ids:
        transcript_feature = db[transcript_id]
        contained_out.write(str(transcript_feature) + "\n")
        for feature in db.children(transcript_id, featuretype='exon', order_by='start'):
            contained_out.write(str(feature) + "\n")
# Write ALL contained gene transcripts to GTF
with open(containedALL_gtf_file, "w") as containedALL_out:
    for transcript_id in gene_transcripts:
        transcript_feature = db[transcript_id]
        containedALL_out.write(str(transcript_feature) + "\n")
        for feature in db.children(transcript_id, featuretype='exon', order_by='start'):
            containedALL_out.write(str(feature) + "\n")

# Write Clean GTF no containgin OPRNs nor OpGenes
with open(clean_gtf_file, "w") as clean_out:
    for transcript in db.features_of_type("transcript"):
        if transcript.id not in operon_ids and transcript.id not in gene_transcripts:
            transcript_feature = db[transcript.id]
            clean_out.write(str(transcript_feature) + "\n")
            for feature in db.children(transcript.id, featuretype='exon', order_by='start'):
               clean_out.write(str(feature) + "\n")

print(f"GTF files saved: \n {operon_gtf_file} (operons) \n {contained_gtf_file} (non-overlaped contained genes) \n {containedALL_gtf_file} ( ALL contained genes) \n {clean_gtf_file} (clean)")
