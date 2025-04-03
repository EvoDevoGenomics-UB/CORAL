import gffutils # type: ignore
import sys
import os
from collections import defaultdict, Counter

# Check if the user provided a GTF file
if len(sys.argv) != 2:
    print("Usage: python detect_operons.py <input.gtf>")
    sys.exit(1)

# Get the input GTF file from the command line
gtf_file = sys.argv[1]

# Ensure the file exists
if not os.path.isfile(gtf_file):
    print(f"Error: File '{gtf_file}' not found.")
    sys.exit(1)

# Create a database from the GTF file (stored in memory)
db_filename = os.path.splitext(gtf_file)[0] + "_annotation.db"
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
output_file = os.path.splitext(gtf_file)[0] + "_operons_found_v3.tsv"

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
        #progress = (idx / len(transcripts)) * 100
        #print(f"Chrom {chrom} Progress: {progress:.2f}%", end="\r")

        for sub_transcript in transcripts:  # Now only compares within the same chromosome
            if transcript.id == sub_transcript.id:
                continue
            if transcript.strand == sub_transcript.strand:
                if transcript.start <= (sub_transcript.start + 25) and transcript.end >= (sub_transcript.end - 25):
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

print(f"Filtered TSV file saved as {output_file}")
