import gffutils # type: ignore
import sys
import os
from collections import defaultdict, Counter

# Check if the user provided a GTF file
if len(sys.argv) != 2:
    print("Usage: python Longest_transcript_filter.py <input.gtf>")
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
output_file = os.path.splitext(gtf_file)[0] + "_longest_trans_only.gtf"

# Get all transcripts grouped by chromosome
chrom_transcripts = defaultdict(list)
for transcript in db.features_of_type("transcript"):
    chrom_transcripts[transcript.chrom].append(transcript)

# Store all transcripts per gene_id per chromosome
gene_transcripts = defaultdict(lambda: defaultdict(list))

# First, organize transcripts per gene per chromosome
for chrom, transcripts in chrom_transcripts.items():
    for transcript in transcripts:
        gene_id = transcript.attributes['gene_id'][0]
        gene_transcripts[chrom][gene_id].append(transcript)

# Now, select the longest transcript for each gene on each chromosome
longest_transcripts = []
for chrom in gene_transcripts:
    print(f"Processing {chrom}...")
    for gene_id, tx_list in gene_transcripts[chrom].items():
        max_len = 0
        longest_tx = None
        for tx in tx_list:
            length = db.children_bp(transcript, child_featuretype='exon')  # sum of exon lengths
            if length > max_len:
                max_len = length
                longest_tx = tx
        if longest_tx:
            longest_transcripts.append(longest_tx.id)

# Write longest transcripts GTF
with open(output_file, "w") as contained_out:
    for transcript_id in longest_transcripts:
        transcript_feature = db[transcript_id]
        contained_out.write(str(transcript_feature) + "\n")
        for feature in db.children(transcript_id, featuretype='exon', order_by='start'):
            contained_out.write(str(feature) + "\n")

print(f"GTF with only longest transcripts saved to {output_file}")
