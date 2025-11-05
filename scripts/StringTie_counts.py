#!/usr/bin/env python3
import subprocess
import sys
import os
import logging
import argparse

# Argument parser setup
parser = argparse.ArgumentParser(
    description="Creates StringTie -e GTF files and the sample list to run prepDE.py3"
)
parser.add_argument(
    "-f", "--file", 
    required=True, 
    help="Path to the input GTF file."
)
parser.add_argument(
    "-b", "--bam",
    nargs="+", 
    required=True, 
    help="List of BAMs to use."
)
parser.add_argument(
    "--outdir",
    help="Prefix of the output directory [default: input GTF base name]."
)
parser.add_argument(
    "-t","--threads",
    type=int,
    default=2,
    help="Number of threads to use [default: 2]"
)
#parser.add_argument(
#    "--log",
#    help="Log file [default: GTF basename + strigntie_counts.log]."
#)

# Parse arguments
args = parser.parse_args()
# Extract values
gtf_file = args.file
threads = args.threads

# If no output was specified, use the input base name
if args.outdir:
    outdir = args.outdir / os.path.splitext(os.path.basename(gtf_file))[0]
else:
    outdir = os.path.splitext(os.path.basename(gtf_file))[0]
#Set log output file
if args.log:
    log_file = args.log
else:
    log_file = os.path.splitext(os.path.basename(gtf_file))[0] + "_stringtie_counts.log"

#Define logger
logging.basicConfig(filename= log_file, 
					format='%(asctime)s %(levelname)s - %(message)s', 
					filemode='w',
                    level=logging.INFO) 

# Ensure the file exists
if not os.path.isfile(gtf_file):
    logging.error(f'File {gtf_file} not found.')
    sys.exit(1)

logging.info(f"Using GTF: {gtf_file}")
logging.info(f"Output directory: {outdir}")

outdir.mkdir(parents=True, exist_ok=True)
list_path = outdir / "ALL_sample_list.txt"

with open(list_path, "w") as f_list:
    for bam in args.bam:
        bam_name = os.path.splitext(os.path.basename(bam))[0]
        logging.info(f"Processing: {bam}.bam")
        sample = bam_name.stem.split("_reads_aln_")[0]
        sample_dir = outdir / sample
        sample_dir.mkdir(exist_ok=True)
        out_gtf = sample_dir / f"{sample}.gtf"

        if not out_gtf.exists():
            logging.info(f"Running stringtie for sample {sample}")
            subprocess.run([
                "stringtie", "-eB", "-G", str(gtf_file),
                "-p", str(threads), "-o", str(out_gtf), str(bam)
            ], check=True)
        else:
            logging.info(f"GTF for {sample} already exists — skipping.")

        # Build sample list for prepDE.py3    
        f_list.write(f"{sample}\t{out_gtf}\n")

logging.info(f"Sample list file created")
print(f"{list_path}")
