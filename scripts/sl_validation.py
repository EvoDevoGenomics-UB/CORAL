#!/usr/bin/env python3
import gffutils # type: ignore
import sys
import os
import argparse
import logging
import re
import numpy as np
from collections import defaultdict, Counter
import pandas as pd
from rpy2 import robjects
from rpy2.robjects import conversion, default_converter, pandas2ri

# Argument parser setup
parser = argparse.ArgumentParser(
    description="Count operons and their included genes from a GTF file."
)
parser.add_argument(
    "-gtf", "--gtf", 
    required=True, 
    help="Path to the GTF file use as reference."
)
parser.add_argument(
    "-sl", "--sl", 
    required=True, 
    help="Path to the input GTF file with SL information."
)
parser.add_argument(
    "-p","--prefix",
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
gtf_file = args.gtf
sl_file = args.sl

# If no output was specified, use the input base name
if args.prefix:
    out_prefix = args.prefix
else:
    out_prefix = os.path.splitext(gtf_file)[0]

#Set log output file
if args.log:
    log_file = args.log
else:
    log_file = os.path.splitext(gtf_file)[0] + "_SLvalidation.log"
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
# Ensure the file exists
if not os.path.isfile(sl_file):
    print(f"Error: File '{sl_file}' not found.")
    logging.error(f'File {sl_file} not found.')
    sys.exit(1)

# Create a database from the GTF file (stored in memory)
db_filename = os.path.splitext(gtf_file)[0] + "_SLvalidation.db"
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

# Create a database from the GTF file (stored in memory)
db_slname = os.path.splitext(sl_file)[0] + "_SLvalidation.db"
db_sl = gffutils.create_db(
    sl_file,
    dbfn=db_slname,
    force=True,
    keep_order=True,
    disable_infer_transcripts=True,
    disable_infer_genes=True
)
print(f"File '{db_slname}' created.")
logging.info(f"File '{db_slname}' created.")

# Define output file
output_file_counts = out_prefix + "_SLvalidation.counts.tsv"

#########################################
# Organize transcripts and SLs by chromosome
chrom_transcripts = defaultdict(list)
chrom_sl_sites = defaultdict(list)

for transcript in db.features_of_type("transcript"):
    chrom_transcripts[transcript.chrom].append(transcript)

for sl_site in db_sl.features_of_type("transcript"):
    chrom_sl_sites[sl_site.chrom].append(sl_site)

#########################################
# Find SL sites matching OpG 5' ends inside operons
contained_pairs = []  # (chrom, strand, operon_id, opg_gene_id, sl_id, sl_pos)
tolerance = 50  # bp window around OpG 5' end

for chrom, transcripts in chrom_transcripts.items():
    if chrom not in chrom_sl_sites:
        continue

    print(f"Checking SL sites on {chrom}...")
    logging.info(f"Checking SL sites on {chrom}...")

    # Build operon groups for this chromosome
    operon_groups = defaultdict(list)
    for t in transcripts:
        if "operon_id" in t.attributes:
            op_id = t.attributes["operon_id"][0]
            operon_groups[op_id].append(t)

    # For each SL site, check which operon it's in
    for sl in chrom_sl_sites[chrom]:
        sl_start, sl_end, sl_strand, sl_id = sl.start, sl.end, sl.strand, sl.id

        for op_id in operon_groups.keys():
            for opg in operon_groups[op_id]:
                if opg.strand != sl_strand:
                    continue
                # Check if SL is within OpG coordinates
                if opg.start <= sl_start <= opg.end:
                    opg_gene_id = opg.attributes["gene_id"][0]

                    if opg.strand == "+":
                        opg_5p = opg.start
                        distance = sl_start - opg_5p
                    if opg.strand == "-":
                        opg_5p = opg.end
                        distance = opg_5p - sl_end
                    # Check SL proximity to OpG 5' end
                    if distance <= tolerance:
                        contained_pairs.append(
                            (chrom, sl_strand, op_id, opg_gene_id, sl_id, sl_start)
                        )
                        logging.info(
                            f"SL {sl_id} at {chrom}:{sl_start} matches OpG {opg_gene_id} (5' of {op_id})"
                        )

# Write match results
output_file = out_prefix + "_SLvalidation.matches.tsv"
with open(output_file, "w") as out:
    out.write("Chrom\tStrand\tOperon\tOpG\tSL_ID\tSL_Position\n")
    for chrom, strand, operon, opg, sl_id, sl_pos in contained_pairs:
        out.write(f"{chrom}\t{strand}\t{operon}\t{opg}\t{sl_id}\t{sl_pos}\n")

print(f"SL–OpG matches saved to {output_file}")
logging.info(f"SL–OpG matches saved to {output_file}")

#########################################
# Summary statistics
logging.info("Generating summary statistics for SL–OpG matches...")

# Build mappings from the GTF itself
operon_to_opgs = defaultdict(set)
for transcript in db.features_of_type("transcript"):
    if "operon_id" in transcript.attributes:
        op_id = transcript.attributes["operon_id"][0]
        gene_id = transcript.attributes["gene_id"][0]
        if op_id != gene_id:  # skip the operon itself
            operon_to_opgs[op_id].add(gene_id)

# Build dictionary of OpGs that have SLs
operon_with_sl = defaultdict(set)
for _, _, op_id, opg_gene_id, _, _ in contained_pairs:
    operon_with_sl[op_id].add(opg_gene_id)

# Compute stats
total_opgs = sum(len(v) for v in operon_to_opgs.values())
total_opgs_with_sl = sum(len(v) for v in operon_with_sl.values())
total_opgs_without_sl = total_opgs - total_opgs_with_sl
percent_opgs_with_sl = (total_opgs_with_sl / total_opgs) * 100

total_oprns = len(operon_to_opgs)
oprns_with_sl = len(operon_with_sl)
oprns_without_sl = total_oprns - oprns_with_sl
percent_oprns_with_sl = ( oprns_with_sl / total_oprns ) * 100

# Average percentage of OpGs with SL per operon
percentages = []
for op_id, opgs in operon_to_opgs.items():
    n_total = len(opgs)
    n_with = len(operon_with_sl.get(op_id, []))
    if n_total > 0:
        percentages.append((op_id, n_total, ((n_with / n_total) * 100)))
#avg_percentage = sum(percentages) / len(percentages) if percentages else 0.0
percentages_file = out_prefix + "_SLvalidation.percentages.tsv"
with open(percentages_file, "w") as out_prcent:
    out_prcent.write(f"OpID\tNumOpGs\tPercentage\n")
    for op_id, n_opgs, percent in percentages:
        out_prcent.write(f"{op_id}\t{n_opgs}\t{percent}\n")
oprn_summary = pd.read_csv(percentages_file, sep = "\t")
# Compute average and standard deviation
mean_percent = oprn_summary["Percentage"].mean()
std_percent = oprn_summary["Percentage"].std()

# Write summary table
summary_file = out_prefix + "_SLvalidation.summary.tsv"
with open(summary_file, "w") as out_sum:
    out_sum.write("Metric\tValue\n")
    out_sum.write(f"Total OpGs:\t{total_opgs}\n")
    out_sum.write(f"OpGs with SL:\t{total_opgs_with_sl} ({percent_opgs_with_sl:.2f}%)\n")
    out_sum.write(f"OpGs without SL:\t{total_opgs_without_sl}\n")
    out_sum.write(f"Total OPRNs:\t{total_oprns}\n")
    out_sum.write(f"OPRNs with SL sites:\t{oprns_with_sl} ({percent_oprns_with_sl:.2f}%)\n")
    out_sum.write(f"OPRNs without SL sites:\t{oprns_without_sl}\n")
    out_sum.write(f"Average % of OpGs with SLs per OPRN: {mean_percent:.2f}% ± {std_percent:.2f}%")

# -----------------------------
# Operons with SL summary by number of genes
# -----------------------------
# Create bins
bins = {
    "2-genes": {"with_sl": 0, "total": 0},
    "3-genes": {"with_sl": 0, "total": 0},
    "4-genes": {"with_sl": 0, "total": 0},
    "5-genes": {"with_sl": 0, "total": 0},
    ">5-genes": {"with_sl": 0, "total": 0}
}

for op_id, genes in operon_to_opgs.items():
    n_genes = len(genes)
    # Determine bin
    if n_genes == 2:
        bin_key = "2-genes"
    elif n_genes == 3:
        bin_key = "3-genes"
    elif n_genes == 4:
        bin_key = "4-genes"
    elif n_genes == 5:
        bin_key = "5-genes"
    elif n_genes > 5:
        bin_key = ">5-genes"
    else:
        continue  # Skip operons with 1 gene (if any)

    bins[bin_key]["total"] += 1
    if op_id in operon_with_sl and len(operon_with_sl[op_id]) > 0:
        bins[bin_key]["with_sl"] += 1

# Append results to summary file
with open(summary_file, "a") as out_sum:
    out_sum.write("\nOperons with SL by number of genes:\n")
    for k in ["2-genes","3-genes","4-genes","5-genes",">5-genes"]:
        out_sum.write(f"{k}: {bins[k]['with_sl']} / {bins[k]['total']}\n")

print("Operons by gene count with SL / total written to summary.tsv")
logging.info("Operons by gene count with SL / total written to summary.tsv")


print("Operons by number of genes with SL added to summary.tsv")
logging.info("Operons by number of genes with SL added to summary.tsv")


print(f"Summary written to {summary_file}")
logging.info(f"Summary written to {summary_file}\n")
print(f"OPRNs with SL sites:\t{oprns_with_sl} ({percent_oprns_with_sl:.2f}%)")
logging.info(f"OPRNs with SL sites:\t{oprns_with_sl} ({percent_oprns_with_sl:.2f}%)")
print(f"Average % of OpGs with SL sites per OPRN: {mean_percent:.2f}% ± {std_percent:.2f}%\n")
logging.info(f"Average % of OpGs with SL sites per OPRN: {mean_percent:.2f}% ± {std_percent:.2f}%")

# Create safe converter for pandas → R
converter = conversion.Converter('pandas2ri', template=default_converter + pandas2ri.converter)
# Convert safely to R data frame
with converter.context():
    robjects.globalenv['a'] = oprn_summary

# Run R code to create the violin plot
robjects.r('''
library(ggplot2)
library(dplyr)

a$split <- ifelse(a$Percentage==0, "WihtOutSL", "WithSL")
a$split <- factor(a$split, levels = c("WithSL","WihtOutSL"))

p1 <- ggplot(a[a$split=="WithSL",], aes(x = split, y = Percentage)) +
  geom_violin(fill = "lightgray", alpha = 1, color = "darkgray") +
  geom_boxplot(width = 0.05, alpha = 0.6, outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 100)) +
  # Median
  stat_summary(fun = mean, geom = "point", color = "darkgreen", size = 0.5) +
  stat_summary(fun = mean, geom = "text",
               aes(label = format(round(..y.., 2), big.mark = ",")),
               vjust = -1, color = "darkgreen", size = 3, fontface = "bold") +
  # Min
  stat_summary(fun = min, geom = "text",
               aes(label = paste0("Min: ", format(round(..y.., 2), big.mark = ","))),
               vjust = 1, color = "blue", size = 3, fontface = "italic") +
  # Max
  stat_summary(fun = max, geom = "text",
               aes(label = paste0("Max: ", format(round(..y.., 2), big.mark = ","))),
               vjust = -0.25, color = "red", size = 3, fontface = "italic") +
  ylab("%OpGs with SL") +
  xlab("NOR_OPRNs") +
  theme_minimal()

ggsave("Violinplot_PercentOpGs_withSL_new.png", p1, width = 3, height = 5, dpi = 300)
ggsave("Violinplot_PercentOpGs_withSL_new.svg", p1, width = 3, height = 5, dpi = 300)

a$group <- "All"
g2 <- a[a$NumOpGs==2 & a$Percentage!=0,]
g2$group <- "2 OpGs"
g3 <- a[a$NumOpGs==3 & a$Percentage!=0,]
g3$group <- "3 OpGs"
g4 <- a[a$NumOpGs==4 & a$Percentage!=0,]
g4$group <- "4 OpGs"
g5 <- a[a$NumOpGs==5 & a$Percentage!=0,]
g5$group <- "5 OpGs" 
over5g <- a[a$NumOpGs>5 & a$Percentage!=0,]
over5g$group <- ">5 OpGs"
b <- bind_rows(g2,g3,g4,g5,over5g,a[a$Percentage!=0,])
b$group <- factor(b$group, levels = c("All","2 OpGs","3 OpGs","4 OpGs","5 OpGs",">5 OpGs"))

p2 <- ggplot(b, aes(x = group, y = Percentage)) +
  geom_violin(fill = "lightgray", alpha = 0.8, color = "darkgray") +
  geom_boxplot(width = 0.05, alpha = 0.6, outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 100)) +
  # Median
  stat_summary(fun = mean, geom = "point", color = "darkgreen", size = 0.5) +
  stat_summary(fun = mean, geom = "text",
               aes(label = format(round(..y.., 2), big.mark = ",")),
               vjust = -1, color = "darkgreen", size = 3, fontface = "bold") +
  # Min
  stat_summary(fun = min, geom = "text",
               aes(label = paste0("Min: ", format(round(..y.., 2), big.mark = ","))),
               vjust = 1, color = "blue", size = 3, fontface = "italic") +
  # Max
  stat_summary(fun = max, geom = "text",
               aes(label = paste0("Max: ", format(round(..y.., 2), big.mark = ","))),
               vjust = -0.25, color = "red", size = 3, fontface = "italic") +
  ylab("%OpG with SL") +
  xlab("NOR_OPRNs") +
  theme_minimal()

ggsave("SL_violin_plot_R.png", p2, width = 7, height = 4, dpi = 300)
ggsave("SL_violin_plot_R.svg", p2, width = 7, height = 4, dpi = 300)
''')

print("Violin plots generated with R (saved as SL_violin_plot_R.png and Violinplot_PercentOpGs_withSL_new.png)")
