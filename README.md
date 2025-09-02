# compact-genome-annotation

This is a set of scripts to <b>annotate compact genomes</b> using <b>long-read RNAseq data</b>.

It uses as input clean (primer-trimmed) pre-processed fastq files and maps them to the given genome using <i>Minimap2</i>.
Then creates non-assembled annotation for each samples using <i>Stringtie (v3.0)</i> and look for the potential operons within those annotations (<i>operon-finder-rust v1.3.1</i>). After identifing operon transcripts, operon-conteined transcripts, and non-operon-related transcriptsit, it generates a consensus annotations for the three sets of transcirpts.

Those sets are merge together to generate the final consensus annotation. From the annotation no containg operons there is selected the longest annotation for each gene (<i>Longest_transcript_filter.py</i>) to generate a fasta file, which is use to asses the quality of the annotation with <i>BUSCO (v5.8)</i>. The transcriptome obtained from the consensus annotations are also assed with <i>BUSCO</i>.

<img width="1795" height="2376" alt="Figure1" src="https://github.com/user-attachments/assets/1e2ed9b4-c782-410e-8d5b-5aadafac6791" />
