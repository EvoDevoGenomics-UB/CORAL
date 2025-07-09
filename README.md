# compact-genome-annotation

This is a set of scripts to <b>annotate compact genomes</b> using <b>long-read RNAseq data</b>.

It uses as input clean (primer-trimmed) pre-processed fastq files and maps it to the given genome using <i>Minimap2</i>.
Then creats non-assembled annotation for each samples using <i>Stringtie (v3.0)</i> and look for the potential operons within those annotations (<i>operon_finder_v9.7.py</i>). After identifing operon transcripts, operon-conteined transcripts, and non-operon-related transcriptsit, it generates a consensus annotations for the three sets of transcirpts.

Those sets are merge together to generate the final consensus annotation. From the annotation no containg operons there is selected the longest annotation for each gene (<i>Longest_transcript_filter.py</i>) to generate a fasta file, which is use to asses the quality of the annotation with <i>BUSCO (v5.8)</i>. The transcriptome obtained from the consensus annotations are also assed with <i>BUSCO</i>.
