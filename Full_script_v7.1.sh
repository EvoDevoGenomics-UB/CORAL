#!/bin/bash
##Note: it is importnat to have an enviroment of conda active that contains:
# - minimap2
# - samtools
# - seqkit
# - Gffread
# - Stringtie (v3.0)

###Setting the local path to the stringtie executable:
path_stringtie=/data2/cristian/ntorres2/LongRead_RNAseq/Stringtie_tests_BP10new/stringtie/stringtie
gffread_path=/data2/cristian/ntorres2/LongRead_RNAseq/Workspaces/BAR_longreads_C05C06/STRINGTIE_annotations_v2/gffread/gffread
#gffcompare_path=/data2/cristian/ntorres2/LongRead_RNAseq/Workspaces/BAR_longreads_C05C06/STRINGTIE_annotations_v2/gffcompare/gffcompare
WORKDIR=/data2/cristian/ntorres2/LongRead_RNAseq/Workspaces/BAR_longreads_C05C06
VERSION=$3
operon_finder=./operon_finder_v6.py 
longest_filter=./Longest_transcript_filter.py
GENOME="/data2/cristian/ntorres2/LongRead_RNAseq/Snakemake_Tutorial/REF_genomes/Bar2_p4.Flye.masked.fa"
THRESHOLD=$4 #Value used in operon_finder
#$1 = samples
#$2 = max intron size

#Minimap2 alingment of the reads
echo "Performing mapping of the reads..."
for i in $1 ; do for x in $2 ; do \
    if [ $i == "C04" ] ; then
        echo "$i SKIPED!"
        continue
    fi
    output_bam="${WORKDIR}/alignments/FAX17881_1_${i}_reads_aln_sorted_v${x}.bam"
    if [ -f "$output_bam" ]; then
        echo "BAM file with -G ${x} for $i already exists. Skipping..."
        continue
    fi
    minimap2 -t 8 -ax splice \
    -G ${x} ${WORKDIR}/index/genome_index.mmi \
    ${WORKDIR}/processed_reads/FAX17881_1_${i}_full_length_reads.clean.fq > ${i}_reads_aln_v${x}.sam ; \
    samtools view ${i}_reads_aln_v${x}.sam -O BAM -o ${i}_reads_aln_v${x}.bam ; \
    seqkit bam -j 8 -q 40 -x - ${i}_reads_aln_v${x}.bam | samtools sort -@ 8 -O BAM -o $output_bam ;
    samtools index $output_bam ; \
    rm ${i}_reads_aln_v${x}.sam ; \
done ; done

#Stringtie v3
echo "Executing non-guide StringTie..."
for i in $1 ; do for x in $2 ; do for cov in 1.5 ; do 
 output_gtf=${i}_no_assembled_v${x}_c${cov}_f0.2.gtf
 cov_single_exon=$cov*2
    if [ -f "$output_gtf" ]; then
        echo "Already done for ${i}_cov${cov}. Skipping..."
        continue
    fi
 $path_stringtie --rf -L -R -p 8 -M 0.75 -j 2 -a 15\
 -c $cov -f 0.2 -s $cov_single_exon \
 -o $output_gtf ${WORKDIR}/alignments/FAX17881_1_${i}_reads_aln_sorted_v${x}.bam ; \
 echo "Stringtie no-guide no-assembly ${i}_${x}_cov${cov} done" ; \
done ; done ; done

#Python script
echo "Performing operon search on each sample..."
for i in $1 ; do for x in $2 ; do
    for cov in 1.5 ; do \
        sample_gtf=${i}_no_assembled_v${x}_c${cov}_f0.2
        operon_file=${sample_gtf}_operons_found_v6.t${THRESHOLD}.tsv
        if [ -f "$operon_file" ]; then
            echo "Operon search already done for ${i}_v${x}_cov${cov}. Skipping..."
            continue
        fi
        python ${operon_finder} -f ${sample_gtf}.gtf --threshold ${THRESHOLD}; \
    done ; \
done ; done

#Make suere format is okay
echo "Sanatizing gtf files..."
for i in $1 ; do for x in $2 ; do for cov in 1.5 ; do
    for type in Operons OperonGenes opCLEAN ; do \
        sample_gtf=${i}_no_assembled_v${x}_c${cov}_f0.2
        gtf_file=${sample_gtf}_${type}_v6.t${THRESHOLD}
        if [ -f "$gtf_file.${VERSION}.clean.gtf" ]; then
            echo "Already done for $type file of ${i}_v${x}_cov${cov}. Skipping..."
            continue
        fi
        awk '{if($4>$5) print $1,$2,$3,$5,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18 ; else print $0}' $gtf_file.gtf > $gtf_file.clean.tmp.gtf; \
        $gffread_path --sort-alpha -F -T -o $gtf_file.${VERSION}.clean.gtf $gtf_file.clean.tmp.gtf ; \
        rm $gtf_file.clean.tmp.gtf ; \
    done ; \
    echo "  ${i}_v${x}_cov${cov} done!" ;\
done ; done ; done

#Creating annoations with potential operons
echo "Creating annoations of potential operons"
for x in $2 ; do for cov in 1.5 ; do \
    output_operons_clean=OdiBAR_LRannot_v${x}_c${cov}_OPRNs.${VERSION}.gtf
    if [ -f "$output_operons_clean" ]; then
        echo "Skipping..."
        continue
    fi
    ls *OperonGenes_v6.t${THRESHOLD}.${VERSION}.clean.gtf > List_merge_OpGs.${VERSION}.txt ; \
    $path_stringtie --merge -l OpG -f 0.2 -F 0 -T 0 -c 0 -g '-50' -o OdiBAR_LRannot_v${x}_c${cov}_OpGs.${VERSION}.gtf List_merge_OpGs.${VERSION}.txt; \
    ls *Operons_v6.t${THRESHOLD}.${VERSION}.clean.gtf > List_merge_OPRNs.${VERSION}.txt ; \
    $path_stringtie --merge -l OPRN -f 0.2 -F 0 -T 0 -c 0 -g 0 -o $output_operons_clean List_merge_OPRNs.${VERSION}.txt ; \
done ; done

#Creating final GTFs
echo 'Creating final GTFs'
for x in $2 ; do for cov in 1.5 ; do
    final_gtf_file=OdiBAR_LRannot_v${x}_c${cov}_StringtieMerge.clean-and-OPRNs.${VERSION}.gtf
    if [ -f "$final_gtf_file" ]; then
        echo "GTF already done. Skipping..."
        continue
    fi
    cat OdiBAR_LRannot_v${x}_c${cov}_OPRNs.${VERSION}.gtf OdiBAR_LRannot_v${x}_c${cov}_OpGs.${VERSION}.gtf > tmp.gtf ; \
    $gffread_path --sort-alpha -F -T -o Merged_OPRNs-and-OpGs_v${x}_c${cov}.${VERSION}.sorted.gtf tmp.gtf ; rm tmp.gtf ; \
    ls *opCLEAN_v6.t${THRESHOLD}.${VERSION}.clean.gtf > List_merge_opCLEAN.${VERSION}.txt ; \
    $path_stringtie --merge List_merge_opCLEAN.${VERSION}.txt \
     -l g -f 0.2 -F 0 -T 0 -c 0 -g '-50' \
     -o OdiBAR_LRannot_v${x}_c${cov}_mergeCLEAN.${VERSION}.gtf ; \
    echo "  Final merge CLEAN done" ; \
    $path_stringtie --merge List_merge_opCLEAN.${VERSION}.txt \
     -G OdiBAR_LRannot_v${x}_c${cov}_OpGs.${VERSION}.gtf \
     -l g -f 0.2 -F 0 -T 0 -c 0 -g '-50' \
     -o OdiBAR_LRannot_v${x}_c${cov}_StringtieMerge.clean-noOPRNs.${VERSION}.gtf ; \
    echo "  Final CLEAN-noOPRNs done" ; \
    $path_stringtie --merge List_merge_opCLEAN.${VERSION}.txt \
     -G Merged_OPRNs-and-OpGs_v${x}_c${cov}.${VERSION}.sorted.gtf \
     -l g -f 0.2 -F 0 -T 0 -c 0 -g '-50' \
     -o $final_gtf_file ; \
    echo "  Final CLEAN-and-OPRNs done" ; \
done ; done

#Obtaining fasta files for BUSCO analysis
echo "Obtianing fasta file with only longest transcript per gene..."
for x in $2 ; do for cov in 1.5 ; do for final in clean-noOPRNs ; do
    file_name=OdiBAR_LRannot_v${x}_c${cov}_StringtieMerge.${final}.${VERSION}
    final_fasta_file=$file_name.fasta
    if [ -f "$final_fasta_file" ]; then
        echo "Fasta already done. Skipping $final ..."
        continue
    fi
    #OdiBAR_LRannot_v${x}_c${cov}_mergeCLEAN.${VERSION}.gtf
    python $longest_filter $file_name.gtf ; \
    $gffread_path -g $GENOME -w $final_fasta_file ${file_name}_longest_trans_only.gtf ; \
done ; done ; done

#Obtaining coverage of final annotation
echo "Obtaining coverage of final annotation..."
SAMPLE_BAMs=$(for i in $1 ; do for x in $2 ; do \
    if [ $i == "C04" ] ; then
        echo "/data2/cristian/ntorres2/LongRead_RNAseq/Workspaces/BAR_longreads_C04/alignments/OdiBAR_C04_pool_reads_aln_v10k.sorted.bam"
        continue
    fi
    output_bam="${WORKDIR}/alignments/FAX17881_1_${i}_reads_aln_sorted_v${x}.bam" ; \
    echo "$output_bam"
done ; done )
for x in $2 ; do for cov in 1.5 ; do
    final_gtf_file_name=OdiBAR_LRannot_v${x}_c${cov}_StringtieMerge.clean-and-OPRNs.${VERSION}
    if [ -f "$final_gtf_file_name.counts.gtf" ]; then
        echo "GTF already done. Skipping..."
        continue
    fi
    $path_stringtie -G $final_gtf_file_name.gtf -e -o $final_gtf_file_name.counts.gtf $SAMPLE_BAMs ; \
done ; done

#Final check of operons in the annotation.
echo "Performing FINAL operon search..."
for x in $2 ; do for cov in 1.5 ; do \
    final_gtf_file_name=OdiBAR_LRannot_v${x}_c${cov}_StringtieMerge.clean-and-OPRNs.${VERSION}.counts
    final_operon_file=${final_gtf_file_name}_operons_found_v6.tsv
    if [ -f "$final_operon_file" ]; then
        echo "Operon search already done. Skipping..."
        continue
    fi
    python ${operon_finder} -f $final_gtf_file_name.gtf --threshold $4 ; \
done ; done

#Moving all files to a organized folder:
for x in $2 ; do for cov in 1.5 ; do \
    final_gtf_file_name=OdiBAR_LRannot_v${x}_c${cov}_StringtieMerge.clean-and-OPRNs.${VERSION}.counts
    final_operon_file=${final_gtf_file_name}_operons_found_v6.tsv
    if [ -f "$final_operon_file" ]; then
#        rm List*.tmp ; rm *.tmp.gtf ; \
        mkdir -p GTFs_${VERSION} ; \
        mv *${VERSION}* ./GTFs_${VERSION} ; \
    fi
done ; done
