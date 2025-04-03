#!/bin/bash
##Note: it is importnat to have an enviroment of conda active that contains:
# - minimap2
# - samtools
# - seqkit
# - Gffcompare

###Setting the local path to the stringtie executable:
path_stringtie=/data2/cristian/ntorres2/LongRead_RNAseq/Stringtie_tests_BP10new/stringtie/stringtie
gffread_path=/data2/cristian/ntorres2/LongRead_RNAseq/Workspaces/BAR_longreads_C05C06/STRINGTIE_annotations_v2/gffread/gffread
gffcompare_path=/data2/cristian/ntorres2/LongRead_RNAseq/Workspaces/BAR_longreads_C05C06/STRINGTIE_annotations_v2/gffcompare/gffcompare
WORKDIR=/data2/cristian/ntorres2/LongRead_RNAseq/Workspaces/BAR_longreads_C05C06
VERSION=$3
operon_finder=./operon_fidner_v3.py
#$1 = samples
#$2 = max intron size

#Minimap2 alingment of the reads
echo "Performing mapping of the reads..."
for i in $1 ; do for x in $2 ; do \
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

#Merging all annotations with gffcompare
for x in $2 ; do for cov in 1.5 ; do
    output_merged_gtf=OdiBAR_LRannot_v${x}_c${cov}_stringtie.no_clean.${VERSION}.combined.gtf
    if [ -f "$output_merged_gtf" ]; then
        echo "$output_merged_gtf already exists. Skipping..."
        continue
    fi
    ls *_no_assembled_v${x}_c${cov}_f0.2.gtf > GTFs_${VERSION}.v${x}_c${cov}.raw.list ; \
    $gffcompare_path -i GTFs_${VERSION}.v${x}_c${cov}.raw.list -o OdiBAR_LRannot_v${x}_c${cov}_stringtie.no_clean.${VERSION} ; \
done ; done

#Python script
echo "Performing operon search..."
for x in $2 ; do \
    for cov in 1.5 ; do \
        operon_file=OdiBAR_LRannot_v${x}_c${cov}_stringtie.no_clean.${VERSION}.combined_operons_found_v3.tsv
        if [ -f "$operon_file" ]; then
            echo "Operon search already done for v${x}_cov${cov}. Skipping..."
            continue
        fi
        python ${operon_finder} OdiBAR_LRannot_v${x}_c${cov}_stringtie.no_clean.${VERSION}.combined.gtf ; \
    done ; \
done

#Creating annoations with potential operons
echo "Creating annoations with potential operons"
for x in $2 ; do for cov in 1.5 ; do \
    output_operon_gtf=OdiBAR_LRannot_v${x}_c${cov}_stringtie.no_clean.${VERSION}.operons.gtf
    output_opgenes_gtf=OdiBAR_LRannot_v${x}_c${cov}_stringtie.no_clean.${VERSION}.operongenes.gtf
    output_operons_clean=OdiBAR_LRannot_v${x}_c${cov}_allGoodOperons.${VERSION}.gtf
     if [ -f "$output_operons_clean" ]; then
        echo "Skipping..."
        continue
    fi
    awk '{print $1}' OdiBAR_LRannot_v${x}_c${cov}_stringtie.no_clean.${VERSION}.combined_operons_found_v3.tsv | sort -u > List.operons.tmp; \
    awk '{print $2}' OdiBAR_LRannot_v${x}_c${cov}_stringtie.no_clean.${VERSION}.combined_operons_found_v3.tsv | sort -u > List.genes.tmp; \
    $gffread_path --ids List.operons.tmp OdiBAR_LRannot_v${x}_c${cov}_stringtie.no_clean.${VERSION}.combined.gtf -F -T -o $output_operon_gtf ; \
    $gffread_path --ids List.genes.tmp OdiBAR_LRannot_v${x}_c${cov}_stringtie.no_clean.${VERSION}.combined.gtf -F -T -o $output_opgenes_gtf ; \
    wc -l List.*.tmp > Counts_${VERSION}_no_assembled_v${x}.txt ; rm *.tmp ; \
    $path_stringtie --merge $output_opgenes_gtf -l OpG -f 0.01 -F 0 -T 0 -c 0 -g '-100' -o OdiBAR_LRannot_v${x}_c${cov}_allGoodOperonGenes.${VERSION}.gtf ; \
    $path_stringtie --merge $output_operon_gtf -l OPRN -c 0 -F 0 -T 0 -f 0.01 -g 0 -o $output_operons_clean ; \
done ; done

#Creating clean version of final GTF
echo 'Creating clean version of final GTF'
for x in $2 ; do for cov in 1.5 ; do \
    output_clean_gtf=OdiBAR_LRannot_v${x}_c${cov}_stringtie_${VERSION}.cleanDEF.gtf
    if [ -f "$output_clean_gtf" ]; then
        echo "Clean GTF v${x}_cov${cov} already exists. Skipping..."
        continue
    fi
    $gffcompare_path -r OdiBAR_LRannot_v${x}_c${cov}_allGoodOperons.${VERSION}.gtf OdiBAR_LRannot_v${x}_c${cov}_stringtie.no_clean.${VERSION}.combined.gtf -o Gffcompare_GoodOperonsVSall_v${x}_c${cov}_${VERSION} ; \
    awk '{if($3!="u" && $3!="x") print $5}' Gffcompare_GoodOperonsVSall_v${x}_c${cov}_${VERSION}*.tmap > GoodOperons_all_v${x}_c${cov}_${VERSION}.loci_exclusion.txt ; \
    $gffread_path --nids GoodOperons_all_v${x}_c${cov}_${VERSION}.loci_exclusion.txt OdiBAR_LRannot_v${x}_c${cov}_stringtie.no_clean.${VERSION}.combined.gtf -F -T -o $output_clean_gtf ; \
done ; done

#Creating final GTFs
echo 'Creating final GTFs'
for x in $2 ; do for cov in 1.5 ; do
    final_gtf_file=OdiBAR_LRannot_v${x}_c${cov}_stringtie.clean-and-OPRNs.${VERSION}.gtf
    if [ -f "$final_gtf_file" ]; then
        echo "GTF already done. Skipping..."
        continue
    fi
    cat OdiBAR_LRannot_v${x}_c${cov}_allGoodOperons.${VERSION}.gtf OdiBAR_LRannot_v${x}_c${cov}_allGoodOperonGenes.${VERSION}.gtf > tmp.gtf ; \
    $gffread_path --sort-alpha -F -T -o Merged_OPRNs-and-OpGs_v${x}_c${cov}.${VERSION}.sorted.gtf tmp.gtf ; rm tmp.gtf ; \
    $path_stringtie --merge OdiBAR_LRannot_v${x}_c${cov}_stringtie_${VERSION}.cleanDEF.gtf -G OdiBAR_LRannot_v${x}_c${cov}_allGoodOperonGenes.${VERSION}.gtf -f 0.2 -F 0 -T 0 -g '-50' -o OdiBAR_LRannot_v${x}_c${cov}_stringtie.clean-noOPRNs.${VERSION}.gtf ; \
    $path_stringtie --merge OdiBAR_LRannot_v${x}_c${cov}_stringtie_${VERSION}.cleanDEF.gtf -G Merged_OPRNs-and-OpGs_v${x}_c${cov}.${VERSION}.sorted.gtf -f 0.2 -F 0 -T 0 -g '-50' -o OdiBAR_LRannot_v${x}_c${cov}_stringtie.clean-and-OPRNs.${VERSION}.gtf ; \
done ; done

#Checking for trasncripts to add:
echo 'Checking for trasncripts to add'
for x in $2 ; do for cov in 1.5 ; do \
    final_gffcompare=OdiBAR_LRannot_v${x}_c${cov}.${VERSION}.Gffcmp_final
    if [ -f "$final_gffcompare.loci" ]; then
        echo "Final gffcompare already done. Skipping..."
        continue
    fi
    $gffcompare_path -r OdiBAR_LRannot_v${x}_c${cov}_stringtie.clean-and-OPRNs.${VERSION}.gtf OdiBAR_LRannot_v${x}_c${cov}_stringtie.no_clean.${VERSION}.combined.gtf -o $final_gffcompare ; \
    (for type in "=" c k m n j e o s x i y p r u ; do \
        count=$(awk -v a="$type" '{if($4==a) print $5,$4}' $final_gffcompare.tracking | wc -l) ; \
        echo "${type} $count"; \
    done ) > $final_gffcompare.transcript_types.txt ; \
done ; done

for x in $2 ; do for cov in 1.5 ; do \
    final_gffcompare=OdiBAR_LRannot_v${x}_c${cov}.${VERSION}.Gffcmp_final
    for type2 in k o x i y u ; do \
        type_list_file=${final_gffcompare}.${type2}.transcript_types.list ;\
        echo $type2 ;\
        if [ -f "$type_list_file" ]; then
            echo 'Already done. Skipping..'
            continue
        fi
        awk -v a="$type2" '{if($3==a) print $5}' ${final_gffcompare}*.tmap > $type_list_file ; \
        echo "Done" ; \
    done ; \
done ; done

for x in $2 ; do for cov in 1.5 ; do \
    first_output_gtf_name=OdiBAR_LRannot_v${x}_c${cov}
    final_gffcompare=$first_output_gtf_name.${VERSION}.Gffcmp_final
    for type3 in o x i y u ; do \
        echo "Extracting $type3 transcripts..."
        type_gtf_file=$first_output_gtf_name.${VERSION}.type_${type3}.gtf
        if [ -f "$type_gtf_file" ]; then
            echo 'Already done. Skipping..'
            continue
        fi
        $gffread_path --ids ${final_gffcompare}.${type3}.transcript_types.list ${first_output_gtf_name}_stringtie.no_clean.${VERSION}.combined.gtf -F -T -o $type_gtf_file ; \
        echo "Done!" ;\
    done ; \
done ; done

#Creating definitive annotation file
echo 'Creating definitive annotation file'
for x in $2 ; do for cov in 1.5 ; do \
    first_output_gtf_name=OdiBAR_LRannot_v${x}_c${cov}
    def_final_gtf_file=${first_output_gtf_name}_stringtie_${VERSION}.clean-and-OPRNs.DEF-FINAL.gtf
    if [ -f "$def_final_gtf_file" ]; then
        echo "Definitive FINAL GTF already done. Skipping..."
        continue
    fi
    $path_stringtie --merge \
        ${first_output_gtf_name}_stringtie_${VERSION}.cleanDEF.gtf \
        $first_output_gtf_name.${VERSION}.type_u.gtf \
        $first_output_gtf_name.${VERSION}.type_y.gtf \
        $first_output_gtf_name.${VERSION}.type_i.gtf \
        $first_output_gtf_name.${VERSION}.type_x.gtf \
        $first_output_gtf_name.${VERSION}.type_o.gtf \
        -G OdiBAR_LRannot_v${x}_c${cov}_allGoodOperonGenes.${VERSION}.gtf \
        -f 0.2 -F 0 -T 0 -g '-50' -o $first_output_gtf_name.${VERSION}.clean-noOPRNs.DEF-FINAL.gtf ; \
    $path_stringtie --merge \
        ${first_output_gtf_name}_stringtie_${VERSION}.cleanDEF.gtf \
        $first_output_gtf_name.${VERSION}.type_u.gtf \
        $first_output_gtf_name.${VERSION}.type_y.gtf \
        $first_output_gtf_name.${VERSION}.type_i.gtf \
        $first_output_gtf_name.${VERSION}.type_x.gtf \
        $first_output_gtf_name.${VERSION}.type_o.gtf \
        -G Merged_OPRNs-and-OpGs_v${x}_c${cov}.${VERSION}.sorted.gtf \
        -f 0.2 -F 0 -T 0 -g '-50' -o $def_final_gtf_file ; \
done ; done 

#Final check of operons in the annotation.
echo "Performing FINAL operon search..."
for x in $2 ; do for cov in 1.5 ; do \
    first_output_gtf_name=OdiBAR_LRannot_v${x}_c${cov}_stringtie_${VERSION}
    def_final_gtf_file=$first_output_gtf_name.clean-and-OPRNs.DEF-FINAL.gtf
    final_operon_file=$first_output_gtf_name.clean-and-OPRNs.DEF-FINAL_operons_found_v3.tsv
    if [ -f "$final_operon_file" ]; then
        echo "Operon search already done. Skipping..."
        continue
    fi
    python ${operon_finder} $def_final_gtf_file ; \
done ; done

for x in $2 ; do for cov in 1.5 ; do \
    echo 'Creating well named annotation DEF-FINAL' ; \
    first_output_gtf_name=OdiBAR_LRannot_v${x}_c${cov}_stringtie_${VERSION} ; \
    def_final_gtf_file=$first_output_gtf_name.clean-and-OPRNs.DEF-FINAL.gtf ; \
    final_operon_file=$first_output_gtf_name.clean-and-OPRNs.DEF-FINAL_operons_found_v3.tsv ; \
    awk '{print $1}' $final_operon_file | sort -u > List.operons.tmp; \
    awk '{print $2}' $final_operon_file | sort -u > List.genes.tmp; \
    $gffread_path --ids List.operons.tmp $def_final_gtf_file -F -T -o OPRNs.${VERSION}.gtf ; \
    $gffread_path --ids List.genes.tmp $def_final_gtf_file -F -T -o OpGs.${VERSION}.gtf ; \
    wc -l List.*.tmp > Counts_DEF-FINAL_${VERSION}_no_assembled_v${x}.txt ; rm *.tmp ; \
    $path_stringtie --merge OpGs.${VERSION}.gtf -l OpG -f 0.01 -F 0 -T 0 -c 0 -g '-100' -o OdiBAR_LRannot_v${x}_c${cov}_allGoodOperonGenes.${VERSION}.DEF-FINAL.gtf ; \
    $path_stringtie --merge OPRNs.${VERSION}.gtf -l OPRN -c 0 -F 0 -T 0 -f 0.01 -g 0 -o OdiBAR_LRannot_v${x}_c${cov}_allGoodOperons.${VERSION}.DEF-FINAL.gtf ; \
done ; done

echo 'Creating definitive annotation file'
for x in $2 ; do for cov in 1.5 ; do \
    first_output_gtf_name=OdiBAR_LRannot_v${x}_c${cov}
    def_final_gtf_file=${first_output_gtf_name}_stringtie_${VERSION}.clean-and-OPRNs.DEF-FINALv2.gtf
    if [ -f "$def_final_gtf_file" ]; then
        echo "Definitive FINAL GTF already done. Skipping..."
        continue
    fi
    $path_stringtie --merge \
        ${first_output_gtf_name}_stringtie_${VERSION}.cleanDEF.gtf \
        $first_output_gtf_name.${VERSION}.type_u.gtf \
        $first_output_gtf_name.${VERSION}.type_y.gtf \
        $first_output_gtf_name.${VERSION}.type_i.gtf \
        $first_output_gtf_name.${VERSION}.type_x.gtf \
        $first_output_gtf_name.${VERSION}.type_o.gtf \
        -G OdiBAR_LRannot_v${x}_c${cov}_allGoodOperonGenes.${VERSION}.DEF-FINAL.gtf \
        -f 0.2 -F 0 -T 0 -g '-50' -o $first_output_gtf_name.${VERSION}.clean-noOPRNs.DEF-FINALv2.gtf ; \
    cat OdiBAR_LRannot_v${x}_c${cov}_allGoodOperonGenes.${VERSION}.DEF-FINAL.gtf  OdiBAR_LRannot_v${x}_c${cov}_allGoodOperons.${VERSION}.DEF-FINAL.gtf > tmp.gtf ; \
    $gffread_path --sort-alpha -F -T -o Merged_OPRNs-and-OpGs_v${x}_c${cov}.${VERSION}.sorted.DEF-FINAL.gtf tmp.gtf ; rm tmp.gtf ; \
    $path_stringtie --merge \
        ${first_output_gtf_name}_stringtie_${VERSION}.cleanDEF.gtf \
        $first_output_gtf_name.${VERSION}.type_u.gtf \
        $first_output_gtf_name.${VERSION}.type_y.gtf \
        $first_output_gtf_name.${VERSION}.type_i.gtf \
        $first_output_gtf_name.${VERSION}.type_x.gtf \
        $first_output_gtf_name.${VERSION}.type_o.gtf \
        -G Merged_OPRNs-and-OpGs_v${x}_c${cov}.${VERSION}.sorted.DEF-FINAL.gtf \
        -f 0.2 -F 0 -T 0 -g '-50' -o $def_final_gtf_file ; \
done ; done 


#Moving all files to a organized folder:
for x in $2 ; do for cov in 1.5 ; do \
    first_output_gtf_name=OdiBAR_LRannot_v${x}_c${cov}_stringtie_${VERSION}
    def_final_gtf_file=$first_output_gtf_name.clean-and-OPRNs.DEF-FINALv2.gtf
    if [ -f "$def_final_gtf_file" ]; then
#        rm List*.tmp ; rm *.tmp.gtf ; \
        mkdir -p GTFs_${VERSION} ; \
        mv *${VERSION}* ./GTFs_${VERSION} ; \
    fi
done ; done 
