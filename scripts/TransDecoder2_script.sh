#!/bin/bash
GENOME=$1

GTF=$2
NAME=$(basename $2 .gtf)

DIR=$3
DIRSCRIPTS=$5

perl -I$DIRSCRIPTS/PerlLib $DIRSCRIPTS/TD2_util/gtf_genome_to_cdna_fasta.pl $GTF $GENOME > ${DIR}$NAME.fasta 

perl -I$DIRSCRIPTS/PerlLib $DIRSCRIPTS/TD2_util/gtf_to_alignment_gff3.pl $GTF > ${DIR}$NAME.gff3

mkdir -p ${DIR}$NAME

TD2.LongOrfs -t ${DIR}$NAME.fasta -O ${DIR}$NAME --complete-orfs-only
#hmmsearch --cpu 8 -E 1e-10 --domtblout ${NAME}.pfam.domtblout ${DIR}/Pfam-A.hmm $NAME/longest_orfs.pep

for i in $4 ; do 
     db_name=$(basename $i)
     echo $db_name
     mmseqs easy-search ${DIR}$NAME/longest_orfs.pep ${i} ${DIR}$NAME/alnRes_${db_name}.m8 tmp -s 7.0
     TD2.Predict -t ${DIR}$NAME.fasta -O ${DIR}$NAME --retain-mmseqs-hits ${DIR}$NAME/alnRes_${db_name}.m8
done
#cat ${DIR}$NAME/alnRes_*.m8 > ${DIR}$NAME/combined_alnRes.m8

perl -I$DIRSCRIPTS/PerlLib $DIRSCRIPTS/TD2_util/cdna_alignment_orf_to_genome_orf.pl \
     $NAME.fasta.TD2.gff3 \
     ${DIR}$NAME.gff3 \
     ${DIR}$NAME.fasta > ${DIR}$NAME.fasta.TD2.genome.gff3

mkdir -p ${DIR}${NAME}_files

mv ${DIR}${NAME}/ ${DIR}${NAME}_files/${NAME}
mv ${NAME}.fasta.TD2.gff3 ${DIR}${NAME}_files
mv ${DIR}${NAME}.gff3 ${DIR}${NAME}_files
mv ${DIR}${NAME}.fasta ${DIR}${NAME}_files
mv ${NAME}.fasta.TD2.bed ${DIR}${NAME}_files
mv ${NAME}.fasta.TD2.cds ${DIR}${NAME}_files
mv ${NAME}.fasta.TD2.pep ${DIR}${NAME}_files
