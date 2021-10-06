#!/bin/bash

TAG=$1
SREF=/nfs/cellgeni/alexp/GRCh38_2020A/STAR
RREF=/nfs/cellgeni/alexp/GRCh38_2020A/GRCh38_v32_rsem
FQDIR=/lustre/scratch117/cellgen/cellgeni/alexp/Bulk/Bulk_for_tic957/fastqs
PAIRED="--paired-end"
CPUS=16

## for single-end reads, change both STAR and RSEM options
## for strand-specific processing, change --forward-prob in rsem-calculate-expression (1 is FR, 0 is RF, 0.5 is non-strand-specific). 

mkdir $TAG && cd $TAG
R1=`find $FQDIR/* | grep $TAG | grep R1`
R2=`find $FQDIR/* | grep $TAG | grep R2`

STAR --runThreadN $CPUS --genomeDir $SREF --readFilesIn $R1 $R2 --readFilesCommand zcat --outFilterMultimapNmax 20 \
     --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
     --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt \
     --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
     --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --sjdbScore 1 --limitBAMsortRAM 30000000000

rsem-calculate-expression $PAIRED -p $CPUS --bam --estimate-rspd --seed 12345 -p 4 --no-bam-output --forward-prob 0.5 Aligned.toTranscriptome.out.bam $RREF $TAG.rsem
