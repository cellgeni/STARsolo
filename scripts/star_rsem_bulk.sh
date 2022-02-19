#!/bin/bash

TAG=$1
SREF=/nfs/cellgeni/STAR/human/2020A/index
RREF=/nfs/cellgeni/STAR/human/2020A/2020A_rsem
FQDIR=/lustre/scratch117/cellgen/cellgeni/TIC-bulkseq/tic-1333/fastqs
PAIRED="--paired-end"
STRAND="--forward-prob 0"
CPUS=16

## for single-end reads, change both STAR and RSEM options
## for strand-specific processing, change --forward-prob in rsem-calculate-expression (1 is FR, 0 is RF, 0.5 is non-strand-specific). 

mkdir $TAG && cd $TAG

R1=""
R2=""
if [[ `find $FQDIR/* | grep $TAG | grep "_1\.fastq"` != "" ]]
then
  R1=`find $FQDIR/* | grep $TAG | grep "_1\.fastq" | sort | tr '\n' ',' | sed "s/,$//g"`
  R2=`find $FQDIR/* | grep $TAG | grep "_2\.fastq" | sort | tr '\n' ',' | sed "s/,$//g"`
elif [[ `find $FQDIR/* | grep $TAG | grep "R1\.fastq"` != "" ]]
then
  R1=`find $FQDIR/* | grep $TAG | grep "R1\.fastq" | sort | tr '\n' ',' | sed "s/,$//g"`
  R2=`find $FQDIR/* | grep $TAG | grep "R2\.fastq" | sort | tr '\n' ',' | sed "s/,$//g"`
elif [[ `find $FQDIR/* | grep $TAG | grep "_R1_.*\.fastq"` != "" ]]
then
  R1=`find $FQDIR/* | grep $TAG | grep "_R1_" | sort | tr '\n' ',' | sed "s/,$//g"`
  R2=`find $FQDIR/* | grep $TAG | grep "_R2_" | sort | tr '\n' ',' | sed "s/,$//g"`
else
  >&2 echo "ERROR: No appropriate fastq files were found! Please check file formatting, and check if you have set the right FQDIR."
  exit 1
fi

GZIP=""
if [[ `find $FQDIR/* | grep $TAG | grep "\.gz$"` != "" ]]
then
  GZIP="--readFilesCommand zcat"
fi

## ENCODE options for RNA-seq alignment
STAR --runThreadN $CPUS --genomeDir $SREF --readFilesIn $R1 $R2 $GZIP --outFilterMultimapNmax 20 \
     --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
     --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt \
     --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
     --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --sjdbScore 1 --limitBAMsortRAM 30000000000

samtools index -@$CPUS Aligned.sortedByCoord.out.bam

rsem-calculate-expression $PAIRED -p $CPUS --bam --estimate-rspd --seed 12345 --no-bam-output $STRAND Aligned.toTranscriptome.out.bam $RREF $TAG.rsem
