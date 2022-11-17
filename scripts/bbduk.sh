#!/bin/bash 

TAG=$1
ADAPTERS=/nfs/users/nfs_a/ap41/bbmap/resources/adapters.fa #need to change to a cellgeni path

bbduk.sh in1=${TAG}_1.fastq.gz in2=${TAG}_2.fastq.gz out1=$TAG.bbduk.R1.fastq out2=$TAG.bbduk.R2.fastq ref=$ADAPTERS trimpolya=10 ktrim=r k=23 mink=11 hdist=1 tpe tbo &> $TAG.bbduk.log
