#!/bin/bash

## newest version of the script uses STAR v2.7.9a with EM multimapper processing in STARsolo 
## velocyto extra processing also became unnecessary 

## based on https://github.com/Teichlab/mapcloud/blob/master/scripts/starsolo/starsolo-wrapper.sh
## this one is V3 3'. UMI=12 for v3, 10 for v2

###### Change these options before running!
###### UMILEN = 10 for v2/v1, 12 for v3 (if R1 = 26 bp that's v2, if 28 bp, that's v3); 
###### STR = Forward for 3', Reverse for 5'; also, ALWAYS use 737K-august-2016.txt for all 5' - at least as of Summer 2021!


TAG=$1
CPUS=16      ## typically bsub this into normal queue with 16 cores and 64 Gb RAM.   
REF=/nfs/cellgeni/alexp/STARsolo-GRCh38-3.0.0/STAR
FQDIR=/lustre/scratch117/cellgen/cellgeni/alexp/10x/STARsolo/Luz_5prime/2new_human/fastqs
#BC=/nfs/users/nfs_a/ap41/whitelists/3M-february-2018.txt  ## 10x v3
BC=/nfs/users/nfs_a/ap41/whitelists/737K-august-2016.txt ## 10x v2 or all 5' 
#UMILEN=12    ## 10x v3
UMILEN=10   ## 10x v2
STR=Reverse  ## 5' 10x
#STR=Forward ## 3' 10x
SORTEDBAM="--outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 2 --limitBAMsortRAM 60000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB GX GN"
NOBAM="--outSAMtype None"

mkdir $TAG && cd $TAG
## for multiple fastq files:
R1=`find $FQDIR/* | grep $TAG | grep R1 | sort | tr '\n' ',' | sed "s/,$//g"`
R2=`find $FQDIR/* | grep $TAG | grep R2 | sort | tr '\n' ',' | sed "s/,$//g"`

STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn $R2 $R1 --runDirPerm All_RWX --readFilesCommand zcat $NOBAM \
     --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloBarcodeReadLength 0 --soloUMIlen $UMILEN --soloStrand $STR \
     --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
     --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
     --soloFeatures Gene GeneFull Velocyto --soloOutFileNames output/ genes.tsv barcodes.tsv matrix.mtx

## gzip all outputs
cd output
for i in Gene/raw Gene/filtered GeneFull/raw GeneFull/filtered Velocyto/raw Velocyto/filtered
do 
  cd $i; for j in *; do gzip $j & done
  cd ../../
done

wait
echo "ALL DONE!"
 
