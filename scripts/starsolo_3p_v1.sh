#!/bin/bash -e 

## newest version of the script uses STAR v2.7.9a with EM multimapper processing in STARsolo (but it's not on by default; turn this on with --soloMultiMappers EM)
## velocyto extra processing also became unnecessary 

TAG=$1
if [[ $TAG == "" ]]
then
  >&2 echo "Usage: ./starsolo_3p_v1.sh <sample_tag>"
  >&2 echo "(make sure you set the correct REF, FQDIR, and SORTEDBAM/NOBAM variables)"
  exit 1
fi

CPUS=16      ## typically bsub this into normal queue with 16 cores and 64 Gb RAM.   
REF=/nfs/cellgeni/STAR/human/2020A/index  ## choose the appropriate reference 
FQDIR=/lustre/scratch117/cellgen/cellgeni/alexp/fastqs  ### change to the directory with fastq files/folders

BC=/nfs/cellgeni/STAR/whitelists/737K-april-2014_rc.txt ## 10x v1 - 14 bp barcodes
#BC=/nfs/cellgeni/STAR/whitelists/3M-february-2018.txt  ## 10x v3
#BC=/nfs/cellgeni/STAR/whitelists/737K-august-2016.txt ## 10x v2 or all 5' 
UMILEN=10   ## 10x v2
#UMILEN=12    ## 10x v3
STR=Forward  ## 3' 10x
#STR=Reverse  ## 5' 10x

## choose one of the two otions, depending on whether you need a BAM file 
SORTEDBAM="--outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 2 --limitBAMsortRAM 60000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB GX GN"
NOBAM="--outSAMtype None"

mkdir $TAG && cd $TAG
## for multiple fastq files; change grep options according to your fastq file format 
R1=`find $FQDIR/* | grep $TAG | grep "R1.fastq.gz" | sort | tr '\n' ',' | sed "s/,$//g"`
R2=`find $FQDIR/* | grep $TAG | grep "R2.fastq.gz" | sort | tr '\n' ',' | sed "s/,$//g"`

if [[ $R1 == "" || $R2 == "" ]]
then
  >&2 echo "No appropriate R1 or R2 read files was found for sample tag $TAG! Make sure you have set the correct FQDIR."
  >&2 echo "Usage: ./starsolo_3p_v1.sh <sample_tag>"
  exit 1
fi

## for v1 10x, we added --soloCBlen 14 --soloUMIstart 15, since barcodes are 14 bp long (default is 16)
STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn $R2 $R1 --runDirPerm All_RWX --readFilesCommand zcat $NOBAM \
     --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloBarcodeReadLength 0 --soloCBlen 14 --soloUMIstart 15 --soloUMIlen $UMILEN --soloStrand $STR \
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
 
