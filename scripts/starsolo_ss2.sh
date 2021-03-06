#!/bin/bash -e 

## STARsolo processing of Smart-seq2 
## you need to have a manifest file named $TAG.manifest.tsv
## with 3 tab-separated columns: 1) full path to read1; 2) full path to read2; 3) cell name

TAG=$1

if [[ ! -s $TAG.manifest.tsv ]]
then 
  >&2 echo "ERROR: File $TAG.manifest.tsv does not exist or is empty! Please provide appropriate manifest file."
  >&2 echo "Usage: ./starsolo_ss2.sh <manifest_tsv>"
  exit 1
fi

CPUS=16      ## typically bsub this into normal queue with 16 cores and 64 Gb RAM.   
REF=/nfs/cellgeni/STAR/human/2020A/index  ## choose the proper reference 
GZIP="--readFilesCommand zcat"
#BAM="--outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 2 --limitBAMsortRAM 60000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM GX GN"
BAM="--outSAMtype None"

mkdir $TAG.solo.SS2 && cd $TAG.solo.SS2

## outFilter* options can be adjusted according to the mapping rate and mapped length
## alignment often works well without clipping; if not, --clip3pAdapterSeq <3' adapter sequence> option can be added, or separate trimming using bbduk.sh works well too
STAR --runThreadN $CPUS --genomeDir $REF --runDirPerm All_RWX $GZIP $BAM \
     --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
     --soloType SmartSeq --readFilesManifest ../$TAG.manifest.tsv --soloUMIdedup Exact --soloStrand Unstranded \
     --soloFeatures Gene GeneFull --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx
