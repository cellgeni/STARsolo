#!/bin/bash -e 

## v3.0 of STARsolo wrappers is set up to guess the chemistry automatically
## newest version of the script uses STAR v2.7.9a with EM multimapper processing in STARsolo (but it's not on by default; turn this on with --soloMultiMappers EM)
## velocyto extra processing also became unnecessary 

TAG=$1
if [[ $TAG == "" ]]
then
  >&2 echo "Usage: ./starsolo_auto.sh <sample_tag>"
  >&2 echo "(make sure you set the correct REF, FQDIR, and SORTEDBAM/NOBAM variables)"
  exit 1
fi

CPUS=16      ## typically bsub this into normal queue with 16 cores and 64 Gb RAM.   
REF=/nfs/cellgeni/STAR/human/2020A/index  ## choose the appropriate reference 
FQDIR=/lustre/scratch117/cellgen/cellgeni/TIC-starsolo/tic-1258/fastq_HCA_Dropseq  ### change to the directory with fastq files/folders
## choose one of the two otions, depending on whether you need a BAM file 
BAM="--outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 2 --limitBAMsortRAM 120000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB GX GN"
#BAM="--outSAMtype None"

###################### DONT CHANGE OPTIONS BELOW THIS LINE ###########################

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

## let's see if the files are archived or not. Gzip is the only one we test for, but bgzip archives should work too since they are gzip-compatible.
GZIP=""
if [[ `find $FQDIR/* | grep $TAG | grep "\.gz$"` != "" ]]
then  
  GZIP="--readFilesCommand zcat"
fi

STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn $R2 $R1 --runDirPerm All_RWX $GZIP $BAM \
     --soloType CB_UMI_Simple --soloCBwhitelist None --soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 --soloBarcodeReadLength 0 \
     --soloFeatures Gene GeneFull --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx

## finally, let's gzip all outputs
cd output
for i in Gene/raw Gene/filtered GeneFull/raw GeneFull/filtered
do 
  cd $i; for j in *; do gzip $j & done
  cd ../../
done

wait
echo "ALL DONE!"
