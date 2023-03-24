#!/bin/bash -e 

## v3.1 of STARsolo wrappers is set up to guess the chemistry automatically
## newest version of the script uses STAR v2.7.10a with EM multimapper processing 
## in STARsolo which on by default; the extra matrix can be found in /raw subdir 

SIF="/nfs/cellgeni/singularity/images/starsolo_2-7-10a-alpha-220818_samtools_1-15-1_seqtk-1-13_bbmap_38-97_RSEM-1-3-3.sif"
CMD="singularity run --nv --bind /nfs,/lustre $SIF"

FQDIR=$1
TAG=$2

if [[ $FQDIR == "" || $TAG == "" ]]
then
  >&2 echo "Usage: ./starsolo_indrops.sh <fastq_dir> <sample_id>"
  >&2 echo "(make sure you set the correct REF, WL, ADAPTER, BC1/BC2, and BAM variables below)"
  exit 1
fi

FQDIR=`readlink -f $FQDIR`
CPUS=16                                                                ## typically bsub this into normal queue with 16 cores and 64 Gb RAM.   
REF=/nfs/cellgeni/STAR/human/2020A/index                               ## choose the appropriate reference 
WL=/nfs/cellgeni/STAR/whitelists                                       ## directory with all barcode whitelists

ADAPTER=GAGTGATTGCTTGTGACGCCTT                                         ## these could be GAGTGATTGCTTGTGACGCCTT or GAGTGATTGCTTGTGACGCCAA 
BC1=$WL/inDrops_Ambrose2_bc1.txt
BC2=$WL/inDrops_Ambrose2_bc2.txt

## choose one of the two otions, depending on whether you need a BAM file 
#BAM="--outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 500 --limitBAMsortRAM 60000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN"
BAM="--outSAMtype None"

###################################################################### DONT CHANGE OPTIONS BELOW THIS LINE ##############################################################################################

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

## increased soloAdapterMismatchesNmax to 3, as per discussions in STAR issues
$CMD STAR  --runThreadN $CPUS --genomeDir $REF --readFilesIn $R2 $R1 --runDirPerm All_RWX $GZIP $BAM \
     --soloType CB_UMI_Complex --soloCBwhitelist $BC1 $BC2 --soloAdapterSequence $ADAPTER  \
     --soloAdapterMismatchesNmax 3 --soloCBmatchWLtype 1MM --soloCBposition 0_0_2_-1 3_1_3_8 --soloUMIposition 3_9_3_14 \
     --soloFeatures Gene GeneFull --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx


## finally, let's gzip all outputs
cd output
for i in Gene/raw Gene/filtered GeneFull/raw GeneFull/filtered
do 
  cd $i; for j in *; do gzip $j & done
  cd ../../
done

## index the BAM file
if [[ -s Aligned.sortedByCoord.out.bam ]]
then
  $CMD samtools index -@16 Aligned.sortedByCoord.out.bam &
fi

wait
echo "ALL DONE!"
