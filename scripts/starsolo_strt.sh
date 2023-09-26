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
  >&2 echo "Usage: ./starsolo_strt.sh <fastq_dir> <sample_id>"
  >&2 echo "(make sure you set the correct REF, WL, and BAM variables below)"
  exit 1
fi

FQDIR=`readlink -f $FQDIR`
CPUS=16                                                                ## typically bsub this into normal queue with 16 cores and 64 Gb RAM.   
REF=/nfs/cellgeni/STAR/human/2020A/index                               ## choose the appropriate reference 
WL=/nfs/cellgeni/STAR/whitelists                                       ## directory with all barcode whitelists
CBLEN=8
UMILEN=8

## choose one of the two otions, depending on whether you need a BAM file 
#BAM="--outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 500 --limitBAMsortRAM 60000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN"
BAM="--outSAMtype None"

###################################################################### DONT CHANGE OPTIONS BELOW THIS LINE ##############################################################################################

BC=$WL/96_barcodes.list

mkdir $TAG && cd $TAG
## for multiple fastq files; change grep options according to your fastq file format 
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
## see if the original fastq files are archived: 
if [[ `find $FQDIR/* | grep $TAG | grep "\.gz$"` != "" ]]
then
  GZIP="--readFilesCommand zcat"
fi

## note the switched R2 and R1 compared to 10x! R1 is biological read in STRT-seq
$CMD STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn $R2 $R1 --runDirPerm All_RWX $GZIP $BAM \
     --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloBarcodeReadLength 0 --soloCBlen $CBLEN --soloUMIstart $((CBLEN+1)) --soloUMIlen $UMILEN --soloStrand $STRAND \
     --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
     --limitOutSJcollapsed 10000000 --soloCellFilter None \
     --soloFeatures Gene GeneFull --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx --outReadsUnmapped Fastx

## max-CR bzip all unmapped reads with multicore pbzip2 
pbzip2 -9 -m2000 -p$CPUS Unmapped.out.mate1
pbzip2 -9 -m2000 -p$CPUS Unmapped.out.mate2

## gzip all outputs
cd output
for i in Gene/raw GeneFull/raw
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
 
