#!/bin/bash -e 

## newest version of the script uses STAR v2.7.9a with EM multimapper processing in STARsolo (but it's not on by default; turn this on with --soloMultiMappers EM)
## velocyto extra processing also became unnecessary 

## this is for STRT - very popular in China!

TAG=$1
if [[ $TAG == "" ]]
then
  >&2 echo "Usage: ./starsolo_strt.sh <sample_tag>"
  >&2 echo "(make sure you set the correct REF, FQDIR, and SORTEDBAM/NOBAM variables)"
  exit 1
fi

CPUS=16                                                                ## typically bsub this into normal queue with 16 cores and 64 Gb RAM.   
REF=/nfs/cellgeni/STAR/human/2020A/index                               ## choose the appropriate reference 
WL=/nfs/cellgeni/STAR/whitelists                                       ## directory with all barcode whitelists
FQDIR=/lustre/scratch117/cellgen/cellgeni/TIC-starsolo/tic-XXX/fastqs  ## directory with your fastq files - can be in subdirs, just make sure tag is unique and greppable (e.g. no Sample1 and Sample 10). 

## choose one of the two otions, depending on whether you need a BAM file 
#BAM="--outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 2 --limitBAMsortRAM 120000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB GX GN"
BAM="--outSAMtype None"

###################################################################### DONT CHANGE OPTIONS BELOW THIS LINE ##############################################################################################

if [[ `which samtools` == "" || `which STAR` == "" ]]
then
  echo "ERROR: Please make sure you have STAR (v2.7.9a or above) and samtools are installed and available in PATH!"
  exit 1
fi

BC=$WL/96_barcodes.list

mkdir $TAG && cd $TAG
## for multiple fastq files; change grep options according to your fastq file format 
R1=`find $FQDIR/* | grep $TAG | grep "_f1.fastq.gz" | sort | tr '\n' ',' | sed "s/,$//g"`
R2=`find $FQDIR/* | grep $TAG | grep "_r2.fastq.gz" | sort | tr '\n' ',' | sed "s/,$//g"`

if [[ $R1 == "" || $R2 == "" ]]
then
  >&2 echo "No appropriate R1 or R2 read files was found for sample tag $TAG! Make sure you have set the correct FQDIR."
  >&2 echo "Usage: ./starsolo_strt.sh <sample_tag>"
  exit 1
fi

## note the switched R2 and R1 compared to 10x! R1 is biological read in STRT-seq
STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn $R1 $R2 --runDirPerm All_RWX --readFilesCommand zcat $NOBAM \
     --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloBarcodeReadLength 0 --soloCBlen 8 --soloUMIstart 9 --soloUMIlen $UMILEN --soloStrand $STR \
     --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
     --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
     --soloFeatures Gene GeneFull Velocyto --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx

## gzip all outputs
cd output
for i in Gene/raw Gene/filtered GeneFull/raw GeneFull/filtered Velocyto/raw Velocyto/filtered
do 
  cd $i; for j in *; do gzip $j & done
  cd ../../
done

## index the BAM file
if [[ -s Aligned.sortedByCoord.out.bam ]]
then
  samtools index -@16 Aligned.sortedByCoord.out.bam
fi

wait
echo "ALL DONE!"
 
