#!/bin/bash -e 

## v3.1 of STARsolo wrappers is set up to guess the chemistry automatically
## newest version of the script uses STAR v2.7.10a with EM multimapper processing where possible 

SIF="/nfs/cellgeni/singularity/images/starsolo_2-7-10a-alpha-220818_samtools_1-15-1_seqtk-1-13_bbmap_38-97_RSEM-1-3-3.sif"
CMD="singularity run --bind /nfs,/lustre $SIF"

## you need to have a manifest file named $TAG.manifest.tsv
## with 3 tab-separated columns: 1) full path to read1; 2) full path to read2; 3) cell name

TSV=$1
TAG=${TSV%%.manifest.tsv} 

if [[ ! -s $TSV ]]
then 
  >&2 echo "Usage: ./starsolo_ss2.sh <manifest_tsv>"
  exit 1
fi

TSV=`readlink -f $TSV` 
CPUS=16                                                                ## typically bsub this into normal queue with 16 cores and 64 Gb RAM.   
REF=/nfs/cellgeni/STAR/human/2020A/index                               ## choose the appropriate reference 

## choose one of the two otions, depending on whether you need a BAM file 
## BAM options are for 10x and not tested with other methods
#BAM="--outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 500 --limitBAMsortRAM 60000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM GX GN RG"
BAM="--outSAMtype None"

###################################################################### DONT CHANGE OPTIONS BELOW THIS LINE ##############################################################################################

mkdir $TAG && cd $TAG

GZIP=""
if [[ `grep -P "\.gz\t" $TSV` != "" ]]
then  
  GZIP="--readFilesCommand zcat"
fi

## outFilter* options can be adjusted according to the mapping rate and mapped length
## alignment often works well without clipping; if not, --clip3pAdapterSeq <3' adapter sequence> option can be added, or separate trimming using bbduk.sh works well too
$CMD STAR --runThreadN $CPUS --genomeDir $REF --runDirPerm All_RWX $GZIP $BAM \
     --soloType SmartSeq --readFilesManifest $TSV --soloUMIdedup Exact --soloStrand Unstranded \
     --limitOutSJcollapsed 10000000 --soloCellFilter None \
     --soloFeatures Gene GeneFull --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx --outReadsUnmapped Fastx

## index the BAM file
if [[ -s Aligned.sortedByCoord.out.bam ]]
then
  $CMD samtools index -@16 Aligned.sortedByCoord.out.bam
fi

## max-CR bzip all unmapped reads with multicore pbzip2 
pbzip2 -9 Unmapped.out.mate1 &
pbzip2 -9 Unmapped.out.mate2 &
wait

cd output
for i in Gene/raw GeneFull/raw
do 
  cd $i; for j in *; do gzip $j & done
  cd ../../
done

wait
echo "ALL DONE!"
