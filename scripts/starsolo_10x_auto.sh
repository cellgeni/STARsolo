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
FQDIR=/lustre/scratch117/cellgen/cellgeni/TIC-starsolo/tic-1328/fastqs  ### change to the directory with fastq files/folders
## choose one of the two otions, depending on whether you need a BAM file 
#BAM="--outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 2 --limitBAMsortRAM 120000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB GX GN"
BAM="--outSAMtype None"

###################### DONT CHANGE OPTIONS BELOW THIS LINE ###########################

mkdir $TAG && cd $TAG

## three popular cases: <sample>_1.fastq/<sample>_2.fastq, <sample>.R1.fastq/<sample>.R2.fastq, and <sample>_L001_R1_S001.fastq/<sample>_L001_R2_S001.fastq
## the command below will generate a comma-separated list for each read
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

## also define one file from R1/R2; we choose the largest one  
R1F=`echo $R1 | tr ',' ' ' | xargs ls -s | tail -n1 | awk '{print $2}'`
R2F=`echo $R2 | tr ',' ' ' | xargs ls -s | tail -n1 | awk '{print $2}'`

## let's see if the files are archived or not. Gzip is the only one we test for, but bgzip archives should work too since they are gzip-compatible.
GZIP=""
BC=""
NBC1=""
NBC2=""
NBC3=""
NBCA=""
R1LEN=""
R2LEN=""
R1DIS=""

## depending on whether the files are archived or not,  
if [[ `find $FQDIR/* | grep $TAG | grep "\.gz$"` != "" ]]
then  
  GZIP="--readFilesCommand zcat"
  NBC1=`zcat $R1F | awk 'NR%4==2' | grep -v N | head -n10000 | grep -F -f /nfs/cellgeni/STAR/whitelists/737K-april-2014_rc.txt | wc -l`
  NBC2=`zcat $R1F | awk 'NR%4==2' | grep -v N | head -n10000 | grep -F -f /nfs/cellgeni/STAR/whitelists/737K-august-2016.txt | wc -l`
  NBC3=`zcat $R1F | awk 'NR%4==2' | grep -v N | head -n10000 | grep -F -f /nfs/cellgeni/STAR/whitelists/3M-february-2018.txt | wc -l`
  NBCA=`zcat $R1F | awk 'NR%4==2' | grep -v N | head -n10000 | grep -F -f /nfs/cellgeni/STAR/whitelists/737K-arc-v1.txt | wc -l`
  R1LEN=`zcat $R1F | awk 'NR%4==2' | head -n1000 | awk '{sum+=length($0)} END {printf "%d\n",sum/NR+0.5}'`
  R2LEN=`zcat $R2F | awk 'NR%4==2' | head -n1000 | awk '{sum+=length($0)} END {printf "%d\n",sum/NR+0.5}'`
  R1DIS=`zcat $R1F | awk 'NR%4==2' | head -n1000 | awk '{print length($0)}' | sort | uniq -c | wc -l`
  zcat $R1F | head -n 100000 > test.R1.fastq
  zcat $R2F | head -n 100000 > test.R2.fastq 
else 
  NBC1=`cat $R1F | awk 'NR%4==2' | grep -v N | head -n10000 | grep -F -f /nfs/cellgeni/STAR/whitelists/737K-april-2014_rc.txt | wc -l`
  NBC2=`cat $R1F | awk 'NR%4==2' | grep -v N | head -n10000 | grep -F -f /nfs/cellgeni/STAR/whitelists/737K-august-2016.txt | wc -l`
  NBC3=`cat $R1F | awk 'NR%4==2' | grep -v N | head -n10000 | grep -F -f /nfs/cellgeni/STAR/whitelists/3M-february-2018.txt | wc -l`
  NBCA=`cat $R1F | awk 'NR%4==2' | grep -v N | head -n10000 | grep -F -f /nfs/cellgeni/STAR/whitelists/737K-arc-v1.txt | wc -l`
  R1LEN=`cat $R1F | awk 'NR%4==2' | head -n1000 | awk '{sum+=length($0)} END {printf "%d\n",sum/NR+0.5}'`
  R2LEN=`cat $R2F | awk 'NR%4==2' | head -n1000 | awk '{sum+=length($0)} END {printf "%d\n",sum/NR+0.5}'`
  R1DIS=`cat $R1F | awk 'NR%4==2' | head -n1000 | awk '{print length($0)}' | sort | uniq -c | wc -l`
  cat $R1F | head -n 100000 > test.R1.fastq
  cat $R2F | head -n 100000 > test.R2.fastq 
fi

## elucidate the right barcode whitelist to use. Grepping out N saves us some trouble. Note the special list for multiome experiments (737K-arc-v1.txt):
if (( $NBC2 > 5000 )) 
then 
  BC=/nfs/cellgeni/STAR/whitelists/737K-august-2016.txt
elif (( $NBC3 > 5000 ))
then
  BC=/nfs/cellgeni/STAR/whitelists/3M-february-2018.txt
elif (( $NBCA > 5000 ))
then
  BC=/nfs/cellgeni/STAR/whitelists/737K-arc-v1.txt
elif (( $NBC1 > 5000 )) 
then
  BC=/nfs/cellgeni/STAR/whitelists/737K-april-2014_rc.txt
else 
  >&2 echo "ERROR: No whitelist has matched first 10000 barcodes!"
  exit 1
fi 

## check read lengths, fail if something funky is going on: 
PAIRED=False
UMILEN=""
CBLEN=""
if (( $R1DIS > 1 && $R1LEN <= 30 ))
then 
  >&2 echo "ERROR: Read 1 (barcode) has varying length; possibly someone thought it's a good idea to quality-trim it. Please check the fastq files."
  exit 1
elif (( $R1LEN < 24 )) 
then
  >&2 echo "ERROR: Read 1 (barcode) is less than 24 bp in length. Please check the fastq files."
  exit 1
elif (( $R2LEN < 50 )) 
then
  >&2 echo "ERROR: Read 2 (biological read) is less than 50 bp in length. Please check the fastq files."
  exit 1
fi

## assign the necessary variables for barcode/UMI length/paired-end processing:  
if (( $R1LEN > 50 )) 
then
  PAIRED=True
  UMILEN=10
  CBLEN=16
elif (( $NBC1 > 5000 )) 
then
  CBLEN=14
  UMILEN=$((R1LEN-14))
else
  CBLEN=16
  UMILEN=$((R1LEN-16)) 
fi 

## finally, see if you have 5' or 3' experiment. I don't know and easier way:  
STRAND=Forward

STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn test.R2.fastq test.R1.fastq --runDirPerm All_RWX --outSAMtype None \
     --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloBarcodeReadLength 0 --soloCBlen $CBLEN --soloUMIstart $((CBLEN+1)) \
     --soloUMIlen $UMILEN --soloStrand Forward \
     --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
     --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
     --soloFeatures Gene --soloOutFileNames test_strand/ features.tsv barcodes.tsv matrix.mtx &> /dev/null 
rm test.R1.fastq test.R2.fastq

GENEPCT=`grep "Reads Mapped to Gene: Unique+Multipe Gene" test_strand/Gene/Summary.csv | awk -F "," '{printf "%d\n",$2*100}'`
if (( $GENEPCT < 10 )) 
then
  STRAND=Reverse
fi

## finally, if paired-end experiment turned out to be 3', process it as single-end: 
if [[ $STRAND == "Forward" && $PAIRED == "True" ]]
then
  PAIRED=False
  if [[ $BC == "/nfs/cellgeni/STAR/whitelists/3M-february-2018.txt" ]] 
  then
    UMILEN=12
  fi
fi

echo "Done setting up the STARsolo run; here are final processing options:"
echo "============================================================================="
echo "Sample: $TAG"
echo "Paired-end mode: $PAIRED"
echo "Strand (Forward = 3', Reverse = 5'): $STRAND"
echo "CB whitelist: $BC"
echo "CB length: $CBLEN"
echo "UMI length: $UMILEN"
echo "GZIP: $GZIP"
echo "-----------------------------------------------------------------------------"
echo "Read 1 files: $R1"
echo "-----------------------------------------------------------------------------"
echo "Read 2 files: $R2" 
echo "-----------------------------------------------------------------------------"

if [[ $PAIRED == "True" ]]
then
  STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn $R1 $R2 --runDirPerm All_RWX $GZIP $BAM --soloBarcodeMate 1 --clip5pNbases 39 0 \
     --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloCBstart 1 --soloCBlen $CBLEN --soloUMIstart $((CBLEN+1)) --soloUMIlen $UMILEN --soloStrand Forward \
     --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
     --soloCellFilter EmptyDrops_CR --outFilterScoreMin 30 \
     --soloFeatures Gene GeneFull Velocyto --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx
else 
  STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn $R2 $R1 --runDirPerm All_RWX $GZIP $BAM \
     --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloBarcodeReadLength 0 --soloCBlen $CBLEN --soloUMIstart $((CBLEN+1)) --soloUMIlen $UMILEN --soloStrand $STRAND \
     --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
     --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
     --soloFeatures Gene GeneFull Velocyto --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx
fi

## index the BAM file
if [[ -s Aligned.sortedByCoord.out.bam ]]
then
  samtools index -@16 Aligned.sortedByCoord.out.bam
fi

## finally, let's gzip all outputs
cd output
for i in Gene/raw Gene/filtered GeneFull/raw GeneFull/filtered Velocyto/raw Velocyto/filtered
do 
  cd $i; for j in *; do gzip $j & done
  cd ../../
done

wait
echo "ALL DONE!"
