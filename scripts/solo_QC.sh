#!/bin/bash 

## check that STAR temporary dir is removed for all samples, and that archived unmapped reads are created
>&2 echo "Checking that all STARsolo jobs went to completion .." 
for i in *
do
  if [[ -d $i && -d $i/output && -s $i/Log.final.out ]]
  then
    if [[ -d $i/_STARtmp ]]
    then
      >&2 echo "WARNING: Sample $i did not run to completion: _STARtmp is still present!" 
    fi 

    if [[ ! -s $i/Unmapped.out.mate1.bz2 || ! -s $i/Unmapped.out.mate2.bz2 ]]
    then
      >&2 echo "WARNING: Unmapped reads (Unmapped.out.mate[1,2].bz2) not found for sample $i!" 
    fi 
  fi
done


## check if samples overlap - this is often a sign that something went quite wrong 
## e.g. starsolo script has chosen wrong reads, or they were swapped during upload, or some lab mix-up
>&2 echo "Checking potential sample cross-contamination .." 
for i in *
do
  for j in *
  do
    if [[ -d $i/output && -s $i/Log.final.out && -d $j/output && -s $j/Log.final.out && $i != $j ]]
    then
      N1=`zcat $i/output/Gene/filtered/barcodes.tsv.gz | wc -l`
      N2=`zcat $j/output/Gene/filtered/barcodes.tsv.gz | wc -l`
      MIN=`(( $N1 <= $N2 )) && echo $N1 || echo $N2`
      COMM=`comm -12 <(zcat $i/output/Gene/filtered/barcodes.tsv.gz) <(zcat $j/output/Gene/filtered/barcodes.tsv.gz) | wc -l`
      PCT=`echo $COMM | awk -v v=$MIN '{printf "%d\n",100*$1/v+0.5}'`
      if (( $PCT >= 20 )) 
      then
        >&2 echo "WARNING: Samples $i ($N1 barcodes) and $j ($N2 barcodes) have $COMM ($PCT%) common barcodes, which is higher than expected by chance! Please investigate .."
      fi 
    fi 
  done
done


## finally, calculate and output STARsolo stats 
>&2 echo "Extracting STARsolo stats .." 
>&2 echo 
echo -e "Sample\tRd_all\tRd_in_cells\tFrc_in_cells\tUMI_in_cells\tCells\tMed_nFeature\tGood_BC\tWL\tSpecies\tPaired\tStrand\tall_u+m\tall_u\texon_u+m\texon_u\tfull_u+m\tfull_u"


for i in *
do
  if [[ -d $i && -d $i/output && -s $i/Log.final.out ]]
  then 
    PAIRED="Single"
    if [[ `grep "clip5pNbases 39 0" $i/Log.out` != "" ]]
    then
      PAIRED="Paired"
    fi

    REF="Other"
    if [[ `grep "^genomeDir" $i/Log.out | tail -n1 | grep "/human/"` != "" ]]
    then
      REF="Human"
    elif [[ `grep "^genomeDir" $i/Log.out | tail -n1 | grep "/mouse/"` != "" ]]
    then
      REF="Mouse"
    fi

    WL="Undef"
    if [[ `grep "^soloCBwhitelist" $i/Log.out | tail -n1 | grep 3M-february-2018.txt` != "" ]]
    then
      WL="v3"
    elif [[ `grep "^soloCBwhitelist" $i/Log.out | tail -n1 | grep 737K-august-2016.txt` != "" ]]
    then
      WL="v2"
    elif [[ `grep "^soloCBwhitelist" $i/Log.out | tail -n1 | grep 737K-april-2014_rc.txt` != "" ]]
    then
      WL="v1"
    elif [[ `grep "^soloCBwhitelist" $i/Log.out | tail -n1 | grep 737K-arc-v1.txt` != "" ]]
    then
      WL="arc"
    fi
    

    R1=`grep "Number of Reads," $i/output/Gene/Summary.csv | awk -F "," '{print $2}'`
    B=`grep "Reads With Valid Barcodes," $i/output/Gene/Summary.csv | awk -F "," '{print $2}'`
    G1=`grep "Reads Mapped to Genome: Unique+Multiple," $i/output/Gene/Summary.csv | awk -F "," '{print $2}'`
    G2=`grep "Reads Mapped to Genome: Unique," $i/output/Gene/Summary.csv | awk -F "," '{print $2}'`
    E1=`grep "Reads Mapped to Gene: Unique+Multip.*e Gene," $i/output/Gene/Summary.csv | awk -F "," '{print $2}'`
    E2=`grep "Reads Mapped to Gene: Unique Gene," $i/output/Gene/Summary.csv | awk -F "," '{print $2}'`
    F1=`grep "Reads Mapped to GeneFull: Unique+Multip.*e GeneFull," $i/output/GeneFull/Summary.csv | awk -F "," '{print $2}'`
    F2=`grep "Reads Mapped to GeneFull: Unique GeneFull," $i/output/GeneFull/Summary.csv | awk -F "," '{print $2}'`
    C=`grep "Estimated Number of Cells," $i/output/GeneFull/Summary.csv | awk -F "," '{print $2}'`
    R2=`grep "Unique Reads in Cells Mapped to GeneFull," $i/output/GeneFull/Summary.csv | awk -F "," '{print $2}'`
    CF=`echo $R1 | awk -v v=$R2 '{printf "%.3f\n",v/$1}'`
    R3=`grep "UMIs in Cells," $i/output/GeneFull/Summary.csv | awk -F "," '{print $2}'`
    GC=`grep "Median GeneFull per Cell," $i/output/GeneFull/Summary.csv | awk -F "," '{print $2}'`
    ST=`grep "^soloStrand" $i/Log.out | grep RE-DEFINED | awk '{print $2}'`

    if [[ $ST == "" ]]
    then
      ST="Undef"
    elif [[ $PAIRED == "Paired" ]]
    then
      ## only 5' experiments can be processed as paired-end; however, actual STAR command has "--soloStrand Forward"
      ## since read order for PE processing is R1 R2 (it's R2 R1 for regular single-end 10X)
      ST="Reverse"
    fi

		## warn about odd mapping stats; should add few more
		F2PCT=`echo $F2 | awk '{printf "%d\n",$1*100+0.5}'`
		if (( $F2PCT <= 20 ))
		then
			>&2 echo "WARNING: Sample $i : GeneFull percentage ($F2PCT) is too low! Please make sure the strand-specificity evaluation worked correctly."
		fi

    echo -e "$i\t$R1\t$R2\t$CF\t$R3\t$C\t$GC\t$B\t$WL\t$REF\t$PAIRED\t$ST\t$G1\t$G2\t$E1\t$E2\t$F1\t$F2"
  fi
done
