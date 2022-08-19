#!/bin/bash 

echo -e "Sample\tRd_all\tRd_in_cells\tFrc_in_cells\tUMI_in_cells\tCells\tMed_nFeature\tGood_BC\tSpecies\tPaired\tStrand\tall_u+m\tall_u\texon_u+m\texon_u\tfull_u+m\tfull_u"


for i in *
do
  if [[ -d $i && -d $i/output ]]
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
    echo -e "$i\t$R1\t$R2\t$CF\t$R3\t$C\t$GC\t$B\t$REF\t$PAIRED\t$ST\t$G1\t$G2\t$E1\t$E2\t$F1\t$F2"
  fi
done
