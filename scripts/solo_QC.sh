#!/bin/bash 

TAG=$1 ## some common element of folder name - e.g. if samples are called SRR124**, use SRR124 or SRR1
if [[ $TAG == "" ]]
then
  >&2 echo "Usage: ./solo_QC.sh <starsolo_ouput_tag>"
  >&2 echo "(starsolo_ouput_tag is some common element of folder name - e.g. if samples are called SRR124**, use SRR124 or SRR1)"
  exit 1
fi

KK=`ls | grep $TAG`

echo -e "Sample\tRd_all\tRd_in_cells\tUMI_in_cells\tCells\tMed_nFeature\tGood_BC\tall_u+m\tall_u\texon_u+m\texon_u\tfull_u+m\tfull_u"

for i in $KK
do
  R1=`grep "Number of Reads," $i/output/Gene/Summary.csv | awk -F "," '{print $2}'`
  B=`grep "Reads With Valid Barcodes," $i/output/Gene/Summary.csv | awk -F "," '{print $2}'`
  G1=`grep "Reads Mapped to Genome: Unique+Multiple," $i/output/Gene/Summary.csv | awk -F "," '{print $2}'`
  G2=`grep "Reads Mapped to Genome: Unique," $i/output/Gene/Summary.csv | awk -F "," '{print $2}'`
  T1=`grep "Reads Mapped to Gene: Unique+Multipe Gene," $i/output/Gene/Summary.csv | awk -F "," '{print $2}'`
  T2=`grep "Reads Mapped to Gene: Unique Gene," $i/output/Gene/Summary.csv | awk -F "," '{print $2}'`
  F1=`grep "Reads Mapped to GeneFull: Unique+Multipe GeneFull," $i/output/GeneFull/Summary.csv | awk -F "," '{print $2}'`
  F2=`grep "Reads Mapped to GeneFull: Unique GeneFull," $i/output/GeneFull/Summary.csv | awk -F "," '{print $2}'`
  C=`grep "Estimated Number of Cells," $i/output/GeneFull/Summary.csv | awk -F "," '{print $2}'`
  R2=`grep "Unique Reads in Cells Mapped to GeneFull," $i/output/GeneFull/Summary.csv | awk -F "," '{print $2}'`
  R3=`grep "UMIs in Cells," $i/output/GeneFull/Summary.csv | awk -F "," '{print $2}'`
  GC=`grep "Median GeneFull per Cell," $i/output/GeneFull/Summary.csv | awk -F "," '{print $2}'`
  echo -e "$i\t$R1\t$R2\t$R3\t$C\t$GC\t$B\t$G1\t$G2\t$T1\t$T2\t$F1\t$F2"
done
