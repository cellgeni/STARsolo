#!/bin/bash 

echo -e "Sample\tN_reads\tSTAR_uniq\tPct_uniq\tSTAR_multi\tPct_multi\tFC_assn\tPct_assn\tFC_unmap\tPct_unmap\tFC_multi\tPct_multi\tFC_nofeat\tPct_nofeat\tFC_ambig\tPct_ambig"

for i in *
do
  if [[ -d $i && -s $i/Log.out ]]
  then 
    N1=`grep "Number of input reads " $i/Log.final.out | awk '{print $6}'`	
    N2=`grep "Uniquely mapped reads number " $i/Log.final.out | awk '{print $6}'`
    N3=`grep "Number of reads mapped to multiple loci " $i/Log.final.out | awk '{print $9}'`

    N4=`grep "^Assigned" $i/$i.feature_counts.tsv.summary | awk '{print $2}'`
    N5=`grep "^Unassigned_Unmapped" $i/$i.feature_counts.tsv.summary | awk '{print $2}'`
    N6=`grep "^Unassigned_MultiMapping" $i/$i.feature_counts.tsv.summary | awk '{print $2}'`
    N7=`grep "^Unassigned_NoFeatures" $i/$i.feature_counts.tsv.summary | awk '{print $2}'`
    N8=`grep "^Unassigned_Ambiguity" $i/$i.feature_counts.tsv.summary | awk '{print $2}'`
    
    P2=`echo $N2 | awk -v v=$N1 '{printf "%.2f\n",100*$1/v}'`
    P3=`echo $N3 | awk -v v=$N1 '{printf "%.2f\n",100*$1/v}'`
    P4=`echo $N4 | awk -v v=$N1 '{printf "%.2f\n",100*$1/v}'`
    P5=`echo $N5 | awk -v v=$N1 '{printf "%.2f\n",100*$1/v}'`
    P6=`echo $N6 | awk -v v=$N1 '{printf "%.2f\n",100*$1/v}'`
    P7=`echo $N7 | awk -v v=$N1 '{printf "%.2f\n",100*$1/v}'`
    P8=`echo $N8 | awk -v v=$N1 '{printf "%.2f\n",100*$1/v}'`
    echo -e "$i\t$N1\t$N2\t$P2\t$N3\t$P3\t$N4\t$P4\t$N5\t$P5\t$N6\t$P6\t$N7\t$P7\t$N8\t$P8"
  fi
done 
