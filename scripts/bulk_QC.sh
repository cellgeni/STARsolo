#!/bin/bash 

echo -e "Sample\tN_reads\tN_uniq\tN_multi\tN_gene\tPct_uniq\tPct_multi\tPct_gene"

for i in *
do
  if [[ -d $i && -s $i/Log.out ]]
  then 
    N1=`grep "Number of input reads " $i/Log.final.out | awk '{print $6}'`	
    N2=`grep "Uniquely mapped reads number " $i/Log.final.out | awk '{print $6}'`
    N3=`grep "Number of reads mapped to multiple loci " $i/Log.final.out | awk '{print $9}'`
    N4=`cat $i/$i.rsem.genes.results | awk '{if (NR>1) {sum+=$5}} END {printf "%d\n",sum}'`
    
    P1=`echo $N2 | awk -v v=$N1 '{printf "%.2f\n",100*$1/v}'`
    P2=`echo $N3 | awk -v v=$N1 '{printf "%.2f\n",100*$1/v}'`
    P3=`echo $N4 | awk -v v=$N1 '{printf "%.2f\n",100*$1/v}'`
    echo -e "$i\t$N1\t$N2\t$N3\t$N4\t$P1\t$P2\t$P3"
  fi
done 
