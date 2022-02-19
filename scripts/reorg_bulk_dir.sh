#!/bin/bash 

## change sample tag for reuse!
## reorganize directories for bulk RNA-seq STAR+rsem processing 

mkdir counts rsem_gene rsem_tran

for i in PR*
do
  echo "Processing sample $i - moving files.."
  cp $i/*genes.results rsem_gene
  cp $i/*isoforms.results rsem_tran
done 

cd rsem_gene
## we assume Ensembl IDs without dots

for i in *genes.results
do
  TAG=${i%%.rsem.genes.results}
  echo "Processing file $TAG - generating count tables.." 
  echo -e "geneName\tgeneEffLength\t$TAG" > $TAG.counts.tsv
  grep -v expected_count $i | grep -v _PAR_ | cut -f 1,4,5 >> $TAG.counts.tsv
done 

mv *.counts.tsv ../counts 

echo "ALL DONE!" 
  
