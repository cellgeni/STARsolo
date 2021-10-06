# Synchronized processing of scRNA-seq and bulk RNA-seq 

These are the scripts used for CellGenIT for synchronized processing of scRNA-seq and bulk RNA-seq. Both use [STAR](https://github.com/alexdobin/STAR) aligner to align reads to the reference genome. 

## Software installation

### STAR and RSEM versions

`STAR` of version 2.7.9a or above is recommended. The newest update includes the ability to correctly process multi-mapping reads, and adds many important options and bug fixes. 

In order to use settings that closely mimic those of `Cell Ranger` v4 or above (see explanations below, particularly `--clipAdapterType CellRanger4` option), `STAR` needs to be re-compiled from source with `make STAR CXXFLAGS_SIMD="-msse4.2"` (see [this issue](https://github.com/alexdobin/STAR/issues/1218) for more info). If you get the "Illegal instruction" error, that's what you need to do. 

There's also Martin Prete's awesome `icpc`-compiled version of `STAR` that's being tested right now - stay tuned for the updates. 

`RSEM` of version 1.3.3 should be installed using bioconda in a separate virtual environment. 

## Reference genome and annotation

The issues of reference genome and annotation used for **mouse** and **human** scRNA-seq experiments are described [here](https://www.singlecellcourse.org/processing-raw-scrna-seq-sequencing-data-from-reads-to-a-count-matrix.html) in some detail. 

Briefly,

  - there are several standard annotation versions used by `Cell Ranger`; They are referred to as versions `1.2.0`, `2.1.0`, `3.0.0`, and `2020-A`; 
  - the references are filtered to remove pseudogenes and small RNAs; exact filtering scripts are available [here](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#header); 
  - number of genes in `Cell Ranger` filtered human reference is about 35k, in full human genome is about 60k.

Depending on the task at hand, bulk RNA-seq can be processed either using the `Cell Ranger` filtered reference, or full annotation with 60k genes. Usually we do the latter.

Once you have downloaded the needed fasta and GTF files, and applied the necessary filtering (see the 10x genomics link above), run the following command using 16 CPUs/64 Gbs of RAM: 

`STAR --runThreadN 16 --runMode genomeGenerate --genomeDir STAR --genomeFastaFiles $FA --sjdbGTFfile $GTF`

## Processing scRNA-seq with STARsolo

### Reprodicing `Cell Ranger` v4 and above (but much faster)

Full script with the latest settings is available in `/scripts/starsolo_newest.sh`. It contains *many* options that frequently change; some of which will be explained below. In general, it's tuned in such way that the results with be very close to those of `Cell Ranger` v4 and above. 

Before running, barcode whitelists need to be downloaded [from here](https://github.com/10XGenomics/cellranger/tree/master/lib/python/cellranger/barcodes). 

Below are the explanations for some of the options (note that 5' experiments **always** use `737K-august-2016.txt` barcode file): 
  - `--runDirPerm All_RWX` allows a directory readable by all users, which becomes an issue when sharing results on Farm; 
  - `--soloCBwhitelist $BC --soloBarcodeReadLength 0 --soloUMIlen $UMILEN --soloStrand $STR` are settings that change with the used 10x chemistry:
    - BC=`737K-april-2014_rc.txt`, UMILEN=10, STR=Forward for 3' v1; 
    - BC=`737K-august-2016.txt`, UMILEN=10, STR=Forward for 3' v2; 
    - BC=`3M-february-2018.txt`, UMILEN=12, STR=Forward for 3' v3 and v3.1; 
    - BC=`737K-august-2016.txt`, UMILEN=10, STR=Reverse, for 5' v2;
    - BC=`737K-august-2016.txt`, UMILEN=12, STR=Reverse, for 5' v3. 
  - `--soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30` are options that define UMI collapsing, barcode collapsing, and read clipping algorithms that are closest to ones used by `Cell Ranger`; 
  - `--soloCellFilter EmptyDrops_CR` specifies the cell filtering algorithm used in [EmptyDrops](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html), which is the default algorithm in later versions of `Cell Ranger`; 
  - `--soloFeatures Gene GeneFull Velocyto` output conventional (exon-only) UMI counts, as well as exon+intron UMI counts (analog of `Cell Ranger` premrna option), as well as matrices preprocessed for `Velocyto`; 
  - `--soloMultiMappers Unique EM` is to count multimappers; 
  - `--readFilesCommand zcat` is used if your input fastq files are gzipped;
  - options grouped as `$SORTEDBAM` should be used if you need a genomic bam file; otherwise, use `$NOBAM`.  

Actual `STAR` command being run (note the order of read files passed to `--readFilesIn`): 

```bash
STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn $R2 $R1 --runDirPerm All_RWX \
     --readFilesCommand zcat $NOBAM --soloType CB_UMI_Simple --soloCBwhitelist $BC \
     --soloBarcodeReadLength 0 --soloUMIlen $UMILEN --soloStrand $STR --soloUMIdedup 1MM_CR \
     --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
     --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
     --soloFeatures Gene GeneFull Velocyto \
     --soloOutFileNames output/ genes.tsv barcodes.tsv matrix.mtx
```

### Counting the multimapping reads

Default approach used by `Cell Ranger` (and `STARsolo` scripts above) is to discard all reads that map to multiple genomic locations with equal mapping quality. This approach creates a bias in gene expression estimation. Pseudocount-based methods correctly quantify multimapping reads, but generate false counts due to pseudo-alignment errors. These issues are described in good detail [here](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1). 

If you would like to process multimappers, add the following options: `--soloMultiMappers Uniform EM`. This will generate an extra matrix in the /raw output folders. There will be non-integer numbers in the matrix because of split reads. If the downstream processing requires integers, you can round with a tool of your liking (e.g. `awk`). 

## Processing bulk RNA-seq with STAR/RSEM

`RSEM` reference files need to be prepared from genome fasta and GTF using the following command: 

`rsem-prepare-reference --gtf $GTF $FA <rsem_ref_name>`

Full script with the latest settings for STAR/RSEM processing of bulk RNA-seq is available in `/scripts/star_rsem_bulk.sh`. Most options for `STAR` alignment are following the **ENCODE** presets (see STAR manual for more details). 

Actual `STAR` command being run: 

```bash
STAR --runThreadN $CPUS --genomeDir $SREF --readFilesIn $R1 $R2 --readFilesCommand zcat \
     --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 \
     --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt \
     --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout \
     --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --sjdbScore 1\
     --quantMode TranscriptomeSAM GeneCounts --limitBAMsortRAM 30000000000
```

`RSEM` command for quantification using the transcriptomic BAM file: 

```bash
rsem-calculate-expression $PAIRED -p $CPUS --bam --estimate-rspd --seed 12345 -p 4 --no-bam-output \
     --forward-prob 0.5 Aligned.toTranscriptome.out.bam $RREF $TAG.rsem
```

`$RREF` here is `RSEM` reference prepared with the `rsem-prepare-reference` shown above. `$PAIRED` is set to "--paired-end" for paired-end experiments and to empty string if the experiment is single-end. For strand-specific processing, `--forward-prob` can be changed to 1 (forward read matches the direction of the gene), or 0 (reverse strand specificity, common for dUTP-based protocols).

`RSEM` output generates both per-gene and per-transcript tables, with raw read counts and TPM and FPKM normalized counts.  
