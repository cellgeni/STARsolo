# Wrapper scripts for using `STARsolo` with various types of single cell RNA-seq

These are the scripts used for CellGenIT for uniform processing of scRNA-seq - both 10X and quite a few other types (see below for supported platforms). Additionally, bulk RNA-seq could be processed using a 10X reference in a matched way - this should improve scRNA-seq-to-bulk mapping. All listed methods use [STAR](https://github.com/alexdobin/STAR) aligner to align reads to the reference genome. 

## Software installation

### STAR and RSEM versions

`STAR` of version 2.7.9a or above is recommended (2.7.10a is the latest and greatest, as of August'22). The newest update includes the ability to correctly process multi-mapping reads, and adds many important options and bug fixes. 

In order to use settings that closely mimic those of `Cell Ranger` v4 or above (see explanations below, particularly `--clipAdapterType CellRanger4` option), `STAR` needs to be re-compiled from source with `make STAR CXXFLAGS_SIMD="-msse4.2"` (see [this issue](https://github.com/alexdobin/STAR/issues/1218) for more info). If you get the _Illegal instruction_ error, that's what you need to do. 

There's also Martin Prete's awesome `icpc`-compiled version of `STAR` that's being tested right now - stay tuned for the updates. 

`RSEM` of version 1.3.3 should be installed using bioconda in a separate virtual environment. 

## Reference genome and annotation

The issues of reference genome and annotation used for **mouse** and **human** scRNA-seq experiments are described [here](https://www.singlecellcourse.org/processing-raw-scrna-seq-sequencing-data-from-reads-to-a-count-matrix.html) in some detail. 

Briefly,

  - there are several standard annotation versions used by `Cell Ranger`; They are referred to as versions `1.2.0`, `2.1.0`, `3.0.0`, and `2020-A`; 
  - the references are filtered to remove pseudogenes and small RNAs; exact filtering scripts are available [here](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#header); 
  - number of genes in `Cell Ranger` filtered human reference is about 35k; in full human annotation is about 60k.

Depending on the task at hand, bulk RNA-seq can be processed either using the `Cell Ranger` filtered reference, or full annotation with 60k genes. Usually we do the latter.

Once you have downloaded the needed fasta and GTF files, and applied the necessary filtering (see the 10x genomics link above), run the following command using 16 CPUs/64 Gbs of RAM: 

`STAR --runThreadN 16 --runMode genomeGenerate --genomeDir STAR --genomeFastaFiles $FA --sjdbGTFfile $GTF`

`STARsolo` reference is not different from a normal `STAR` reference - `STARsolo` workflow is invoked using the options listed below. Thus, it can be run on both `Cell Ranger` filtered and un-filtered index. We use the filtered reference to match the `Cell Ranger` output as closely as possible. 

### Pre-made reference location

All **CellGenIT** pre-made `STAR` references are located in `/nfs/cellgeni/STAR/`. Barcode whitelist files are located in `/nfs/cellgeni/STAR/whitelists`. 

## Processing scRNA-seq with STARsolo

### 10X: reproducing `Cell Ranger` v4 and above (but much faster)

Full scripts with the latest settings are available in `/scripts` (there are several scripts according to 10x chemistry version; e.g. `starsolo_3p_v3.sh` should be used for v3 of 3' 10x, while `starsolo_5p_v2.sh` should be used for v2 of 5'. The scripts contain *many* options that frequently change; some of which will be explained below. In general, commands are tuned in such way that the results will be very close to those of `Cell Ranger` v4 and above. 

Before running, barcode whitelists need to be downloaded [from here](https://github.com/10XGenomics/cellranger/tree/master/lib/python/cellranger/barcodes). 

Below are the explanations for some of the options (note that 5' experiments **always** use `737K-august-2016.txt` barcode file): 
  - `--runDirPerm All_RWX` allows a directory readable by all users, which becomes an issue when sharing results on Farm; 
  - `--soloCBwhitelist $BC --soloBarcodeReadLength 0 --soloUMIlen $UMILEN --soloStrand $STR` are settings that change with the used 10x chemistry:

<div align="center">

| 10X VERSION | BC | UMILEN | STR |
|:-:|:-:|:-:|:-:|
| 3' v1 | 737K-april-2014_rc.txt | 10 | Forward |
| 3' v2 | 737K-august-2016.txt | 10 | Forward |
| 3' v3, v3.1 | 3M-february-2018.txt | 12 | Forward |
| 5' v1.1, v2 | 737K-august-2016.txt | 10 | Reverse |
| 5' v3 | 737K-august-2016.txt | 12 | Reverse |
| multiome | 737K-arc-v1.txt | 12 | Forward |

</div>

  - `--soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30` are options that define UMI collapsing, barcode collapsing, and read clipping algorithms that are closest to ones used by `Cell Ranger`; 
  - `--soloCellFilter EmptyDrops_CR` specifies the cell filtering algorithm used in [EmptyDrops](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html), which is the default algorithm in later versions of `Cell Ranger`; 
  - `--soloFeatures Gene GeneFull Velocyto` output conventional (exon-only) UMI counts, as well as exon+intron UMI counts (analog of `Cell Ranger` premrna option), as well as matrices preprocessed for `Velocyto`; 
  - `--soloMultiMappers Unique EM` is to count multimappers (on by default in v3.0+ of these scripts; does not influence the main output, but creates an additional matrix in `/raw` subdir of `Gene` and `GeneFull`); 
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
     --soloOutFileNames output/ genes.tsv barcodes.tsv matrix.mtx --soloMultiMappers EM
```

### Using STARsolo for SMART-seq/SMART-seq2

For plate-based methods that don't use UMIs (such as [SMART-Seq and SMART-Seq2](https://teichlab.github.io/scg_lib_structs/methods_html/SMART-seq_family.html)), `STARsolo` can be used as well. Fastq files for these methods usually come as separate, paired-end files; all of these should be listed in a *manifest* file - plain text, tab-separated file containing three columns per line: 1) full path to R1; 2) full path to R2; 3) cell name or ID. 

Example of a script used to process Smart-seq2 data can be found in `/scripts/starsolo_ss2.sh`. Key parameters that could be adjusted are `--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3`; the higher they are, the less permissive is the alignment. Lower values can help you "rescue" a large proportion of reads with high adapter content (see below for adapter trimming). Actual `STAR` command being run:

```bash
STAR --runThreadN $CPUS --genomeDir $REF --runDirPerm All_RWX --readFilesCommand zcat $SORTEDBAM \
     --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
     --soloType SmartSeq --readFilesManifest ../$TAG.manifest.tsv \
     --soloUMIdedup Exact --soloStrand Unstranded \
     --soloFeatures Gene GeneFull --soloOutFileNames output/ genes.tsv barcodes.tsv matrix.mtx
```

Often, SMART-seq2 reads can benefit from trimming adapters, which can be turned on using `--clip3pAdapterSeq <3' adapter sequence>` option. Alternatively, `/scripts/bbduk.sh` can be used to trim adapters from reads prior to the alignment and quantification (this is the preferred option).  

### Counting the multimapping reads

Default approach used by `Cell Ranger` (and `STARsolo` scripts above) is to discard all reads that map to multiple genomic locations with equal mapping quality. This approach creates a bias in gene expression estimation. Pseudocount-based methods correctly quantify multimapping reads, but generate false counts due to pseudo-alignment errors. These issues are described in good detail [here](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1). 

If you would like to process multimappers, add the following options: `--soloMultiMappers Uniform EM` (or simply `--soloMultiMappers EM`; this is on by default in v3.0+ of these scripts). This will generate an extra matrix in the `/raw` output folders. There will be non-integer numbers in the matrix because of split reads. If the downstream processing requires integers, you can round with a tool of your liking (e.g. `awk`). 

As of `STAR` v2.7.10a, **multimapper counting still does not work for SMART-seq2 or bulk RNA-seq processing**. 

### Running STARsolo on other scRNA-seq platforms 

STARsolo is very flexible and can be used with almost any scRNA-seq method, provided you know the library structure - i.e. where cell barcodes, UMIs, and biological parts of the read are located in the sequencing fragment or reads. A great source of information about scRNA-seq library structures is [this page](https://teichlab.github.io/scg_lib_structs/).

Currently, our scripts directory provides dedicated scripts for
  - Drop-seq;
  - inDrops; 
  - STRT-seq. 

Please contact `CellGenIT` if you need to process an unusual dataset. 

## Quick evaluation of multiple STARsolo runs

If you've used these scripts to process multiple 10X samples, you can get a quick look at the results by copying `solo_QC.sh` script from this repo to the directory with `STARsolo` output folders, and running

```bash
./solo_QC.sh | column -t 
```

The **output_tag** argument could be any common part of the folder name - e.g. if samples are called SRR124444 .. SRR124534, you can use **SRR124** or **SRR**.

The script is designed for 10X or other droplet-based methods; the output will make a lot less sense for SMART-seq2. 

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
