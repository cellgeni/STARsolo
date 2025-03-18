# Wrapper scripts for `STARsolo` scRNA-seq pipeline

These are the scripts used for CellGenIT for uniform processing of scRNA-seq - both 10X and quite a few other types (see below for supported platforms). All listed methods use [STAR](https://github.com/alexdobin/STAR) aligner to align reads to the reference genome. 

## Software installation

### STAR version

`STAR` of version 2.7.9a or above is recommended (2.7.10a is the latest and greatest, as of August'22). The newest update includes the ability to correctly process multi-mapping reads, and adds many important options and bug fixes. 

In order to use settings that closely mimic those of `Cell Ranger` v4 or above (see explanations below, particularly `--clipAdapterType CellRanger4` option), `STAR` needs to be re-compiled from source with `make STAR CXXFLAGS_SIMD="-msse4.2"` (see [this issue](https://github.com/alexdobin/STAR/issues/1218) for more info). If you get the _Illegal instruction_ error, that's what you need to do. 

### Docker

We also have a Dockerfile which builds with specific versions of all the tools used. Once the container is built you can do `cat /versions.txt` and it will show the versions of the tools in the container. 

## Reference genome and annotation

The issues of reference genome and annotation used for **mouse** and **human** scRNA-seq experiments are described [here](https://www.singlecellcourse.org/processing-raw-scrna-seq-sequencing-data-from-reads-to-a-count-matrix.html) in some detail. 

Briefly,

  - there are several standard annotation versions used by `Cell Ranger`; They are referred to as versions `1.2.0`, `2.1.0`, `3.0.0`, and `2020-A`; 
  - the references are filtered to remove pseudogenes and small RNAs; exact filtering scripts are available [here](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#header); 
  - number of genes in `Cell Ranger` filtered human reference is about 35k; in full human annotation there are about 60k.

Once you have downloaded the needed fasta and GTF files, and applied the necessary filtering (see the 10x genomics link above), run the following command using 16 CPUs/64Gb of RAM: 

`STAR --runThreadN 16 --runMode genomeGenerate --genomeDir STAR --genomeFastaFiles $FA --sjdbGTFfile $GTF`

`STARsolo` reference is not different from a normal `STAR` reference - `STARsolo` workflow is invoked using the options listed below. Thus, it can be run on both `Cell Ranger` filtered and un-filtered index. We use the filtered reference to match the `Cell Ranger` output as closely as possible. 

### Pre-made reference location

All **CellGenIT** pre-made `STAR` references are located in `/nfs/cellgeni/STAR/`. Barcode whitelist files are located in `/nfs/cellgeni/STAR/whitelists`. 
## Resource requirements 

By default, all processing is done using 16 CPUs and 128Gb of RAM. Sanger Farm typically has 8 Gb of RAM per core, so 8 CPUs/64Gb RAM, or 16 CPUs/128Gb RAM is probably optimal. Unfortunately, many datasets fail with OOM error when 64Gb is requested, so we reverted our defaults to 128Gb. 

## Processing scRNA-seq with STARsolo

### 10X: reproducing `Cell Ranger` v4 and above (but much faster)

The current relevant wrapper script for all 10X scRNA-seq processing is called `starsolo_10x_auto.sh`. The scripts contain *many* options that sometimes change; some of which will be explained below. In general, commands are tuned in such way that the results will be very close to those of `Cell Ranger` v4 and above. 

Before running, barcode whitelists need to be downloaded [from here](https://github.com/10XGenomics/cellranger/tree/master/lib/python/cellranger/barcodes). 

Below are the explanations for some of the options (note that 5' experiments of versions v1.1/v2/v3 all use `737K-august-2016.txt` barcode file): 
  - `--runDirPerm All_RWX` allows a directory readable by all users, which becomes an issue when sharing results on Farm; 
  - `--soloCBwhitelist $BC --soloBarcodeReadLength 0 --soloUMIlen $UMILEN --soloStrand $STR` are settings that change with the used 10x chemistry:

<div align="center">

| 10X VERSION | BC | UMILEN | STR |
|:-:|:-:|:-:|:-:|
| 3' v1 | 737K-april-2014_rc.txt | 10 | Forward |
| 3' v2 | 737K-august-2016.txt | 10 | Forward |
| 3' v3, v3.1 | 3M-february-2018.txt | 12 | Forward |
| 3' v4 | 3M-3pgex-may-2023.txt | 12 | Forward |
| 5' v1.1, v2 | 737K-august-2016.txt | 10 | Reverse |
| 5' v3 | 737K-august-2016.txt | 12 | Reverse |
| 5' v4 | 3M-5pgex-jan-2023.txt | 12 | Reverse |
| multiome | 737K-arc-v1.txt | 12 | Forward |

</div>

  - `--soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30` are options that define UMI collapsing, barcode collapsing, and read clipping algorithms that are closest to ones used by `Cell Ranger` v4+; 
  - `--soloCellFilter EmptyDrops_CR` specifies the cell filtering algorithm used in [EmptyDrops](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html), which is the default algorithm in later versions of `Cell Ranger`; 
  - `--soloFeatures Gene GeneFull Velocyto` output conventional (exon-only) UMI counts, as well as exon+intron UMI counts (analog of `Cell Ranger` premrna option), as well as matrices preprocessed for `Velocyto`; 
  - `--soloMultiMappers EM` is to count multimappers (on by default in v3.0+ of our wrapper scripts; does not influence the main output, but creates an additional matrix in `/raw` subdir of `Gene` and `GeneFull`); 
  - `--readFilesCommand zcat` is used if your input fastq files are gzipped (auto-deteted);
  - `--outReadsUnmapped Fastx` to output the unmapped reads that can be used for metagenomic analysis; 
  - options grouped as `$SORTEDBAM` should be used if you need a genomic bam file; otherwise, use `$NOBAM`.

Actual `STAR` command being run (note the order of read files passed to `--readFilesIn`): 

```bash
STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn $R2 $R1 --runDirPerm All_RWX $GZIP $BAM \
   --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloBarcodeReadLength 0 \
   --soloCBlen $CBLEN --soloUMIstart $((CBLEN+1)) --soloUMIlen $UMILEN --soloStrand $STRAND \
   --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
   --soloCellFilter EmptyDrops_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
   --soloFeatures Gene GeneFull Velocyto --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx \
   --soloMultiMappers EM --outReadsUnmapped Fastx
```
### Processing 10X datasets in paired-end mode

Although some 3' 10X datasets are sequenced with a "long" R1, usually they are aligned as single-end, since it would take a very long R1 to accommodate for barcode, UMI, TSO + adapter, and a poly-A tail. On the other hand, 5' 10X datasets, especially ones sequenced with 2x150 bp paired-end reads, could be trimmed and aligned as paired-end (this is what `Cell Ranger` does, too). Effective `STARsolo` command for paired-end processing is as follows (note the 5' clipping flag, read input order, and flipped `soloStrand` because of the changed read order):

```bash
STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn $R1 $R2 --runDirPerm All_RWX $GZIP $BAM \
   --soloBarcodeMate 1 --clip5pNbases 39 0 \
   --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloCBstart 1 --soloCBlen $CBLEN \
   --soloUMIstart $((CBLEN+1)) --soloUMIlen $UMILEN --soloStrand Forward \
   --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
   --soloCellFilter EmptyDrops_CR --outFilterScoreMin 30 \
   --soloFeatures Gene GeneFull Velocyto --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx \
   --soloMultiMappers EM --outReadsUnmapped Fastx
```

### Using STARsolo for SMART-seq/SMART-seq2

For plate-based methods that don't use UMIs (such as [SMART-Seq and SMART-Seq2](https://teichlab.github.io/scg_lib_structs/methods_html/SMART-seq_family.html)), `STARsolo` can be used as well. Fastq files for these methods usually come as separate, paired-end files; all of these should be listed in a *manifest* file - plain text, tab-separated file containing three columns per line: 1) full path to R1; 2) full path to R2; 3) cell name or ID. 

Example of a script used to process Smart-seq2 data can be found in `/scripts/starsolo_ss2.sh`. Key parameters that could be adjusted are `--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3`; the higher they are, the less permissive is the alignment (we use default settings usually). Lower values can help you "rescue" a larger proportion of reads with high adapter content (see below for adapter trimming). Actual `STAR` command being run:

```bash
STAR --runThreadN $CPUS --genomeDir $REF --runDirPerm All_RWX $GZIP $BAM \
     --soloType SmartSeq --readFilesManifest ../$TAG.manifest.tsv \
     --soloUMIdedup Exact --soloStrand Unstranded \
     --limitOutSJcollapsed 10000000 --soloCellFilter None \
     --soloFeatures Gene GeneFull --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx
```

Often, SMART-seq2 reads can benefit from trimming adapters, which can be turned on using `--clip3pAdapterSeq <3' adapter sequence>` option. Alternatively, `/scripts/bbduk.sh` can be used to trim adapters from reads prior to the alignment and quantification (this is the preferred option).  

### Counting the multimapping reads

Default approach used by `Cell Ranger` (and `STARsolo` scripts above) is to discard all reads that map to multiple genomic locations with equal mapping quality. This approach creates a bias in gene expression estimation. Pseudocount-based methods correctly quantify multimapping reads, but generate false counts due to pseudo-alignment errors. These issues are described in good detail [here](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1). 

If you would like to process multimappers, add the following options: `--soloMultiMappers Uniform EM` (or simply `--soloMultiMappers EM`; this is on by default in v3.0+ of these scripts). This will generate an extra matrix in the `/raw` output folders. There will be non-integer numbers in the matrix because of split reads. If the downstream processing requires integers, you can round with a tool of your liking (e.g. `awk`). 

As of `STAR` v2.7.10a, **multimapper counting still does not work for SMART-seq2 or bulk RNA-seq processing**. 

### Running STARsolo on other scRNA-seq platforms 

STARsolo is very flexible and can be used with almost any scRNA-seq method, provided you know the library structure - i.e. where cell barcodes, UMIs, and biological parts of the read are located in the sequencing fragment or reads. A great source of information about scRNA-seq library structures is [this page](https://teichlab.github.io/scg_lib_structs/).

Currently, our scripts directory provides dedicated scripts for:
  - SMART-seq/SMART-seq2;
  - Drop-seq;
  - inDrops;
  - Microwell-seq;
  - BD Rhapsody;
  - CEL-seq;
  - STRT-seq. 

Please contact `CellGenIT` if you need to process an unusual dataset. 

## Quick evaluation of multiple STARsolo runs

If you've used these scripts to process multiple 10X samples, you can get a quick look at the results by copying `solo_QC.sh` script from this repo to the directory with `STARsolo` output folders, and running

```bash
./solo_QC.sh | column -t 
```

The script is designed for 10X or other droplet-based methods. For SMART-seq2, if processing was done per plate, it will also produce an informative output. 
