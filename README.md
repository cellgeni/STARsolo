# `starsolo` — unified CLI for STARsolo scRNA-seq processing

A single command-line tool that wraps [STAR](https://github.com/alexdobin/STAR) in `STARsolo` mode for uniform processing of scRNA-seq data across multiple platforms.

## Supported platforms

| Subcommand | Platform | Barcode type | Chemistry auto-detection |
|:--|:--|:--|:-:|
| `10x` | 10x Genomics (v1 – v4, multiome) | CB_UMI_Simple | ✅ |
| `smartseq` | Smart-seq / Smart-seq2 | SmartSeq (plate-based, no UMIs) | — |
| `dropseq` | Drop-seq | CB_UMI_Simple (no whitelist) | — |
| `rhapsody` | BD Rhapsody | CB_UMI_Complex (3 segments) | — |
| `indrops` | inDrops | CB_UMI_Complex (2 segments + adapter) | — |
| `strt` | STRT-seq | CB_UMI_Simple (96-barcode list) | — |
| `qc` | Aggregate QC stats across runs | — | — |

## Installation

### Prerequisites

The following tools **must be available in `$PATH`** before running `starsolo`:

| Tool | Tested version | Purpose |
|:--|:--|:--|
| [STAR](https://github.com/alexdobin/STAR) | 2.7.10a | Alignment & quantification |
| [samtools](https://github.com/samtools/samtools) | 1.15.1 | BAM indexing |
| [seqtk](https://github.com/lh3/seqtk) | 1.3+ | Read subsampling (10x chemistry detection) |
| [pbzip2](https://launchpad.net/pbzip2) | 1.1+ | Compressing unmapped reads |
| [BBMap](https://sourceforge.net/projects/bbmap/) | 38.97 | Adapter trimming (optional, via `bbduk.sh`) |

> **Tip:** The included [Dockerfile](Dockerfile) builds all of the above into a single container image.

STAR must be compiled with `make STAR CXXFLAGS_SIMD="-msse4.2"` to support `--clipAdapterType CellRanger4` (see [STAR#1218](https://github.com/alexdobin/STAR/issues/1218)).

### Install `starsolo` CLI

```bash
git clone <this-repo> && cd STARsolo

# Default: symlinks into ~/.local/bin
./install.sh

# Or specify a directory:
./install.sh /usr/local/bin

# Verify
starsolo --version
```

The installer creates a single `starsolo` symlink. All library code is resolved relative to the repo, so you can `git pull` to update in-place.

## Reference genome

Use a standard STAR genome index. [Cell Ranger-filtered](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#header) references are recommended for comparability with Cell Ranger output.

```bash
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir STAR --genomeFastaFiles $FA --sjdbGTFfile $GTF
```

By default, `starsolo` resolves references via `--species`:

```
$STARSOLO_REF_BASE/<species>/2020A/index
```

Override with `--ref /your/path/to/index` on any subcommand, or change `STARSOLO_REF_BASE` in [etc/defaults.conf](etc/defaults.conf).

### Barcode whitelists

Download 10x barcode whitelists [from Cell Ranger](https://github.com/10XGenomics/cellranger/tree/master/lib/python/cellranger/barcodes) and place them in the whitelist directory (default: `/nfs/cellgeni/STAR/whitelists`, configurable via `--whitelist-dir` or `STARSOLO_WL_DIR`).

## Usage

```
starsolo <platform> <args…> [options]
```

### Global options

| Flag | Description | Default |
|:--|:--|:--|
| `-s, --species <name>` | Species name (`human`, `mouse`, …) — resolves reference automatically | — |
| `-r, --ref <path>` | Explicit STAR index path (overrides `--species`) | — |
| `-w, --whitelist-dir <dir>` | Barcode whitelist directory | from config |
| `-c, --cpus <N>` | Number of threads | `16` |
| `--bam` | Output coordinate-sorted BAM | off |
| `--no-bam` | Suppress BAM output (default) | ✅ |
| `-h, --help` | Show help (global or per-platform) | — |
| `--version` | Print version | — |

### 10x Genomics

Automatically detects chemistry (v1–v4, multiome), strand specificity, and paired-end mode.

```bash
starsolo 10x /data/fastqs SAMPLE1 --species human
starsolo 10x /data/fastqs SAMPLE1 --ref /path/to/index --bam
```

#### 10x chemistry detection

The script subsamples 200,000 reads and matches barcodes against all known whitelists:

| 10x version | Whitelist | CB len | UMI len | Strand |
|:-:|:-:|:-:|:-:|:-:|
| 3' v1 | 737K-april-2014_rc.txt | 14 | 10 | Forward |
| 3' v2 | 737K-august-2016.txt | 16 | 10 | Forward |
| 3' v3/v3.1 | 3M-february-2018.txt | 16 | 12 | Forward |
| 3' v4 | 3M-3pgex-may-2023.txt | 16 | 12 | Forward |
| 5' v1.1/v2 | 737K-august-2016.txt | 16 | 10 | Reverse |
| 5' v3 | 737K-august-2016.txt | 16 | 12 | Reverse |
| 5' v4 | 3M-5pgex-jan-2023.txt | 16 | 12 | Reverse |
| multiome | 737K-arc-v1.txt | 16 | 12 | Forward |

#### Paired-end processing

5' libraries with long R1 reads (>50 bp) are automatically processed in paired-end mode with `--soloBarcodeMate 1 --clip5pNbases 39 0`. 3' libraries are always single-end.

### Smart-seq / Smart-seq2

Requires a **manifest file** — a tab-separated file with three columns: `R1_path`, `R2_path`, `cell_name`.

```bash
starsolo smartseq manifest.tsv --species mouse
```

Aliases: `smart-seq`, `ss2`

### Drop-seq

12 bp cell barcode, 8 bp UMI, no whitelist.

```bash
starsolo dropseq /data/fastqs SAMPLE1 --species human
```

Alias: `drop-seq`

### BD Rhapsody

3 barcode segments + 8 bp UMI. Requires `Rhapsody_bc1/2/3.txt` in the whitelist directory.

```bash
starsolo rhapsody /data/fastqs SAMPLE1 --species human
```

Aliases: `bd-rhapsody`, `bd`

### inDrops

2 barcode segments + adapter + 6 bp UMI. Requires `inDrops_Ambrose2_bc1/2.txt` in the whitelist directory.

```bash
starsolo indrops /data/fastqs SAMPLE1 --species human
starsolo indrops /data/fastqs SAMPLE1 --species human --adapter GAGTGATTGCTTGTGACGCCAA
```

Default adapter: `GAGTGATTGCTTGTGACGCCTT`

### STRT-seq

96-barcode whitelist, 8 bp cell barcode, 8 bp UMI. Requires `96_barcodes.list` in the whitelist directory.

```bash
starsolo strt /data/fastqs SAMPLE1 --species human --strand Forward
```

Extra option: `--strand <Forward|Reverse>` (default: `Forward`)

Alias: `strt-seq`

### QC aggregation

Run from the parent directory containing per-sample output folders:

```bash
starsolo qc | column -t
```

Checks for:
- Incomplete runs (`_STARtmp` still present)
- Barcode cross-contamination between samples
- Low mapping percentages

Outputs a tab-separated table with read counts, mapping rates, cell counts, and configuration for each sample.

## Configuration

Edit [etc/defaults.conf](etc/defaults.conf) to change default paths and parameters. All values can also be set as environment variables:

```bash
export STARSOLO_REF_BASE=/my/references
export STARSOLO_WL_DIR=/my/whitelists
export STARSOLO_CPUS=8
export STARSOLO_BAM_MODE=bam
```

CLI flags always take precedence over environment variables, which take precedence over config file defaults.

## Helper scripts

| Script | Purpose |
|:--|:--|
| [scripts/bbduk_trim.sh](scripts/bbduk_trim.sh) | Adapter/polyA trimming via BBMap |
| [scripts/bsub_submit.sh](scripts/bsub_submit.sh) | Submit any starsolo command as an LSF job |

### LSF submission

```bash
./scripts/bsub_submit.sh starsolo 10x /data/fastqs SAMPLE1 --species human
```

This submits with 16 CPUs / 64 GB RAM to the `normal` queue.

## Resource requirements

Default: **16 CPUs, 64–128 GB RAM**. Some datasets require 128 GB; if jobs fail with OOM errors, increase RAM in your job submission.

## Docker

The included [Dockerfile](Dockerfile) builds all required tools:

```bash
docker build -t starsolo .
docker run -v /data:/data starsolo starsolo 10x /data/fastqs SAMPLE1 --species human
```

## Project structure

```
STARsolo/
├── bin/
│   └── starsolo              # CLI entrypoint (add to PATH)
├── lib/
│   ├── common.sh             # Shared functions (FASTQ discovery, compression, post-processing)
│   ├── platform_10x.sh       # 10x Genomics logic
│   ├── platform_smartseq.sh  # Smart-seq2 logic
│   ├── platform_dropseq.sh   # Drop-seq logic
│   ├── platform_rhapsody.sh  # BD Rhapsody logic
│   ├── platform_indrops.sh   # inDrops logic
│   └── platform_strt.sh      # STRT-seq logic
├── etc/
│   └── defaults.conf          # Default configuration
├── scripts/
│   ├── solo_qc.sh             # QC aggregation
│   ├── bbduk_trim.sh          # Adapter trimming helper
│   └── bsub_submit.sh         # LSF job submission helper
├── Dockerfile
├── install.sh
├── LICENSE
└── README.md
```

## License

See [LICENSE](LICENSE).
