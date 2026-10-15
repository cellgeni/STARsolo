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

The following tools **must be available in `$PATH`** before running `starsolo` (which itself needs bash ≥ 4.4):

| Tool | Tested version | Purpose |
|:--|:--|:--|
| [STAR](https://github.com/alexdobin/STAR) | 2.7.10a_alpha_220818 (image: `_velofix` build, see [Velocyto](#velocyto)) | Alignment & quantification |
| [samtools](https://github.com/samtools/samtools) | 1.21 | BAM indexing |
| [seqtk](https://github.com/lh3/seqtk) | 1.4 | Read subsampling (10x chemistry detection) |
| [pbzip2](https://launchpad.net/pbzip2) | 1.1.13 | Compressing unmapped reads |
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
| `--wl <chemistry>` | **10x only** — use this chemistry instead of detecting it | auto-detect |
| `--skip-length-checks` | **10x only** — warn instead of aborting on unexpected read lengths | off |
| `--no-velocyto` | **10x only** — do not count Velocyto matrices (see [Velocyto](#velocyto)) | off |
| `-c, --cpus <N>` | Number of threads | `16` |
| `--bam` | Output coordinate-sorted BAM | off |
| `--no-bam` | Suppress BAM output (default) | ✅ |
| `-h, --help` | Show help (global or per-platform) | — |
| `--version` | Print version | — |
| `--` | Pass all following arguments to `STAR`, replacing built-in flags of the same name | — |

### Passing extra options to STAR

Any arguments after a `--` sentinel are passed to the underlying `STAR`
command. This works on every platform and lets you override defaults or set
options the CLI doesn't expose:

```bash
starsolo 10x /data/fastqs SAMPLE1 --species human -- --outFilterScoreMin 20 --soloFeatures Gene GeneFull GeneFull_Ex50pAS Velocyto
```

The `--` must come last, after all `starsolo` options.

**Overriding built-in flags.** STAR rejects a parameter given twice, so
`starsolo` removes each built-in flag (with its values) that the passthrough
sets again, and appends the passthrough. `-- --outFilterScoreMin 20` replaces
the default `--outFilterScoreMin 30`. For multi-value parameters this is a full
replacement, not a merge: `-- --soloFeatures Gene GeneFull SJ` replaces the
built-in feature list entirely. For example, to add `GeneFull_Ex50pAS` (STAR's
intronic counting closest to Cell Ranger's "include introns") to the 10x
defaults, pass `-- --soloFeatures Gene GeneFull GeneFull_Ex50pAS Velocyto`.

A few flags are set by `starsolo` and cannot be passed after `--`:
`--runThreadN` (use `--cpus`), `--genomeDir` (use `--ref`/`--species`),
`--readFilesIn`, `--readFilesManifest`, `--outFileNamePrefix` and
`--soloOutFileNames`.

### 10x Genomics

Automatically detects chemistry (v1–v4, multiome), strand specificity, and paired-end mode. Chemistry detection can be bypassed with [`--wl`](#skipping-chemistry-detection---wl) when it is already known.

```bash
starsolo 10x /data/fastqs SAMPLE1 --species human
starsolo 10x /data/fastqs SAMPLE1 --ref /path/to/index --bam
starsolo 10x /data/fastqs SAMPLE1 --species human --wl gex_3pv3_family
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
| 5' v3 (GEM-X) | 3M-5pgex-jan-2023.txt | 16 | 12 | Reverse |
| multiome | 737K-arc-v1.txt | 16 | 12 | Forward |

#### Skipping chemistry detection (`--wl`)

If the chemistry is already known, `--wl <chemistry>` selects the whitelist and
CB/UMI geometry directly and skips the barcode-matching step. Names match the
chemistry ids emitted by `infer_10x_run.py`, so a value taken from that script's
JSON report can be passed straight through:

| `--wl` value | Alias | Whitelist | CB len | UMI len |
|:--|:-:|:--|:-:|:-:|
| `gex_3pv1` | `v1` | 737K-april-2014_rc.txt | 14 | 10 |
| `gex_3pv2_or_5pv1v2` | `v2` | 737K-august-2016.txt | 16 | 10 |
| `gex_3pv3_family` | `v3` | 3M-february-2018.txt | 16 | 12 |
| `gex_3pv4_gemx` | `v4` | 3M-3pgex-may-2023.txt | 16 | 12 |
| `gex_5pv3_gemx` | `5pv3` | 3M-5pgex-jan-2023.txt | 16 | 12 |
| `gex_multiome_arc_v1` | `multiome`, `arc` | gex_737K-arc-v1.txt | 16 | 12 |

```bash
starsolo 10x /data/fastqs SAMPLE1 --species human --wl gex_3pv3_family
starsolo 10x /data/fastqs SAMPLE1 --species human --wl v3
```

Values are case-insensitive. The multiome entry falls back to `737K-arc-v1.txt`
if the `gex_`-prefixed copy is absent. ATAC chemistry ids are rejected —
`starsolo 10x` handles gene expression only.

Only the whitelist matching is skipped: reads are still subsampled, and strand
specificity is still determined by test alignment (it cannot be inferred from
the whitelist). The run's `strand.txt` reports `NA` for the per-whitelist match
counts in this mode.

#### Read-length validation (`--skip-length-checks`)

Before alignment the 10x workflow aborts on read lengths that indicate a problem
with the input: a variable-length R1, an R1 under 24 bp, or an R2 under 40 bp.
`--skip-length-checks` reports these as warnings instead and continues:

```bash
starsolo 10x /data/fastqs SAMPLE1 --species human --wl v3 --skip-length-checks
```

CB/UMI reconciliation still applies — if R1 is too short to hold the full
barcode+UMI, the UMI length is trimmed to fit and a warning is logged. An R1
shorter than the cell barcode itself remains a fatal error, since it leaves no
UMI at all.

#### Paired-end processing

5' libraries with long R1 reads (>50 bp) are automatically processed in paired-end mode with `--soloBarcodeMate 1 --clip5pNbases <CB+UMI+13> 0` (39 for 5' v1/v2, 41 for GEM-X 5' with its 12 bp UMI). 3' libraries are always single-end.

Only the two whitelists used by 5' chemistries (`737K-august-2016.txt` and
`3M-5pgex-jan-2023.txt`) can be processed as 5'. With any other (3'-only)
whitelist the strand is always `Forward`: if the strand test favours Reverse
(this happens on some single-nucleus 3' data with a long R1), a warning with
both percentages is logged and the run uses Forward, single-end. `strand.txt`
notes the override.

#### Velocyto

`starsolo 10x` counts `Gene GeneFull Velocyto` by default. Stock STAR
(2.7.10a to 2.7.11b) breaks in Velocyto in two cases
([#27](https://github.com/cellgeni/STARsolo/issues/27)):

1. **Paired-end reads** whose mate lies outside the transcript make
   Velocyto's transcript check read out of bounds. With the 3M GEM-X 5'
   whitelist this segfaults while mapping; with 737K it reads garbage for a
   handful of reads.
2. **2³² (about 4.3 × 10⁹) or more input reads**: the read index STAR keeps
   for Velocyto wraps, so STAR segfaults in Velocyto counting after Gene and
   GeneFull are written, or finishes with wrong Velocyto counts and BAM
   `CB`/`UB` tags ([alexdobin/STAR#2710](https://github.com/alexdobin/STAR/issues/2710)).

The Docker image therefore uses **STAR `2.7.10a_alpha_220818_velofix`**: the
same 2.7.10a_alpha_220818 with only these two fixes (the second is a backport
of [alexdobin/STAR#2711](https://github.com/alexdobin/STAR/pull/2711)). On
single-end 10x 3' data its output is byte-identical to stock 2.7.10a
(matrices, `Log.final.out`, BAM), so new batches stay comparable with earlier
ones; paired-end Velocyto changes only for the reads that were read from
garbage. The Dockerfile's `STAR_REPO`/`STAR_REF` select the source.

The wrapper also copes with stock STAR, so that a Velocyto failure costs the
Velocyto output, not the sample:

- `--no-velocyto` drops Velocyto from the feature list.
- With a STAR that is not a `velofix` build, paired-end GEM-X 5' runs
  (`3M-5pgex-jan-2023.txt`) drop Velocyto automatically.
- If STAR crashes in Velocyto counting after the other features are written,
  those features are kept and compressed, the partial `output/Velocyto` is
  removed, and the run exits 0. If a BAM was requested, it is missing (STAR
  sorts it after Solo counting), so the run exits with STAR's code instead.
- If STAR crashes while mapping with Velocyto on, the run is repeated once
  without Velocyto. The first attempt's log is kept as `Log.attempt1.out`.
- With 2³² or more input reads, `output/Velocyto` is removed after the run
  and a warning about BAM `CB`/`UB` tags is logged, unless `Log.out` shows
  that STAR stored the read index in 64 bits (the `velofix` build does).

Each of these is recorded in `strand.txt` (`Velocyto:` line), and whenever
Velocyto output is missing for a reason other than `--no-velocyto`, a
`velocyto_failed.txt` marker (`key=value` lines: `star_exit_code`, `stage`,
`action`, …) is written next to it and reported by `starsolo qc`. After any
STAR failure `_STARtmp` is deleted, and the genome that the strand test loaded
into shared memory is released.

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

```bash
starsolo qc <dir1> [dir2 … dirN] [options]
```

Pass one or more STARsolo output directories as positional arguments:

```bash
starsolo qc sample1/ sample2/ sample3/
starsolo qc results/GSM* | column -t
starsolo qc sample1/ sample2/ --species human
starsolo qc sample1/ sample2/ --no-contamination-check
```

#### QC options

| Flag | Description | Default |
|:--|:--|:--|
| `-s, --species <name>` | Override species label in output (e.g. `human`, `mouse`) | guessed from `genomeDir` in `Log.out` |
| `--no-contamination-check` | Skip cross-sample barcode overlap check | off (check runs by default) |
| `-h, --help` | Show help | — |

#### Checks performed

- **Completion** — warns if `_STARtmp` is still present or `Log.final.out` is missing (run did not finish)
- **Velocyto** — warns if `starsolo 10x` dropped the Velocyto output (`velocyto_failed.txt`)
- **Cross-contamination** — compares filtered barcodes between every pair of samples; warns if overlap exceeds 20%
- **Low mapping** — warns if `GeneFull` unique mapping rate is ≤ 20%

#### Output columns

Outputs a tab-separated table to stdout (pipe through `column -t` for aligned display):

| Column | Description |
|:--|:--|
| `Sample` | Directory name |
| `Rd_all` | Total reads |
| `Rd_in_cells` | Reads in cells (GeneFull) |
| `Frc_in_cells` | Fraction of reads in cells |
| `UMI_in_cells` | UMIs in cells |
| `Cells` | Estimated number of cells |
| `Med_nFeature` | Median genes per cell |
| `Good_BC` | Fraction of reads with valid barcodes |
| `WL` | Detected whitelist (`v1`, `v2`, `v3`, `v4-3p` = GEM-X 3' v4, `v4-5p` = GEM-X 5' v3, `arc`) |
| `Species` | Species label |
| `Paired` | `Single` or `Paired` |
| `Strand` | Strand specificity |
| `all_u+m` / `all_u` | Genome mapping rate (unique+multi / unique) |
| `exon_u+m` / `exon_u` | Exonic (Gene) mapping rate |
| `full_u+m` / `full_u` | Full-gene (GeneFull) mapping rate |

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
| [scripts/bsub_submit.sh](scripts/bsub_submit.sh) | Submit any `starsolo` command as an LSF job |

### LSF submission

```bash
./scripts/bsub_submit.sh starsolo 10x /data/fastqs SAMPLE1 --species human
```

This submits with 16 CPUs / 64 GB RAM to the `normal` queue.

## Resource requirements

Default: **16 CPUs, 64–128 GB RAM**. Some datasets require 128 GB; if jobs fail with OOM errors, increase RAM in your job submission.

## Testing

```bash
bash tests/test_lib.sh                  # library functions; needs only bash and pbzip2
bash tests/test_recovery.sh             # 10x crash recovery; needs STAR (run in the image)
bash tests/test_platforms.sh [case …]   # every platform on data/test.tar.gz; needs STAR
```

`data/test.tar.gz` is generated, not edited: change `scripts/generate_test_data.py`
and run `scripts/build_test_data.sh` inside the image.

## Docker

The included [Dockerfile](Dockerfile) builds all required tools:

```bash
docker build -t starsolo .
docker run -v /data:/data starsolo starsolo 10x /data/fastqs SAMPLE1 --species human
```

STAR is compiled in a separate build stage from the GitHub source selected by
`STAR_REPO` and `STAR_REF` (a tag, branch or commit), and only the binary is
copied into the image; `/versions.txt` records which source was used:

```bash
docker build --build-arg STAR_REPO=alexdobin/STAR --build-arg STAR_REF=2.7.10a_alpha_220818 -t starsolo .
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
│   ├── platform_strt.sh      # STRT-seq logic
│   └── qc.sh                 # QC aggregation logic
├── etc/
│   └── defaults.conf          # Default configuration
├── tests/
│   ├── test_lib.sh            # Unit tests of the library functions (no STAR needed)
│   ├── test_recovery.sh       # 10x Velocyto/crash-recovery tests (stub STAR)
│   └── test_platforms.sh      # Every platform on the test data, checked against data/test/truth
├── scripts/
│   ├── bbduk_trim.sh          # Adapter trimming helper
│   ├── bsub_submit.sh         # LSF job submission helper
│   ├── generate_test_data.py  # Synthetic test reads, reference and whitelists
│   ├── build_test_data.sh     # Rebuilds data/test.tar.gz (generator + STAR index)
│   └── legacy/                # The original per-platform scripts
├── data/
│   ├── test.tar.gz            # Synthetic test dataset (see data/test/README.md inside)
│   └── whitelists.tar.gz      # Barcode whitelists (unpacked into the Docker image)
├── Dockerfile
├── install.sh
├── AGENTS.md                  # Notes for contributors and coding agents
├── LICENSE
└── README.md
```

## License

See [LICENSE](LICENSE).
