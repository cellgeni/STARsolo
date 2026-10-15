# AGENTS.md

Notes for anyone (human or coding agent) changing this repository. User
documentation is in [README.md](README.md).

## What this is

`starsolo` is a bash CLI around [STAR](https://github.com/alexdobin/STAR) in
STARsolo mode. It picks the barcode/UMI layout for each scRNA-seq platform,
detects 10x chemistry and strand, runs STAR, and post-processes the output.
There is no build step: `bin/starsolo` sources everything from `lib/`.

| Path | Role |
|:--|:--|
| `bin/starsolo` | Entry point: option parsing (`parse_droplet_args`), dispatch, version |
| `lib/common.sh` | Shared helpers: FASTQ discovery, STAR argument handling, `run_star`, `process_output_files` |
| `lib/platform_<name>.sh` | One file per platform: `run_<name>` and `show_<name>_help` |
| `lib/platform_10x.sh` | Chemistry/strand detection, paired-end routing, Velocyto policy and crash recovery |
| `lib/qc.sh` | `starsolo qc` report over finished runs |
| `etc/defaults.conf` | Defaults, overridable by `STARSOLO_*` environment variables |
| `tests/` | `test_lib.sh`, `test_recovery.sh`, `test_platforms.sh` (see below) |
| `scripts/` | Test-data generator/builder, LSF and bbduk helpers, `legacy/` original scripts |
| `data/` | `whitelists.tar.gz` (shipped in the image), `test.tar.gz` (generated test data) |

## Rules that are easy to break

- **Never call `STAR` directly from a platform.** Build the arguments with
  `merge_star_args <built-in args…>` and run them with `run_star`.
  STAR 2.7.10a rejects a parameter given twice, so `merge_star_args` drops
  every built-in flag that the user's `-- …` passthrough sets again; `run_star`
  cleans up `_STARtmp` when STAR fails. Pass `$GZIP $BAM` unquoted (they
  hold several words).
- **Options that only make sense for one platform are rejected elsewhere**
  in `parse_droplet_args` (see `--wl`, `--skip-length-checks`,
  `--no-velocyto`), and every option that takes a value calls `_need_value`.
- **Flags starsolo owns cannot be passed after `--`** (`--runThreadN`,
  `--genomeDir`, `--readFilesIn`, `--readFilesManifest`, `--outFileNamePrefix`,
  `--soloOutFileNames`): recovery and `starsolo qc` rely on where outputs are.
- **Velocyto (10x):** `_10x_velocyto_policy` decides whether Velocyto stays,
  `_10x_align` recovers from a crash in it, `_10x_check_read_index_limit`
  handles runs with >= 2^32 reads. A run that loses its Velocyto output writes
  `velocyto_failed.txt`, which `starsolo qc` reports. Background: issue #27,
  alexdobin/STAR#2710.
- **STAR version.** The image uses STAR 2.7.10a_alpha_220818 so that new
  batches match earlier ones. The Dockerfile builds it from `STAR_REPO` /
  `STAR_REF` (pinned to a commit). The cellgeni build
  `2.7.10a_alpha_220818_velofix` (branch
  [`2.7.10a_alpha_220818-velocyto`](https://github.com/cellgeni/STAR/tree/2.7.10a_alpha_220818-velocyto)
  of `cellgeni/STAR`) adds only the two
  Velocyto fixes from #27; on single-end data it was verified byte-identical
  to stock. The wrapper detects it (`STAR --version` contains `velofix`,
  `Log.out` has the 64-bit read-index line) and still copes with stock STAR.
  Do not upgrade STAR casually: compare outputs on real data first.

## Bash pitfalls (everything runs under `set -euo pipefail`)

- `x=$(f)` runs `f` in a subshell: a function that sets globals
  (`merge_star_args`, `_10x_velocyto_policy`) must be called as a statement.
  `die` inside `$(…)` only ends the subshell; the failed assignment then stops
  the script.
- `local x=$(f)` hides `f`'s exit status. Declare first: `local x; x=$(f)`.
- A function whose last command is `[[ cond ]] && action` returns 1 when the
  condition is false, and stops the caller. End with `if … fi` or `return 0`.
- `grep pattern file | awk …` fails the pipeline when nothing matches. Use
  awk alone or the `_qc_csv`/`_qc_stat` helpers.
- Inside `if`, `&&`, `||` or `$(…)` in a condition, `set -e` is off, so a
  test that calls a function that way can pass while the CLI would fail.
  Test through `bin/starsolo` when the exit status matters.
- A bare `wait` returns 0 whatever the background jobs did; wait on PIDs.
- `pbzip2` takes `-p<N>` with no space; `-p N` aborts.
- The image's awk is mawk: it prints large numbers in exponent form, so
  read counts use bash arithmetic.

## Tests

```bash
bash tests/test_lib.sh                  # library functions; bash + pbzip2 only
bash tests/test_recovery.sh             # 10x crash recovery with a stub STAR; needs STAR
bash tests/test_platforms.sh [case …]   # every platform on data/test.tar.gz vs data/test/truth
```

`test_recovery.sh` and `test_platforms.sh` need STAR, seqtk, samtools and
pbzip2: run them inside the image. They use `$STARSOLO_TEST_DATA`, else
`data/test`, else unpack `data/test.tar.gz`; a local `data/test` from an older
tarball is used as is, so re-extract it after the tarball changes.
`KEEP_WORK=1` keeps their temp directories.

On the Sanger farm (no STAR on the head node), run them on a compute node:

```bash
IMG=/nfs/cellgeni/singularity/images/quay.io-cellgeni-starsolo-v4.3.img
bsub -q normal -n 4 -M 16000 -R "select[mem>16000] rusage[mem=16000] span[hosts=1]" -K \
    singularity exec -B /lustre,/nfs "$IMG" bash tests/test_platforms.sh
```

Inside the image `/opt/STARsolo/bin` comes first on `PATH`: call this
checkout's `bin/starsolo` by path (the tests already do).

CI (`.github/workflows/test.yml`) runs `test_platforms.sh` per case, plus
`test_lib.sh`, `test_recovery.sh` and a 10x run with BAM output, in
`quay.io/cellgeni/starsolo:base`.

## Test data

`data/test.tar.gz` is generated, never edited by hand. Change
`scripts/generate_test_data.py` (standard library only, Python 3.8,
deterministic for a seed) and rebuild inside the image:

```bash
bash scripts/build_test_data.sh          # generator + STAR index + test_platforms.sh + tarball
```

The genome is tiny, so the index is built with `--genomeSAindexNbases 4`.
Each droplet sample needs at least 31 real cells, otherwise STAR's cell
filter keeps every barcode; `data/test/truth/` lists the expected cells.

## Releases

- Bump `STARSOLO_VERSION` (and the header comment) in `bin/starsolo`.
- Pushing a tag builds and pushes `quay.io/cellgeni/starsolo:<tag>` (and
  `latest` for tags without a hyphen).
- Every push to `main` opens an issue in `cellgeni/nf-starsolo`.
- Work happens on `dev`; `main` gets it through pull requests.

## Conventions

- Match the existing style: `local` variables, `UPPERCASE` for values passed
  between steps, `log_info`/`log_warn`/`die` for messages, one comment line
  above each helper describing its arguments and output.
- Keep `README.md` and the `--help` texts in step with behaviour changes.
- One logical change per commit, with a message that says why.
