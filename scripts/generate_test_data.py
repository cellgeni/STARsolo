#!/usr/bin/env python3
"""Generate a tiny synthetic dataset for the STARsolo CLI wrapper.

- Creates a minimal synthetic reference (FASTA + GTF).
- Creates minimal whitelist files with the exact filenames expected by the wrapper.
- Creates gzipped FASTQs for each supported technology:
    10x, dropseq, rhapsody, indrops, strt, smartseq.

This script uses ONLY the Python standard library.
"""

from __future__ import annotations

import argparse
import gzip
import tarfile
from pathlib import Path
from typing import Iterable, Tuple

DNA_ALPHABET = "ACGT"


def base4_dna(n: int, length: int) -> str:
    """Deterministic pseudo-random-ish DNA string from an integer."""
    if length <= 0:
        return ""
    chars = []
    x = n
    for _ in range(length):
        chars.append(DNA_ALPHABET[x & 3])
        x >>= 2
        if x == 0:
            # mix a bit so we don't get long runs of 'A' for small n
            x = (n * 1103515245 + 12345) & 0xFFFFFFFF
    return "".join(chars)


def repeat_to_length(motif: str, length: int) -> str:
    if not motif:
        raise ValueError("motif must be non-empty")
    return (motif * ((length // len(motif)) + 1))[:length]


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_fastq_gz(path: Path, records: Iterable[Tuple[str, str]]) -> None:
    """Write gzipped FASTQ from (name, sequence) pairs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", compresslevel=6) as fh:
        for name, seq in records:
            qual = "I" * len(seq)
            fh.write(f"@{name}\n{seq}\n+\n{qual}\n")


def gen_records(prefix: str, mate: int, n: int, seq_fn) -> Iterable[Tuple[str, str]]:
    for i in range(n):
        yield (f"{prefix}_{i}/{mate}", seq_fn(i))


def make_reference(outdir: Path) -> None:
    """Create a synthetic reference with two genes on chrTest."""
    ref_dir = outdir / "reference"
    ref_dir.mkdir(parents=True, exist_ok=True)

    # Build a 1500bp contig. Embed two easy-to-hit regions for mapping.
    # geneA: 101-500 (1-based, inclusive) contains lots of ACGT repeats
    # geneB: 1001-1400 contains lots of TGCATGCA repeats
    contig_len = 1500
    seq = list(
        repeat_to_length(
            "CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTAC",
            contig_len,
        )
    )

    geneA = repeat_to_length("ACGT", 400)
    geneB = repeat_to_length("TGCATGCA", 400)

    # Insert (convert to 0-based slicing)
    seq[100:500] = list(geneA)  # 101..500
    seq[1000:1400] = list(geneB)  # 1001..1400

    fasta = ">chrTest\n" + "".join(seq) + "\n"
    write_text(ref_dir / "test.fa", fasta)

    gtf = "\n".join(
        [
            'chrTest\tsynth\tgene\t101\t500\t.\t+\t.\tgene_id "geneA"; gene_name "geneA";',
            'chrTest\tsynth\ttranscript\t101\t500\t.\t+\t.\tgene_id "geneA"; transcript_id "txA"; gene_name "geneA";',
            'chrTest\tsynth\texon\t101\t500\t.\t+\t.\tgene_id "geneA"; transcript_id "txA"; gene_name "geneA"; exon_number "1";',
            'chrTest\tsynth\tgene\t1001\t1400\t.\t+\t.\tgene_id "geneB"; gene_name "geneB";',
            'chrTest\tsynth\ttranscript\t1001\t1400\t.\t+\t.\tgene_id "geneB"; transcript_id "txB"; gene_name "geneB";',
            'chrTest\tsynth\texon\t1001\t1400\t.\t+\t.\tgene_id "geneB"; transcript_id "txB"; gene_name "geneB"; exon_number "1";',
            "",
        ]
    )
    write_text(ref_dir / "test.gtf", gtf)


def make_whitelists(outdir: Path) -> None:
    wl = outdir / "whitelists"
    wl.mkdir(parents=True, exist_ok=True)

    # 10x whitelists: minimal, but filenames MUST match the wrapper.
    write_text(wl / "3M-february-2018.txt", "ACGTTGCAACGTTGCA\n")
    write_text(wl / "737K-august-2016.txt", "TTTTTTTTTTTTTTTT\n")
    write_text(wl / "737K-arc-v1.txt", "TGCATGCATGCATGCA\n")
    write_text(wl / "737K-april-2014_rc.txt", "CCCCCCCCCCCCCC\n")  # 14bp
    write_text(wl / "3M-3pgex-may-2023.txt", "GGGGGGGGGGGGGGGG\n")
    write_text(wl / "3M-5pgex-jan-2023.txt", "AAAAAAAAAAAAAAAA\n")

    # BD Rhapsody
    write_text(wl / "Rhapsody_bc1.txt", "ACGTACGTA\n")
    write_text(wl / "Rhapsody_bc2.txt", "TGCATGCAT\n")
    write_text(wl / "Rhapsody_bc3.txt", "GATCGATCG\n")

    # inDrops
    write_text(wl / "inDrops_Ambrose2_bc1.txt", "ACGTACGTACG\n")
    write_text(wl / "inDrops_Ambrose2_bc2.txt", "TGCATGCA\n")

    # STRT-seq
    write_text(wl / "96_barcodes.list", "AACCGGTT\nTTGGAACC\n")


def make_fastqs(outdir: Path, tenx_reads: int) -> None:
    fq_root = outdir / "fastqs"

    # Shared cDNA read used across droplet tests (maps into geneA region)
    cdna_50 = repeat_to_length("ACGT", 50)
    cdna_75 = repeat_to_length("ACGT", 75)

    # 10x
    bc_10x = "ACGTTGCAACGTTGCA"  # must match 3M-february-2018.txt

    def tenx_r1(i: int) -> str:
        return bc_10x + base4_dna(i, 12)

    write_fastq_gz(
        fq_root / "10x" / "TENX1_1.fastq.gz",
        gen_records("10x", 1, tenx_reads, tenx_r1),
    )
    write_fastq_gz(
        fq_root / "10x" / "TENX1_2.fastq.gz",
        gen_records("10x", 2, tenx_reads, lambda i: cdna_50),
    )

    # Drop-seq (12bp CB + 8bp UMI)
    bc_drop = "AACCGGTTAACC"

    def drop_r1(i: int) -> str:
        return bc_drop + base4_dna(i, 8)

    write_fastq_gz(
        fq_root / "dropseq" / "DROP1_1.fastq.gz",
        gen_records("drop", 1, 2000, drop_r1),
    )
    write_fastq_gz(
        fq_root / "dropseq" / "DROP1_2.fastq.gz",
        gen_records("drop", 2, 2000, lambda i: cdna_50),
    )

    # BD Rhapsody (3 barcode segments + 8bp UMI)
    bc1, bc2, bc3 = "ACGTACGTA", "TGCATGCAT", "GATCGATCG"
    filler1, filler2 = "A" * 12, "A" * 13

    def rhap_r1(i: int) -> str:
        return bc1 + filler1 + bc2 + filler2 + bc3 + base4_dna(i, 8)

    write_fastq_gz(
        fq_root / "rhapsody" / "RHAP1_1.fastq.gz",
        gen_records("rhap", 1, 1000, rhap_r1),
    )
    write_fastq_gz(
        fq_root / "rhapsody" / "RHAP1_2.fastq.gz",
        gen_records("rhap", 2, 1000, lambda i: cdna_50),
    )

    # inDrops (BC1 + adapter + BC2 + 6bp UMI)
    bc1_i, bc2_i = "ACGTACGTACG", "TGCATGCA"
    adapter = "GAGTGATTGCTTGTGACGCCTT"  # wrapper default

    def ind_r1(i: int) -> str:
        return bc1_i + adapter + bc2_i + base4_dna(i, 6)

    write_fastq_gz(
        fq_root / "indrops" / "IND1_1.fastq.gz",
        gen_records("ind", 1, 1000, ind_r1),
    )
    write_fastq_gz(
        fq_root / "indrops" / "IND1_2.fastq.gz",
        gen_records("ind", 2, 1000, lambda i: cdna_50),
    )

    # STRT-seq (8bp CB + 8bp UMI) — in this wrapper, barcode is in *_1.fastq
    bc_strt = "AACCGGTT"  # one of 96_barcodes.list

    def strt_r1(i: int) -> str:
        return bc_strt + base4_dna(i, 8)

    write_fastq_gz(
        fq_root / "strt" / "STRT1_1.fastq.gz",
        gen_records("strt", 1, 500, strt_r1),
    )
    write_fastq_gz(
        fq_root / "strt" / "STRT1_2.fastq.gz",
        gen_records("strt", 2, 500, lambda i: cdna_50),
    )

    # Smart-seq (manifest with absolute paths)
    smart_dir = fq_root / "smartseq"
    smart_dir.mkdir(parents=True, exist_ok=True)

    for cell in ("CELL1", "CELL2"):
        write_fastq_gz(
            smart_dir / f"{cell}_R1.fastq.gz",
            gen_records(cell, 1, 200, lambda i: cdna_75),
        )
        write_fastq_gz(
            smart_dir / f"{cell}_R2.fastq.gz",
            gen_records(cell, 2, 200, lambda i: cdna_75),
        )

    manifest_path = smart_dir / "smartseq_example.manifest.tsv"
    cell1_r1 = (smart_dir / "CELL1_R1.fastq.gz").resolve()
    cell1_r2 = (smart_dir / "CELL1_R2.fastq.gz").resolve()
    cell2_r1 = (smart_dir / "CELL2_R1.fastq.gz").resolve()
    cell2_r2 = (smart_dir / "CELL2_R2.fastq.gz").resolve()

    manifest = "\n".join(
        [
            f"{cell1_r1}\t{cell1_r2}\tCELL1",
            f"{cell2_r1}\t{cell2_r2}\tCELL2",
            "",
        ]
    )
    write_text(manifest_path, manifest)


def make_readme(outdir: Path, tenx_reads: int) -> None:
    text = f"""# STARsolo synthetic smoke-test dataset

**Note on 10x:** the wrapper’s chemistry auto-detection requires >50,000 barcode matches in the
subsampled reads. This dataset generates {tenx_reads:,} reads for the 10x sample.

## Build a STAR index

```bash
mkdir -p reference/index
STAR --runThreadN 4 --runMode genomeGenerate \\
  --genomeDir reference/index \\
  --genomeFastaFiles reference/test.fa \\
  --sjdbGTFfile reference/test.gtf

  ## Example runs
  WL=$(pwd)/whitelists
REF=$(pwd)/reference/index
```
## Example runs
```bash
starsolo 10x fastqs/10x TENX1 --ref "$REF" --whitelist-dir "$WL"
starsolo dropseq fastqs/dropseq DROP1 --ref "$REF"
starsolo rhapsody fastqs/rhapsody RHAP1 --ref "$REF" --whitelist-dir "$WL"
starsolo indrops fastqs/indrops IND1 --ref "$REF" --whitelist-dir "$WL"
starsolo strt fastqs/strt STRT1 --ref "$REF" --whitelist-dir "$WL"
starsolo smartseq fastqs/smartseq/smartseq_example.manifest.tsv --ref "$REF"
```
"""
    write_text(outdir / "README.md", text)


def maybe_tar(outdir: Path, tar_path: Path) -> None:
    tar_path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "w:gz") as tf:
        tf.add(outdir, arcname=outdir.name)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Generate synthetic STARsolo wrapper test data"
    )
    ap.add_argument("--outdir", default="starsolo_synth", help="Output directory")
    ap.add_argument(
        "--tenx-reads", type=int, default=60000, help="Number of 10x reads (R1/R2)"
    )
    ap.add_argument(
        "--tar", action="store_true", help="Also create <outdir>.tar.gz next to outdir"
    )
    args = ap.parse_args()

    outdir = Path(args.outdir).resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit(f"Refusing to write into non-empty directory: {outdir}")

    outdir.mkdir(parents=True, exist_ok=True)

    make_reference(outdir)
    make_whitelists(outdir)
    make_fastqs(outdir, tenx_reads=args.tenx_reads)
    make_readme(outdir, tenx_reads=args.tenx_reads)

    if args.tar:
        tar_path = outdir.with_suffix(".tar.gz")
        maybe_tar(outdir, tar_path)
        print(f"Wrote {outdir} and {tar_path}")
    else:
        print(f"Wrote {outdir}")


if __name__ == "__main__":
    main()
