#!/usr/bin/env python3
"""Generate the synthetic test dataset used by the starsolo tests and CI.

Reads are sampled from a small random genome with three genes, so every
platform maps reads, counts genes, calls cells and (for 10x) counts Velocyto
spliced/unspliced reads:

- geneA (+) and geneC (-) have two exons and an intron, geneB (+) one exon;
- 80% of genic molecules come from mature mRNA, 20% from pre-mRNA (intronic),
  plus 5% intergenic and 2% unmappable molecules;
- each droplet sample has 40-50 cells with many reads and hundreds of empty
  barcodes with a few reads, so STAR's cell filtering keeps exactly the cells.
  (With fewer than 31 cells every barcode passes STAR's filter, and with
  <= 45,000 barcodes EmptyDrops_CR falls back to the CellRanger2.2 filter.)

Samples (fastqs/<sample>/, the directory name is the sample id used by CI):
  10x       10x 3' v3, single-end, 3M-february-2018 whitelist (auto-detected)
  10x5p     10x 5' v2, paired-end, 737K-august-2016 (run with --wl v2)
  dropseq, rhapsody, indrops, strt   the layouts the wrapper expects
  smartseq  paired-end, manifest smartseq_example.manifest.tsv

truth/<sample>.cells.txt lists the barcodes of the true cells, in the form
STAR writes them. The STAR index is built by scripts/build_test_data.sh.
Standard library only (Python >= 3.8); deterministic for a given --seed.
"""

import argparse
import gzip
import random
from pathlib import Path

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

GENOME_LEN = 6000
# name: (strand, exons as 1-based inclusive genome coordinates)
GENES = {
    "geneA": ("+", [(501, 900), (1401, 1800)]),
    "geneB": ("+", [(2501, 3100)]),
    "geneC": ("-", [(3901, 4200), (4701, 5000)]),
}
INTERGENIC = [(1, 500), (1801, 2500), (3101, 3900), (5001, 6000)]
# Two cell types with different expression, and the mix for empty droplets
CELL_TYPES = [(0.6, 0.3, 0.1), (0.1, 0.3, 0.6)]
AMBIENT = (0.35, 0.3, 0.35)

TSO_TAIL = "TTTCTTATATGGG"                 # last 13 nt of the 10x 5' TSO
INDROPS_ADAPTER = "GAGTGATTGCTTGTGACGCCTT"
RHAP_LINKER1 = "ACTGGCCTGCGA"
RHAP_LINKER2 = "GGTAGCGGTGACA"


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def random_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


# ---------- reference ----------------------------------------------------------


def genome_ok(seq):
    """No homopolymer >= 8, every 20-mer unique on both strands, no adapters."""
    for base in "ACGT":
        if base * 8 in seq:
            return False
    seen = set()
    for s in (seq, revcomp(seq)):
        for i in range(len(s) - 19):
            kmer = s[i:i + 20]
            if kmer in seen:
                return False
            seen.add(kmer)
    for motif in (TSO_TAIL, INDROPS_ADAPTER, RHAP_LINKER1, RHAP_LINKER2, "A" * 7):
        if motif in seq or revcomp(motif) in seq:
            return False
    return True


def make_genome(rng):
    while True:
        seq = [rng.choice("ACGT") for _ in range(GENOME_LEN)]
        for strand, exons in GENES.values():
            for (_, end1), (start2, _) in zip(exons, exons[1:]):
                first, last = end1 + 1, start2 - 1       # intron, 1-based
                donor, acceptor = ("GT", "AG") if strand == "+" else ("CT", "AC")
                seq[first - 1:first + 1] = donor
                seq[last - 2:last] = acceptor
        genome = "".join(seq)
        if genome_ok(genome):
            return genome


def make_gtf():
    lines = []
    for name, (strand, exons) in GENES.items():
        attrs = f'gene_id "{name}"; gene_name "{name}"; gene_biotype "protein_coding";'
        tx = f'gene_id "{name}"; transcript_id "tx{name[-1]}"; gene_name "{name}"; gene_biotype "protein_coding";'
        start, end = exons[0][0], exons[-1][1]
        lines.append(f"chrTest\tsynth\tgene\t{start}\t{end}\t.\t{strand}\t.\t{attrs}")
        lines.append(f"chrTest\tsynth\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t{tx}")
        order = exons if strand == "+" else exons[::-1]
        for i, (s, e) in enumerate(order, 1):
            lines.append(f'chrTest\tsynth\texon\t{s}\t{e}\t.\t{strand}\t.\t{tx} exon_number "{i}";')
    return "\n".join(lines) + "\n"


def make_transcripts(genome):
    """name -> (mature mRNA, pre-mRNA), both in transcript (sense) orientation."""
    tx = {}
    for name, (strand, exons) in GENES.items():
        mrna = "".join(genome[s - 1:e] for s, e in exons)
        pre = genome[exons[0][0] - 1:exons[-1][1]]
        if strand == "-":
            mrna, pre = revcomp(mrna), revcomp(pre)
        tx[name] = (mrna, pre)
    return tx


# ---------- molecules --------------------------------------------------------


class Molecules:
    """Draws cDNA molecules (sense orientation) from the genome."""

    def __init__(self, rng, genome, tx):
        self.rng, self.genome, self.tx = rng, genome, tx

    def draw(self, weights, length, anchor):
        """A sense cDNA fragment of `length` nt.

        anchor: "3p" (3' end of the transcript, droplet 3' protocols),
                "5p" (5' end) or "any" (uniform, Smart-seq).
        """
        rng = self.rng
        r = rng.random()
        if r < 0.02:                                   # unmappable
            return random_seq(rng, length)
        if r < 0.07:                                   # intergenic, either strand
            s, e = rng.choice(INTERGENIC)
            start = rng.randint(s - 1, e - length)
            frag = self.genome[start:start + length]
            return frag if rng.random() < 0.5 else revcomp(frag)
        gene = rng.choices(list(GENES), weights=weights)[0]
        mrna, pre = self.tx[gene]
        mature = rng.random() < 0.8
        src = mrna if mature else pre
        n = len(src)
        if anchor == "3p":
            end = n - rng.randint(0, min(350, n - length)) if mature else rng.randint(length, n)
            return src[end - length:end]
        if anchor == "5p" and mature:
            start = rng.randint(0, min(20, n - length))
        else:
            start = rng.randint(0, n - length)
        return src[start:start + length]


def split_reads(rng, n_reads, n_molecules):
    """Reads per molecule: at least one each, the rest are PCR duplicates."""
    counts = [1] * n_molecules
    for _ in range(n_reads - n_molecules):
        counts[rng.randrange(n_molecules)] += 1
    return counts


def unique_umis(rng, n, length):
    umis, seen = [], set()           # list keeps the order independent of str hashing
    while len(umis) < n:
        umi = random_seq(rng, length)
        if umi not in seen:
            seen.add(umi)
            umis.append(umi)
    return umis


def barcode_reads(rng, n_reads, dup):
    """(number of molecules, reads per molecule) for one barcode."""
    n_mol = max(1, round(n_reads / dup))
    return n_mol, split_reads(rng, n_reads, n_mol)


# ---------- barcodes ---------------------------------------------------------


def neighbours(bc):
    for i, base in enumerate(bc):
        for b in "ACGT":
            if b != base:
                yield bc[:i] + b + bc[i + 1:]


class BarcodePool:
    """Random barcodes with pairwise Hamming distance >= 3 within the pool."""

    def __init__(self, rng):
        self.rng = rng
        self.near = set()          # barcodes and their 1-mismatch neighbours

    def draw(self, n, length):
        out = []
        while len(out) < n:
            bc = random_seq(self.rng, length)
            if "AAAAA" in bc or "TTTTT" in bc:
                continue
            if bc in self.near or any(nb in self.near for nb in neighbours(bc)):
                continue
            out.append(bc)
            self.near.add(bc)
            self.near.update(neighbours(bc))
        return out

    def invalid(self, length):
        """A barcode that is not within one mismatch of any pooled barcode."""
        while True:
            bc = random_seq(self.rng, length)
            if bc not in self.near:
                return bc


def mutate(rng, seq, pos_range):
    """One substitution at a random position in pos_range."""
    i = rng.randrange(*pos_range)
    return seq[:i] + rng.choice([b for b in "ACGT" if b != seq[i]]) + seq[i + 1:]


# ---------- output -----------------------------------------------------------


def write_text(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_lines(path, items):
    write_text(path, "".join(f"{x}\n" for x in sorted(items)))


def write_fastq(path, records):
    """records: (name, seq); gzip without a timestamp so output is reproducible."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        for name, seq in records:
            qual = "".join("#" if c == "N" else "F" for c in seq)
            gz.write(f"@{name}\n{seq}\n+\n{qual}\n".encode())


def write_pairs(outdir, prefix_r1, prefix_r2, sample, pairs):
    write_fastq(outdir / prefix_r1, ((f"{sample}:{i}", r1) for i, (r1, _) in enumerate(pairs)))
    write_fastq(outdir / prefix_r2, ((f"{sample}:{i}", r2) for i, (_, r2) in enumerate(pairs)))


# ---------- samples ----------------------------------------------------------


def droplet_barcodes(rng, n_cells, cell_reads, n_empty, empty_reads):
    """[(is_cell, weights, n_reads)] for the cells and empty droplets."""
    out = [(True, CELL_TYPES[i % 2], rng.randint(*cell_reads)) for i in range(n_cells)]
    out += [(False, AMBIENT, rng.randint(*empty_reads)) for _ in range(n_empty)]
    return out


def sample_10x(rng, mol, pool, n_reads_target):
    """10x 3' v3: R1 = CB16 + UMI12 (28 nt), R2 = 90 nt sense cDNA."""
    n_cells, n_empty = 50, 2000
    empty_reads = (1, 4)
    per_cell = (n_reads_target - n_empty * 2.5) / n_cells
    plan = droplet_barcodes(rng, n_cells, (int(per_cell * 0.9), int(per_cell * 1.1)), n_empty, empty_reads)
    cbs = pool.draw(len(plan), 16)
    pairs, cells = [], []
    for (is_cell, weights, n_reads), cb in zip(plan, cbs):
        if is_cell:
            cells.append(cb)
        n_mol, reads = barcode_reads(rng, n_reads, 1.3)
        for umi, n in zip(unique_umis(rng, n_mol, 12), reads):
            r = rng.random()
            if r < 0.03:                               # polyA tail for CellRanger4 clipping
                gene = rng.choices(list(GENES), weights=weights)[0]
                cdna = mol.tx[gene][0][-60:] + "A" * 30
            else:
                cdna = mol.draw(weights, 90, "3p")
            for _ in range(n):
                q = rng.random()
                cb_read = pool.invalid(16) if q < 0.005 else mutate(rng, cb, (0, 16)) if q < 0.015 else cb
                pairs.append((cb_read + umi, cdna))
    rng.shuffle(pairs)
    return pairs, cells, cbs


def sample_10x5p(rng, mol, pool):
    """10x 5' v2: R1 = CB16 + UMI10 + TSO tail 13 + 51 nt sense cDNA (90 nt),
    R2 = 90 nt antisense from the other end of the same fragment."""
    plan = droplet_barcodes(rng, 40, (100, 160), 400, (1, 3))
    cbs = pool.draw(len(plan), 16)
    pairs, cells = [], []
    for (is_cell, weights, n_reads), cb in zip(plan, cbs):
        if is_cell:
            cells.append(cb)
        n_mol, reads = barcode_reads(rng, n_reads, 1.2)
        for umi, n in zip(unique_umis(rng, n_mol, 10), reads):
            frag = mol.draw(weights, rng.randint(200, 450), "5p")
            r1_cdna, r2 = frag[:51], revcomp(frag[-90:])
            for _ in range(n):
                q = rng.random()
                cb_read = pool.invalid(16) if q < 0.005 else mutate(rng, cb, (0, 16)) if q < 0.015 else cb
                pairs.append((cb_read + umi + TSO_TAIL + r1_cdna, r2))
    rng.shuffle(pairs)
    return pairs, cells, cbs


def sample_dropseq(rng, mol, pool):
    """Drop-seq: R1 = CB12 + UMI8, R2 = 60 nt sense cDNA; no whitelist."""
    plan = droplet_barcodes(rng, 40, (150, 250), 400, (1, 3))
    cbs = pool.draw(len(plan), 12)
    pairs, cells = [], []
    for (is_cell, weights, n_reads), cb in zip(plan, cbs):
        if is_cell:
            cells.append(cb)
        n_mol, reads = barcode_reads(rng, n_reads, 1.2)
        for umi, n in zip(unique_umis(rng, n_mol, 8), reads):
            cdna = mol.draw(weights, 60, "3p")
            pairs.extend([(cb + umi, cdna)] * n)
    rng.shuffle(pairs)
    return pairs, cells


def sample_rhapsody(rng, mol, wl):
    """BD Rhapsody: R1 = CL1 9 + linker 12 + CL2 9 + linker 13 + CL3 9 + UMI 8 + T15
    (75 nt), R2 = 75 nt sense cDNA."""
    combos = set()
    while len(combos) < 440:
        combos.add(tuple(rng.choice(w) for w in wl))
    combos = sorted(combos)
    rng.shuffle(combos)
    plan = droplet_barcodes(rng, 40, (150, 250), 400, (1, 3))
    pairs, cells = [], []
    for (is_cell, weights, n_reads), (c1, c2, c3) in zip(plan, combos):
        if is_cell:
            cells.append(f"{c1}_{c2}_{c3}")
        n_mol, reads = barcode_reads(rng, n_reads, 1.2)
        for umi, n in zip(unique_umis(rng, n_mol, 8), reads):
            cdna = mol.draw(weights, 75, "3p")
            for _ in range(n):
                cl1 = mutate(rng, c1, (0, 9)) if rng.random() < 0.01 else c1
                pairs.append((cl1 + RHAP_LINKER1 + c2 + RHAP_LINKER2 + c3 + umi + "T" * 15, cdna))
    rng.shuffle(pairs)
    return pairs, cells


def sample_indrops(rng, mol, wl1, wl2):
    """inDrops: R1 = CB1 (8-11) + adapter 22 + CB2 8 + UMI 6, padded with T to
    60 nt; R2 = 60 nt sense cDNA."""
    combos = set()
    while len(combos) < 440:
        combos.add((rng.choice(wl1), rng.choice(wl2)))
    combos = sorted(combos)
    rng.shuffle(combos)
    plan = droplet_barcodes(rng, 40, (150, 250), 400, (1, 3))
    pairs, cells = [], []
    for (is_cell, weights, n_reads), (c1, c2) in zip(plan, combos):
        if is_cell:
            cells.append(f"{c1}_{c2}")
        n_mol, reads = barcode_reads(rng, n_reads, 1.2)
        for umi, n in zip(unique_umis(rng, n_mol, 6), reads):
            cdna = mol.draw(weights, 60, "3p")
            for _ in range(n):
                adapter = mutate(rng, INDROPS_ADAPTER, (0, 22)) if rng.random() < 0.01 else INDROPS_ADAPTER
                r1 = c1 + adapter + c2 + umi
                pairs.append((r1 + "T" * (60 - len(r1)), cdna))
    rng.shuffle(pairs)
    return pairs, cells


def sample_strt(rng, mol, wl):
    """STRT-seq as the wrapper runs it: *_1 = CB8 + UMI8 (barcode read),
    *_2 = 75 nt sense cDNA from the 5' end."""
    wells = wl[:48]
    pairs = []
    for i, cb in enumerate(wells):
        n_mol, reads = barcode_reads(rng, rng.randint(60, 120), 1.2)
        for umi, n in zip(unique_umis(rng, n_mol, 8), reads):
            cdna = mol.draw(CELL_TYPES[i % 2], 75, "5p")
            pairs.extend([(cb + umi, cdna)] * n)
    rng.shuffle(pairs)
    return pairs, wells


def sample_smartseq(rng, mol, n_cells=8):
    """Smart-seq2: unstranded 75+75 nt pairs from full-length fragments."""
    out = {}
    for i in range(n_cells):
        pairs = []
        for _ in range(rng.randint(200, 300)):
            frag = mol.draw(CELL_TYPES[i % 2], rng.randint(150, 400), "any")
            if rng.random() < 0.5:
                frag = revcomp(frag)
            pairs.append((frag[:75], revcomp(frag[-75:])))
        out[f"CELL{i + 1:02d}"] = pairs
    return out


# ---------- main -------------------------------------------------------------


README = """# STARsolo synthetic test dataset

Generated by `scripts/generate_test_data.py` (seed {seed}); packaged with its
STAR index by `scripts/build_test_data.sh`. Reads are sampled from a 6 kb
random genome (`reference/`) with three genes: geneA (+) and geneC (-) with an
intron, geneB (+) without. Every sample maps, counts genes and calls cells;
`truth/<sample>.cells.txt` lists the true cells as STAR writes them.

| Sample | Platform | Run with |
|:--|:--|:--|
{rows}

Run from the repository root (the Smart-seq manifest uses paths relative to it):

```bash
tar -xzf data/test.tar.gz
bin/starsolo 10x data/test/fastqs/10x 10x --ref data/test/reference/index --whitelist-dir data/test/whitelists
bin/starsolo 10x data/test/fastqs/10x5p 10x5p --wl v2 --ref data/test/reference/index --whitelist-dir data/test/whitelists
bin/starsolo smartseq data/test/fastqs/smartseq/smartseq_example.manifest.tsv --ref data/test/reference/index
```
"""


def main():
    ap = argparse.ArgumentParser(description="Generate the synthetic starsolo test dataset")
    ap.add_argument("--outdir", default="data/test", help="Output directory (must be empty)")
    ap.add_argument("--seed", default="20261015", help="Random seed")
    ap.add_argument("--tenx-reads", type=int, default=60000,
                    help="Approximate 10x 3' reads (chemistry detection needs > 50,000)")
    ap.add_argument("--manifest-prefix", default="data/test/fastqs/smartseq",
                    help="Directory written into the Smart-seq manifest")
    args = ap.parse_args()

    out = Path(args.outdir)
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"Refusing to write into non-empty directory: {out}")
    rng = lambda name: random.Random(f"{args.seed}/{name}")

    genome = make_genome(rng("genome"))
    write_text(out / "reference" / "test.fa",
               ">chrTest\n" + "\n".join(genome[i:i + 60] for i in range(0, len(genome), 60)) + "\n")
    write_text(out / "reference" / "test.gtf", make_gtf())
    tx = make_transcripts(genome)

    # Whitelists. 10x barcodes come from one pool, so they differ across files too.
    wl, fq, truth = out / "whitelists", out / "fastqs", out / "truth"
    pool10x = BarcodePool(rng("10x-barcodes"))
    pairs, cells, used = sample_10x(rng("10x"), Molecules(rng("10x-mol"), genome, tx), pool10x, args.tenx_reads)
    write_lines(wl / "3M-february-2018.txt", used + pool10x.draw(1000, 16))
    write_pairs(fq / "10x", "TENX1_1.fastq.gz", "TENX1_2.fastq.gz", "10x", pairs)
    write_lines(truth / "10x.cells.txt", cells)
    summary = [("10x", "10x", "auto-detected (3' v3)", len(pairs), len(cells))]
    exact = sum(r1[:16] in set(used) for r1, _ in pairs)
    if exact <= 50000:
        raise SystemExit(f"10x: only {exact} reads have a whitelisted barcode; chemistry detection needs > 50,000")

    pairs, cells, used = sample_10x5p(rng("10x5p"), Molecules(rng("10x5p-mol"), genome, tx), pool10x)
    write_lines(wl / "737K-august-2016.txt", used + pool10x.draw(1000, 16))
    half = len(pairs) // 2
    write_pairs(fq / "10x5p", "10x5p_S1_L001_R1_001.fastq.gz", "10x5p_S1_L001_R2_001.fastq.gz", "10x5p:L1", pairs[:half])
    write_pairs(fq / "10x5p", "10x5p_S1_L002_R1_001.fastq.gz", "10x5p_S1_L002_R2_001.fastq.gz", "10x5p:L2", pairs[half:])
    write_lines(truth / "10x5p.cells.txt", cells)
    summary.append(("10x5p", "10x", "`--wl v2` (5' v2, paired-end)", len(pairs), len(cells)))

    for name in ("3M-3pgex-may-2023.txt", "3M-5pgex-jan-2023.txt", "737K-arc-v1.txt"):
        write_lines(wl / name, pool10x.draw(500, 16))
    decoys14 = BarcodePool(rng("10x-v1"))
    decoys14.near |= {bc[:14] for bc in pool10x.near if len(bc) == 16}
    write_lines(wl / "737K-april-2014_rc.txt", decoys14.draw(500, 14))

    pairs, cells = sample_dropseq(rng("dropseq"), Molecules(rng("dropseq-mol"), genome, tx), BarcodePool(rng("dropseq-bc")))
    write_pairs(fq / "dropseq", "DROP1_1.fastq.gz", "DROP1_2.fastq.gz", "dropseq", pairs)
    write_lines(truth / "dropseq.cells.txt", cells)
    summary.append(("dropseq", "dropseq", "", len(pairs), len(cells)))

    rpool = BarcodePool(rng("rhapsody-bc"))
    rwl = [rpool.draw(24, 9) for _ in range(3)]
    for i, w in enumerate(rwl, 1):
        write_lines(wl / f"Rhapsody_bc{i}.txt", w)
    pairs, cells = sample_rhapsody(rng("rhapsody"), Molecules(rng("rhapsody-mol"), genome, tx), rwl)
    write_pairs(fq / "rhapsody", "RHAP1.R1.fastq.gz", "RHAP1.R2.fastq.gz", "rhapsody", pairs)
    write_lines(truth / "rhapsody.cells.txt", cells)
    summary.append(("rhapsody", "rhapsody", "", len(pairs), len(cells)))

    # inDrops bc1 must contain every length from the shortest to the longest
    ipool = BarcodePool(rng("indrops-bc"))
    iwl1 = [bc for n in (8, 9, 10, 11) for bc in ipool.draw(12, n)]
    iwl2 = ipool.draw(48, 8)
    write_lines(wl / "inDrops_Ambrose2_bc1.txt", iwl1)
    write_lines(wl / "inDrops_Ambrose2_bc2.txt", iwl2)
    pairs, cells = sample_indrops(rng("indrops"), Molecules(rng("indrops-mol"), genome, tx), iwl1, iwl2)
    write_pairs(fq / "indrops", "IND1_1.fastq.gz", "IND1_2.fastq.gz", "indrops", pairs)
    write_lines(truth / "indrops.cells.txt", cells)
    summary.append(("indrops", "indrops", "", len(pairs), len(cells)))

    swl = BarcodePool(rng("strt-bc")).draw(96, 8)
    write_lines(wl / "96_barcodes.list", swl)
    pairs, cells = sample_strt(rng("strt"), Molecules(rng("strt-mol"), genome, tx), swl)
    write_pairs(fq / "strt", "STRT1_1.fastq.gz", "STRT1_2.fastq.gz", "strt", pairs)
    write_lines(truth / "strt.cells.txt", cells)
    summary.append(("strt", "strt", "", len(pairs), len(cells)))

    cells = sample_smartseq(rng("smartseq"), Molecules(rng("smartseq-mol"), genome, tx))
    manifest = []
    for cell, pairs in cells.items():
        write_pairs(fq / "smartseq", f"{cell}_R1.fastq.gz", f"{cell}_R2.fastq.gz", cell, pairs)
        manifest.append(f"{args.manifest_prefix}/{cell}_R1.fastq.gz\t{args.manifest_prefix}/{cell}_R2.fastq.gz\t{cell}")
    write_text(fq / "smartseq" / "smartseq_example.manifest.tsv", "\n".join(manifest) + "\n")
    write_lines(truth / "smartseq.cells.txt", cells)
    summary.append(("smartseq", "smartseq", "manifest", sum(len(p) for p in cells.values()), len(cells)))

    write_text(truth / "summary.tsv", "sample\tplatform\treads\tcells\n"
               + "".join(f"{s}\t{p}\t{n}\t{c}\n" for s, p, _, n, c in summary))
    rows = "\n".join(f"| `{s}` | {p} | {how or '—'} |" for s, p, how, _, _ in summary)
    write_text(out / "README.md", README.format(seed=args.seed, rows=rows))
    print("\n".join(f"{s}\t{n} reads\t{c} cells" for s, _, _, n, c in summary))


if __name__ == "__main__":
    main()
