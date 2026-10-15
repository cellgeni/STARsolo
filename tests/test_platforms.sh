#!/bin/bash
# =============================================================================
# Runs starsolo on each sample of the synthetic test dataset and checks the
# results against data/test/truth.
#   bash tests/test_platforms.sh [case …]          (default: all cases)
# Cases: 10x 10x5p dropseq dropseq_bam rhapsody indrops strt smartseq
# Needs STAR, seqtk, samtools and pbzip2 on PATH (e.g. inside the starsolo
# image). Data: $STARSOLO_TEST_DATA, data/test, or data/test.tar.gz unpacked
# into a temp directory. KEEP_WORK=1 keeps the run directories.
# =============================================================================

set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$ROOT/lib/common.sh"
source "$ROOT/lib/qc.sh"            # _qc_csv

WORK=$(mktemp -d)
if [[ -n ${KEEP_WORK:-} ]]; then
    echo "# work directory: $WORK"
else
    trap 'rm -rf -- "$WORK"' EXIT
fi

DATA=${STARSOLO_TEST_DATA:-$ROOT/data/test}
if [[ ! -d $DATA ]]; then
    tar -xzf "$ROOT/data/test.tar.gz" -C "$WORK"
    DATA=$WORK/data/test
fi
DATA=$(readlink -f "$DATA")

PASS=0
FAIL=0
ok()   { PASS=$((PASS + 1)); echo "ok   - $1"; }
fail() { FAIL=$((FAIL + 1)); echo "FAIL - $1"; }
expect() {   # expect NAME CONDITION…
    local name=$1; shift
    if "$@"; then ok "$name"; else fail "$name"; fi
}
nnz()      { zcat "$1" 2>/dev/null | awk '!/^%/ { print $3; exit }'; }    # non-zero entries
ncols()    { zcat "$1" 2>/dev/null | awk '!/^%/ && n++ { print $2 }' | sort -u | wc -l; }
ge()       { awk -v a="${1:-0}" -v b="$2" 'BEGIN { exit !(a + 0 >= b) }'; }
same_set() { cmp -s <(zcat "$1" 2>/dev/null | sort) <(sort "$2"); }
has()      { grep -qF -- "$1" "$2" 2>/dev/null; }

# case → "platform|sample|extra starsolo args"
case_spec() {
    case "$1" in
        10x)         echo "10x|10x|" ;;
        10x5p)       echo "10x|10x5p|--wl v2" ;;
        dropseq)     echo "dropseq|dropseq|" ;;
        dropseq_bam) echo "dropseq|dropseq|--bam -- --limitBAMsortRAM 2000000000" ;;
        rhapsody)    echo "rhapsody|rhapsody|" ;;
        indrops)     echo "indrops|indrops|" ;;
        strt)        echo "strt|strt|" ;;
        smartseq)    echo "smartseq|smartseq|" ;;
        *)           return 1 ;;
    esac
}

run_case() {
    local name=$1 spec platform sample extra rc=0 O N
    spec=$(case_spec "$name") || { fail "$name: unknown case"; return; }
    IFS='|' read -r platform sample extra <<< "$spec"

    # Run from a directory where data/test is the dataset: the Smart-seq
    # manifest uses paths relative to the repository root.
    mkdir -p "$WORK/$name/data"
    ln -sfn "$DATA" "$WORK/$name/data/test"
    cd "$WORK/$name"
    local -a ARGS=(--ref data/test/reference/index --whitelist-dir data/test/whitelists --cpus 2)
    if [[ $platform == smartseq ]]; then
        O=smartseq_example
        "$ROOT/bin/starsolo" smartseq data/test/fastqs/smartseq/smartseq_example.manifest.tsv \
            --ref data/test/reference/index --cpus 2 > run.log 2>&1 || rc=$?
    else
        O=$sample
        # shellcheck disable=SC2086  # extra holds whitespace-separated options
        "$ROOT/bin/starsolo" "$platform" "data/test/fastqs/$sample" "$sample" "${ARGS[@]}" $extra \
            > run.log 2>&1 || rc=$?
    fi
    N=$(wc -l < "$DATA/truth/$sample.cells.txt")
    local GS=$O/output/Gene/Summary.csv FS=$O/output/GeneFull/Summary.csv
    echo "#    $name: exit $rc; genome unique $(_qc_csv "$GS" "Reads Mapped to Genome: Unique"); GeneFull unique $(_qc_csv "$FS" "Reads Mapped to GeneFull: Unique GeneFull"); cells $(_qc_csv "$FS" "Estimated Number of Cells") (truth $N)"

    expect "$name: exit 0" [ "$rc" = 0 ]
    expect "$name: no _STARtmp left" [ ! -e "$O/_STARtmp" ]
    expect "$name: >= 85% of reads map uniquely" ge "$(_qc_csv "$GS" "Reads Mapped to Genome: Unique")" 0.85
    expect "$name: >= 70% of reads map uniquely to GeneFull" ge "$(_qc_csv "$FS" "Reads Mapped to GeneFull: Unique GeneFull")" 0.70
    expect "$name: Gene and GeneFull matrices have counts" \
        ge "$(( $(nnz "$O/output/Gene/raw/matrix.mtx.gz") > 0 && $(nnz "$O/output/GeneFull/raw/matrix.mtx.gz") > 0 ))" 1

    case "$platform" in
        strt|smartseq)
            expect "$name: counts for every true cell ($N)" [ "$(ncols "$O/output/Gene/raw/matrix.mtx.gz")" = "$N" ] ;;
        *)
            expect "$name: Gene filtered cells = truth" same_set "$O/output/Gene/filtered/barcodes.tsv.gz" "$DATA/truth/$sample.cells.txt"
            expect "$name: GeneFull filtered cells = truth" same_set "$O/output/GeneFull/filtered/barcodes.tsv.gz" "$DATA/truth/$sample.cells.txt" ;;
    esac

    local QC qcrc=0
    QC=$("$ROOT/bin/starsolo" qc "$O" --no-contamination-check 2> qc.err) || qcrc=$?
    expect "$name: starsolo qc exits 0 with one row" [ "$qcrc" = 0 -a "$(wc -l <<< "$QC")" = 2 ]
    if [[ $platform != strt && $platform != smartseq ]]; then
        expect "$name: qc reports $N cells" [ "$(awk -F'\t' 'NR == 2 { print $6 }' <<< "$QC")" = "$N" ]
    fi

    case "$name" in
        10x)
            expect "10x: single-end, Forward, 3M-february-2018" \
                eval 'has "Paired-end mode: False" 10x/strand.txt && has "5'"'"'): Forward" 10x/strand.txt && has "3M-february-2018.txt" 10x/strand.txt'
            expect "10x: Velocyto spliced and unspliced counts" \
                ge "$(( $(nnz 10x/output/Velocyto/raw/spliced.mtx.gz) > 0 && $(nnz 10x/output/Velocyto/raw/unspliced.mtx.gz) > 0 ))" 1
            expect "10x: no velocyto_failed.txt" [ ! -e 10x/velocyto_failed.txt ] ;;
        10x5p)
            expect "10x5p: paired-end, Reverse, 737K-august-2016" \
                eval 'has "Paired-end mode: True" 10x5p/strand.txt && has "5'"'"'): Reverse" 10x5p/strand.txt && has "737K-august-2016.txt" 10x5p/strand.txt'
            expect "10x5p: --clip5pNbases 39 0 (CB 16 + UMI 10 + TSO 13)" has "--clip5pNbases 39 0" 10x5p/Log.out
            expect "10x5p: qc reports Paired" [ "$(awk -F'\t' 'NR == 2 { print $11 }' <<< "$QC")" = Paired ]
            expect "10x5p: Velocyto counts" ge "$(nnz 10x5p/output/Velocyto/raw/spliced.mtx.gz)" 1 ;;
        dropseq_bam)
            expect "dropseq_bam: indexed BAM" [ -s dropseq/Aligned.sortedByCoord.out.bam.bai ]
            expect "dropseq_bam: CB tags in the BAM" \
                ge "$(samtools view dropseq/Aligned.sortedByCoord.out.bam | grep -c 'CB:Z:[ACGT]' || true)" 1 ;;
    esac
    cd "$WORK"
}

CASES=("$@")
if (( ${#CASES[@]} == 0 )); then
    CASES=(10x 10x5p dropseq dropseq_bam rhapsody indrops strt smartseq)
fi
for c in "${CASES[@]}"; do
    run_case "$c"
done

echo
echo "$PASS passed, $FAIL failed"
if (( FAIL > 0 )); then
    for c in "${CASES[@]}"; do echo "===== $c"; tail -n 15 "$WORK/$c/run.log" 2>/dev/null || true; done
fi
(( FAIL == 0 ))
