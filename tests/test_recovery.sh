#!/bin/bash
# =============================================================================
# End-to-end tests of 'starsolo 10x' Velocyto handling and crash recovery.
#   bash tests/test_recovery.sh
# Needs STAR, seqtk, pbzip2 and samtools on PATH (e.g. inside the starsolo
# image) and the test dataset: data/test, or data/test.tar.gz which is
# unpacked into a temp directory. A stub STAR placed first on PATH runs the
# real STAR for the strand test and fakes crashes in the main run.
# =============================================================================

set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
REAL_STAR=$(command -v STAR) || { echo "STAR not on PATH" >&2; exit 1; }
export REAL_STAR

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

# ---------- stub STAR ---------------------------------------------------------

mkdir -p "$WORK/bin"
cat > "$WORK/bin/STAR" <<'EOF'
#!/bin/bash
# Strand test and genome removal: the real thing
case " $* " in
    *" test_forward/ "*|*" test_reverse/ "*|*" --genomeLoad Remove "*) exec "$REAL_STAR" "$@" ;;
esac

echo "$*" >> "$STUB_LOG"
ulimit -c 0
VELO=False
if [[ " $* " == *" Velocyto "* ]]; then VELO=True; fi

crash() { kill -SEGV $$; }
genome_loaded() { echo "Finished loading the genome: now" > Log.out; mkdir -p output _STARtmp; }

case "$STUB_MODE" in
    ok)
        exec "$REAL_STAR" "$@" ;;
    huge)       # finish, but report more than 2^32 input reads
        "$REAL_STAR" "$@" || exit $?
        sed -i 's/\(Number of input reads |\s*\)[0-9]*/\15696594413/' Log.final.out ;;
    segv_velocyto)  # write everything, then die in Velocyto counting
        # A real crash leaves the shared genome loaded: keep it, as STAR would
        args=("$@")
        for i in "${!args[@]}"; do
            if [[ ${args[i]} == LoadAndRemove ]]; then args[i]=LoadAndKeep; fi
        done
        "$REAL_STAR" "${args[@]}" || exit $?
        if [[ $VELO == True ]]; then
            sed -i '/Starting Solo post-map for Velocyto/q' Log.out
            echo "... Velocyto counting: allocated arrays" >> Log.out
            rm -rf output/Velocyto/filtered Log.final.out Aligned.sortedByCoord.out.bam
            mkdir -p _STARtmp && echo junk > _STARtmp/junk
            crash
        fi ;;
    segv_mapping)   # die while mapping, but only with Velocyto
        if [[ $VELO == True ]]; then genome_loaded; crash; fi
        exec "$REAL_STAR" "$@" ;;
    segv_always)
        genome_loaded; crash ;;
    segv_gene)      # die in Solo counting before Velocyto
        genome_loaded
        echo "... Starting Solo post-map for GeneFull" >> Log.out
        crash ;;
    fatal_param)
        echo "EXITING: FATAL INPUT ERROR: stub" >&2
        exit 102 ;;
esac
EOF
chmod +x "$WORK/bin/STAR"

# SysV segments owned by this user (ipcs cuts user names to 10 characters)
shm_count() {
    if [[ ! -r /proc/sysvipc/shm ]]; then echo 0; return 0; fi
    awk -v u="$(id -u)" 'NR > 1 && $8 == u { n++ } END { print n + 0 }' /proc/sysvipc/shm
}
SHM0=$(shm_count)

# run NAME MODE [starsolo 10x options…] → runs in $WORK/NAME; sets RC, D, CALLS
run() {
    local name=$1 mode=$2
    shift 2
    D=$WORK/$name/10x
    mkdir -p "$WORK/$name"
    export STUB_MODE=$mode STUB_LOG=$WORK/$name/stub.log
    : > "$STUB_LOG"
    RC=0
    ( cd "$WORK/$name" && PATH="$WORK/bin:$PATH" "$ROOT/bin/starsolo" 10x "$DATA/fastqs/10x" 10x \
        --ref "$DATA/reference/index" --whitelist-dir "$DATA/whitelists" --cpus 2 "$@" \
        > "$WORK/$name/run.log" 2>&1 ) || RC=$?
    CALLS=$(wc -l < "$STUB_LOG")
    echo "#    $name: exit $RC, $CALLS main STAR call(s)"
    if (( $(shm_count) > SHM0 )); then
        fail "$name: shared-memory genome released"
        "$REAL_STAR" --genomeDir "$DATA/reference/index" --genomeLoad Remove \
            --outFileNamePrefix "$WORK/$name/rm_" >/dev/null 2>&1 || true
    fi
}
has() { grep -q -- "$1" "$2" 2>/dev/null; }
marker() { grep -qx -- "$1" "$D/velocyto_failed.txt" 2>/dev/null; }

# ---------- scenarios ------------------------------------------------------------

run normal ok
expect "normal: exit 0, Velocyto counted" [ "$RC" = 0 -a -s "$D/output/Velocyto/raw/spliced.mtx.gz" ]
expect "normal: no marker, no _STARtmp" [ ! -e "$D/velocyto_failed.txt" -a ! -e "$D/_STARtmp" ]
expect "normal: strand.txt lists features" has "Solo features: Gene GeneFull Velocyto" "$D/strand.txt"

run override ok -- --soloFeatures Gene GeneFull GeneFull_Ex50pAS Velocyto
expect "override: exit 0, GeneFull_Ex50pAS written" [ "$RC" = 0 -a -d "$D/output/GeneFull_Ex50pAS" ]
expect "override: a single --soloFeatures reaches STAR" \
    [ "$(grep -o -- "--soloFeatures" "$WORK/override/stub.log" | wc -l)" = 1 ]

run novelo ok --no-velocyto
expect "--no-velocyto: exit 0, no Velocyto" [ "$RC" = 0 -a ! -e "$D/output/Velocyto" ]
expect "--no-velocyto: STAR not asked for Velocyto" [ "$(grep -c Velocyto "$WORK/novelo/stub.log")" = 0 ]
expect "--no-velocyto: reason in strand.txt" has "Velocyto: off (--no-velocyto)" "$D/strand.txt"

run conflict ok --no-velocyto -- --soloFeatures Gene Velocyto
expect "--no-velocyto + passthrough Velocyto: exit 1 before any work" [ "$RC" = 1 -a ! -e "$D" ]

run segv_velocyto segv_velocyto
expect "segv_velocyto: exit 0" [ "$RC" = 0 ]
expect "segv_velocyto: Gene kept, Velocyto removed" [ -s "$D/output/Gene/raw/matrix.mtx.gz" -a ! -e "$D/output/Velocyto" ]
expect "segv_velocyto: marker" marker "star_exit_code=139"
expect "segv_velocyto: marker stage" marker "stage=solo_postmap"
expect "segv_velocyto: _STARtmp removed" [ ! -e "$D/_STARtmp" ]
expect "segv_velocyto: strand.txt says FAILED" has "Velocyto: FAILED" "$D/strand.txt"
( cd "$WORK/segv_velocyto" && "$ROOT/bin/starsolo" qc 10x > qc.tsv 2> qc.err ) || true
expect "segv_velocyto: starsolo qc reports it" has "10x has no Velocyto output" "$WORK/segv_velocyto/qc.err"

run segv_velocyto_bam segv_velocyto --bam -- --limitBAMsortRAM 2000000000
expect "segv_velocyto --bam: exit 139 (BAM missing)" [ "$RC" = 139 ]
expect "segv_velocyto --bam: matrices kept and compressed" [ -s "$D/output/GeneFull/filtered/matrix.mtx.gz" ]
expect "segv_velocyto --bam: marker bam=missing" marker "bam=missing"

run segv_mapping segv_mapping
expect "segv_mapping: exit 0 after one rerun" [ "$RC" = 0 -a "$CALLS" = 2 ]
expect "segv_mapping: rerun without Velocyto" [ "$(sed -n 2p "$WORK/segv_mapping/stub.log" | grep -c Velocyto)" = 0 ]
expect "segv_mapping: first log kept, marker" [ -s "$D/Log.attempt1.out" ]
expect "segv_mapping: marker stage" marker "stage=mapping"
expect "segv_mapping: Gene written" [ -s "$D/output/Gene/raw/matrix.mtx.gz" ]

run segv_always segv_always
expect "segv_always: exit 139 after one rerun" [ "$RC" = 139 -a "$CALLS" = 2 ]
expect "segv_always: _STARtmp removed" [ ! -e "$D/_STARtmp" ]

run segv_always_novelo segv_always --no-velocyto
expect "segv_always --no-velocyto: exit 139, no rerun" [ "$RC" = 139 -a "$CALLS" = 1 ]

run segv_gene segv_gene
expect "segv_gene: exit 139, no rerun, no marker" [ "$RC" = 139 -a "$CALLS" = 1 -a ! -e "$D/velocyto_failed.txt" ]

run fatal_param fatal_param
expect "fatal_param: exit 102, no rerun" [ "$RC" = 102 -a "$CALLS" = 1 ]

run huge huge
expect "huge: exit 0, Velocyto removed" [ "$RC" = 0 -a ! -e "$D/output/Velocyto" -a -s "$D/output/Gene/raw/matrix.mtx.gz" ]
expect "huge: marker stage" marker "stage=read_index_limit"

# ---------- summary -------------------------------------------------------------

echo
echo "$PASS passed, $FAIL failed"
if (( FAIL > 0 )); then
    for f in "$WORK"/*/run.log; do echo "===== $f"; tail -n 15 "$f"; done
fi
(( FAIL == 0 ))
