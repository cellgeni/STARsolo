#!/bin/bash
# =============================================================================
# Unit tests for the starsolo library functions (no alignment needed).
#   bash tests/test_lib.sh
# pbzip2 is needed for the process_output_files tests (skipped otherwise);
# STAR, if on PATH, is used to check that merged arguments have no duplicates.
# =============================================================================

set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$ROOT/etc/defaults.conf"
source "$ROOT/lib/common.sh"
source "$ROOT/lib/platform_10x.sh"
source "$ROOT/lib/qc.sh"

PASS=0
FAIL=0
ok()   { PASS=$((PASS + 1)); echo "ok   - $1"; }
fail() { FAIL=$((FAIL + 1)); echo "FAIL - $1"; }
# check NAME EXPECTED ACTUAL
check() {
    if [[ $2 == "$3" ]]; then
        ok "$1"
    else
        fail "$1"
        echo "       expected: $2"
        echo "       actual:   $3"
    fi
}
args() { printf '%q ' "${STAR_ARGS[@]}"; }

WORK=$(mktemp -d)
trap 'rm -rf -- "$WORK"' EXIT
cd "$WORK"

# ---------- merge_star_args --------------------------------------------------

BUILTIN=(--runThreadN 4 --readFilesIn "my reads.fq" bc.fq --soloCBposition 0_0_2_-1 3_1_3_8
         --outFilterScoreMin 30 --soloFeatures Gene GeneFull Velocyto
         --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx)

unset STAR_EXTRA_ARGS
merge_star_args "${BUILTIN[@]}"
check "merge: unset passthrough keeps built-ins" "$(printf '%q ' "${BUILTIN[@]}")" "$(args)"

STAR_EXTRA_ARGS=()
merge_star_args "${BUILTIN[@]}"
check "merge: empty passthrough keeps built-ins" "$(printf '%q ' "${BUILTIN[@]}")" "$(args)"

STAR_EXTRA_ARGS=(--outFilterScoreMin 20 --winAnchorMultimapNmax 50)
merge_star_args "${BUILTIN[@]}"
check "merge: scalar override replaces the built-in value" \
    "--runThreadN 4 --readFilesIn my\ reads.fq bc.fq --soloCBposition 0_0_2_-1 3_1_3_8 --soloFeatures Gene GeneFull Velocyto --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx --outFilterScoreMin 20 --winAnchorMultimapNmax 50 " \
    "$(args)"

STAR_EXTRA_ARGS=(--soloFeatures Gene GeneFull GeneFull_Ex50pAS Velocyto)
merge_star_args "${BUILTIN[@]}"
check "merge: multi-value override replaces the whole list" \
    "--runThreadN 4 --readFilesIn my\ reads.fq bc.fq --soloCBposition 0_0_2_-1 3_1_3_8 --outFilterScoreMin 30 --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx --soloFeatures Gene GeneFull GeneFull_Ex50pAS Velocyto " \
    "$(args)"

STAR_EXTRA_ARGS=()
GZIP="--readFilesCommand zcat"
BAM="--outSAMtype None"
merge_star_args --runThreadN 1 $GZIP $BAM
check "merge: unquoted GZIP/BAM strings are split" "--runThreadN 1 --readFilesCommand zcat --outSAMtype None " "$(args)"
GZIP=""
merge_star_args --runThreadN 1 $GZIP $BAM
check "merge: empty GZIP adds no empty argument" "--runThreadN 1 --outSAMtype None " "$(args)"

if command -v STAR >/dev/null; then
    STAR_EXTRA_ARGS=(--outFilterScoreMin 20 --soloFeatures Gene GeneFull)
    merge_star_args --genomeDir /nonexistent --outFilterScoreMin 30 --soloFeatures Gene GeneFull Velocyto \
        --outFileNamePrefix "$WORK/dup/"
    mkdir -p "$WORK/dup"
    STAR "${STAR_ARGS[@]}" > "$WORK/dup.out" 2>&1 || true
    if grep -qi "duplicate parameter" "$WORK/dup.out" "$WORK/dup/Log.out" 2>/dev/null; then
        fail "merge: STAR accepts the merged arguments (no duplicate parameter)"
    else
        ok "merge: STAR accepts the merged arguments (no duplicate parameter)"
    fi
else
    echo "skip - STAR not on PATH: duplicate-parameter check"
fi

# ---------- feature helpers -------------------------------------------------

STAR_EXTRA_ARGS=()
merge_star_args --soloFeatures Gene Velocyto GeneFull --soloOutFileNames output/ a b c
if has_solo_feature Velocyto "${STAR_ARGS[@]}"; then ok "has_solo_feature: found"; else fail "has_solo_feature: found"; fi
if has_solo_feature output/ "${STAR_ARGS[@]}"; then fail "has_solo_feature: ignores other params"; else ok "has_solo_feature: ignores other params"; fi
if has_solo_feature Velocyto --outSAMtype None; then fail "has_solo_feature: no --soloFeatures"; else ok "has_solo_feature: no --soloFeatures"; fi

star_args_drop_feature Velocyto
check "drop_feature: removes from the middle" "--soloFeatures Gene GeneFull --soloOutFileNames output/ a b c " "$(args)"
check "star_arg_values" "Gene GeneFull" "$(star_arg_values soloFeatures)"
check "star_arg_first: present" "output/" "$(star_arg_first soloOutFileNames -)"
check "star_arg_first: default" "-" "$(star_arg_first outTmpDir -)"

merge_star_args --soloFeatures Velocyto
if ( star_args_drop_feature Velocyto ) 2>/dev/null; then
    fail "drop_feature: dies when no feature would be left"
else
    ok "drop_feature: dies when no feature would be left"
fi

# ---------- remove_star_tmp -----------------------------------------------------

mkdir -p _STARtmp/x custom_tmp/y
merge_star_args --runThreadN 1
remove_star_tmp
if [[ ! -e _STARtmp ]]; then ok "remove_star_tmp: default _STARtmp"; else fail "remove_star_tmp: default _STARtmp"; fi
STAR_EXTRA_ARGS=(--outTmpDir custom_tmp)
merge_star_args --runThreadN 1
remove_star_tmp
if [[ ! -e custom_tmp ]]; then ok "remove_star_tmp: --outTmpDir"; else fail "remove_star_tmp: --outTmpDir"; fi
STAR_EXTRA_ARGS=()

# ---------- process_output_files --------------------------------------------

if command -v pbzip2 >/dev/null; then
    REAL_PBZIP2=$(command -v pbzip2)
    mkdir -p "$WORK/bin"
    printf '#!/bin/bash\necho "$*" >> "$PBLOG"\nexec %q "$@"\n' "$REAL_PBZIP2" > "$WORK/bin/pbzip2"
    chmod +x "$WORK/bin/pbzip2"

    for c in 1 4 5; do
        mkdir -p "$WORK/pof$c/output/Gene/raw" "$WORK/pof$c/output/SJ/raw"
        (
            cd "$WORK/pof$c"
            seq 1 50000 > Unmapped.out.mate1
            seq 1 50000 > Unmapped.out.mate2
            echo g > output/Gene/raw/features.tsv
            echo m > output/Gene/raw/matrix.mtx
            echo sj > SJ.out.tab
            ln -s ../../../SJ.out.tab output/SJ/raw/features.tsv
            export PBLOG="$WORK/pof$c/pb.log"
            PATH="$WORK/bin:$PATH" process_output_files output/ "$c" 2>/dev/null
            PATH="$WORK/bin:$PATH" process_output_files output/ "$c" 2>/dev/null    # rerun
        ) || true
        case $c in
            1) want="-f -9 -p1 Unmapped.out.mate1|-f -9 -p1 Unmapped.out.mate2" ;;
            4) want="-f -9 -p2 Unmapped.out.mate1|-f -9 -p2 Unmapped.out.mate2" ;;
            5) want="-f -9 -p3 Unmapped.out.mate1|-f -9 -p2 Unmapped.out.mate2" ;;
        esac
        check "process_output_files: pbzip2 threads with $c CPUs" "$want" \
            "$(sort -k4 "$WORK/pof$c/pb.log" 2>/dev/null | paste -sd'|')"
        if bzip2 -t "$WORK/pof$c/Unmapped.out.mate1.bz2" 2>/dev/null \
            && [[ -s $WORK/pof$c/output/Gene/raw/matrix.mtx.gz && -L $WORK/pof$c/output/SJ/raw/features.tsv ]]; then
            ok "process_output_files: outputs compressed, SJ symlink left alone, rerun ok ($c CPUs)"
        else
            fail "process_output_files: outputs compressed, SJ symlink left alone, rerun ok ($c CPUs)"
        fi
    done

    printf '#!/bin/bash\nexit 255\n' > "$WORK/bin/pbzip2"
    mkdir -p "$WORK/poffail/output"
    if ( cd "$WORK/poffail" && seq 1 10 > Unmapped.out.mate1 && seq 1 10 > Unmapped.out.mate2 \
           && PATH="$WORK/bin:$PATH" process_output_files output/ 4 ) 2>/dev/null; then
        fail "process_output_files: a pbzip2 failure is fatal"
    else
        ok "process_output_files: a pbzip2 failure is fatal"
    fi
else
    echo "skip - pbzip2 not on PATH: process_output_files tests"
fi

# ---------- 10x layout and Velocyto policy --------------------------------------

for wl in 737K-april-2014_rc.txt 3M-february-2018.txt 3M-3pgex-may-2023.txt 737K-arc-v1.txt gex_737K-arc-v1.txt; do
    check "layout: 3'-only $wl forces Forward, single-end" \
        "Forward|False|forced Forward: 3'-only whitelist" \
        "$(_10x_resolve_layout "/wl/$wl" Reverse True 30 40 2>/dev/null)"
done
check "layout: 737K Reverse + long R1 is paired-end" "Reverse|True|" "$(_10x_resolve_layout /wl/737K-august-2016.txt Reverse True 10 80)"
check "layout: 3M-5pgex Reverse + long R1 is paired-end" "Reverse|True|" "$(_10x_resolve_layout /wl/3M-5pgex-jan-2023.txt Reverse True 10 80)"
check "layout: 737K Reverse + short R1 is single-end" "Reverse|False|" "$(_10x_resolve_layout /wl/737K-august-2016.txt Reverse False 10 80)"
check "layout: Forward is single-end" "Forward|False|" "$(_10x_resolve_layout /wl/3M-february-2018.txt Forward True 80 10)"

policy() {   # policy NOVELO PAIRED BC → "<features>|<status>"
    merge_star_args --soloFeatures Gene GeneFull Velocyto
    _10x_velocyto_policy "$@" 2>/dev/null
    echo "$(star_arg_values soloFeatures)|$VELO_STATUS"
}
check "policy: default keeps Velocyto" "Gene GeneFull Velocyto|on" "$(policy False False /wl/3M-february-2018.txt)"
check "policy: --no-velocyto drops it" "Gene GeneFull|off (--no-velocyto)" "$(policy True False /wl/3M-february-2018.txt)"
check "policy: paired-end 737K keeps it" "Gene GeneFull Velocyto|on" "$(policy False True /wl/737K-august-2016.txt)"
check "policy: paired-end 3M-5pgex drops it" \
    "Gene GeneFull|off (paired-end GEM-X 5': Velocyto crashes STAR 2.7.10a)" \
    "$(policy False True /wl/3M-5pgex-jan-2023.txt)"
merge_star_args --soloFeatures Gene GeneFull
_10x_velocyto_policy False False /wl/3M-february-2018.txt
check "policy: Velocyto not requested" "off (not in --soloFeatures)" "$VELO_STATUS"

# ---------- 10x recovery helpers ----------------------------------------------

mkdir -p "$WORK/rec" && cd "$WORK/rec"
merge_star_args --soloFeatures Gene GeneFull Velocyto --soloCellFilter EmptyDrops_CR
for f in Gene GeneFull; do
    mkdir -p "output/$f/raw" "output/$f/filtered"
    echo s > "output/$f/Summary.csv"; echo m > "output/$f/raw/matrix.mtx"; echo m > "output/$f/filtered/matrix.mtx"
done
if _10x_solo_outputs_complete output/; then ok "outputs_complete: Gene + GeneFull"; else fail "outputs_complete: Gene + GeneFull"; fi
: > output/GeneFull/filtered/matrix.mtx
if _10x_solo_outputs_complete output/; then fail "outputs_complete: empty filtered matrix"; else ok "outputs_complete: empty filtered matrix"; fi

printf 'Finished loading the genome: x\n' > Log.out
if _10x_died_while_mapping 139; then ok "died_while_mapping: SIGSEGV while mapping"; else fail "died_while_mapping: SIGSEGV while mapping"; fi
if _10x_died_while_mapping 102; then fail "died_while_mapping: parameter error"; else ok "died_while_mapping: parameter error"; fi
printf '... Starting Solo post-map for GeneFull\n' >> Log.out
if _10x_died_while_mapping 139; then fail "died_while_mapping: crash after Solo started"; else ok "died_while_mapping: crash after Solo started"; fi

printf 'x ... Allocated and initialized readInfo array, nReadsInput = 4823088474\n' >> Log.out
check "n_reads_input: from Log.out" "4823088473" "$(_10x_n_reads_input)"
printf '                          Number of input reads |\t5696594413\n' > Log.final.out
check "n_reads_input: from Log.final.out" "5696594413" "$(_10x_n_reads_input)"

mkdir -p output/Velocyto/raw
: > strand.txt
_10x_check_read_index_limit output/ 2>/dev/null
if [[ ! -e output/Velocyto ]] && grep -q "^stage=read_index_limit" velocyto_failed.txt; then
    ok "read_index_limit: Velocyto removed above 2^32 reads"
else
    fail "read_index_limit: Velocyto removed above 2^32 reads"
fi
rm -f velocyto_failed.txt
mkdir -p output/Velocyto/raw
printf '                          Number of input reads |\t4294967295\n' > Log.final.out
_10x_check_read_index_limit output/
if [[ -d output/Velocyto && ! -e velocyto_failed.txt ]]; then
    ok "read_index_limit: Velocyto kept below 2^32 reads"
else
    fail "read_index_limit: Velocyto kept below 2^32 reads"
fi
cd "$WORK"

# ---------- qc --------------------------------------------------------------

make_qc_sample() {   # make_qc_sample DIR LOG_LINE
    mkdir -p "$1/output/Gene/filtered" "$1/output/GeneFull"
    printf '%s\n' "Number of Reads,1000" "Reads With Valid Barcodes,0.9" \
        "Reads Mapped to Genome: Unique+Multiple,0.9" "Reads Mapped to Genome: Unique,0.8" \
        "Reads Mapped to Gene: Unique+Multipe Gene,0.6" "Reads Mapped to Gene: Unique Gene,0.5" \
        > "$1/output/Gene/Summary.csv"
    printf '%s\n' "Reads Mapped to GeneFull: Unique+Multipe GeneFull,0.8" "Reads Mapped to GeneFull: Unique GeneFull,0.7" \
        "Estimated Number of Cells,10" "Unique Reads in Cells Mapped to GeneFull,500" \
        "UMIs in Cells,400" "Median GeneFull per Cell,30" \
        > "$1/output/GeneFull/Summary.csv"
    printf '%s\n' "genomeDir                         /ref/human/index" \
        "soloStrand                        Forward     ~RE-DEFINED" "$2" > "$1/Log.out"
    echo done > "$1/Log.final.out"
}
make_qc_sample qc_pe  "soloBarcodeMate                   1"
make_qc_sample qc_old "soloBarcodeMate                   1     ~RE-DEFINED"
make_qc_sample qc_se  "soloBarcodeMate                   0"
printf 'star_exit_code=139\nstage=solo_postmap\naction=velocyto_output_removed\n' > qc_se/velocyto_failed.txt
QC_OUT=$(run_qc qc_pe qc_old qc_se --no-contamination-check 2> qc.err)
check "qc: paired-end from soloBarcodeMate 1" "Paired Paired" "$(awk -F'\t' '$1 == "qc_pe" || $1 == "qc_old" { print $11 }' <<< "$QC_OUT" | paste -sd' ')"
check "qc: single-end otherwise" "Single" "$(awk -F'\t' '$1 == "qc_se" { print $11 }' <<< "$QC_OUT")"
if grep -q "qc_se has no Velocyto output: STAR crashed in Velocyto counting (exit 139)" qc.err; then
    ok "qc: reports velocyto_failed.txt"
else
    fail "qc: reports velocyto_failed.txt"
fi

# ---------- summary -------------------------------------------------------------

echo
echo "$PASS passed, $FAIL failed"
(( FAIL == 0 ))
