#!/bin/bash
# =============================================================================
# STARsolo CLI - QC aggregation
# =============================================================================
# Aggregate QC statistics for one or more STARsolo output directories.
# =============================================================================

# ---------- Helpers ----------------------------------------------------------

_qc_guess_species() {
    local logfile="$1"
    local genomedir
    genomedir=$(grep "^genomeDir" "$logfile" | tail -n1)
    if echo "$genomedir" | grep -q "/human/"; then
        echo "Human"
    elif echo "$genomedir" | grep -q "/mouse/"; then
        echo "Mouse"
    else
        echo "Other"
    fi
}

_qc_guess_wl() {
    local logfile="$1"
    local wlfile
    wlfile=$(grep "^soloCBwhitelist" "$logfile" | tail -n1 | awk '{print $NF}')
    wlfile=$(basename "$wlfile")
    case "$wlfile" in
        3M-february-2018.txt)    echo "v3"    ;;
        3M-3pgex-may-2023.txt)   echo "v4-3p" ;;
        3M-5pgex-jan-2023.txt)   echo "v4-5p" ;;
        737K-august-2016.txt)    echo "v2"    ;;
        737K-april-2014_rc.txt)  echo "v1"    ;;
        737K-arc-v1.txt|gex_737K-arc-v1.txt) echo "arc" ;;
        *)                       echo "Undef" ;;
    esac
}

# Prints the value of the first line of a STARsolo Summary.csv whose name
# matches REGEX (anchored at the start); empty if there is no such line.
#   $1: Summary.csv  $2: REGEX
_qc_csv() {
    { grep -m1 -- "^$2," "$1" 2>/dev/null || true; } | cut -d, -f2
}

# Prints the count for KEY in a STARsolo Features.stats file; empty if absent.
_qc_stat() {
    awk -v k="$2" '$1 == k { print $2; exit }' "$1" 2>/dev/null || true
}

# Prints A/B with 3 decimals, or NA if either is missing or B is 0.
_qc_ratio() {
    awk -v a="$1" -v b="$2" 'BEGIN { if (a == "" || b + 0 == 0) print "NA"; else printf "%.3f\n", a / b }'
}

# Prints a fraction (0.83) as a rounded percentage (83); 0 if not a number.
_qc_pct() {
    awk -v x="$1" 'BEGIN { printf "%d\n", (x ~ /^[0-9.eE+-]+$/ ? x * 100 + 0.5 : 0) }'
}

# Describes why Velocyto output is missing, from a velocyto_failed.txt marker.
_qc_velocyto_reason() {
    local marker="$1" stage rc
    stage=$(sed -n 's/^stage=//p' "$marker")
    rc=$(sed -n 's/^star_exit_code=//p' "$marker")
    case "$stage" in
        solo_postmap)     echo "STAR crashed in Velocyto counting (exit $rc)" ;;
        mapping)          echo "STAR crashed while mapping with Velocyto (exit $rc); rerun without Velocyto" ;;
        read_index_limit) echo "removed, >= 2^32 input reads make Velocyto counts wrong" ;;
        *)                echo "see $marker" ;;
    esac
}

# ---------- Main QC workflow -------------------------------------------------

run_qc() {
    local FORCE_SPECIES="" CHECK_CONTAMINATION=1
    local -a DIRS=()

    # --- Argument parsing ---
    while (( $# > 0 )); do
        case "$1" in
            -h|--help)
                show_qc_help; exit 0 ;;
            -s|--species|--specie)
                [[ -n ${2:-} && $2 != -* ]] || die "Option $1 needs a value. Run 'starsolo qc --help'."
                FORCE_SPECIES="$2"; shift 2 ;;
            --no-contamination-check)
                CHECK_CONTAMINATION=0; shift ;;
            -*)
                die "Unknown option: $1. Run 'starsolo qc --help'." ;;
            *)
                DIRS+=("$1"); shift ;;
        esac
    done

    if (( ${#DIRS[@]} == 0 )); then
        show_qc_help
        exit 0
    fi

    # Strip trailing slashes
    local -a CLEAN_DIRS=()
    local d
    for d in "${DIRS[@]}"; do
        CLEAN_DIRS+=("${d%/}")
    done

    # --- Check completion status ---
    log_info "Checking that all STARsolo jobs went to completion …"
    local i
    local MARK
    for i in "${CLEAN_DIRS[@]}"; do
        [[ -d "$i/output" ]] || continue
        MARK="$i/velocyto_failed.txt"
        if [[ -d "$i/_STARtmp" ]]; then
            log_warn "Sample $i did not run to completion: _STARtmp is still present!"
        elif [[ -s "$MARK" ]]; then
            log_warn "Sample $i has no Velocyto output: $(_qc_velocyto_reason "$MARK"). Other features were kept."
            if grep -q "^bam=missing" "$MARK"; then
                log_warn "Sample $i: the requested BAM was not written."
            fi
        elif [[ -s "$i/Log.out" && ! -s "$i/Log.final.out" ]]; then
            log_warn "Sample $i did not run to completion: Log.final.out is missing!"
        fi
    done

    # --- Cross-contamination check ---
    if (( CHECK_CONTAMINATION )); then
        log_info "Checking potential sample cross-contamination …"
        local TMPDIR_SORT
        TMPDIR_SORT=$(mktemp -d)

        local j local_i local_j N1 N2 MIN COMM PCT
        for i in "${CLEAN_DIRS[@]}"; do
            if [[ -s "$i/output/Gene/filtered/barcodes.tsv.gz" ]]; then
                for j in "${CLEAN_DIRS[@]}"; do
                    if [[ -s "$j/output/Gene/filtered/barcodes.tsv.gz" && "$i" != "$j" ]]; then
                        local_i="$TMPDIR_SORT/$(basename "$i").barcodes.tsv"
                        local_j="$TMPDIR_SORT/$(basename "$j").barcodes.tsv"
                        if [[ ! -s "$local_i" ]]; then
                            zcat "$i/output/Gene/filtered/barcodes.tsv.gz" | sort > "$local_i" &
                        fi
                        if [[ ! -s "$local_j" ]]; then
                            zcat "$j/output/Gene/filtered/barcodes.tsv.gz" | sort > "$local_j" &
                        fi
                        wait
                        N1=$(wc -l < "$local_i")
                        N2=$(wc -l < "$local_j")
                        MIN=$(( N1 <= N2 ? N1 : N2 ))
                        COMM=$(comm -12 "$local_i" "$local_j" | wc -l)
                        PCT=0
                        if (( MIN > 0 )); then PCT=$(( (200 * COMM / MIN + 1) / 2 )); fi
                        if (( PCT >= 20 )); then
                            log_warn "Samples $i ($N1 barcodes) and $j ($N2 barcodes) share $COMM ($PCT%) barcodes — higher than expected!"
                        fi
                    fi
                done
            elif [[ -d "$i/output" && -s "$i/Log.out" ]]; then
                log_info "Sample $i: no filtered output detected (probably plate-based)."
            fi
        done
        rm -rf "$TMPDIR_SORT"
    else
        log_info "Skipping cross-contamination check (--no-contamination-check)."
    fi

    # --- Extract and output statistics ---
    log_info "Extracting STARsolo stats …"
    >&2 echo

    echo -e "Sample\tRd_all\tRd_in_cells\tFrc_in_cells\tUMI_in_cells\tCells\tMed_nFeature\tGood_BC\tWL\tSpecies\tPaired\tStrand\tall_u+m\tall_u\texon_u+m\texon_u\tfull_u+m\tfull_u"

    local REF WL R1 B G1 G2 E1 E2 F1 F2 C R2 CF R3 GC ST PAIRED F2PCT GS FS
    for i in "${CLEAN_DIRS[@]}"; do

        if [[ -d "$i/output/Gene/filtered/" && -s "$i/Log.out" ]]; then
            # --- Droplet-based ---
            PAIRED="Single"
            if grep -qE '^soloBarcodeMate[[:space:]]+1([[:space:]]|$)' "$i/Log.out"; then
                PAIRED="Paired"
            fi

            if [[ -n "$FORCE_SPECIES" ]]; then
                REF="$FORCE_SPECIES"
            else
                REF=$(_qc_guess_species "$i/Log.out")
            fi

            WL=$(_qc_guess_wl "$i/Log.out")

            GS="$i/output/Gene/Summary.csv"
            FS="$i/output/GeneFull/Summary.csv"
            R1=$(_qc_csv "$GS" "Number of Reads")
            B=$( _qc_csv "$GS" "Reads With Valid Barcodes")
            G1=$(_qc_csv "$GS" "Reads Mapped to Genome: Unique+Multiple")
            G2=$(_qc_csv "$GS" "Reads Mapped to Genome: Unique")
            E1=$(_qc_csv "$GS" "Reads Mapped to Gene: Unique+Multip.*e Gene")
            E2=$(_qc_csv "$GS" "Reads Mapped to Gene: Unique Gene")
            F1=$(_qc_csv "$FS" "Reads Mapped to GeneFull: Unique+Multip.*e GeneFull")
            F2=$(_qc_csv "$FS" "Reads Mapped to GeneFull: Unique GeneFull")
            C=$( _qc_csv "$FS" "Estimated Number of Cells")
            R2=$(_qc_csv "$FS" "Unique Reads in Cells Mapped to GeneFull")
            CF=$(_qc_ratio "$R2" "$R1")
            R3=$(_qc_csv "$FS" "UMIs in Cells")
            GC=$(_qc_csv "$FS" "Median GeneFull per Cell")
            ST=$(awk '$1 == "soloStrand" && /RE-DEFINED/ { print $2; exit }' "$i/Log.out")

            if [[ -z "$ST" ]]; then
                ST="Undef"
            elif [[ "$PAIRED" == "Paired" ]]; then
                ST="Reverse"
            fi

            F2PCT=$(_qc_pct "$F2")
            if (( F2PCT <= 20 )); then
                log_warn "Sample $i: GeneFull percentage ($F2PCT%) is too low!"
            fi

            echo -e "$i\t$R1\t$R2\t$CF\t$R3\t$C\t$GC\t$B\t$WL\t$REF\t$PAIRED\t$ST\t$G1\t$G2\t$E1\t$E2\t$F1\t$F2"

        elif [[ -d "$i/output" && -s "$i/Log.out" ]]; then
            # --- Plate-based (Smart-seq2 etc.) ---
            PAIRED="Single"
            if grep -q "mate 2" "$i/Log.out"; then
                PAIRED="Paired"
            fi

            if [[ -n "$FORCE_SPECIES" ]]; then
                REF="$FORCE_SPECIES"
            else
                REF=$(_qc_guess_species "$i/Log.out")
            fi

            WL="Undef"
            if grep -q "manifest" "$i/Log.out"; then
                WL="None(Smart-seq)"
            fi

            GS="$i/output/Gene/Summary.csv"
            FS="$i/output/GeneFull/Summary.csv"
            R1=$(_qc_csv "$GS" "Number of Reads")
            B=$( _qc_csv "$GS" "Reads With Valid Barcodes")
            G1=$(_qc_csv "$GS" "Reads Mapped to Genome: Unique+Multiple")
            G2=$(_qc_csv "$GS" "Reads Mapped to Genome: Unique")
            E1=$(_qc_csv "$GS" "Reads Mapped to Gene: Unique+Multip.*e Gene")
            E2=$(_qc_csv "$GS" "Reads Mapped to Gene: Unique Gene")
            F1=$(_qc_csv "$FS" "Reads Mapped to GeneFull: Unique+Multip.*e GeneFull")
            F2=$(_qc_csv "$FS" "Reads Mapped to GeneFull: Unique GeneFull")
            C=$( zcat "$i/output/Gene/raw/barcodes.tsv.gz" | wc -l)
            R2=$(_qc_stat "$i/output/GeneFull/Features.stats" yessubWLmatch_UniqueFeature)
            CF=$(_qc_ratio "$R2" "$R1")
            R3=$(_qc_stat "$i/output/GeneFull/Features.stats" yesUMIs)
            GC=$(zcat "$i/output/GeneFull/raw/matrix.mtx.gz" | awk 'NR>3 {print $2}' | uniq -c | awk '{print $1}' | sort -n | awk '{a[i++]=$1} END {print a[int(i/2)]}')
            ST=$(awk '$1 == "soloStrand" && /RE-DEFINED/ { print $2; exit }' "$i/Log.out")

            [[ -z "$ST" ]] && ST="Undef"

            F2PCT=$(_qc_pct "$F2")
            if (( F2PCT <= 20 )); then
                log_warn "Sample $i: GeneFull percentage ($F2PCT%) is too low!"
            fi

            echo -e "$i\t$R1\t$R2\t$CF\t$R3\t$C\t$GC\t$B\t$WL\t$REF\t$PAIRED\t$ST\t$G1\t$G2\t$E1\t$E2\t$F1\t$F2"
        fi
    done
}

# ---------- Help -------------------------------------------------------------

show_qc_help() {
    cat <<'EOF'
Usage: starsolo qc <dir1> [dir2 … dirN] [options]

Aggregate QC statistics across one or more STARsolo output directories.

Arguments:
  <dir>                     One or more STARsolo output directories.

Checks: unfinished runs (_STARtmp left or no Log.final.out), Velocyto output
dropped by 'starsolo 10x' (velocyto_failed.txt), cross-sample barcode overlap
and low GeneFull mapping.

Options:
  -s, --species <name>      Override species label in output (e.g. human, mouse).
                            By default guessed from genomeDir in Log.out.
  --no-contamination-check  Skip cross-sample barcode overlap check.
  -h, --help                Show this help and exit.

Examples:
  starsolo qc sample1/ sample2/ | column -t
  starsolo qc results/GSM* --species human
  starsolo qc sample1/ sample2/ --no-contamination-check
EOF
}
