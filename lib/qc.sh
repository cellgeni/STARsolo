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
        737K-arc-v1.txt)         echo "arc"   ;;
        *)                       echo "Undef" ;;
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
    for i in "${CLEAN_DIRS[@]}"; do
        if [[ -d "$i/output" && -s "$i/Log.final.out" ]]; then
            if [[ -d "$i/_STARtmp" ]]; then
                log_warn "Sample $i did not run to completion: _STARtmp is still present!"
            fi
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
                        PCT=$(echo "$COMM" | awk -v v="$MIN" '{printf "%d\n",100*$1/v+0.5}')
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

    local REF WL R1 B G1 G2 E1 E2 F1 F2 C R2 CF R3 GC ST PAIRED F2PCT
    for i in "${CLEAN_DIRS[@]}"; do

        if [[ -d "$i/output/Gene/filtered/" && -s "$i/Log.out" ]]; then
            # --- Droplet-based ---
            PAIRED="Single"
            if grep -q "clip5pNbases 39 0" "$i/Log.out"; then
                PAIRED="Paired"
            fi

            if [[ -n "$FORCE_SPECIES" ]]; then
                REF="$FORCE_SPECIES"
            else
                REF=$(_qc_guess_species "$i/Log.out")
            fi

            WL=$(_qc_guess_wl "$i/Log.out")

            R1=$(grep "Number of Reads,"                                   "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
            B=$( grep "Reads With Valid Barcodes,"                         "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
            G1=$(grep "Reads Mapped to Genome: Unique+Multiple,"           "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
            G2=$(grep "Reads Mapped to Genome: Unique,"                    "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
            E1=$(grep "Reads Mapped to Gene: Unique+Multip.*e Gene,"       "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
            E2=$(grep "Reads Mapped to Gene: Unique Gene,"                 "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
            F1=$(grep "Reads Mapped to GeneFull: Unique+Multip.*e GeneFull," "$i/output/GeneFull/Summary.csv" | awk -F "," '{print $2}')
            F2=$(grep "Reads Mapped to GeneFull: Unique GeneFull,"         "$i/output/GeneFull/Summary.csv" | awk -F "," '{print $2}')
            C=$( grep "Estimated Number of Cells,"                         "$i/output/GeneFull/Summary.csv" | awk -F "," '{print $2}')
            R2=$(grep "Unique Reads in Cells Mapped to GeneFull,"          "$i/output/GeneFull/Summary.csv" | awk -F "," '{print $2}')
            CF=$(echo "$R1" | awk -v v="$R2" '{printf "%.3f\n",v/$1}')
            R3=$(grep "UMIs in Cells,"                                     "$i/output/GeneFull/Summary.csv" | awk -F "," '{print $2}')
            GC=$(grep "Median GeneFull per Cell,"                          "$i/output/GeneFull/Summary.csv" | awk -F "," '{print $2}')
            ST=$(grep "^soloStrand" "$i/Log.out" | grep RE-DEFINED | awk '{print $2}')

            if [[ -z "$ST" ]]; then
                ST="Undef"
            elif [[ "$PAIRED" == "Paired" ]]; then
                ST="Reverse"
            fi

            F2PCT=$(echo "$F2" | awk '{printf "%d\n",$1*100+0.5}')
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

            R1=$(grep "Number of Reads,"                                   "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
            B=$( grep "Reads With Valid Barcodes,"                         "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
            G1=$(grep "Reads Mapped to Genome: Unique+Multiple,"           "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
            G2=$(grep "Reads Mapped to Genome: Unique,"                    "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
            E1=$(grep "Reads Mapped to Gene: Unique+Multip.*e Gene,"       "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
            E2=$(grep "Reads Mapped to Gene: Unique Gene,"                 "$i/output/Gene/Summary.csv"     | awk -F "," '{print $2}')
            F1=$(grep "Reads Mapped to GeneFull: Unique+Multip.*e GeneFull," "$i/output/GeneFull/Summary.csv" | awk -F "," '{print $2}')
            F2=$(grep "Reads Mapped to GeneFull: Unique GeneFull,"         "$i/output/GeneFull/Summary.csv" | awk -F "," '{print $2}')
            C=$( zcat "$i/output/Gene/raw/barcodes.tsv.gz" | wc -l)
            R2=$(grep "yessubWLmatch_UniqueFeature" "$i/output/GeneFull/Features.stats" | awk '{print $2}')
            CF=$(echo "$R1" | awk -v v="$R2" '{printf "%.3f\n",v/$1}')
            R3=$(grep "yesUMIs" "$i/output/GeneFull/Features.stats" | awk '{print $2}')
            GC=$(zcat "$i/output/GeneFull/raw/matrix.mtx.gz" | awk 'NR>3 {print $2}' | uniq -c | awk '{print $1}' | sort -n | awk '{a[i++]=$1} END {print a[int(i/2)]}')
            ST=$(grep "^soloStrand" "$i/Log.out" | grep RE-DEFINED | awk '{print $2}')

            [[ -z "$ST" ]] && ST="Undef"

            F2PCT=$(echo "$F2" | awk '{printf "%d\n",$1*100+0.5}')
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
