#!/bin/bash
# =============================================================================
# STARsolo CLI - 10x Genomics platform
# =============================================================================
# Automatic chemistry detection (v1 / v2 / v3 / v4 / multiome),
# strand-specificity detection, and paired-end support.
# =============================================================================

# ---------- 10x-specific helpers -------------------------------------------

# Extracts a random subsample of reads for chemistry testing.
#   $1: R1 (comma-sep)   $2: R2 (comma-sep)   $3: ZCMD
_10x_extract_test_reads() {
    local R1=$1 R2=$2 ZCMD=$3
    local COUNT=0

    for i in $(echo "$R1" | tr ',' ' '); do
        $ZCMD "$i" | head -4000000 > "$COUNT.R1_head" &
        COUNT=$((COUNT + 1))
    done
    wait

    COUNT=0
    for i in $(echo "$R2" | tr ',' ' '); do
        $ZCMD "$i" | head -4000000 > "$COUNT.R2_head" &
        COUNT=$((COUNT + 1))
    done
    wait

    cat *.R1_head | seqtk sample -s100 - 200000 > test.R1.fastq &
    cat *.R2_head | seqtk sample -s100 - 200000 > test.R2.fastq &
    wait
    rm -f *.R1_head *.R2_head
}

# Resolves a chemistry name to its whitelist file and CB/UMI geometry.
# Names match the chemistry ids emitted by infer_10x_run.py, plus short
# aliases for interactive use.
#   $1: chemistry name  $2: WL directory
# Prints: "BC|CBLEN|UMILEN"
_10x_resolve_chemistry() {
    local NAME=$1 WL=$2
    local FILE CBLEN UMILEN

    case "${NAME,,}" in
        gex_3pv1|v1|3pv1)
            FILE=737K-april-2014_rc.txt;  CBLEN=14; UMILEN=10 ;;
        gex_3pv2_or_5pv1v2|v2|3pv2|5pv1|5pv2)
            FILE=737K-august-2016.txt;    CBLEN=16; UMILEN=10 ;;
        gex_3pv3_family|v3|3pv3)
            FILE=3M-february-2018.txt;    CBLEN=16; UMILEN=12 ;;
        gex_3pv4_gemx|v4|3pv4)
            FILE=3M-3pgex-may-2023.txt;   CBLEN=16; UMILEN=12 ;;
        gex_5pv3_gemx|5pv3)
            FILE=3M-5pgex-jan-2023.txt;   CBLEN=16; UMILEN=12 ;;
        gex_multiome_arc_v1|multiome|arc)
            # infer_10x_run.py prefers the gex_-prefixed copy; fall back to the
            # unprefixed name for older whitelist directories.
            FILE=gex_737K-arc-v1.txt
            [[ -s "$WL/$FILE" ]] || FILE=737K-arc-v1.txt
            CBLEN=16; UMILEN=12 ;;
        atac_*|atac)
            die "Chemistry '$NAME' is ATAC; STARsolo 10x handles gene expression only." ;;
        *)
            die "Unknown chemistry '$NAME'. Valid: gex_3pv1, gex_3pv2_or_5pv1v2, gex_3pv3_family, gex_3pv4_gemx, gex_5pv3_gemx, gex_multiome_arc_v1 (aliases: v1, v2, v3, v4, 5pv3, multiome)." ;;
    esac

    [[ -s "$WL/$FILE" ]] || die "Whitelist for '$NAME' not found: $WL/$FILE"
    echo "$WL/$FILE|$CBLEN|$UMILEN"
}

# Measures read lengths from the subsampled test reads.
# Prints: "R1LEN|R2LEN|R1DIS"
_10x_read_metrics() {
    local R1LEN R2LEN R1DIS
    R1LEN=$(awk 'NR%4==2' test.R1.fastq | awk '{sum+=length($0)} END {printf "%d\n",sum/NR+0.5}')
    R2LEN=$(awk 'NR%4==2' test.R2.fastq | awk '{sum+=length($0)} END {printf "%d\n",sum/NR+0.5}')
    R1DIS=$(awk 'NR%4==2' test.R1.fastq | awk '{print length($0)}' | sort | uniq -c | wc -l)

    echo "$R1LEN|$R2LEN|$R1DIS"
}

# Determines the correct barcode whitelist.
#   $1: WL directory
# Prints: "BC|NBC1|NBC2|NBC3|NBC4|NBC5|NBCA"
_10x_determine_whitelist() {
    local WL=$1

    local NBC1 NBC2 NBC3 NBC4 NBC5 NBCA
    NBC1=$(awk 'NR%4==2' test.R1.fastq | cut -c-14 | grep -F -f "$WL/737K-april-2014_rc.txt" | wc -l)
    NBC2=$(awk 'NR%4==2' test.R1.fastq | cut -c-16 | grep -F -f "$WL/737K-august-2016.txt"   | wc -l)
    NBC3=$(awk 'NR%4==2' test.R1.fastq | cut -c-16 | grep -F -f "$WL/3M-february-2018.txt"   | wc -l)
    NBC4=$(awk 'NR%4==2' test.R1.fastq | cut -c-16 | grep -F -f "$WL/3M-3pgex-may-2023.txt"  | wc -l)
    NBC5=$(awk 'NR%4==2' test.R1.fastq | cut -c-16 | grep -F -f "$WL/3M-5pgex-jan-2023.txt"  | wc -l)
    NBCA=$(awk 'NR%4==2' test.R1.fastq | cut -c-16 | grep -F -f "$WL/737K-arc-v1.txt"        | wc -l)

    local BC=""
    if   (( NBC3 > 50000 )); then BC="$WL/3M-february-2018.txt"
    elif (( NBC2 > 50000 )); then BC="$WL/737K-august-2016.txt"
    elif (( NBCA > 50000 )); then BC="$WL/737K-arc-v1.txt"
    elif (( NBC1 > 50000 )); then BC="$WL/737K-april-2014_rc.txt"
    elif (( NBC4 > 50000 )); then BC="$WL/3M-3pgex-may-2023.txt"
    elif (( NBC5 > 50000 )); then BC="$WL/3M-5pgex-jan-2023.txt"
    else
        die "No whitelist matched 200,000 random barcodes! Counts: v1=$NBC1 v2=$NBC2 v3=$NBC3 v4-3p=$NBC4 v4-5p=$NBC5 multiome=$NBCA"
    fi

    echo "$BC|$NBC1|$NBC2|$NBC3|$NBC4|$NBC5|$NBCA"
}

# Maps a whitelist file name to its CB/UMI geometry.
#   $1: BC (whitelist path)
# Prints: "CBLEN|UMILEN"
_10x_whitelist_geometry() {
    local BC=$1 CBLEN UMILEN

    case "$BC" in
        *3M-february-2018.txt|*737K-arc-v1.txt|*3M-3pgex-may-2023.txt|*3M-5pgex-jan-2023.txt)
            CBLEN=16; UMILEN=12 ;;
        *737K-august-2016.txt)
            CBLEN=16; UMILEN=10 ;;
        *737K-april-2014_rc.txt)
            CBLEN=14; UMILEN=10 ;;
        *)
            die "Cannot determine CB/UMI geometry for whitelist: $BC" ;;
    esac

    echo "$CBLEN|$UMILEN"
}

# Checks read lengths and reconciles barcode/UMI parameters against R1.
#   $1: R1DIS  $2: R1LEN  $3: R2LEN  $4: CBLEN  $5: UMILEN
#   $6: SKIPLEN – "True" downgrades the fatal sanity checks to warnings
# Prints: "PAIRED|CBLEN|UMILEN"
_10x_check_read_lengths() {
    local R1DIS=$1 R1LEN=$2 R2LEN=$3 CBLEN=$4 UMILEN=$5 SKIPLEN=${6:-False}
    local PAIRED=False

    # With --skip-length-checks the same problems are reported but not fatal.
    local ONBAD=die
    [[ $SKIPLEN == "True" ]] && ONBAD=log_warn

    if (( R1DIS > 1 && R1LEN <= 30 )); then
        $ONBAD "Read 1 (barcode) has varying length; possibly quality-trimmed."
    elif (( R1LEN < 24 )); then
        $ONBAD "Read 1 (barcode) is less than 24 bp. Check FASTQ files."
    elif (( R2LEN < 40 )); then
        $ONBAD "Read 2 (biological read) is less than 40 bp. Check FASTQ files."
    fi

    (( R1LEN > 50 )) && PAIRED=True

    local BCUMI=$((CBLEN + UMILEN))
    if (( BCUMI > R1LEN )); then
        local NEWUMI=$((R1LEN - CBLEN))
        # Always fatal: STAR cannot run without at least a 1 bp UMI.
        (( NEWUMI > 0 )) || die "R1 length ($R1LEN) leaves no room for a UMI after the ${CBLEN} bp barcode."
        log_warn "R1 length ($R1LEN) < barcode+UMI ($BCUMI). Setting UMI length to $NEWUMI."
        UMILEN=$NEWUMI
    elif (( BCUMI < R1LEN )); then
        log_warn "R1 length ($R1LEN) > barcode+UMI ($BCUMI)."
    fi

    echo "$PAIRED|$CBLEN|$UMILEN"
}

# Succeeds if the whitelist is used by a 5' chemistry. Every other
# whitelist (v1, v3, v4-3p, multiome) is 3'-only.
#   $1: BC (whitelist path)
_10x_wl_can_be_5p() {
    case "$(basename "$1")" in
        737K-august-2016.txt|3M-5pgex-jan-2023.txt) return 0 ;;
        *)                                          return 1 ;;
    esac
}

# Determines strand specificity via test alignments.
#   $1: REF  $2: CBLEN  $3: UMILEN  $4: CPUS  $5: BC
# Prints: "STRAND|PCTFWD|PCTREV"
_10x_determine_strand() {
    local REF=$1 CBLEN=$2 UMILEN=$3 CPUS=$4 BC=$5
    local STRAND=Forward

    log_info "Running test alignment (Forward) …"
    STAR --runThreadN "$CPUS" --genomeDir "$REF" --readFilesIn test.R2.fastq test.R1.fastq \
        --runDirPerm All_RWX --outSAMtype None \
        --soloType CB_UMI_Simple --soloCBwhitelist "$BC" --soloBarcodeReadLength 0 \
        --soloCBlen "$CBLEN" --soloUMIstart $((CBLEN + 1)) --soloUMIlen "$UMILEN" \
        --soloStrand Forward --genomeLoad LoadAndKeep \
        --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIfiltering MultiGeneUMI_CR --soloCellFilter EmptyDrops_CR \
        --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
        --soloFeatures Gene GeneFull \
        --soloOutFileNames test_forward/ features.tsv barcodes.tsv matrix.mtx &>/dev/null

    log_info "Running test alignment (Reverse) …"
    STAR --runThreadN "$CPUS" --genomeDir "$REF" --readFilesIn test.R2.fastq test.R1.fastq \
        --runDirPerm All_RWX --outSAMtype None \
        --soloType CB_UMI_Simple --soloCBwhitelist "$BC" --soloBarcodeReadLength 0 \
        --soloCBlen "$CBLEN" --soloUMIstart $((CBLEN + 1)) --soloUMIlen "$UMILEN" \
        --soloStrand Reverse --genomeLoad LoadAndKeep \
        --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIfiltering MultiGeneUMI_CR --soloCellFilter EmptyDrops_CR \
        --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
        --soloFeatures Gene GeneFull \
        --soloOutFileNames test_reverse/ features.tsv barcodes.tsv matrix.mtx &>/dev/null

    local PCTFWD PCTREV
    PCTFWD=$(grep "Reads Mapped to GeneFull: Unique GeneFull" test_forward/GeneFull/Summary.csv | awk -F "," '{printf "%d\n",$2*100+0.5}')
    PCTREV=$(grep "Reads Mapped to GeneFull: Unique GeneFull" test_reverse/GeneFull/Summary.csv | awk -F "," '{printf "%d\n",$2*100+0.5}')

    (( PCTREV > PCTFWD )) && STRAND=Reverse

    if (( PCTREV < 50 && PCTFWD < 50 )); then
        log_warn "Low GeneFull mapping: forward=${PCTFWD}%, reverse=${PCTREV}%"
    fi

    echo "$STRAND|$PCTFWD|$PCTREV"
}

# Reconciles the strand test with the whitelist and decides single/paired-end.
# 3'-only whitelists are Forward by definition, so a Reverse call there (seen
# on single-nucleus 3' v3 data with a long R1) is overridden. Only 5' data
# (Reverse) is processed paired-end.
#   $1: BC  $2: STRAND  $3: PAIRED (R1 long enough)  $4: PCTFWD  $5: PCTREV
# Prints: "STRAND|PAIRED|NOTE"
_10x_resolve_layout() {
    local BC=$1 STRAND=$2 PAIRED=$3 PCTFWD=$4 PCTREV=$5 NOTE=""

    if [[ $STRAND == "Reverse" ]] && ! _10x_wl_can_be_5p "$BC"; then
        log_warn "Strand test favours Reverse (fwd=${PCTFWD}%, rev=${PCTREV}%), but $(basename "$BC") is a 3'-only whitelist; using Forward, single-end."
        STRAND=Forward
        NOTE="forced Forward: 3'-only whitelist"
    fi

    # 3' experiments are processed as single-end even with a long R1
    if [[ $STRAND == "Forward" ]]; then PAIRED=False; fi

    echo "$STRAND|$PAIRED|$NOTE"
}

# ---------- Velocyto handling -----------------------------------------------
# Stock STAR (2.7.10a-2.7.11b) breaks in Velocyto in two known cases (#27):
#  - paired-end reads with a mate outside the transcript make
#    alignToTranscriptMinOverlap read out of bounds; with the 3M GEM-X 5'
#    whitelist this segfaults while mapping;
#  - above 2^32 input reads the read index STAR keeps for Velocyto wraps, so
#    Velocyto segfaults after the other features are written, or finishes
#    with wrong counts (alexdobin/STAR#2710).
# The cellgeni build 2.7.10a_alpha_220818_velofix fixes both. With stock STAR
# the helpers below avoid the first case, keep the finished matrices when
# Velocyto crashes, and discard Velocyto output from runs above 2^32 reads.

# Succeeds if STAR is the build with the paired-end Velocyto fix.
_10x_star_has_pe_velocyto_fix() {
    [[ $(STAR --version 2>/dev/null) == *velofix* ]]
}

# Decides whether Velocyto stays in STAR_ARGS. Edits STAR_ARGS and sets the
# global VELO_STATUS, so never call it inside $(...).
#   $1: NOVELO  $2: PAIRED  $3: BC
_10x_velocyto_policy() {
    local NOVELO=$1 PAIRED=$2 BC=$3

    if ! has_solo_feature Velocyto "${STAR_ARGS[@]}"; then
        VELO_STATUS="off (not in --soloFeatures)"
    elif [[ $NOVELO == "True" ]]; then
        star_args_drop_feature Velocyto
        VELO_STATUS="off (--no-velocyto)"
        log_info "Velocyto disabled (--no-velocyto)."
    elif [[ $PAIRED == "True" && $(basename "$BC") == "3M-5pgex-jan-2023.txt" ]] \
            && ! _10x_star_has_pe_velocyto_fix; then
        star_args_drop_feature Velocyto
        VELO_STATUS="off (paired-end GEM-X 5': stock STAR crashes in Velocyto; fixed in 2.7.10a_alpha_220818_velofix)"
        log_warn "Velocyto disabled: this STAR ($(STAR --version 2>/dev/null)) crashes on paired-end data with the 3M-5pgex whitelist."
    else
        VELO_STATUS="on"
    fi
}

# Prints the number of input reads: from Log.final.out, or from the Solo line
# in Log.out (nReadsInput = reads + 1) if STAR stopped before Log.final.out.
_10x_n_reads_input() {
    local N=""
    if [[ -s Log.final.out ]]; then
        N=$(awk -F'|' '/Number of input reads/ { gsub(/[^0-9]/, "", $2); print $2; exit }' Log.final.out)
    fi
    if [[ -z $N && -s Log.out ]]; then
        N=$(grep -m1 -o "nReadsInput = [0-9]*" Log.out | grep -o "[0-9]*$") || N=""
        if [[ -n $N ]]; then N=$((N - 1)); fi
    fi
    echo "$N"
}

# Succeeds if every requested feature other than Velocyto has a non-empty raw
# matrix, plus Summary.csv and (when cells are filtered) a filtered matrix.
#   $1: SOLO_OUT
_10x_solo_outputs_complete() {
    local OUT=${1%/} FILTER f n=0
    FILTER=$(star_arg_first soloCellFilter CellRanger2.2)
    for f in $(star_arg_values soloFeatures); do
        case "$f" in Velocyto|VelocytoSimple) continue ;; esac
        [[ -s $OUT/$f/raw/matrix.mtx ]] || return 1
        if [[ $f != "SJ" ]]; then
            [[ -s $OUT/$f/Summary.csv ]] || return 1
            if [[ $FILTER != "None" ]]; then
                [[ -s $OUT/$f/filtered/matrix.mtx ]] || return 1
            fi
        fi
        n=$((n + 1))
    done
    (( n > 0 ))
}

# Succeeds if STAR was killed by a crash signal while mapping: the genome was
# loaded but Solo counting had not started.
#   $1: RC (STAR exit code)
_10x_died_while_mapping() {
    case "$1" in
        134|135|139) ;;     # SIGABRT, SIGBUS, SIGSEGV
        *) return 1 ;;
    esac
    [[ -s Log.out ]] || return 1
    grep -qF "Finished loading the genome" Log.out || return 1
    ! grep -qF "Starting Solo post-map" Log.out
}

# Writes velocyto_failed.txt (key=value lines), read by 'starsolo qc'.
#   $1: RC  $2: STAGE  $3: ACTION  $4…: extra key=value lines
_10x_velocyto_marker() {
    local RC=$1 STAGE=$2 ACTION=$3 l
    shift 3
    {
        echo "star_exit_code=$RC"
        echo "stage=$STAGE"
        echo "action=$ACTION"
        for l in "$@"; do echo "$l"; done
        echo "date=$(date '+%Y-%m-%dT%H:%M:%S')"
    } > velocyto_failed.txt
}

# Frees the strand test's shared-memory genome unless the run asked to keep
# it (--genomeLoad LoadAndKeep). Also used as the EXIT trap of run_10x.
#   $1: REF
_10x_release_genome() {
    if [[ $(star_arg_first genomeLoad LoadAndRemove) != "LoadAndKeep" ]]; then
        release_shared_genome "$1"
    fi
}

# Runs the main alignment and recovers from a Velocyto crash:
#  - crash in Velocyto after the other features were written: keep them,
#    remove the partial Velocyto output, write velocyto_failed.txt;
#  - crash while mapping with Velocyto on: rerun once without Velocyto.
# Anything else exits with STAR's code. Sets SALVAGE_RC to STAR's code if a
# requested BAM is missing after a salvage (STAR sorts the BAM after Solo).
#   $1: REF  $2: SOLO_OUT
_10x_align() {
    local REF=$1 OUT=${2%/} rc=0 VELO=False
    local -a EXTRA=()
    SALVAGE_RC=0
    if has_solo_feature Velocyto "${STAR_ARGS[@]}"; then VELO=True; fi

    # Logs left by the strand test must not be mistaken for this run's
    rm -f Log.out Log.progress.out Log.final.out SJ.out.tab

    run_star || rc=$?
    if (( rc == 0 )); then return 0; fi
    if [[ $VELO != "True" ]]; then exit "$rc"; fi

    if grep -qF "Starting Solo post-map for Velocyto" Log.out 2>/dev/null; then
        if ! _10x_solo_outputs_complete "$OUT"; then
            log_error "STAR crashed in Velocyto (exit $rc) and the other features are incomplete."
            exit "$rc"
        fi
        rm -rf -- "$OUT/Velocyto"
        EXTRA=("n_reads_input=$(_10x_n_reads_input)")
        if [[ $(star_arg_first outSAMtype None) == "BAM" && ! -s Aligned.sortedByCoord.out.bam && ! -s Aligned.out.bam ]]; then
            EXTRA+=("bam=missing")
            SALVAGE_RC=$rc
        fi
        _10x_velocyto_marker "$rc" solo_postmap velocyto_output_removed "${EXTRA[@]}"
        echo "Velocyto: FAILED (STAR exit $rc in Velocyto counting); other features kept, see velocyto_failed.txt" >> strand.txt
        log_warn "STAR crashed in Velocyto (exit $rc) after the other features were written; keeping them and removing the Velocyto output (see velocyto_failed.txt)."
        if (( SALVAGE_RC != 0 )); then
            log_error "The requested BAM was not written (STAR sorts it after Solo counting)."
        fi
        _10x_release_genome "$REF"
        return 0
    fi

    if _10x_died_while_mapping "$rc"; then
        log_warn "STAR crashed while mapping with Velocyto on (exit $rc); rerunning once without Velocyto."
        mv -f Log.out Log.attempt1.out
        rm -rf -- "$OUT" Log.progress.out SJ.out.tab Unmapped.out.mate1 Unmapped.out.mate2 Aligned.*
        star_args_drop_feature Velocyto
        _10x_velocyto_marker "$rc" mapping rerun_without_velocyto
        echo "Velocyto: FAILED (STAR exit $rc while mapping); rerun without Velocyto, see velocyto_failed.txt" >> strand.txt
        rc=0
        run_star || rc=$?
        if (( rc == 0 )); then return 0; fi
        log_error "STAR failed again without Velocyto (exit $rc)."
    fi
    exit "$rc"
}

# Above 2^32 input reads stock STAR's Solo read index wraps
# (alexdobin/STAR#2710): Velocyto counts are wrong even when STAR finishes,
# so the output is removed; BAM CB/UB tags come from the same index and are
# flagged. STAR builds with the fix (alexdobin/STAR#2711, and the image's
# patched 2.7.10a) store the index in 64 bits and say so in Log.out.
#   $1: SOLO_OUT
_10x_check_read_index_limit() {
    local OUT=${1%/} N MSG
    N=$(_10x_n_reads_input)
    if ! [[ $N =~ ^[0-9]+$ ]] || (( N < 4294967296 )); then return 0; fi
    if grep -qF "Solo read index will be stored in 64 bits" Log.out 2>/dev/null; then
        log_info "$N input reads >= 2^32: STAR stored the Solo read index in 64 bits; Velocyto output kept."
        return 0
    fi

    if [[ -d $OUT/Velocyto ]]; then
        rm -rf -- "$OUT/Velocyto"
        _10x_velocyto_marker 0 read_index_limit velocyto_output_removed "n_reads_input=$N"
        MSG="$N input reads >= 2^32: STAR's read index wraps, so Velocyto counts are wrong; Velocyto output removed (see velocyto_failed.txt)."
        log_warn "$MSG"
        echo "Velocyto: REMOVED ($MSG)" >> strand.txt
    fi
    if [[ -s Aligned.sortedByCoord.out.bam ]]; then
        MSG="$N input reads >= 2^32: STAR's read index wraps, so BAM CB/UB tags are unreliable."
        log_warn "$MSG"
        echo "WARNING: $MSG" >> strand.txt
    fi
    return 0
}

# Writes a human-readable configuration summary.
_10x_write_config() {
    local TAG=$1 PAIRED=$2 STRAND=$3 PCTFWD=$4 PCTREV=$5
    local BC=$6 NBC1=$7 NBC2=$8 NBC3=$9 NBCA=${10} NBC4=${11} NBC5=${12}
    local CBLEN=${13} UMILEN=${14} GZIP=${15} R1=${16} R2=${17} STRAND_NOTE=${18:-}
    local FEATURES=${19:-} VELO=${20:-}

    log_info "STARsolo 10x run configuration:"
    {
        echo "============================================================================="
        echo "Sample: $TAG"
        echo "Paired-end mode: $PAIRED"
        echo "Strand (Forward=3', Reverse=5'): $STRAND  (%fwd=$PCTFWD, %rev=$PCTREV)${STRAND_NOTE:+  [$STRAND_NOTE]}"
        echo "CB whitelist: $BC"
        echo "  Matches /200k: v3=$NBC3 v2=$NBC2 v1=$NBC1 v4-3p=$NBC4 v4-5p=$NBC5 multiome=$NBCA (NA = chemistry given via --wl)"
        echo "CB length: $CBLEN    UMI length: $UMILEN"
        echo "Compression: ${GZIP:-none}"
        echo "Solo features: $FEATURES"
        echo "Velocyto: $VELO"
        echo "-----------------------------------------------------------------------------"
        echo "R1: $R1"
        echo "R2: $R2"
        echo "============================================================================="
    } | tee strand.txt
}

# ---------- Main 10x workflow -----------------------------------------------

# Entry point called by bin/starsolo.
#   Arguments come from the CLI parser.
run_10x() {
    local FQDIR=$1 TAG=$2 CPUS=$3 REF=$4 WL=$5 BAM=$6 CHEM=${7:-} SKIPLEN=${8:-False} NOVELO=${9:-False}

    local FQDIR_ABS
    FQDIR_ABS=$(validate_fqdir "$FQDIR") || exit 1
    TAG=$(validate_tag "$TAG") || exit 1

    # Convert paths to absolute before cd
    REF=$(readlink -f "$REF")
    WL=$(readlink -f "$WL")

    # 1. Discover FASTQs
    log_info "Discovering FASTQ files for sample '$TAG' …"
    local fq_info R1 R2
    fq_info=$(find_fastq_files "$FQDIR_ABS" "$TAG")
    IFS='|' read -r R1 R2 <<< "$fq_info"

    # 2. Compression
    local comp_info GZIP ZCMD
    comp_info=$(check_compression "$FQDIR_ABS" "$TAG")
    IFS='|' read -r GZIP ZCMD <<< "$comp_info"

    # Create output directory (only once the FASTQs are found) and enter it
    mkdir -p "$TAG" && cd "$TAG" || exit

    # 3. Chemistry: subsample reads, then either use the given chemistry or detect it
    log_info "Subsampling reads …"
    _10x_extract_test_reads "$R1" "$R2" "$ZCMD"

    local rm_info R1LEN R2LEN R1DIS
    rm_info=$(_10x_read_metrics)
    IFS='|' read -r R1LEN R2LEN R1DIS <<< "$rm_info"

    local BC NBC1 NBC2 NBC3 NBC4 NBC5 NBCA CBLEN UMILEN
    if [[ -n $CHEM ]]; then
        local chem_info
        chem_info=$(_10x_resolve_chemistry "$CHEM" "$WL")
        IFS='|' read -r BC CBLEN UMILEN <<< "$chem_info"
        NBC1=NA; NBC2=NA; NBC3=NA; NBC4=NA; NBC5=NA; NBCA=NA
        log_info "Using chemistry '$CHEM' → $BC (skipping whitelist matching)"
    else
        log_info "Determining barcode whitelist …"
        local wl_info geom_info
        wl_info=$(_10x_determine_whitelist "$WL")
        IFS='|' read -r BC NBC1 NBC2 NBC3 NBC4 NBC5 NBCA <<< "$wl_info"
        geom_info=$(_10x_whitelist_geometry "$BC")
        IFS='|' read -r CBLEN UMILEN <<< "$geom_info"
    fi

    # 4. Read length checks
    local rl_info PAIRED
    [[ $SKIPLEN == "True" ]] && log_warn "Read-length validation disabled (--skip-length-checks)."
    rl_info=$(_10x_check_read_lengths "$R1DIS" "$R1LEN" "$R2LEN" "$CBLEN" "$UMILEN" "$SKIPLEN")
    IFS='|' read -r PAIRED CBLEN UMILEN <<< "$rl_info"

    # 5. Strand specificity. The test alignments leave the genome in shared
    #    memory (LoadAndKeep) for the main run; free it if we exit before
    #    a LoadAndRemove run has done so.
    STAR_ARGS=()
    trap "_10x_release_genome $(printf %q "$REF")" EXIT
    log_info "Determining strand specificity …"
    local strand_info layout_info STRAND PCTFWD PCTREV STRAND_NOTE
    strand_info=$(_10x_determine_strand "$REF" "$CBLEN" "$UMILEN" "$CPUS" "$BC")
    IFS='|' read -r STRAND PCTFWD PCTREV <<< "$strand_info"
    rm -rf test.R?.fastq test_forward test_reverse
    layout_info=$(_10x_resolve_layout "$BC" "$STRAND" "$PAIRED" "$PCTFWD" "$PCTREV")
    IFS='|' read -r STRAND PAIRED STRAND_NOTE <<< "$layout_info"

    # 6. STAR arguments
    local SOLOFILENAMES="output/"

    if [[ $PAIRED == "True" ]]; then
        # 5' R1 carries CB + UMI + 13 bp TSO before the cDNA
        merge_star_args --runThreadN "$CPUS" --genomeDir "$REF" \
            --readFilesIn "$R1" "$R2" --runDirPerm All_RWX $GZIP $BAM \
            --soloBarcodeMate 1 --clip5pNbases $((CBLEN + UMILEN + 13)) 0 \
            --soloType CB_UMI_Simple --soloCBwhitelist "$BC" --soloCBstart 1 \
            --soloCBlen "$CBLEN" --soloUMIstart $((CBLEN + 1)) --soloUMIlen "$UMILEN" \
            --soloStrand Forward \
            --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
            --soloUMIfiltering MultiGeneUMI_CR --soloCellFilter EmptyDrops_CR \
            --outFilterScoreMin 30 --genomeLoad LoadAndRemove \
            --soloFeatures Gene GeneFull Velocyto \
            --soloOutFileNames "$SOLOFILENAMES" features.tsv barcodes.tsv matrix.mtx \
            --soloMultiMappers EM
    else
        merge_star_args --runThreadN "$CPUS" --genomeDir "$REF" \
            --readFilesIn "$R2" "$R1" --runDirPerm All_RWX $GZIP $BAM \
            --soloType CB_UMI_Simple --soloCBwhitelist "$BC" --soloBarcodeReadLength 0 \
            --soloCBlen "$CBLEN" --soloUMIstart $((CBLEN + 1)) --soloUMIlen "$UMILEN" \
            --soloStrand "$STRAND" \
            --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
            --soloUMIfiltering MultiGeneUMI_CR --soloCellFilter EmptyDrops_CR \
            --clipAdapterType CellRanger4 --outFilterScoreMin 30 --genomeLoad LoadAndRemove \
            --soloFeatures Gene GeneFull Velocyto \
            --soloOutFileNames "$SOLOFILENAMES" features.tsv barcodes.tsv matrix.mtx \
            --soloMultiMappers EM
    fi
    _10x_velocyto_policy "$NOVELO" "$PAIRED" "$BC"

    # 7. Log configuration
    _10x_write_config "$TAG" "$PAIRED" "$STRAND" "$PCTFWD" "$PCTREV" \
        "$BC" "$NBC1" "$NBC2" "$NBC3" "$NBCA" "$NBC4" "$NBC5" \
        "$CBLEN" "$UMILEN" "$GZIP" "$R1" "$R2" "$STRAND_NOTE" \
        "$(star_arg_values soloFeatures)" "$VELO_STATUS"

    # 8. Run STAR (recovers from a Velocyto crash)
    log_info "Running STARsolo alignment …"
    _10x_align "$REF" "$SOLOFILENAMES"
    case "$(star_arg_first genomeLoad LoadAndRemove)" in
        LoadAndRemove|LoadAndKeep) trap - EXIT ;;   # STAR handled the shared genome
    esac

    # 9. Checks & compress
    _10x_check_read_index_limit "$SOLOFILENAMES"
    process_output_files "$SOLOFILENAMES" "$CPUS"
    if (( SALVAGE_RC != 0 )); then exit "$SALVAGE_RC"; fi
}

# ---------- Subcommand help -------------------------------------------------

show_10x_help() {
    cat <<'EOF'
Usage: starsolo 10x <fastq_dir> <sample_id> [options]

Automatically detects 10x chemistry (v1–v4, multiome), strand specificity,
and paired-end vs single-end mode.

Required:
  <fastq_dir>               Directory containing FASTQ files
  <sample_id>               Sample identifier (used to match FASTQ filenames)

Options:
  -s, --species <name>      Species name (e.g. human, mouse).
                             Resolves reference via $REF_BASE/<species>/2020A/index
  -r, --ref <path>          Explicit STAR reference index (overrides --species)
  -w, --whitelist-dir <dir> Barcode whitelist directory  [default: from config]
      --wl <chemistry>      Use this chemistry instead of detecting it. Names
                             chemistry ids:
                               gex_3pv1              (alias v1)
                               gex_3pv2_or_5pv1v2    (alias v2)
                               gex_3pv3_family       (alias v3)
                               gex_3pv4_gemx         (alias v4)
                               gex_5pv3_gemx         (alias 5pv3)
                               gex_multiome_arc_v1   (alias multiome, arc)
      --skip-length-checks  Do not abort on unexpected read lengths (short R1,
                             short R2, or variable-length R1); the same problems
                             are reported as warnings and the run continues.
                             A UMI-less R1 remains a fatal error.
      --no-velocyto         Do not count Velocyto (spliced/unspliced) matrices.
                             With stock STAR (not the image's ..._velofix
                             build) Velocyto is also dropped for paired-end
                             GEM-X 5' data, and removed from runs with
                             >= 2^32 reads. If STAR crashes in Velocyto anyway,
                             Gene/GeneFull are kept and velocyto_failed.txt is
                             written.
  -c, --cpus <N>            Number of threads             [default: 16]
  --bam                     Output sorted BAM file
  --no-bam                  Do not output BAM (default)
  --                        Pass the rest to STAR; replaces built-in flags
  -h, --help                Show this help
EOF
}
