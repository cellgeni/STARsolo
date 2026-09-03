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

# Determines strand specificity via test alignments.
#   $1: REF  $2: CBLEN  $3: UMILEN  $4: CPUS  $5: BC  $6: PAIRED
# Prints: "STRAND|PCTFWD|PCTREV|PAIRED"
_10x_determine_strand() {
    local REF=$1 CBLEN=$2 UMILEN=$3 CPUS=$4 BC=$5 PAIRED=$6
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

    # 3' paired-end experiments should still be processed as single-end
    [[ $STRAND == "Forward" && $PAIRED == "True" ]] && PAIRED=False

    echo "$STRAND|$PCTFWD|$PCTREV|$PAIRED"
}

# Writes a human-readable configuration summary.
_10x_write_config() {
    local TAG=$1 PAIRED=$2 STRAND=$3 PCTFWD=$4 PCTREV=$5
    local BC=$6 NBC1=$7 NBC2=$8 NBC3=$9 NBCA=${10} NBC4=${11} NBC5=${12}
    local CBLEN=${13} UMILEN=${14} GZIP=${15} R1=${16} R2=${17}

    log_info "STARsolo 10x run configuration:"
    {
        echo "============================================================================="
        echo "Sample: $TAG"
        echo "Paired-end mode: $PAIRED"
        echo "Strand (Forward=3', Reverse=5'): $STRAND  (%fwd=$PCTFWD, %rev=$PCTREV)"
        echo "CB whitelist: $BC"
        echo "  Matches /200k: v3=$NBC3 v2=$NBC2 v1=$NBC1 v4-3p=$NBC4 v4-5p=$NBC5 multiome=$NBCA (NA = chemistry given via --wl)"
        echo "CB length: $CBLEN    UMI length: $UMILEN"
        echo "Compression: ${GZIP:-none}"
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
    local FQDIR=$1 TAG=$2 CPUS=$3 REF=$4 WL=$5 BAM=$6 CHEM=${7:-} SKIPLEN=${8:-False}

    local FQDIR_ABS
    FQDIR_ABS=$(validate_fqdir "$FQDIR") || exit 1
    TAG=$(validate_tag "$TAG") || exit 1

    # Convert paths to absolute before cd
    REF=$(readlink -f "$REF")
    WL=$(readlink -f "$WL")

    # Create output directory and enter it
    mkdir -p "$TAG" && cd "$TAG" || exit

    # 1. Discover FASTQs
    log_info "Discovering FASTQ files for sample '$TAG' …"
    local fq_info R1 R2
    fq_info=$(find_fastq_files "$FQDIR_ABS" "$TAG")
    IFS='|' read -r R1 R2 <<< "$fq_info"

    # 2. Compression
    local comp_info GZIP ZCMD
    comp_info=$(check_compression "$FQDIR_ABS" "$TAG")
    IFS='|' read -r GZIP ZCMD <<< "$comp_info"

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

    # 5. Strand specificity
    log_info "Determining strand specificity …"
    local strand_info STRAND PCTFWD PCTREV
    strand_info=$(_10x_determine_strand "$REF" "$CBLEN" "$UMILEN" "$CPUS" "$BC" "$PAIRED")
    IFS='|' read -r STRAND PCTFWD PCTREV PAIRED <<< "$strand_info"

    # 6. Log configuration
    _10x_write_config "$TAG" "$PAIRED" "$STRAND" "$PCTFWD" "$PCTREV" \
        "$BC" "$NBC1" "$NBC2" "$NBC3" "$NBCA" "$NBC4" "$NBC5" \
        "$CBLEN" "$UMILEN" "$GZIP" "$R1" "$R2"

    # 7. Run STAR
    local SOLOFILENAMES="output/"
    log_info "Running STARsolo alignment …"

    if [[ $PAIRED == "True" ]]; then
        STAR --runThreadN "$CPUS" --genomeDir "$REF" \
            --readFilesIn "$R1" "$R2" --runDirPerm All_RWX $GZIP $BAM \
            --soloBarcodeMate 1 --clip5pNbases 39 0 \
            --soloType CB_UMI_Simple --soloCBwhitelist "$BC" --soloCBstart 1 \
            --soloCBlen "$CBLEN" --soloUMIstart $((CBLEN + 1)) --soloUMIlen "$UMILEN" \
            --soloStrand Forward \
            --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
            --soloUMIfiltering MultiGeneUMI_CR --soloCellFilter EmptyDrops_CR \
            --outFilterScoreMin 30 --genomeLoad LoadAndRemove \
            --soloFeatures Gene GeneFull Velocyto \
            --soloOutFileNames "$SOLOFILENAMES" features.tsv barcodes.tsv matrix.mtx \
            --soloMultiMappers EM \
            "${STAR_EXTRA_ARGS[@]}"
    else
        STAR --runThreadN "$CPUS" --genomeDir "$REF" \
            --readFilesIn "$R2" "$R1" --runDirPerm All_RWX $GZIP $BAM \
            --soloType CB_UMI_Simple --soloCBwhitelist "$BC" --soloBarcodeReadLength 0 \
            --soloCBlen "$CBLEN" --soloUMIstart $((CBLEN + 1)) --soloUMIlen "$UMILEN" \
            --soloStrand "$STRAND" \
            --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
            --soloUMIfiltering MultiGeneUMI_CR --soloCellFilter EmptyDrops_CR \
            --clipAdapterType CellRanger4 --outFilterScoreMin 30 --genomeLoad LoadAndRemove \
            --soloFeatures Gene GeneFull Velocyto \
            --soloOutFileNames "$SOLOFILENAMES" features.tsv barcodes.tsv matrix.mtx \
            --soloMultiMappers EM \
            "${STAR_EXTRA_ARGS[@]}"
    fi

    # 8. Cleanup & compress
    rm -rf test.R?.fastq test_forward test_reverse
    process_output_files "$SOLOFILENAMES"
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
  -c, --cpus <N>            Number of threads             [default: 16]
  --bam                     Output sorted BAM file
  --no-bam                  Do not output BAM (default)
  --                        Pass all following arguments verbatim to STAR
  -h, --help                Show this help
EOF
}
