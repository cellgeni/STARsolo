#!/bin/bash
# =============================================================================
# STARsolo CLI - Shared library functions
# =============================================================================
# Sourced by bin/starsolo and all platform scripts. Never executed directly.
# =============================================================================

# ---------- Logging --------------------------------------------------------

log_info()  { echo "[INFO]  $*" >&2; }
log_warn()  { echo "[WARN]  $*" >&2; }
log_error() { echo "[ERROR] $*" >&2; }
die()       { log_error "$@"; exit 1; }

# ---------- FASTQ file discovery -------------------------------------------

# Finds paired-end FASTQ files based on common naming conventions.
#   $1 : FQDIR  – directory containing FASTQ files
#   $2 : TAG    – sample identifier
# Prints: "R1_FILES|R2_FILES" (comma-separated within each group)
find_fastq_files() {
    local FQDIR=$1 TAG=$2
    local R1="" R2=""

    if [[ $(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_1\\.f.*q") != "" ]]; then
        R1=$(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_1\\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g")
        R2=$(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_2\\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g")
    elif [[ $(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "R1\\.f.*q") != "" ]]; then
        R1=$(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "R1\\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g")
        R2=$(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "R2\\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g")
    elif [[ $(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_R1_.*\\.f.*q") != "" ]]; then
        R1=$(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_R1_.*\\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g")
        R2=$(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_R2_.*\\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g")
    elif [[ $(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_f1\\.f.*q") != "" ]]; then
        R1=$(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_f1\\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g")
        R2=$(find "$FQDIR"/* | grep -P "/$TAG[\\/\\._]" | grep "_r2\\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g")
    else
        die "No appropriate FASTQ files found for sample '$TAG' in '$FQDIR'. Check file names and directory."
    fi

    echo "$R1|$R2"
}

# ---------- Compression detection ------------------------------------------

# Detects GZIP/BZIP2 compression for a sample's FASTQ files.
#   $1 : FQDIR
#   $2 : TAG
# Prints: "STAR_FLAG|DECOMPRESS_CMD"
check_compression() {
    local FQDIR=$1 TAG=$2
    local STAR_FLAG="" ZCMD="cat"

    if [[ $(find "$FQDIR"/* | grep -P "$TAG[\\/\\._]" | grep "\\.gz$") != "" ]]; then
        STAR_FLAG="--readFilesCommand zcat"
        ZCMD="zcat"
    elif [[ $(find "$FQDIR"/* | grep -P "$TAG[\\/\\._]" | grep "\\.bz2$") != "" ]]; then
        STAR_FLAG="--readFilesCommand bzcat"
        ZCMD="bzcat"
    fi

    echo "$STAR_FLAG|$ZCMD"
}

# ---------- STAR arguments --------------------------------------------------

# Builds the global STAR_ARGS array: the built-in arguments ("$@") minus every
# --param that the passthrough (STAR_EXTRA_ARGS, everything after '--') sets
# again, together with that param's values, followed by the passthrough.
# STAR rejects a repeated parameter, so this is what makes '-- --param value'
# replace a built-in flag. STAR values never start with "--".
merge_star_args() {
    local -A override=()
    local a skip=0
    for a in ${STAR_EXTRA_ARGS[@]+"${STAR_EXTRA_ARGS[@]}"}; do
        if [[ $a == --* ]]; then override[$a]=1; fi
    done

    STAR_ARGS=()
    for a in "$@"; do
        if [[ $a == --* ]]; then
            skip=0
            if [[ -n ${override[$a]:-} ]]; then skip=1; fi
        fi
        if (( skip == 0 )); then STAR_ARGS+=("$a"); fi
    done
    STAR_ARGS+=(${STAR_EXTRA_ARGS[@]+"${STAR_EXTRA_ARGS[@]}"})
}

# Succeeds if FEATURE is one of the --soloFeatures values in the given args.
#   $1 : FEATURE (e.g. Velocyto)   $2… : STAR arguments
has_solo_feature() {
    local FEAT=$1 a in=0
    shift
    for a in "$@"; do
        if [[ $a == --* ]]; then
            in=0
            if [[ $a == --soloFeatures ]]; then in=1; fi
        elif (( in )) && [[ $a == "$FEAT" ]]; then
            return 0
        fi
    done
    return 1
}

# Removes FEATURE from the --soloFeatures values in STAR_ARGS.
#   $1 : FEATURE
star_args_drop_feature() {
    local FEAT=$1 a in=0 left=0
    local -a out=()
    for a in "${STAR_ARGS[@]}"; do
        if [[ $a == --* ]]; then
            in=0
            if [[ $a == --soloFeatures ]]; then in=1; fi
        elif (( in )); then
            if [[ $a == "$FEAT" ]]; then continue; fi
            left=$((left + 1))
        fi
        out+=("$a")
    done
    (( left > 0 )) || die "Removing $FEAT would leave --soloFeatures empty."
    STAR_ARGS=("${out[@]}")
}

# Prints the first value of --NAME in STAR_ARGS, or DEFAULT if it is absent.
#   $1 : NAME (without dashes)   $2 : DEFAULT
star_arg_first() {
    local NAME=--$1 DEFAULT=${2:-} i
    local -a A=(${STAR_ARGS[@]+"${STAR_ARGS[@]}"})
    for (( i = 0; i + 1 < ${#A[@]}; i++ )); do
        if [[ ${A[i]} == "$NAME" ]]; then
            echo "${A[i+1]}"
            return 0
        fi
    done
    echo "$DEFAULT"
}

# Prints all values of --NAME in STAR_ARGS, space-separated.
#   $1 : NAME (without dashes)
star_arg_values() {
    local NAME=--$1 a in=0
    local -a vals=()
    for a in ${STAR_ARGS[@]+"${STAR_ARGS[@]}"}; do
        if [[ $a == --* ]]; then
            in=0
            if [[ $a == "$NAME" ]]; then in=1; fi
        elif (( in )); then
            vals+=("$a")
        fi
    done
    echo "${vals[*]}"
}

# Removes STAR's temporary directory (_STARtmp, or --outTmpDir if given).
# After a crash it can hold hundreds of GB.
remove_star_tmp() {
    local TMP
    TMP=$(star_arg_first outTmpDir -)
    if [[ $TMP == - ]]; then TMP=_STARtmp; fi
    case "${TMP%/}" in
        ""|/|.|..) log_warn "Not removing suspicious temp dir '$TMP'." ;;
        *)         rm -rf -- "$TMP" ;;
    esac
}

# Runs STAR with STAR_ARGS. On failure removes STAR's temp directory and
# returns STAR's exit code (128+N if it was killed by signal N).
run_star() {
    local rc=0
    STAR "${STAR_ARGS[@]}" || rc=$?
    if (( rc != 0 )); then
        log_error "STAR exited with code $rc."
        remove_star_tmp
    fi
    return "$rc"
}

# Removes a genome left in shared memory by --genomeLoad LoadAndKeep when no
# later LoadAndRemove run got to do it (e.g. the main run crashed). Processes
# still attached keep working; the segment is freed when they detach.
#   $1 : REF
release_shared_genome() {
    local REF=$1 D="_genome_remove.$$"
    mkdir -p "$D"
    STAR --genomeDir "$REF" --genomeLoad Remove --outFileNamePrefix "$D/" &>/dev/null || true
    rm -rf -- "$D"
}

# ---------- Post-alignment processing --------------------------------------

# Compresses outputs and optionally indexes BAM.
#   $1 : OUTPUT_DIR   – STARsolo output directory (e.g. "output/")
#   $2 : CPUS         – thread budget shared by samtools / pbzip2
#   No singularity; tools expected in PATH.
process_output_files() {
    local OUTPUT_DIR=${1:-output} CPUS=${2:-1}
    [[ $CPUS =~ ^[1-9][0-9]*$ ]] || die "process_output_files: invalid CPU count '$CPUS'."

    # Index BAM if it exists (-@ counts threads in addition to the main one)
    if [[ -s "Aligned.sortedByCoord.out.bam" ]]; then
        log_info "Indexing BAM file …"
        samtools index -@ $(( CPUS > 1 ? CPUS - 1 : 0 )) Aligned.sortedByCoord.out.bam
    fi

    # Compress unmapped reads. pbzip2 takes "-p<N>" with no space, and
    # otherwise uses every core on the node, so split CPUS between the mates.
    local -a MATES=()
    local m
    for m in Unmapped.out.mate1 Unmapped.out.mate2; do
        if [[ -s "$m" ]]; then MATES+=("$m"); fi
    done

    if (( ${#MATES[@]} == 1 )); then
        pbzip2 -f -9 -p"$CPUS" "${MATES[0]}" || die "pbzip2 failed on ${MATES[0]}."
    elif (( ${#MATES[@]} == 2 && CPUS == 1 )); then
        pbzip2 -f -9 -p1 "${MATES[0]}" || die "pbzip2 failed on ${MATES[0]}."
        pbzip2 -f -9 -p1 "${MATES[1]}" || die "pbzip2 failed on ${MATES[1]}."
    elif (( ${#MATES[@]} == 2 )); then
        local pid1 pid2 rc=0
        pbzip2 -f -9 -p"$(( (CPUS + 1) / 2 ))" "${MATES[0]}" & pid1=$!
        pbzip2 -f -9 -p"$(( CPUS / 2 ))" "${MATES[1]}" & pid2=$!
        wait "$pid1" || rc=$?
        wait "$pid2" || rc=$?
        (( rc == 0 )) || die "pbzip2 failed (exit $rc) compressing unmapped reads."
    fi

    # Gzip all matrix / barcode / feature files. -type f skips the SJ feature's
    # features.tsv symlink (gzip refuses symlinks); -f allows a rerun.
    find "$OUTPUT_DIR" -type f \( -name "*.tsv" -o -name "*.mtx" \) -exec gzip -f {} +

    log_info "ALL DONE!"
}

# ---------- BAM flag builder -----------------------------------------------

# Resolves the BAM flags for STAR depending on --bam / --no-bam.
#   $1 : "bam" or "nobam"
# Prints the STAR flags string.
resolve_bam_flags() {
    local MODE=$1
    if [[ $MODE == "bam" ]]; then
        echo "--outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 500 --limitBAMsortRAM 60000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN"
    else
        echo "--outSAMtype None"
    fi
}

# ---------- Reference path resolver -----------------------------------------

# Resolves a STAR reference index path from species name.
#   $1 : species (e.g. "human", "mouse")
#   $2 : REF_BASE directory
# Prints absolute path and validates it exists.
resolve_reference() {
    local SPECIES=$1 REF_BASE=$2
    local REF="${REF_BASE}/${SPECIES}/2020A/index"

    if [[ ! -d "$REF" ]]; then
        log_error "Reference directory not found: $REF"
        return 1
    fi
    echo "$REF"
}

# ---------- Input validation ------------------------------------------------

validate_fqdir() {
    local FQDIR=$1
    if [[ ! -d "$FQDIR" ]]; then
        log_error "FASTQ directory does not exist: $FQDIR"
        return 1
    fi
    readlink -f "$FQDIR"
}

validate_tag() {
    local TAG=$1
    if [[ -z "$TAG" ]]; then
        log_error "Sample TAG must not be empty."
        return 1
    fi
    echo "$TAG"
}
