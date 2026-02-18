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

# ---------- Post-alignment processing --------------------------------------

# Compresses outputs and optionally indexes BAM.
#   $1 : OUTPUT_DIR   – STARsolo output directory (e.g. "output/")
#   No singularity; tools expected in PATH.
process_output_files() {
    local OUTPUT_DIR=${1:-output}

    # Index BAM if it exists
    if [[ -s "Aligned.sortedByCoord.out.bam" ]]; then
        log_info "Indexing BAM file …"
        samtools index -@16 Aligned.sortedByCoord.out.bam
    fi

    # Compress unmapped reads
    if [[ -s "Unmapped.out.mate1" ]]; then
        pbzip2 -9 Unmapped.out.mate1 &
    fi
    if [[ -s "Unmapped.out.mate2" ]]; then
        pbzip2 -9 Unmapped.out.mate2 &
    fi
    wait

    # Gzip all matrix / barcode / feature files
    find "$OUTPUT_DIR" -name "*.tsv" -exec gzip {} + 2>/dev/null
    find "$OUTPUT_DIR" -name "*.mtx" -exec gzip {} + 2>/dev/null

    wait
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
