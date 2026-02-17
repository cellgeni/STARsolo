#!/bin/bash
# =============================================================================
# STARsolo CLI - BD Rhapsody platform
# =============================================================================
# CB_UMI_Complex with 3 barcode segments and 8 bp UMI.
# =============================================================================

run_rhapsody() {
    local FQDIR=$1 TAG=$2 CPUS=$3 REF=$4 WL=$5 BAM=$6

    local FQDIR_ABS
    FQDIR_ABS=$(validate_fqdir "$FQDIR") || exit 1
    TAG=$(validate_tag "$TAG") || exit 1

    local BC1="$WL/Rhapsody_bc1.txt"
    local BC2="$WL/Rhapsody_bc2.txt"
    local BC3="$WL/Rhapsody_bc3.txt"

    for f in "$BC1" "$BC2" "$BC3"; do
        [[ -f "$f" ]] || die "Rhapsody barcode file not found: $f"
    done

    mkdir -p "$TAG" && cd "$TAG" || exit

    # Discover FASTQs
    local fq_info R1 R2
    fq_info=$(find_fastq_files "$FQDIR_ABS" "$TAG")
    IFS='|' read -r R1 R2 <<< "$fq_info"

    # Compression
    local comp_info GZIP ZCMD
    comp_info=$(check_compression "$FQDIR_ABS" "$TAG")
    IFS='|' read -r GZIP ZCMD <<< "$comp_info"

    log_info "Running STARsolo (BD Rhapsody) for '$TAG' …"

    STAR --runThreadN "$CPUS" --genomeDir "$REF" \
        --readFilesIn "$R2" "$R1" --runDirPerm All_RWX $GZIP $BAM \
        --soloType CB_UMI_Complex --soloCBwhitelist "$BC1" "$BC2" "$BC3" \
        --soloUMIlen 8 \
        --soloCBmatchWLtype 1MM \
        --soloCBposition 0_0_0_8 0_21_0_29 0_43_0_51 \
        --soloUMIposition 0_52_0_59 \
        --soloFeatures Gene GeneFull \
        --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx \
        --outReadsUnmapped Fastx

    process_output_files "output"
}

show_rhapsody_help() {
    cat <<'EOF'
Usage: starsolo rhapsody <fastq_dir> <sample_id> [options]

Process BD Rhapsody data (3 barcode segments, 8 bp UMI).

Requires Rhapsody barcode whitelist files (Rhapsody_bc1/2/3.txt)
in the whitelist directory.

Required:
  <fastq_dir>               Directory containing FASTQ files
  <sample_id>               Sample identifier

Options:
  -s, --species <name>      Species name (resolves reference automatically)
  -r, --ref <path>          Explicit STAR reference index (overrides --species)
  -w, --whitelist-dir <dir> Barcode whitelist directory  [default: from config]
  -c, --cpus <N>            Number of threads             [default: 16]
  --bam                     Output sorted BAM file
  --no-bam                  Do not output BAM (default)
  -h, --help                Show this help
EOF
}
