#!/bin/bash
# =============================================================================
# STARsolo CLI - Drop-seq platform
# =============================================================================
# CB_UMI_Simple with no whitelist (12 bp CB, 8 bp UMI).
# =============================================================================

run_dropseq() {
    local FQDIR=$1 TAG=$2 CPUS=$3 REF=$4 BAM=$5

    local FQDIR_ABS
    FQDIR_ABS=$(validate_fqdir "$FQDIR") || exit 1
    TAG=$(validate_tag "$TAG") || exit 1

    # Convert to absolute path before cd
    REF=$(readlink -f "$REF")

    mkdir -p "$TAG" && cd "$TAG" || exit

    # Discover FASTQs
    local fq_info R1 R2
    fq_info=$(find_fastq_files "$FQDIR_ABS" "$TAG")
    IFS='|' read -r R1 R2 <<< "$fq_info"

    # Compression
    local comp_info GZIP ZCMD
    comp_info=$(check_compression "$FQDIR_ABS" "$TAG")
    IFS='|' read -r GZIP ZCMD <<< "$comp_info"

    log_info "Running STARsolo (Drop-seq) for '$TAG' …"

    STAR --runThreadN "$CPUS" --genomeDir "$REF" \
        --readFilesIn "$R2" "$R1" --runDirPerm All_RWX $GZIP $BAM \
        --soloType CB_UMI_Simple --soloCBwhitelist None \
        --soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 \
        --soloBarcodeReadLength 0 \
        --soloFeatures Gene GeneFull \
        --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx \
        --outReadsUnmapped Fastx

    process_output_files "output"
}

show_dropseq_help() {
    cat <<'EOF'
Usage: starsolo dropseq <fastq_dir> <sample_id> [options]

Process Drop-seq data (12 bp cell barcode, 8 bp UMI, no whitelist).

Required:
  <fastq_dir>               Directory containing FASTQ files
  <sample_id>               Sample identifier

Options:
  -s, --species <name>      Species name (resolves reference automatically)
  -r, --ref <path>          Explicit STAR reference index (overrides --species)
  -c, --cpus <N>            Number of threads             [default: 16]
  --bam                     Output sorted BAM file
  --no-bam                  Do not output BAM (default)
  -h, --help                Show this help
EOF
}
