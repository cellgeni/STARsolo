#!/bin/bash
# =============================================================================
# STARsolo CLI - STRT-seq platform
# =============================================================================
# CB_UMI_Simple with a 96-barcode whitelist (8 bp CB, 8 bp UMI).
# NOTE: R1 is the biological read; R2 carries the barcode in STRT-seq.
# =============================================================================

run_strt() {
    local FQDIR=$1 TAG=$2 CPUS=$3 REF=$4 WL=$5 BAM=$6 STRAND=$7

    local FQDIR_ABS
    FQDIR_ABS=$(validate_fqdir "$FQDIR") || exit 1
    TAG=$(validate_tag "$TAG") || exit 1

    local BC="$WL/96_barcodes.list"
    [[ -f "$BC" ]] || die "STRT-seq barcode file not found: $BC"

    # Convert to absolute paths before cd
    REF=$(readlink -f "$REF")
    BC=$(readlink -f "$BC")

    local CBLEN=8
    local UMILEN=8

    mkdir -p "$TAG" && cd "$TAG" || exit

    # Discover FASTQs
    local fq_info R1 R2
    fq_info=$(find_fastq_files "$FQDIR_ABS" "$TAG")
    IFS='|' read -r R1 R2 <<< "$fq_info"

    # Compression
    local comp_info GZIP ZCMD
    comp_info=$(check_compression "$FQDIR_ABS" "$TAG")
    IFS='|' read -r GZIP ZCMD <<< "$comp_info"

    log_info "Running STARsolo (STRT-seq) for '$TAG' (strand=$STRAND) …"

    # NOTE: R2 R1 order — R1 is biological read, R2 carries barcode
    STAR --runThreadN "$CPUS" --genomeDir "$REF" \
        --readFilesIn "$R2" "$R1" --runDirPerm All_RWX $GZIP $BAM \
        --soloType CB_UMI_Simple --soloCBwhitelist "$BC" \
        --soloBarcodeReadLength 0 \
        --soloCBlen "$CBLEN" --soloUMIstart $((CBLEN + 1)) --soloUMIlen "$UMILEN" \
        --soloStrand "$STRAND" \
        --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIfiltering MultiGeneUMI_CR \
        --limitOutSJcollapsed 10000000 --soloCellFilter None \
        --soloFeatures Gene GeneFull \
        --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx \
        --outReadsUnmapped Fastx

    process_output_files "output"
}

show_strt_help() {
    cat <<'EOF'
Usage: starsolo strt <fastq_dir> <sample_id> [options]

Process STRT-seq data (96-barcode whitelist, 8 bp CB, 8 bp UMI).

Requires 96_barcodes.list in the whitelist directory.

Required:
  <fastq_dir>               Directory containing FASTQ files
  <sample_id>               Sample identifier

Options:
  -s, --species <name>      Species name (resolves reference automatically)
  -r, --ref <path>          Explicit STAR reference index (overrides --species)
  -w, --whitelist-dir <dir> Barcode whitelist directory  [default: from config]
  -c, --cpus <N>            Number of threads             [default: 16]
  --strand <Forward|Reverse> Strand specificity           [default: Forward]
  --bam                     Output sorted BAM file
  --no-bam                  Do not output BAM (default)
  -h, --help                Show this help
EOF
}
