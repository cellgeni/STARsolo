#!/bin/bash
# =============================================================================
# STARsolo CLI - inDrops platform
# =============================================================================
# CB_UMI_Complex with 2 barcode segments, adapter sequence, and 6 bp UMI.
# =============================================================================

run_indrops() {
    local FQDIR=$1 TAG=$2 CPUS=$3 REF=$4 WL=$5 BAM=$6 ADAPTER=$7

    local FQDIR_ABS
    FQDIR_ABS=$(validate_fqdir "$FQDIR") || exit 1
    TAG=$(validate_tag "$TAG") || exit 1

    local BC1="$WL/inDrops_Ambrose2_bc1.txt"
    local BC2="$WL/inDrops_Ambrose2_bc2.txt"

    for f in "$BC1" "$BC2"; do
        [[ -f "$f" ]] || die "inDrops barcode file not found: $f"
    done

    # Convert to absolute paths before cd
    REF=$(readlink -f "$REF")
    BC1=$(readlink -f "$BC1")
    BC2=$(readlink -f "$BC2")

    mkdir -p "$TAG" && cd "$TAG" || exit

    # Discover FASTQs
    local fq_info R1 R2
    fq_info=$(find_fastq_files "$FQDIR_ABS" "$TAG")
    IFS='|' read -r R1 R2 <<< "$fq_info"

    # Compression
    local comp_info GZIP ZCMD
    comp_info=$(check_compression "$FQDIR_ABS" "$TAG")
    IFS='|' read -r GZIP ZCMD <<< "$comp_info"

    log_info "Running STARsolo (inDrops) for '$TAG' (adapter=$ADAPTER) …"

    STAR --runThreadN "$CPUS" --genomeDir "$REF" \
        --readFilesIn "$R2" "$R1" --runDirPerm All_RWX $GZIP $BAM \
        --soloType CB_UMI_Complex --soloCBwhitelist "$BC1" "$BC2" \
        --soloAdapterSequence "$ADAPTER" --soloAdapterMismatchesNmax 3 \
        --soloCBmatchWLtype 1MM \
        --soloCBposition 0_0_2_-1 3_1_3_8 \
        --soloUMIposition 3_9_3_14 \
        --soloFeatures Gene GeneFull \
        --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx \
        --outReadsUnmapped Fastx

    process_output_files "output"
}

show_indrops_help() {
    cat <<'EOF'
Usage: starsolo indrops <fastq_dir> <sample_id> [options]

Process inDrops data (2 barcode segments + adapter + 6 bp UMI).

Requires inDrops barcode whitelist files (inDrops_Ambrose2_bc1/2.txt)
in the whitelist directory.

Required:
  <fastq_dir>               Directory containing FASTQ files
  <sample_id>               Sample identifier

Options:
  -s, --species <name>      Species name (resolves reference automatically)
  -r, --ref <path>          Explicit STAR reference index (overrides --species)
  -w, --whitelist-dir <dir> Barcode whitelist directory  [default: from config]
  -c, --cpus <N>            Number of threads             [default: 16]
  --adapter <seq>           Adapter sequence
                             [default: GAGTGATTGCTTGTGACGCCTT]
  --bam                     Output sorted BAM file
  --no-bam                  Do not output BAM (default)
  -h, --help                Show this help
EOF
}
