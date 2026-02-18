#!/bin/bash
# =============================================================================
# STARsolo CLI - Smart-seq / Smart-seq2 platform
# =============================================================================
# Plate-based, no UMIs. Requires a manifest TSV file.
# =============================================================================

run_smartseq() {
    local TSV=$1 CPUS=$2 REF=$3 BAM=$4

    [[ -s "$TSV" ]] || die "Manifest file not found or empty: $TSV"

    local TAG
    TAG=$(basename "${TSV%%.manifest.tsv}")
    TAG=$(basename "${TAG%%.tsv}")        # fallback if named differently
    [[ -n "$TAG" ]] || die "Could not derive sample name from manifest filename."

    # Convert to absolute paths before cd
    TSV=$(readlink -f "$TSV")
    REF=$(readlink -f "$REF")

    # Create output directory
    mkdir -p "$TAG" || exit

    # Convert manifest paths to absolute and write to output directory
    local NEW_MANIFEST="$TAG/manifest.tsv"
    while IFS=$'\t' read -r r1 r2 cell; do
        r1=$(readlink -f "$r1")
        r2=$(readlink -f "$r2")
        echo -e "$r1\t$r2\t$cell"
    done < "$TSV" > "$NEW_MANIFEST"

    cd "$TAG" || exit
    TSV="manifest.tsv"

    # Compression detection from manifest content
    local GZIP=""
    if grep -qP "\.gz\t" "$TSV"; then
        GZIP="--readFilesCommand zcat"
    fi

    log_info "Running STARsolo (Smart-seq) for '$TAG' …"

    STAR --runThreadN "$CPUS" --genomeDir "$REF" --runDirPerm All_RWX $GZIP $BAM \
        --soloType SmartSeq --readFilesManifest "$TSV" \
        --soloUMIdedup Exact --soloStrand Unstranded \
        --limitOutSJcollapsed 10000000 --soloCellFilter None \
        --soloFeatures Gene GeneFull \
        --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx \
        --outReadsUnmapped Fastx

    process_output_files "output"
}

show_smartseq_help() {
    cat <<'EOF'
Usage: starsolo smartseq <manifest.tsv> [options]

Process Smart-seq / Smart-seq2 data using a manifest file.
The manifest is a tab-separated file with 3 columns:
  1) Full path to R1   2) Full path to R2   3) Cell name/ID

The sample name is derived from the manifest filename (strips .manifest.tsv).

Required:
  <manifest.tsv>            Path to the manifest TSV file

Options:
  -s, --species <name>      Species name (resolves reference automatically)
  -r, --ref <path>          Explicit STAR reference index (overrides --species)
  -c, --cpus <N>            Number of threads             [default: 16]
  --bam                     Output sorted BAM file
  --no-bam                  Do not output BAM (default)
  -h, --help                Show this help
EOF
}
