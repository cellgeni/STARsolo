#!/bin/bash
# =============================================================================
# bbduk.sh — adapter/polyA trimming helper
# =============================================================================
# Usage: bbduk.sh <sample_id> [adapters_fasta]
#
# Requires BBMap (bbduk.sh) in PATH.
# =============================================================================

set -euo pipefail

TAG=${1:?"Usage: bbduk.sh <sample_id> [adapters_fasta]"}
ADAPTERS=${2:-/software/cellgen/cellgeni/bbmap/resources/adapters.fa}

bbduk.sh \
    in1="${TAG}_1.fastq.gz" \
    in2="${TAG}_2.fastq.gz" \
    out1="${TAG}.bbduk.R1.fastq" \
    out2="${TAG}.bbduk.R2.fastq" \
    ref="$ADAPTERS" \
    trimpolya=10 ktrim=r k=23 mink=11 hdist=1 tpe tbo \
    &> "${TAG}.bbduk.log"

echo "[INFO] Trimming complete. Output: ${TAG}.bbduk.R1.fastq / ${TAG}.bbduk.R2.fastq"
