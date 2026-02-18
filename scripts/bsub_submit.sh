#!/bin/bash
# =============================================================================
# bsub.sh — LSF job submission helper
# =============================================================================
# Wraps any starsolo command in an LSF bsub job.
#
# Usage: bsub.sh <starsolo_command_and_args…>
#
# Example:
#   ./scripts/bsub.sh starsolo 10x /data/fastqs SAMPLE1 --species human
# =============================================================================

set -euo pipefail

if (( $# < 1 )); then
    echo "Usage: bsub.sh <command> [args…]" >&2
    echo "Example: bsub.sh starsolo 10x /data/fastqs SAMPLE1 --species human" >&2
    exit 1
fi

SCRIPT="$*"
AA=("$@")
# Use the last argument as the job tag (basename only)
TAG=$(basename "${AA[-1]}")

GROUP=$(bugroup -w | grep "\b${USER}\b" | cut -d" " -f1)
CPUS=16
RAM=64000
QUE="normal"
WDIR=$(pwd)

bsub -G "$GROUP" \
    -n"$CPUS" \
    -R"span[hosts=1] select[mem>$RAM] rusage[mem=$RAM]" \
    -M"$RAM" \
    -o "$WDIR/$TAG.%J.bsub.log" \
    -e "$WDIR/$TAG.%J.bsub.err" \
    -q "$QUE" \
    $SCRIPT
