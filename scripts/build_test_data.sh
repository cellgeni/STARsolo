#!/bin/bash
# =============================================================================
# Rebuilds data/test.tar.gz: synthetic reads, reference and STAR index.
#
# Run from anywhere inside the starsolo image (needs python3 and STAR), e.g.
#   singularity exec -B "$PWD" docker://quay.io/cellgeni/starsolo:latest \
#       bash scripts/build_test_data.sh
#
# Options:
#   --seed S      generator seed                  [default: 20261015]
#   --out PATH    output tarball                  [default: data/test.tar.gz]
#   --no-verify   skip tests/test_platforms.sh on the new data
# =============================================================================

set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SEED=20261015
OUT=$ROOT/data/test.tar.gz
VERIFY=1

while (( $# > 0 )); do
    case "$1" in
        --seed)      SEED=$2; shift 2 ;;
        --out)       OUT=$2; shift 2 ;;
        --no-verify) VERIFY=0; shift ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

for tool in python3 STAR; do
    command -v "$tool" >/dev/null || { echo "$tool not on PATH; run inside the starsolo image" >&2; exit 1; }
done

WORK=$(mktemp -d)
trap 'rm -rf -- "$WORK"' EXIT

python3 -I "$ROOT/scripts/generate_test_data.py" --outdir "$WORK/data/test" --seed "$SEED"

# Tiny genome: genomeSAindexNbases must be about log2(length)/2 - 1 (here 4);
# the default 14 builds a 1.5 GB SAindex and nothing maps.
(
    cd "$WORK/data/test/reference"
    mkdir index
    STAR --runMode genomeGenerate --runThreadN 1 --genomeDir index \
        --genomeFastaFiles test.fa --sjdbGTFfile test.gtf \
        --sjdbOverhang 100 --genomeSAindexNbases 4 \
        --outFileNamePrefix "$WORK/genome." > /dev/null
    rm -f index/Log.out          # run log with timestamps; STAR does not need it
)

if (( VERIFY )); then
    STARSOLO_TEST_DATA=$WORK/data/test bash "$ROOT/tests/test_platforms.sh"
fi

# Reproducible archive: fixed order, times and owners
tar -C "$WORK" --sort=name --mtime=@0 --owner=0 --group=0 --numeric-owner -cf - data/test \
    | gzip -n -9 > "$OUT.tmp"
mv "$OUT.tmp" "$OUT"
echo "Wrote $OUT"
