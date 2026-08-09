#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 CORE_SNPS_FASTA OUTPUT_PREFIX SEED" >&2
  exit 2
fi

ALIGNMENT="$1"
PREFIX="$2"
SEED="$3"
IQTREE_BIN="${IQTREE_BIN:-iqtree2}"

[[ -s "$ALIGNMENT" ]] || { echo "Alignment not found: $ALIGNMENT" >&2; exit 2; }
mkdir -p "$(dirname "$PREFIX")"

"$IQTREE_BIN" \
  -s "$ALIGNMENT" \
  -pre "$PREFIX" \
  -m GTR+ASC \
  -seed "$SEED" \
  -nt AUTO \
  -redo

[[ -s "${PREFIX}.treefile" ]] || { echo "IQ-TREE did not produce ${PREFIX}.treefile" >&2; exit 1; }
echo "[done] inferred core tree: ${PREFIX}.treefile"
