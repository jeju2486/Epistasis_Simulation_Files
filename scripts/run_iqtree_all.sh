#!/usr/bin/env bash
set -euo pipefail

MANIFEST="${MANIFEST:-manifests/cases.tsv}"
JOBS="${JOBS:-5}"
THREADS_PER_CASE="${THREADS_PER_CASE:-1}"
IQTREE_MODEL="${IQTREE_MODEL:-GTR+ASC}"
IQTREE_BIN="${IQTREE_BIN:-iqtree}"

python3 scripts/run_iqtree_manifest.py --manifest "$MANIFEST" \
  --jobs "$JOBS" --threads-per-case "$THREADS_PER_CASE" \
  --iqtree "$IQTREE_BIN" --model "$IQTREE_MODEL"
