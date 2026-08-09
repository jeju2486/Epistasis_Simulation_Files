#!/usr/bin/env bash
set -euo pipefail

MANIFEST="${MANIFEST:-manifests/cases.tsv}"
JOBS="${JOBS:-1}"
THREADS_PER_CASE="${THREADS_PER_CASE:-4}"
TREE_THREADS="${TREE_THREADS:-1}"
MAX_PAIRS="${MAX_PAIRS:-0}"
TREE_MODE="${TREE_MODE:-core}"
MIN_MAF="${MIN_MAF:-0.05}"
MIN_CELL_COUNT="${MIN_CELL_COUNT:-5}"
SPA_MODE="${SPA_MODE:-off}"
FULL_REFIT_P="${FULL_REFIT_P:-0}"
WORKER_CHUNK_SIZE="${WORKER_CHUNK_SIZE:-1}"
PREDICTOR_BATCH_SIZE="${PREDICTOR_BATCH_SIZE:-256}"

python3 scripts/run_kovar_manifest.py \
  --manifest "$MANIFEST" \
  --jobs "$JOBS" \
  --threads-per-case "$THREADS_PER_CASE" \
  --tree-threads "$TREE_THREADS" \
  --tree-mode "$TREE_MODE" \
  --min-maf "$MIN_MAF" \
  --min-cell-count "$MIN_CELL_COUNT" \
  --spa-mode "$SPA_MODE" \
  --full-refit-p "$FULL_REFIT_P" \
  --worker-chunk-size "$WORKER_CHUNK_SIZE" \
  --predictor-batch-size "$PREDICTOR_BATCH_SIZE" \
  --max-pairs "$MAX_PAIRS"
