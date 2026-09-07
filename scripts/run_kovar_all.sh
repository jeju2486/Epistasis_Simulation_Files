#!/usr/bin/env bash
set -euo pipefail

MANIFEST="${MANIFEST:-manifests/cases.tsv}"
JOBS="${JOBS:-5}"
THREADS_PER_CASE="${THREADS_PER_CASE:-1}"
MIN_MAF="${MIN_MAF:-0.05}"
MIN_CELL_COUNT="${MIN_CELL_COUNT:-0}"
SPA_MODES="${SPA_MODES:-off auto}"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
for SPA_MODE in $SPA_MODES; do
  case "$SPA_MODE" in off|auto) ;; *) echo "Invalid SPA mode: $SPA_MODE" >&2; exit 2 ;; esac
  python3 scripts/run_kovar_manifest.py --manifest "$MANIFEST" \
    --jobs "$JOBS" --threads-per-case "$THREADS_PER_CASE" \
    --min-maf "$MIN_MAF" --min-cell-count "$MIN_CELL_COUNT" --spa-mode "$SPA_MODE"
  python3 scripts/plot_kovar.py --manifest "$MANIFEST" --spa-mode "$SPA_MODE"
  python3 scripts/evaluate_lineage_confounding.py --manifest "$MANIFEST" \
    --spa-mode "$SPA_MODE"
done
