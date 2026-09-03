#!/usr/bin/env bash
set -euo pipefail

MANIFEST="${MANIFEST:-manifests/cases.tsv}"
JOBS="${JOBS:-5}"
THREADS_PER_CASE="${THREADS_PER_CASE:-1}"
MIN_MAF="${MIN_MAF:-0.05}"

for SAMPLE_REWEIGHTING in default none; do
  python3 scripts/run_spydrpick_manifest.py --manifest "$MANIFEST" \
    --jobs "$JOBS" --threads-per-case "$THREADS_PER_CASE" --min-maf "$MIN_MAF" \
    --sample-reweighting "$SAMPLE_REWEIGHTING"
  python3 scripts/plot_spydrpick.py --manifest "$MANIFEST" \
    --sample-reweighting "$SAMPLE_REWEIGHTING"
done
