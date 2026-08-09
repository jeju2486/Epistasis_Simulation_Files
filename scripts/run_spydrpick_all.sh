#!/usr/bin/env bash
set -euo pipefail

MANIFEST="${MANIFEST:-manifests/cases.tsv}"
JOBS="${JOBS:-1}"
THREADS_PER_CASE="${THREADS_PER_CASE:-1}"

python3 scripts/run_spydrpick_manifest.py \
  --manifest "$MANIFEST" \
  --jobs "$JOBS" \
  --threads-per-case "$THREADS_PER_CASE"

python3 scripts/plot_spydrpick.py --manifest "$MANIFEST"
