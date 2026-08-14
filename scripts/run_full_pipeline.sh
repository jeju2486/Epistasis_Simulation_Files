#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="config/five_replicates.toml"
CHECKPOINT_JOBS="${CHECKPOINT_JOBS:-2}"
CASE_JOBS="${CASE_JOBS:-6}"
SPYDRPICK_JOBS="${SPYDRPICK_JOBS:-2}"
SPYDRPICK_THREADS_PER_CASE="${SPYDRPICK_THREADS_PER_CASE:-4}"
KOVAR_JOBS="${KOVAR_JOBS:-1}"
KOVAR_THREADS_PER_CASE="${KOVAR_THREADS_PER_CASE:-8}"
MIN_MAF="${MIN_MAF:-0.05}"

cd "$REPO_ROOT"

for command in python3 slim SpydrPick ko-variation; do
  command -v "$command" >/dev/null 2>&1 || {
    echo "Missing command in PATH: $command" >&2
    exit 2
  }
done

KOVAR_VERSION="$(ko-variation --version)"
[[ "$KOVAR_VERSION" == "KO-Variation 0.8.1" ]] || {
  echo "Expected KO-Variation 0.8.1; found: $KOVAR_VERSION" >&2
  exit 2
}

python3 scripts/prepare_inputs.py --config "$CONFIG"
python3 scripts/build_manifest.py --config "$CONFIG"

python3 scripts/run_manifest.py \
  --manifest manifests/checkpoints.tsv \
  --stage checkpoint \
  --jobs "$CHECKPOINT_JOBS"

python3 scripts/run_manifest.py \
  --manifest manifests/cases.tsv \
  --stage case \
  --jobs "$CASE_JOBS"

python3 scripts/show_status.py manifests/checkpoints.tsv manifests/cases.tsv

python3 scripts/run_spydrpick_manifest.py \
  --manifest manifests/cases.tsv \
  --jobs "$SPYDRPICK_JOBS" \
  --threads-per-case "$SPYDRPICK_THREADS_PER_CASE" \
  --min-maf "$MIN_MAF"

python3 scripts/plot_spydrpick.py --manifest manifests/cases.tsv

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

python3 scripts/run_kovar_manifest.py \
  --manifest manifests/cases.tsv \
  --jobs "$KOVAR_JOBS" \
  --threads-per-case "$KOVAR_THREADS_PER_CASE" \
  --tree-mode oracle \
  --tree-threads 1 \
  --min-maf "$MIN_MAF" \
  --min-cell-count 5 \
  --spa-mode off \
  --full-refit-p 0 \
  --max-pairs 0

python3 scripts/plot_kovar.py \
  --manifest manifests/cases.tsv \
  --tree-mode oracle

echo "[done] simulation -> exhaustive SpydrPick -> exhaustive KOVAR 0.8.1 -> distance plots"
