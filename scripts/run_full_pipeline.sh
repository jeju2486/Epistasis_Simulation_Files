#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${CONFIG:-config/five_replicates.toml}"
CHECKPOINT_JOBS="${CHECKPOINT_JOBS:-5}"
CASE_JOBS="${CASE_JOBS:-5}"
SPYDRPICK_JOBS="${SPYDRPICK_JOBS:-5}"
SPYDRPICK_THREADS_PER_CASE="${SPYDRPICK_THREADS_PER_CASE:-1}"
KOVAR_JOBS="${KOVAR_JOBS:-5}"
KOVAR_THREADS_PER_CASE="${KOVAR_THREADS_PER_CASE:-1}"
MIN_MAF="${MIN_MAF:-0.05}"
MIN_CELL_COUNT="${MIN_CELL_COUNT:-5}"
SPA_MODE="${SPA_MODE:-auto}"

cd "$REPO_ROOT"
for command in python3 slim SpydrPick ko-variation; do
  command -v "$command" >/dev/null 2>&1 || {
    echo "Missing command in PATH: $command" >&2
    exit 2
  }
done

MANIFEST_ROOT="$(python3 -c 'import sys,tomllib; print(tomllib.load(open(sys.argv[1], "rb"))["paths"]["manifest_root"])' "$CONFIG")"
CHECKPOINT_MANIFEST="$MANIFEST_ROOT/checkpoints.tsv"
CASE_MANIFEST="$MANIFEST_ROOT/cases.tsv"

KOVAR_VERSION="$(ko-variation --version)"
[[ "$KOVAR_VERSION" == "KO-Variation 0.8.3" ]] || {
  echo "Expected KO-Variation 0.8.3; found: $KOVAR_VERSION" >&2
  exit 2
}

python3 scripts/prepare_inputs.py --config "$CONFIG"
python3 scripts/build_manifest.py --config "$CONFIG"
python3 scripts/run_manifest.py --manifest "$CHECKPOINT_MANIFEST" \
  --stage checkpoint --jobs "$CHECKPOINT_JOBS"
python3 scripts/run_manifest.py --manifest "$CASE_MANIFEST" \
  --stage case --jobs "$CASE_JOBS"
python3 scripts/show_status.py "$CHECKPOINT_MANIFEST" "$CASE_MANIFEST"

python3 scripts/run_spydrpick_manifest.py --manifest "$CASE_MANIFEST" \
  --jobs "$SPYDRPICK_JOBS" --threads-per-case "$SPYDRPICK_THREADS_PER_CASE" \
  --min-maf "$MIN_MAF"
python3 scripts/plot_spydrpick.py --manifest "$CASE_MANIFEST"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
python3 scripts/run_kovar_manifest.py --manifest "$CASE_MANIFEST" \
  --jobs "$KOVAR_JOBS" --threads-per-case "$KOVAR_THREADS_PER_CASE" \
  --min-maf "$MIN_MAF" --min-cell-count "$MIN_CELL_COUNT" --spa-mode "$SPA_MODE"
python3 scripts/plot_kovar.py --manifest "$CASE_MANIFEST"
python3 scripts/evaluate_lineage_confounding.py --manifest "$CASE_MANIFEST"

echo "[done] five replicates x three modes x three cross-HGT rates: SLiM -> SpydrPick -> KOVAR 0.8.3 -> lineage evaluation"
