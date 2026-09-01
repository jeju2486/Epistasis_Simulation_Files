#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${CONFIG:-config/five_replicates.toml}"
CHECKPOINT_JOBS="${CHECKPOINT_JOBS:-5}"
CASE_JOBS="${CASE_JOBS:-5}"

cd "$REPO_ROOT"
command -v slim >/dev/null 2>&1 || {
  echo "SLiM is not available. Activate the epistasis-sim environment first." >&2
  exit 2
}
python3 scripts/prepare_inputs.py --config "$CONFIG"
python3 scripts/build_manifest.py --config "$CONFIG"
MANIFEST_ROOT="$(python3 -c 'import sys,tomllib; print(tomllib.load(open(sys.argv[1], "rb"))["paths"]["manifest_root"])' "$CONFIG")"
CHECKPOINT_MANIFEST="$MANIFEST_ROOT/checkpoints.tsv"
CASE_MANIFEST="$MANIFEST_ROOT/cases.tsv"
echo "Running neutral checkpoints with $CHECKPOINT_JOBS parallel job(s)..."
python3 scripts/run_manifest.py --manifest "$CHECKPOINT_MANIFEST" \
  --stage checkpoint --jobs "$CHECKPOINT_JOBS"
echo "Running mode-by-cross-HGT cases with $CASE_JOBS parallel job(s)..."
python3 scripts/run_manifest.py --manifest "$CASE_MANIFEST" \
  --stage case --jobs "$CASE_JOBS"
python3 scripts/show_status.py "$CHECKPOINT_MANIFEST" "$CASE_MANIFEST"
