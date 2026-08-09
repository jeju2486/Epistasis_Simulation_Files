#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="$REPO_ROOT/config/five_replicates.toml"
CHECKPOINT_JOBS="${CHECKPOINT_JOBS:-2}"
CASE_JOBS="${CASE_JOBS:-6}"

cd "$REPO_ROOT"

command -v slim >/dev/null 2>&1 || {
  echo "SLiM is not available. Activate the epistasis-sim environment first." >&2
  exit 2
}

python scripts/prepare_inputs.py --config "$CONFIG"
python scripts/build_manifest.py --config "$CONFIG"

echo "Running 5 neutral checkpoints with $CHECKPOINT_JOBS parallel job(s)..."
python scripts/run_manifest.py \
  --manifest manifests/checkpoints.tsv \
  --stage checkpoint \
  --jobs "$CHECKPOINT_JOBS"

echo "Running 30 paired continuations with $CASE_JOBS parallel job(s)..."
python scripts/run_manifest.py \
  --manifest manifests/cases.tsv \
  --stage case \
  --jobs "$CASE_JOBS"

python scripts/show_status.py manifests/checkpoints.tsv manifests/cases.tsv
