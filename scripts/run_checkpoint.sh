#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SLIM_BIN="${SLIM_BIN:-slim}"

OUTDIR="" REF_FILE="" SEED=""
GENOME_LENGTH="" MU="" TRACT_LENGTH="" HGT_WITHIN=""
N_ANCESTRAL="" N_CLADE="" N_TERMINAL=""
ANCESTRAL_GENERATIONS="" DEEP_GENERATIONS="" TERMINAL_GENERATIONS=""
GLOBAL_FREQ_MIN="" GLOBAL_FREQ_MAX="" LINEAGE_FREQ_MIN="" LINEAGE_FREQ_MAX=""
MINIMUM_DISTANCE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    --reference) REF_FILE="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --genome-length) GENOME_LENGTH="$2"; shift 2 ;;
    --mu) MU="$2"; shift 2 ;;
    --tract-length) TRACT_LENGTH="$2"; shift 2 ;;
    --within-hgt) HGT_WITHIN="$2"; shift 2 ;;
    --ancestral-size) N_ANCESTRAL="$2"; shift 2 ;;
    --clade-size) N_CLADE="$2"; shift 2 ;;
    --terminal-size) N_TERMINAL="$2"; shift 2 ;;
    --ancestral-generations) ANCESTRAL_GENERATIONS="$2"; shift 2 ;;
    --deep-generations) DEEP_GENERATIONS="$2"; shift 2 ;;
    --terminal-generations) TERMINAL_GENERATIONS="$2"; shift 2 ;;
    --global-freq-min) GLOBAL_FREQ_MIN="$2"; shift 2 ;;
    --global-freq-max) GLOBAL_FREQ_MAX="$2"; shift 2 ;;
    --lineage-freq-min) LINEAGE_FREQ_MIN="$2"; shift 2 ;;
    --lineage-freq-max) LINEAGE_FREQ_MAX="$2"; shift 2 ;;
    --minimum-distance) MINIMUM_DISTANCE="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

required=(OUTDIR REF_FILE SEED GENOME_LENGTH MU TRACT_LENGTH HGT_WITHIN N_ANCESTRAL N_CLADE N_TERMINAL ANCESTRAL_GENERATIONS DEEP_GENERATIONS TERMINAL_GENERATIONS GLOBAL_FREQ_MIN GLOBAL_FREQ_MAX LINEAGE_FREQ_MIN LINEAGE_FREQ_MAX MINIMUM_DISTANCE)
for name in "${required[@]}"; do
  [[ -n "${!name}" ]] || { echo "Missing required value: $name" >&2; exit 2; }
done

[[ "$OUTDIR" = /* ]] || OUTDIR="$REPO_ROOT/$OUTDIR"
[[ "$REF_FILE" = /* ]] || REF_FILE="$REPO_ROOT/$REF_FILE"
[[ -s "$REF_FILE" ]] || { echo "Reference not found: $REF_FILE" >&2; exit 2; }

if [[ -f "$OUTDIR/_SUCCESS" ]]; then
  echo "[skip] checkpoint complete: $OUTDIR"
  exit 0
fi

mkdir -p "$(dirname "$OUTDIR")"
LOCKDIR="${OUTDIR}.lock"
if ! mkdir "$LOCKDIR" 2>/dev/null; then
  echo "Checkpoint is already locked: $OUTDIR" >&2
  exit 3
fi

TMPDIR_RUN="$(mktemp -d "${OUTDIR}.tmp.XXXXXX")"
cleanup() {
  rc=$?
  rmdir "$LOCKDIR" 2>/dev/null || true
  if [[ $rc -ne 0 && -d "$TMPDIR_RUN" ]]; then
    failed="${OUTDIR}.failed.$(date -u +%Y%m%dT%H%M%SZ)"
    mv "$TMPDIR_RUN" "$failed"
    echo "[failed] retained attempt at $failed" >&2
    if [[ -s "$failed/slim.log" ]]; then
      echo "[failed] SLiM log follows:" >&2
      sed 's/^/[slim] /' "$failed/slim.log" >&2
    fi
  fi
}
trap cleanup EXIT

if [[ -e "$OUTDIR" ]]; then
  stale="${OUTDIR}.incomplete.$(date -u +%Y%m%dT%H%M%SZ)"
  mv "$OUTDIR" "$stale"
  echo "[archive] moved incomplete output to $stale"
fi

cat > "$TMPDIR_RUN/params.tsv" <<EOF
parameter	value
seed	$SEED
reference	$REF_FILE
genome_length	$GENOME_LENGTH
mutation_rate	$MU
tract_length	$TRACT_LENGTH
within_hgt_probability	$HGT_WITHIN
ancestral_size	$N_ANCESTRAL
clade_size	$N_CLADE
terminal_size	$N_TERMINAL
ancestral_generations	$ANCESTRAL_GENERATIONS
deep_generations	$DEEP_GENERATIONS
terminal_generations	$TERMINAL_GENERATIONS
global_frequency_min	$GLOBAL_FREQ_MIN
global_frequency_max	$GLOBAL_FREQ_MAX
lineage_frequency_min	$LINEAGE_FREQ_MIN
lineage_frequency_max	$LINEAGE_FREQ_MAX
minimum_distance	$MINIMUM_DISTANCE
EOF

"$SLIM_BIN" -s "$SEED" \
  -d "REF_FILE=\"$REF_FILE\"" \
  -d "CHECKPOINT_FILE=\"$TMPDIR_RUN/checkpoint.trees\"" \
  -d "LOCI_FILE=\"$TMPDIR_RUN/selected_loci.tsv\"" \
  -d "METRICS_FILE=\"$TMPDIR_RUN/checkpoint_metrics.tsv\"" \
  -d "GENOME_LENGTH=$GENOME_LENGTH" \
  -d "MU=$MU" \
  -d "TRACT_LENGTH=$TRACT_LENGTH" \
  -d "HGT_P_WITHIN=$HGT_WITHIN" \
  -d "N_ANCESTRAL=$N_ANCESTRAL" \
  -d "N_CLADE=$N_CLADE" \
  -d "N_TERMINAL=$N_TERMINAL" \
  -d "ANCESTRAL_GENERATIONS=$ANCESTRAL_GENERATIONS" \
  -d "DEEP_GENERATIONS=$DEEP_GENERATIONS" \
  -d "TERMINAL_GENERATIONS=$TERMINAL_GENERATIONS" \
  -d "GLOBAL_FREQ_MIN=$GLOBAL_FREQ_MIN" \
  -d "GLOBAL_FREQ_MAX=$GLOBAL_FREQ_MAX" \
  -d "LINEAGE_FREQ_MIN=$LINEAGE_FREQ_MIN" \
  -d "LINEAGE_FREQ_MAX=$LINEAGE_FREQ_MAX" \
  -d "MIN_DISTANCE=$MINIMUM_DISTANCE" \
  "$REPO_ROOT/slim/build_checkpoint.slim" \
  > "$TMPDIR_RUN/slim.log" 2>&1

[[ -s "$TMPDIR_RUN/checkpoint.trees" ]] || { echo "Missing checkpoint.trees" >&2; exit 1; }
python "$REPO_ROOT/scripts/validate_selected_loci.py" "$TMPDIR_RUN/selected_loci.tsv"
[[ -s "$TMPDIR_RUN/checkpoint_metrics.tsv" ]] || { echo "Missing checkpoint metrics" >&2; exit 1; }

touch "$TMPDIR_RUN/_SUCCESS"
mv "$TMPDIR_RUN" "$OUTDIR"
echo "[done] checkpoint: $OUTDIR"
