#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SLIM_BIN="${SLIM_BIN:-slim}"

OUTDIR="" REF_FILE="" SEED=""
GENOME_LENGTH="" MU="" TRACT_LENGTH="" HGT_WITHIN="" HGT_CROSS=""
N_ANCESTRAL="" N_CLADE="" N_TERMINAL=""
ANCESTRAL_GENERATIONS="" DEEP_GENERATIONS="" TERMINAL_GENERATIONS=""
A_POSITION="" B_POSITION="" C_POSITION="" D_POSITION=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    --reference) REF_FILE="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --genome-length) GENOME_LENGTH="$2"; shift 2 ;;
    --mu) MU="$2"; shift 2 ;;
    --tract-length) TRACT_LENGTH="$2"; shift 2 ;;
    --within-hgt) HGT_WITHIN="$2"; shift 2 ;;
    --cross-hgt) HGT_CROSS="$2"; shift 2 ;;
    --ancestral-size) N_ANCESTRAL="$2"; shift 2 ;;
    --clade-size) N_CLADE="$2"; shift 2 ;;
    --terminal-size) N_TERMINAL="$2"; shift 2 ;;
    --ancestral-generations) ANCESTRAL_GENERATIONS="$2"; shift 2 ;;
    --deep-generations) DEEP_GENERATIONS="$2"; shift 2 ;;
    --terminal-generations) TERMINAL_GENERATIONS="$2"; shift 2 ;;
    --a-position) A_POSITION="$2"; shift 2 ;;
    --b-position) B_POSITION="$2"; shift 2 ;;
    --c-position) C_POSITION="$2"; shift 2 ;;
    --d-position) D_POSITION="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

required=(OUTDIR REF_FILE SEED GENOME_LENGTH MU TRACT_LENGTH HGT_WITHIN HGT_CROSS N_ANCESTRAL N_CLADE N_TERMINAL ANCESTRAL_GENERATIONS DEEP_GENERATIONS TERMINAL_GENERATIONS A_POSITION B_POSITION C_POSITION D_POSITION)
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
cross_hgt_probability	$HGT_CROSS
ancestral_size	$N_ANCESTRAL
clade_size	$N_CLADE
terminal_size	$N_TERMINAL
ancestral_generations	$ANCESTRAL_GENERATIONS
deep_generations	$DEEP_GENERATIONS
terminal_generations	$TERMINAL_GENERATIONS
A_position	$A_POSITION
B_position	$B_POSITION
C_position	$C_POSITION
D_position	$D_POSITION
AB_distance	$((B_POSITION - A_POSITION))
CD_distance	$((D_POSITION - C_POSITION))
focal_seeding	none_until_experimental_generation_0
EOF

"$SLIM_BIN" -s "$SEED" \
  -d "REF_FILE=\"$REF_FILE\"" \
  -d "CHECKPOINT_FILE=\"$TMPDIR_RUN/checkpoint.trees\"" \
  -d "METRICS_FILE=\"$TMPDIR_RUN/checkpoint_metrics.tsv\"" \
  -d "GENOME_LENGTH=$GENOME_LENGTH" \
  -d "MU=$MU" \
  -d "TRACT_LENGTH=$TRACT_LENGTH" \
  -d "HGT_P_WITHIN=$HGT_WITHIN" \
  -d "HGT_P_CROSS=$HGT_CROSS" \
  -d "N_ANCESTRAL=$N_ANCESTRAL" \
  -d "N_CLADE=$N_CLADE" \
  -d "N_TERMINAL=$N_TERMINAL" \
  -d "ANCESTRAL_GENERATIONS=$ANCESTRAL_GENERATIONS" \
  -d "DEEP_GENERATIONS=$DEEP_GENERATIONS" \
  -d "TERMINAL_GENERATIONS=$TERMINAL_GENERATIONS" \
  -d "A_POSITION=$A_POSITION" \
  -d "B_POSITION=$B_POSITION" \
  -d "C_POSITION=$C_POSITION" \
  -d "D_POSITION=$D_POSITION" \
  "$REPO_ROOT/slim/build_checkpoint.slim" \
  > "$TMPDIR_RUN/slim.log" 2>&1

[[ -s "$TMPDIR_RUN/checkpoint.trees" ]] || { echo "Missing checkpoint.trees" >&2; exit 1; }
[[ -s "$TMPDIR_RUN/checkpoint_metrics.tsv" ]] || { echo "Missing checkpoint metrics" >&2; exit 1; }

touch "$TMPDIR_RUN/_SUCCESS"
mv "$TMPDIR_RUN" "$OUTDIR"
echo "[done] checkpoint: $OUTDIR"
