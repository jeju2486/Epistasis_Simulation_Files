#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SLIM_BIN="${SLIM_BIN:-slim}"

OUTDIR="" CHECKPOINT_DIR="" REF_FILE="" SEED="" MODE="" CROSS_HGT=""
GENOME_LENGTH="" MU="" TRACT_LENGTH="" HGT_WITHIN="" N_TERMINAL=""
SAMPLE_PER_TERMINAL="" EXPERIMENT_GENERATIONS="" S_AB="" S_CD="" MONITOR_EVERY=""
ANCESTRAL_NE="" ORACLE_TREE_POSITION=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    --checkpoint-dir) CHECKPOINT_DIR="$2"; shift 2 ;;
    --reference) REF_FILE="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --cross-hgt) CROSS_HGT="$2"; shift 2 ;;
    --genome-length) GENOME_LENGTH="$2"; shift 2 ;;
    --mu) MU="$2"; shift 2 ;;
    --tract-length) TRACT_LENGTH="$2"; shift 2 ;;
    --within-hgt) HGT_WITHIN="$2"; shift 2 ;;
    --terminal-size) N_TERMINAL="$2"; shift 2 ;;
    --sample-per-terminal) SAMPLE_PER_TERMINAL="$2"; shift 2 ;;
    --experiment-generations) EXPERIMENT_GENERATIONS="$2"; shift 2 ;;
    --s-ab) S_AB="$2"; shift 2 ;;
    --s-cd) S_CD="$2"; shift 2 ;;
    --monitor-every) MONITOR_EVERY="$2"; shift 2 ;;
    --ancestral-ne) ANCESTRAL_NE="$2"; shift 2 ;;
    --oracle-tree-position) ORACLE_TREE_POSITION="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

required=(OUTDIR CHECKPOINT_DIR REF_FILE SEED MODE CROSS_HGT GENOME_LENGTH MU TRACT_LENGTH HGT_WITHIN N_TERMINAL SAMPLE_PER_TERMINAL EXPERIMENT_GENERATIONS S_AB S_CD MONITOR_EVERY ANCESTRAL_NE ORACLE_TREE_POSITION)
for name in "${required[@]}"; do
  [[ -n "${!name}" ]] || { echo "Missing required value: $name" >&2; exit 2; }
done

[[ "$OUTDIR" = /* ]] || OUTDIR="$REPO_ROOT/$OUTDIR"
[[ "$CHECKPOINT_DIR" = /* ]] || CHECKPOINT_DIR="$REPO_ROOT/$CHECKPOINT_DIR"
[[ "$REF_FILE" = /* ]] || REF_FILE="$REPO_ROOT/$REF_FILE"

[[ -f "$CHECKPOINT_DIR/_SUCCESS" ]] || { echo "Checkpoint is incomplete: $CHECKPOINT_DIR" >&2; exit 2; }
[[ -s "$REF_FILE" ]] || { echo "Reference not found: $REF_FILE" >&2; exit 2; }

if [[ -f "$OUTDIR/_SUCCESS" ]]; then
  echo "[skip] case complete: $OUTDIR"
  exit 0
fi

LOCI_FILE="$CHECKPOINT_DIR/selected_loci.tsv"
A_ID="$(awk -F '\t' '$1=="A" {print $2}' "$LOCI_FILE")"
B_ID="$(awk -F '\t' '$1=="B" {print $2}' "$LOCI_FILE")"
C_ID="$(awk -F '\t' '$1=="C" {print $2}' "$LOCI_FILE")"
D_ID="$(awk -F '\t' '$1=="D" {print $2}' "$LOCI_FILE")"
[[ -n "$A_ID" && -n "$B_ID" && -n "$C_ID" && -n "$D_ID" ]] || { echo "Could not parse A-D mutation IDs" >&2; exit 2; }

mkdir -p "$(dirname "$OUTDIR")"
LOCKDIR="${OUTDIR}.lock"
if ! mkdir "$LOCKDIR" 2>/dev/null; then
  echo "Case is already locked: $OUTDIR" >&2
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
  fi
}
trap cleanup EXIT

if [[ -e "$OUTDIR" ]]; then
  stale="${OUTDIR}.incomplete.$(date -u +%Y%m%dT%H%M%SZ)"
  mv "$OUTDIR" "$stale"
  echo "[archive] moved incomplete output to $stale"
fi

cp "$LOCI_FILE" "$TMPDIR_RUN/selected_loci.tsv"
cp "$CHECKPOINT_DIR/checkpoint_metrics.tsv" "$TMPDIR_RUN/checkpoint_metrics.tsv"

cat > "$TMPDIR_RUN/params.tsv" <<EOF
parameter	value
seed	$SEED
mode	$MODE
cross_hgt_probability	$CROSS_HGT
checkpoint	$CHECKPOINT_DIR/checkpoint.trees
A_mutation_id	$A_ID
B_mutation_id	$B_ID
C_mutation_id	$C_ID
D_mutation_id	$D_ID
experimental_generations	$EXPERIMENT_GENERATIONS
s_ab	$S_AB
s_cd	$S_CD
EOF

"$SLIM_BIN" -s "$SEED" \
  -d "REF_FILE=\"$REF_FILE\"" \
  -d "CHECKPOINT_FILE=\"$CHECKPOINT_DIR/checkpoint.trees\"" \
  -d "TREE_OUT=\"$TMPDIR_RUN/out.trees\"" \
  -d "VCF_OUT=\"$TMPDIR_RUN/out.vcf\"" \
  -d "SAMPLE_FILE=\"$TMPDIR_RUN/sample_names.tsv\"" \
  -d "MONITOR_FILE=\"$TMPDIR_RUN/monitor.tsv\"" \
  -d "GENOME_LENGTH=$GENOME_LENGTH" \
  -d "MU=$MU" \
  -d "TRACT_LENGTH=$TRACT_LENGTH" \
  -d "HGT_P_WITHIN=$HGT_WITHIN" \
  -d "HGT_P_CROSS=$CROSS_HGT" \
  -d "N_TERMINAL=$N_TERMINAL" \
  -d "SAMPLE_PER_TERMINAL=$SAMPLE_PER_TERMINAL" \
  -d "EXPERIMENT_GENERATIONS=$EXPERIMENT_GENERATIONS" \
  -d "MONITOR_EVERY=$MONITOR_EVERY" \
  -d "MODE=$MODE" \
  -d "S_AB=$S_AB" \
  -d "S_CD=$S_CD" \
  -d "A_ID=$A_ID" \
  -d "B_ID=$B_ID" \
  -d "C_ID=$C_ID" \
  -d "D_ID=$D_ID" \
  "$REPO_ROOT/slim/continue_experiment.slim" \
  > "$TMPDIR_RUN/slim.log" 2>&1

python "$REPO_ROOT/scripts/vcf_to_fasta.py" \
  --vcf "$TMPDIR_RUN/out.vcf" \
  --sample-map "$TMPDIR_RUN/sample_names.tsv" \
  --output "$TMPDIR_RUN/all_snps.fa" \
  --positions-output "$TMPDIR_RUN/all_snps.positions.tsv"

python "$REPO_ROOT/scripts/vcf_to_fasta.py" \
  --vcf "$TMPDIR_RUN/out.vcf" \
  --sample-map "$TMPDIR_RUN/sample_names.tsv" \
  --exclude-loci "$TMPDIR_RUN/selected_loci.tsv" \
  --output "$TMPDIR_RUN/core_snps.fa" \
  --positions-output "$TMPDIR_RUN/core_snps.positions.tsv"

python "$REPO_ROOT/scripts/export_oracle_tree.py" \
  --trees "$TMPDIR_RUN/out.trees" \
  --sample-map "$TMPDIR_RUN/sample_names.tsv" \
  --output "$TMPDIR_RUN/oracle_local_tree.nwk" \
  --position "$ORACLE_TREE_POSITION" \
  --ancestral-ne "$ANCESTRAL_NE" \
  --seed "$SEED"

for required_file in out.trees out.vcf sample_names.tsv monitor.tsv all_snps.fa core_snps.fa oracle_local_tree.nwk; do
  [[ -s "$TMPDIR_RUN/$required_file" ]] || { echo "Missing required output: $required_file" >&2; exit 1; }
done

touch "$TMPDIR_RUN/_SUCCESS"
mv "$TMPDIR_RUN" "$OUTDIR"
echo "[done] case: $OUTDIR"
