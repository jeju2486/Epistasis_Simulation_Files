#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SLIM_BIN="${SLIM_BIN:-slim}"

OUTDIR="" CHECKPOINT_DIR="" REF_FILE="" SEED="" MODE="" CROSS_HGT=""
GENOME_LENGTH="" MU="" TRACT_LENGTH="" HGT_WITHIN="" N_TERMINAL=""
SAMPLE_PER_TERMINAL="" EXPERIMENT_GENERATIONS="" S_AB="" S_CD="" MONITOR_EVERY=""
ANCESTRAL_NE="" ORACLE_TREE_POSITION=""
A_POSITION="" B_POSITION="" C_POSITION="" D_POSITION=""

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
    --a-position) A_POSITION="$2"; shift 2 ;;
    --b-position) B_POSITION="$2"; shift 2 ;;
    --c-position) C_POSITION="$2"; shift 2 ;;
    --d-position) D_POSITION="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

required=(OUTDIR CHECKPOINT_DIR REF_FILE SEED MODE CROSS_HGT GENOME_LENGTH MU TRACT_LENGTH HGT_WITHIN N_TERMINAL SAMPLE_PER_TERMINAL EXPERIMENT_GENERATIONS S_AB S_CD MONITOR_EVERY ANCESTRAL_NE ORACLE_TREE_POSITION A_POSITION B_POSITION C_POSITION D_POSITION)
for name in "${required[@]}"; do
  [[ -n "${!name}" ]] || { echo "Missing required value: $name" >&2; exit 2; }
done

case "$MODE" in
  0) ACTIVE_PAIR="AB_and_CD_neutral"; SEEDING_DESIGN="balanced_16_haplotype_cycle" ;;
  1) ACTIVE_PAIR="AB"; SEEDING_DESIGN="balanced_four_haplotype_cycle" ;;
  2) ACTIVE_PAIR="CD"; SEEDING_DESIGN="balanced_four_haplotype_cycle" ;;
  *) echo "Mode must be 0, 1 or 2: $MODE" >&2; exit 2 ;;
esac

[[ "$OUTDIR" = /* ]] || OUTDIR="$REPO_ROOT/$OUTDIR"
[[ "$CHECKPOINT_DIR" = /* ]] || CHECKPOINT_DIR="$REPO_ROOT/$CHECKPOINT_DIR"
[[ "$REF_FILE" = /* ]] || REF_FILE="$REPO_ROOT/$REF_FILE"

[[ -f "$CHECKPOINT_DIR/_SUCCESS" ]] || { echo "Checkpoint is incomplete: $CHECKPOINT_DIR" >&2; exit 2; }
[[ -s "$REF_FILE" ]] || { echo "Reference not found: $REF_FILE" >&2; exit 2; }

if [[ -f "$OUTDIR/_SUCCESS" ]]; then
  echo "[skip] case complete: $OUTDIR"
  exit 0
fi

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

cp "$CHECKPOINT_DIR/checkpoint_metrics.tsv" "$TMPDIR_RUN/checkpoint_metrics.tsv"

cat > "$TMPDIR_RUN/params.tsv" <<EOF
parameter	value
seed	$SEED
mode	$MODE
active_pair	$ACTIVE_PAIR
cross_hgt_probability	$CROSS_HGT
checkpoint	$CHECKPOINT_DIR/checkpoint.trees
focal_seeding	experimental_generation_0_mode_specific
seeding_design	$SEEDING_DESIGN
A_position	$A_POSITION
B_position	$B_POSITION
C_position	$C_POSITION
D_position	$D_POSITION
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
  -d "LOCI_FILE=\"$TMPDIR_RUN/selected_loci.tsv\"" \
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
  -d "A_POSITION=$A_POSITION" \
  -d "B_POSITION=$B_POSITION" \
  -d "C_POSITION=$C_POSITION" \
  -d "D_POSITION=$D_POSITION" \
  "$REPO_ROOT/slim/continue_experiment.slim" \
  > "$TMPDIR_RUN/slim.log" 2>&1

python "$REPO_ROOT/scripts/validate_selected_loci.py" \
  "$TMPDIR_RUN/selected_loci.tsv" \
  --mode "$MODE" \
  --a-position "$A_POSITION" \
  --b-position "$B_POSITION" \
  --c-position "$C_POSITION" \
  --d-position "$D_POSITION"

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

for required_file in out.trees out.vcf sample_names.tsv selected_loci.tsv monitor.tsv all_snps.fa core_snps.fa oracle_local_tree.nwk; do
  [[ -s "$TMPDIR_RUN/$required_file" ]] || { echo "Missing required output: $required_file" >&2; exit 1; }
done

touch "$TMPDIR_RUN/_SUCCESS"
mv "$TMPDIR_RUN" "$OUTDIR"
echo "[done] case: $OUTDIR"
