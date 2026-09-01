#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SLIM_BIN="${SLIM_BIN:-slim}"
OUTDIR="" CHECKPOINT_DIR="" REF_FILE="" SEED="" MODE=""
GENOME_LENGTH="" MU="" TRACT_LENGTH="" HGT_WITHIN="" N_TERMINAL=""
SAMPLE_PER_TERMINAL="" END_TICK="" FDS_STRENGTH="" FDS_EPSILON=""
ANCESTRAL_NE="" TREE_POSITION="" A_POSITION="" B_POSITION=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    --checkpoint-dir) CHECKPOINT_DIR="$2"; shift 2 ;;
    --reference) REF_FILE="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --genome-length) GENOME_LENGTH="$2"; shift 2 ;;
    --mu) MU="$2"; shift 2 ;;
    --tract-length) TRACT_LENGTH="$2"; shift 2 ;;
    --within-hgt) HGT_WITHIN="$2"; shift 2 ;;
    --terminal-size) N_TERMINAL="$2"; shift 2 ;;
    --sample-per-terminal) SAMPLE_PER_TERMINAL="$2"; shift 2 ;;
    --end-tick) END_TICK="$2"; shift 2 ;;
    --fds-strength) FDS_STRENGTH="$2"; shift 2 ;;
    --fds-epsilon) FDS_EPSILON="$2"; shift 2 ;;
    --ancestral-ne) ANCESTRAL_NE="$2"; shift 2 ;;
    --tree-position) TREE_POSITION="$2"; shift 2 ;;
    --a-position) A_POSITION="$2"; shift 2 ;;
    --b-position) B_POSITION="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

required=(OUTDIR CHECKPOINT_DIR REF_FILE SEED MODE GENOME_LENGTH MU TRACT_LENGTH HGT_WITHIN N_TERMINAL SAMPLE_PER_TERMINAL END_TICK FDS_STRENGTH FDS_EPSILON ANCESTRAL_NE TREE_POSITION A_POSITION B_POSITION)
for name in "${required[@]}"; do
  [[ -n "${!name}" ]] || { echo "Missing required value: $name" >&2; exit 2; }
done
case "$MODE" in 1|2) ;; *) echo "Mode must be 1 or 2" >&2; exit 2 ;; esac
[[ "$OUTDIR" = /* ]] || OUTDIR="$REPO_ROOT/$OUTDIR"
[[ "$CHECKPOINT_DIR" = /* ]] || CHECKPOINT_DIR="$REPO_ROOT/$CHECKPOINT_DIR"
[[ "$REF_FILE" = /* ]] || REF_FILE="$REPO_ROOT/$REF_FILE"
[[ -f "$CHECKPOINT_DIR/_SUCCESS" ]] || { echo "Incomplete checkpoint: $CHECKPOINT_DIR" >&2; exit 2; }
[[ -s "$REF_FILE" ]] || { echo "Reference not found: $REF_FILE" >&2; exit 2; }
if [[ -f "$OUTDIR/_SUCCESS" ]]; then echo "[skip] case complete: $OUTDIR"; exit 0; fi

mkdir -p "$(dirname "$OUTDIR")"
LOCKDIR="${OUTDIR}.lock"
mkdir "$LOCKDIR" 2>/dev/null || { echo "Case is locked: $OUTDIR" >&2; exit 3; }
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
if [[ -e "$OUTDIR" ]]; then mv "$OUTDIR" "${OUTDIR}.incomplete.$(date -u +%Y%m%dT%H%M%SZ)"; fi
cp "$CHECKPOINT_DIR/checkpoint_metrics.tsv" "$TMPDIR_RUN/checkpoint_metrics.tsv"

cat > "$TMPDIR_RUN/params.tsv" <<EOF
parameter	value
seed	$SEED
mode	$MODE
checkpoint	$CHECKPOINT_DIR/checkpoint.trees
end_tick	$END_TICK
fds_strength	$FDS_STRENGTH
fds_epsilon	$FDS_EPSILON
within_hgt_probability	$HGT_WITHIN
cross_hgt_probability	0
A_position	$A_POSITION
B_position	$B_POSITION
tree_position	$TREE_POSITION
EOF

"$SLIM_BIN" -s "$SEED" \
  -d "REF_FILE=\"$REF_FILE\"" -d "CHECKPOINT_FILE=\"$CHECKPOINT_DIR/checkpoint.trees\"" \
  -d "TREE_OUT=\"$TMPDIR_RUN/out.trees\"" -d "VCF_OUT=\"$TMPDIR_RUN/out.vcf\"" \
  -d "SAMPLE_FILE=\"$TMPDIR_RUN/sample_names.tsv\"" \
  -d "LOCI_FILE=\"$TMPDIR_RUN/selected_loci.tsv\"" \
  -d "ENDPOINT_FILE=\"$TMPDIR_RUN/focal_endpoint.tsv\"" \
  -d "GENOME_LENGTH=$GENOME_LENGTH" -d "MU=$MU" \
  -d "TRACT_LENGTH=$TRACT_LENGTH" -d "HGT_P_WITHIN=$HGT_WITHIN" \
  -d "N_TERMINAL=$N_TERMINAL" -d "SAMPLE_PER_TERMINAL=$SAMPLE_PER_TERMINAL" \
  -d "END_TICK=$END_TICK" -d "MODE=$MODE" \
  -d "FDS_STRENGTH=$FDS_STRENGTH" -d "FDS_EPSILON=$FDS_EPSILON" \
  -d "A_POSITION=$A_POSITION" -d "B_POSITION=$B_POSITION" \
  "$REPO_ROOT/slim/continue_experiment.slim" > "$TMPDIR_RUN/slim.log" 2>&1

python "$REPO_ROOT/scripts/validate_selected_loci.py" "$TMPDIR_RUN/selected_loci.tsv" \
  --mode "$MODE" --a-position "$A_POSITION" --b-position "$B_POSITION"
python "$REPO_ROOT/scripts/vcf_to_fasta.py" --vcf "$TMPDIR_RUN/out.vcf" \
  --sample-map "$TMPDIR_RUN/sample_names.tsv" --output "$TMPDIR_RUN/all_snps.fa" \
  --positions-output "$TMPDIR_RUN/all_snps.positions.tsv"
python "$REPO_ROOT/scripts/export_oracle_tree.py" --trees "$TMPDIR_RUN/out.trees" \
  --sample-map "$TMPDIR_RUN/sample_names.tsv" --output "$TMPDIR_RUN/simulation_tree.nwk" \
  --position "$TREE_POSITION" --ancestral-ne "$ANCESTRAL_NE" --seed "$SEED"

for file in out.trees out.vcf sample_names.tsv selected_loci.tsv focal_endpoint.tsv all_snps.fa all_snps.positions.tsv simulation_tree.nwk; do
  [[ -s "$TMPDIR_RUN/$file" ]] || { echo "Missing $file" >&2; exit 1; }
done
touch "$TMPDIR_RUN/_SUCCESS"
mv "$TMPDIR_RUN" "$OUTDIR"
echo "[done] case: $OUTDIR"
