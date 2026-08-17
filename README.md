# Epistasis simulation files

Forward-time bacterial simulations for evaluating KOVAR against SpydrPick/PAN-GWES.
The primary benchmark separates repeated within-lineage covariation from marginal
association caused by stable phylogenetic structure.

This repository is validation infrastructure; it is not part of the KOVAR runtime.

## Scientific contrast

Each replicate first builds one neutral, structured population and saves a checkpoint.
The checkpoint is then continued through two matched arms:

| Mode | Fitness model | Expected result |
|---|---|---|
| 0 | A, B, C and D remain neutral | No reproducible unlinked peak, apart from stochastic/physical LD |
| 1 | A-B compatibility plus independent clade-aligned effects for C and D | SpydrPick detects A-B and C-D; KOVAR should retain A-B and demote C-D |

In Mode 1 the individual log-fitness contribution is

```text
-s_ab * I(A != B) + z_lineage * s_cd * (C + D)
```

where `z_lineage=+1` for p1/p2 and `-1` for p3/p4. Therefore A-B has a
non-zero interaction contrast, while C-D has no pairwise interaction:

```text
epsilon_AB = 2 * s_ab
epsilon_CD = 0
```

C-D is a false positive only relative to the target "repeated association beyond
shared ancestry". It is genuinely marginally dependent after lineage-aligned
selection.

## Population design

```text
                 ancestral p0
                  /       \
                pL         pR
               /  \       /  \
             p1    p2   p3    p4
```

Default pilot values are defined in [`config/pilot.toml`](config/pilot.toml):

- ancestral census: 10,000;
- two internal clades: 5,000 each;
- four terminal populations: 2,500 each;
- 250 sampled genomes per terminal population (`n=1,000`);
- 100 kb haploid nucleotide genome;
- nucleotide mutation rate `2e-8` per site per generation;
- 500 bp mean HGT tract;
- within-lineage HGT probability 0.02 per offspring;
- 5,000 experimental generations;
- mode-1 coefficients `S_AB=0.003` and `S_CD=0.003`;
- cross-lineage HGT probabilities 0, 0.002 and 0.02.

Cross-lineage values are probabilities, not `rho` values multiplied by genome
length. This removes the ambiguity in the older runner.

## Seeded standing variants A-D

After demographic burn-in, each checkpoint receives four neutral nucleotide mutations
at positions fixed before analysis:

| Locus | Zero-based position |
|---|---:|
| A | 10,000 |
| B | 20,000 |
| C | 60,000 |
| D | 85,000 |

Consequently A-B spans 10 kb and C-D spans 25 kb. Their unequal distances make the
two truth pairs visually distinguishable, while both distances remain much larger than
the 500 bp mean HGT tract.

Within every terminal population, individuals are randomly ordered and assigned a
cycle through all 16 four-locus haplotypes. With 2,500 individuals, every haplotype
occurs 156 or 157 times per population, giving focal prevalences of approximately 0.5
and negligible starting pairwise LD. Background mutation is disabled at the four focal
positions, but HGT can transfer the focal alleles. The same seeded checkpoint supplies
all matched Mode 0 and Mode 1 continuations. Neither SpydrPick nor KOVAR receives the
truth labels or a restricted truth-pair list.

## Workflow

The workflow is manifest-driven and does not require one sweep script or one
`screen` session per parameter combination.

### 1. Create the software environment

```bash
micromamba create -f environment.yml
micromamba activate epistasis-sim
```

SLiM 5.2 is the validated target. The scripts require Python 3.11+ and Bash.

### 2. Generate the deterministic reference and manifests

```bash
python scripts/prepare_inputs.py --config config/pilot.toml
python scripts/build_manifest.py --config config/pilot.toml
```

This creates:

```text
inputs/reference_100kb.fa
manifests/checkpoints.tsv
manifests/cases.tsv
```

### 3. Build neutral checkpoints

```bash
python scripts/run_manifest.py \
  --manifest manifests/checkpoints.tsv \
  --stage checkpoint \
  --jobs 4
```

### 4. Run paired continuations

```bash
python scripts/run_manifest.py \
  --manifest manifests/cases.tsv \
  --stage case \
  --jobs 8
```

For 50 replicates and three cross-lineage HGT values, the design runs 50 demographic
checkpoints and 300 short continuations rather than 300 complete demographic simulations.

### 5. Inspect status

```bash
python scripts/show_status.py manifests/checkpoints.tsv manifests/cases.tsv
```

A case is complete only when its `_SUCCESS` marker exists. Partial outputs are never
treated as successful runs.

### Five replicates per setting

The convenience launcher uses the full pilot parameters but limits the design to five
matched replicates. With two modes and three cross-lineage HGT values, it creates five
checkpoints and 30 continuations:

```bash
micromamba activate epistasis-sim
bash scripts/run_five_replicates.sh
```

Parallelism can be adjusted without editing the script:

```bash
CHECKPOINT_JOBS=1 CASE_JOBS=4 bash scripts/run_five_replicates.sh
```

The launcher is restart-safe: completed output directories containing `_SUCCESS` are
skipped. Its settings are in [`config/five_replicates.toml`](config/five_replicates.toml).

## Outputs

Each checkpoint contains:

```text
checkpoint.trees       resumable population and ancestry
selected_loci.tsv      A-D mutation IDs, positions and frequencies
checkpoint_metrics.tsv population sizes and differentiation diagnostics
params.tsv              exact parameters and seed
```

Each continuation contains:

```text
out.trees               final recorded ancestry
out.vcf                 sampled haploid genotypes
all_snps.fa             sampled SNP alignment
core_snps.fa            alignment excluding A-D
sample_names.tsv        VCF label, population and pedigree mapping
monitor.tsv             time series for A-B and C-D
params.tsv              exact parameters and seed
oracle_local_tree.nwk   local true-ancestry tree (diagnostic)
```

`oracle_local_tree.nwk` is the simulation-produced local genealogy at position 50 kb.
This validation workflow uses it as oracle covariance and must label it accordingly.
IQ-TREE is not run by the primary pipeline.

## All-pair comparison

The simulation does not use a top-N SpydrPick filter. Calculate MI for every eligible
pair, omit ARACNE, and supply the same pair universe to KOVAR in both directions.
Physical distance is retained so that short-range LD peaks can be classified rather
than mistaken for phylogenetic or epistatic effects.

Install the updated environment, then run every completed case and make the plots:

```bash
micromamba env update -n epistasis-sim -f environment.yml
micromamba activate epistasis-sim
JOBS=2 THREADS_PER_CASE=4 bash scripts/run_spydrpick_all.sh
```

This example runs two cases concurrently and gives each SpydrPick process four
threads, for at most eight compute threads. Before pair construction, the shared
binary matrix is filtered at `MAF >= 0.05`. The runner then uses `--mi-threshold=0`,
`--no-aracne`, and `--no-filter-alignment` on that already-filtered matrix and
verifies that the output contains exactly `L(L-1)/2` pairs. Default SpydrPick sample reweighting is retained for the
primary comparison. SpydrPick and KOVAR receive the exact same binary A/C matrix:
reference state is absence and any non-reference state is presence. Because the
packaged SpydrPick 1.2.0 executable does not support
`--mappings-list`, the runner converts its one-based alignment-column output to
zero-based columns and physical distances using `all_snps.positions.tsv`. The
optional unweighted sensitivity analysis is:

```bash
python scripts/run_spydrpick_manifest.py \
  --manifest manifests/cases.tsv \
  --jobs 2 \
  --threads-per-case 4 \
  --unweighted
python scripts/plot_spydrpick.py --manifest manifests/cases.tsv --unweighted
```

Each case receives `spydrpick_all_pairs/spydrpick.edges.gz`, a small
`truth_pairs.tsv` with the MI and all-pair ranks of A-B and C-D, the exact command
metadata, log, and `mi_vs_distance.png`. The PAN-GWES-style figure plots physical
SNP-pair distance in kb against MI. Ordinary eligible pairs are grey; AB, AC, AD,
BC, BD and CD use stable focal-pair colours. At most 250,000 ordinary points are
drawn to keep the raster figure tractable, but this is display-only thinning: all
pairs remain in `spydrpick.edges.gz`, and focal pairs are never thinned. The rank
curve and aggregate rank boxplot are not generated, while the diagnostic ranks
remain in `results/spydrpick_all_pairs/truth_pair_ranks.tsv`. If a selected locus
is lost or fixed before sampling, its truth pair is retained with
`candidate_status=locus_absent` rather than being misreported as a statistical
failure. A sampled truth locus below the predeclared 5% MAF threshold is recorded
separately as `maf_filtered`.

## KOVAR 0.8.1 analysis

KOVAR must already be installed in the active environment. Confirm the exact
experimental version before launching:

```bash
ko-variation --version
# KO-Variation 0.8.1
```

After the exhaustive SpydrPick stage completes, run KOVAR across the same candidate
universe:

```bash
JOBS=1 THREADS_PER_CASE=8 bash scripts/run_kovar_all.sh
```

The primary default uses `oracle_local_tree.nwk`, tests both directions, uses 5% MAF
and a minimum four-cell count of five, disables experimental
SPA, and disables full alternative-model refits. The score tests and directional
p-values are still produced; disabling refits avoids potentially tens of thousands
of expensive effect-estimation fits in an exhaustive scan. Override settings without
editing the runner, for example:

```bash
JOBS=1 THREADS_PER_CASE=8 TREE_MODE=oracle \
MIN_MAF=0.05 MIN_CELL_COUNT=5 SPA_MODE=off FULL_REFIT_P=0 \
bash scripts/run_kovar_all.sh
```

The oracle local tree is deliberate for this simulation validation. It evaluates KOVAR
under known covariance rather than adding tree-inference error to the primary contrast.

KOVAR 0.8.1 cannot read the compressed SpydrPick edge file directly. Each case runner
therefore materializes a temporary plain pair file, runs KOVAR, and removes that copy
after success. It reuses the exact binary FASTA analyzed by SpydrPick rather than
performing an independent conversion. `run_kovar_all.sh` fixes `--max-pairs 0`, so
every MAF-eligible SpydrPick pair is tested by KOVAR in both directions.

Results are written under each case's `kovar_v081_oracle/` directory. The exact commands,
tree mode, candidate-universe label, and pair count are recorded in
`run_metadata.json`; KOVAR's results are in `kovar_v081_oracle/results/`.
Each case also receives `kovar_p_vs_distance.png`, with the same distance axis and
focal-pair colour mapping as the SpydrPick figure. Its y-axis is directional
`-log10(p_primary)`; upward and downward triangles distinguish `u_predicts_v` and
`v_predicts_u`. A dashed line shows the directional Bonferroni 0.05 threshold.
Display-only thinning is applied to ordinary pairs, never to focal pairs or the
underlying exhaustive KOVAR result table.

## Complete five-replicate workflow

Run the full simulation, exhaustive SpydrPick scan, plots, and exhaustive KOVAR scan
with one command:

```bash
conda activate epistasis-sim
bash scripts/run_full_pipeline.sh
```

Parallelism may be changed without altering statistical scope:

```bash
CHECKPOINT_JOBS=2 \
CASE_JOBS=6 \
SPYDRPICK_JOBS=2 \
SPYDRPICK_THREADS_PER_CASE=4 \
KOVAR_JOBS=1 \
KOVAR_THREADS_PER_CASE=8 \
bash scripts/run_full_pipeline.sh
```

The launcher defaults to `config/five_replicates.toml`, the 5% shared MAF filter,
all eligible pairs, both KOVAR directions, and the simulation-produced local tree.
It is restart-safe because successful stage directories are skipped.

## Cross-lineage HGT sensitivity

All checkpoints are built with zero cross-lineage HGT. Values 0, 0.002 and 0.02 are
applied only during the matched 5,000-generation continuation. This isolates the
effect of subsequent mixing while holding the initial population structure fixed.

A distinct future experiment would allow cross-lineage HGT throughout population
formation; that design requires separate checkpoints for each HGT value.

## Validation

Run local Python tests with:

```bash
python -m unittest discover -s tests -v
```

On a machine with SLiM and ShellCheck:

```bash
slim -c slim/build_checkpoint.slim
slim -c slim/continue_experiment.slim
shellcheck scripts/*.sh
```

The GitHub Actions workflow performs these checks in the pinned micromamba environment.
