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
- 500 bp mean HGT tract;
- within-lineage HGT probability 0.02 per offspring;
- 1,000 experimental generations;
- cross-lineage HGT probabilities 0, 0.002 and 0.02.

Cross-lineage values are probabilities, not `rho` values multiplied by genome
length. This removes the ambiguity in the older runner.

## Standing variants A-D

At the checkpoint, the model randomly selects four naturally segregating neutral
mutations subject to predeclared criteria:

- global derived-state frequency 0.30-0.70;
- frequency 0.20-0.80 in every terminal population;
- minimum physical separation 20 kb;
- initial pairwise `r^2 <= 0.02`.

The checkpoint fails visibly if no qualifying quartet exists. It does not silently
seed or choose loci using future association results. The same checkpoint and the
same A-D mutations are used for all paired continuations.

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
inputs/reference.fa
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

`oracle_local_tree.nwk` gives KOVAR true-tree information and must be labelled as an
oracle diagnostic. For the realistic primary analysis, infer a core-genome tree from
`core_snps.fa`:

```bash
bash scripts/infer_core_tree.sh runs/.../core_snps.fa runs/.../core_tree 12345
```

## All-pair comparison

The simulation does not use a top-N SpydrPick filter. Calculate MI for every eligible
pair, omit ARACNE, and supply the same pair universe to KOVAR in both directions.
Physical distance is retained so that short-range LD peaks can be classified rather
than mistaken for phylogenetic or epistatic effects.

## Cross-lineage HGT sensitivity

All checkpoints are built with zero cross-lineage HGT. Values 0, 0.002 and 0.02 are
applied only during the matched 1,000-generation continuation. This isolates the
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
