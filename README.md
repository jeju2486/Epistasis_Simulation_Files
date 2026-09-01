# Frequency-dependent SLiM–SpydrPick–KOVAR benchmark

This repository is a validation benchmark for two distinct questions:

- SpydrPick asks whether two binary loci are marginally dependent in the pooled sample.
- KOVAR 0.8.3 asks whether one unordered pair retains evidence of covariation after adjustment for shared ancestry.

Neither result proves fitness epistasis or causality. The benchmark deliberately controls the joint frequency of two seeded loci so that the statistical distinction is stable across replicates without requiring mutation-selection balance at the focal sites.

## Design

The haploid nonWF genome is 100 kb with mutation rate `2e-8`, mean homologous-transfer tract length 500 bp, and within-population HGT probability 0.02 per offspring. There is no cross-population HGT and no migration.

The absolute SLiM timeline is:

| Tick | Event |
|---:|---|
| 1 | Create ancestral population, N=10,000 |
| 5,000 | Split into two clades, N=5,000 each |
| 10,000 | Split into p1-p4, N=2,500 each; save neutral checkpoint |
| 10,000 | Seed A=10 kb and B=50 kb as 00, 01, 10, and 11 at 25% each in every terminal population |
| 10,001-30,000 | Apply frequency-dependent regulation while mutation and within-population HGT continue |
| 30,000 | Sample 250 isolates per terminal population |

There is no separate “experimental generation” coordinate. All reports use absolute ticks.

For state `h` in population `p`, the juvenile log-weight correction is

```text
K * log((target[p,h] + epsilon) / (observed[p,h] + epsilon))
```

with `K=0.25` and `epsilon=0.0002`. This is a controlled equilibrium assay, not a claim that real bacterial fitness is literally frequency dependent in this form.

### Mode 1: within-population independence, pooled dependence

| Populations | 00 | 01 | 10 | 11 |
|---|---:|---:|---:|---:|
| p1, p2 | 0.0225 | 0.1275 | 0.1275 | 0.7225 |
| p3, p4 | 0.7225 | 0.1275 | 0.1275 | 0.0225 |

A and B are independent within every population, but their clade-aligned marginal frequencies create strong pooled dependence. Expected result: high SpydrPick MI and attenuation by KOVAR if the tree adjustment works.

### Mode 2: within-population dependence, pooled independence

| Populations | 00 | 01 | 10 | 11 |
|---|---:|---:|---:|---:|
| p1, p2 | 0.25 | 0.05 | 0.45 | 0.25 |
| p3, p4 | 0.25 | 0.45 | 0.05 | 0.25 |

Each population has the same positive within-population odds ratio, but opposite A/B marginal-frequency imbalance makes the pooled table exactly uniform. Expected result: low pooled SpydrPick MI and a stronger KOVAR signal if adjustment recovers the replicated within-lineage association.

These are expectations to test, not guaranteed claims about either method. Calibration, convergence, filtering, and replicate variability remain part of the result.

## Analysis contract

- All loci with MAF at least 0.05 enter one shared A/C binary alignment.
- SpydrPick runs once with its default sample weighting, no ARACNE, zero MI threshold, and every eligible pair.
- KOVAR receives every unordered eligible pair as a two-column `u/v` file. Truth labels are never candidate inputs.
- KOVAR 0.8.3 uses `--min-maf 0.05`, explicit `--min-cell-count 5`, and `--spa-mode auto`.
- The rooted tree is the simulation genealogy at the predeclared neutral position 90 kb. No inferred IQ-TREE/ClonalFrameML stage is run.
- Because HGT gives different local genealogies along the genome, the 90-kb tree is exact locally but remains a single-tree approximation for genome-wide adjustment.
- Periodic monitoring is disabled. SLiM writes only checkpoint metadata and the four endpoint A-B tables.

## Installation

Create the simulation environment:

```bash
conda env create -f environment.yml
conda activate epistasis-sim
```

Install KOVAR 0.8.3 in the same environment:

```bash
git clone https://github.com/jeju2486/Pair-logistic-mds.git
cd Pair-logistic-mds
git switch agent/replace-with-pair-lmm-gwes
python -m pip install .
ko-variation --version
```

The final command must print `KO-Variation 0.8.3`. Then return to this repository.

## Run five jobs at a time

The following command runs the whole pipeline. Each stage permits five concurrent cases; stages remain ordered because downstream files depend on upstream completion.

```bash
CHECKPOINT_JOBS=5 \
CASE_JOBS=5 \
SPYDRPICK_JOBS=5 \
SPYDRPICK_THREADS_PER_CASE=1 \
KOVAR_JOBS=5 \
KOVAR_THREADS_PER_CASE=1 \
bash scripts/run_full_pipeline.sh
```

This uses about five CPU cores during each compute stage. On a scheduler, run it inside an allocation; invoking the Bash script directly does not activate `#SBATCH` resources.

If memory is limited during KOVAR, reduce `KOVAR_JOBS` and give each case more threads, for example:

```bash
KOVAR_JOBS=1 KOVAR_THREADS_PER_CASE=5 bash scripts/run_kovar_all.sh
```

## Stage-by-stage and restart commands

Simulation only:

```bash
CHECKPOINT_JOBS=5 CASE_JOBS=5 bash scripts/run_five_replicates.sh
```

SpydrPick only:

```bash
JOBS=5 THREADS_PER_CASE=1 bash scripts/run_spydrpick_all.sh
```

KOVAR only:

```bash
JOBS=5 THREADS_PER_CASE=1 bash scripts/run_kovar_all.sh
```

Completed simulation and SpydrPick directories are skipped by their `_SUCCESS` markers. An interrupted KOVAR score scan is resumed automatically when its `.kovar_checkpoint` is present. Re-running `run_full_pipeline.sh` is therefore the normal restart procedure.

For a quick wiring check, use the smoke configuration:

```bash
CONFIG=config/smoke.toml \
CHECKPOINT_JOBS=1 CASE_JOBS=1 SPYDRPICK_JOBS=1 KOVAR_JOBS=1 \
bash scripts/run_full_pipeline.sh
```

The smoke run validates plumbing only; its population sizes and duration are not scientifically interpretable.

## Outputs

Five neutral checkpoints are written below `checkpoints_100kb_neutral10k_mu2e8/`. Ten cases are written below `runs_100kb_frequency_dependent_mu2e8/rep_####/mode_#/`.

Important per-case outputs are:

- `focal_endpoint.tsv`: target and observed A-B state frequencies in p1-p4.
- `selected_loci.tsv`: A/B positions and mutation identifiers.
- `simulation_tree.nwk`: rooted simulation tree at 90 kb.
- `spydrpick_all_pairs/mi_vs_distance.png` and `mi_vs_distance_0_5kb.png`.
- `spydrpick_all_pairs/mi_distance_quantiles.tsv`: 1-kb-bin median, 95th, and 99th MI quantiles.
- `kovar_v083_simulation_tree/kovar_p_vs_distance.png`.
- `kovar_v083_simulation_tree/kovar_qq_by_distance.png`.
- `kovar_v083_simulation_tree/response_models.tsv` and `execution_metadata.tsv`: convergence and runtime diagnostics.

Across-case focal tables and `focal_ab_by_mode.png` replicate plots are written to `results/spydrpick_all_pairs/` and `results/kovar_v083_simulation_tree/`.
