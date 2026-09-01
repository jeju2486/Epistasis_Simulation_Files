# Global-frequency SLiM–SpydrPick–KOVAR benchmark

This repository tests two separate claims: SpydrPick and KOVAR 0.8.3 should both
detect a planted global two-locus association, while KOVAR should attenuate
nonfocal covariance generated primarily by shared lineage. The focal A-B pair is
a global control; it is not itself the lineage-confounded comparison.

Neither method demonstrates biological fitness epistasis or causality by itself.
The simulation directly regulates the joint frequency of two seeded loci, A and
B, so that the statistical truth remains stable while cross-population HGT changes
the strength of the genome-wide lineage structure.

## Factorial design

Each neutral replicate produces four terminal populations, p1-p4. The same neutral
checkpoint is reused for all nine continuations in that replicate:

- modes: `0`, `1`, and `2`;
- cross-population HGT probabilities per offspring: `0`, `0.002`, and `0.02`;
- five independent replicates, giving 45 cases in total.

The haploid nonWF genome is 100 kb, with mutation rate `2e-8`, mean HGT tract
length 500 bp, and within-population HGT probability `0.02` per offspring. A
cross-HGT event chooses a donor uniformly from the other three populations. The
within- and cross-HGT probabilities are unconditional and mutually exclusive, so
their sum is the total chance of an HGT event in an offspring.

Every population has the same focal target table. Thus p1-p4 can differ throughout
the rest of the genome, but lineage membership never defines the planted A-B truth.

| Mode | Interpretation | 00 | 01 | 10 | 11 | A frequency | B frequency | Odds ratio |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 0 | balanced independent negative control | 0.25 | 0.25 | 0.25 | 0.25 | 0.50 | 0.50 | 1 |
| 1 | high-1-frequency independent control | 0.0625 | 0.1875 | 0.1875 | 0.5625 | 0.75 | 0.75 | 1 |
| 2 | high-1-frequency dependent signal | 0.1875 | 0.0625 | 0.0625 | 0.6875 | 0.75 | 0.75 | 33 |

Modes 1 and 2 deliberately share the same single-locus frequencies. Their focal
difference is therefore association, rather than a difference in A or B frequency.
Mode 0 provides a conventional balanced negative control. Expected pooled A-B MI
is zero in modes 0 and 1 and approximately 0.298 bits in mode 2.

The intended focal result is:

- modes 0 and 1: no systematic focal A-B signal after calibration;
- mode 2: strong focal A-B signal;
- increasing cross-HGT: weaker genome-wide population structure, without changing
  the target focal truth.

The intended nonfocal result is that low cross-HGT creates more lineage-driven
background covariance. SpydrPick measures that marginal covariance, whereas KOVAR
should attenuate the component explained by the phylogeny. The difference between
the methods should narrow as cross-HGT erodes lineage structure.

These are predeclared expectations to test, not guaranteed claims about either
analysis method.

## Timeline and equilibrium sampling

The neutral history is shared within each replicate:

| Absolute tick | Event |
|---:|---|
| 1 | create ancestral population, N=10,000 |
| 5,000 | split into two clades, N=5,000 each |
| 10,000 | split into p1-p4, N=2,500 each and save the neutral checkpoint |
| 10,000 | start each continuation and seed A=10 kb and B=50 kb as 00, 01, 10, 11 at 25% each |

After seeding, SLiM applies a frequency-dependent correction to juveniles. For
state `h` in population `p`, the log-weight correction is

```text
K * log((target[p,h] + epsilon) / (observed[p,h] + epsilon))
```

with `K=0.25` and `epsilon=0.0002`. This mechanism is a controlled frequency
equilibrium assay; it is not presented as a literal biological fitness model.

Every continuation runs for the same fixed duration and is sampled at absolute
tick 30,000, giving exactly 20,000 post-seeding generations. Starting 1,000 ticks
after seeding, all four population tables are checked every 100 ticks. The first
tick at which every cell is within 0.03 of its target for five consecutive checks
is recorded, but it never terminates the simulation. A case that has not reached
that criterion by sampling is retained and marked `not_reached_by_sampling` for
quality-control review.

## Analysis contract

- All loci with MAF at least 0.05 enter one shared A/C binary alignment.
- SpydrPick uses its default sample weighting, no ARACNE, zero MI threshold, and
  returns every eligible unordered pair.
- KOVAR receives every eligible unordered pair as a two-column `u/v` file. Truth
  labels are never candidate inputs.
- KOVAR 0.8.3 uses `--min-maf 0.05`, explicit `--min-cell-count 5`, and
  `--spa-mode auto`.
- The rooted covariate tree is the simulation genealogy at the predeclared neutral
  position 90 kb. It is exact locally but remains a single-tree approximation for
  a genome whose local genealogies are altered by HGT.
- Lineage labels and covariance-component categories are used only after both
  methods finish. They are never supplied as candidate or model inputs.

## Lineage-confounding evaluation

For every tested pair, the post-analysis evaluator decomposes binary covariance as

```text
total covariance
  = weighted within-population covariance
  + covariance of population allele frequencies.
```

Pairs are evaluated separately as focal A-B, focal-proximal, short-distance,
lineage-driven, within-population, or other distant pairs. The predeclared
lineage-driven category requires distance greater than 5 kb, no locus within 5 kb
of A or B, absolute pooled covariance at least 0.01, and at least 80% of the
absolute component magnitude attributable to between-population covariance.

For an equal discovery budget, the evaluator compares the enrichment of these
lineage-driven pairs in each method's top 1%. KOVAR Bonferroni counts are also
reported independently. These categories are evaluation labels, not simulation
truth supplied to either method.

## Installation

Create and activate the simulation environment:

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

This command runs simulation, supplementary SpydrPick, primary KOVAR, and plots.
Stages remain ordered, while up to five cases run concurrently inside each stage.

```bash
CHECKPOINT_JOBS=5 \
CASE_JOBS=5 \
SPYDRPICK_JOBS=5 \
SPYDRPICK_THREADS_PER_CASE=1 \
KOVAR_JOBS=5 \
KOVAR_THREADS_PER_CASE=1 \
bash scripts/run_full_pipeline.sh
```

If KOVAR memory use is limiting, run fewer cases and give a case more threads:

```bash
KOVAR_JOBS=1 KOVAR_THREADS_PER_CASE=5 bash scripts/run_kovar_all.sh
```

Stage-specific commands are:

```bash
CHECKPOINT_JOBS=5 CASE_JOBS=5 bash scripts/run_five_replicates.sh
JOBS=5 THREADS_PER_CASE=1 bash scripts/run_spydrpick_all.sh
JOBS=5 THREADS_PER_CASE=1 bash scripts/run_kovar_all.sh
```

Completed simulation and SpydrPick directories are skipped using `_SUCCESS`.
Interrupted KOVAR scans resume from `.kovar_checkpoint`. Re-running the full
pipeline is therefore the normal restart procedure.

For a wiring-only smoke test:

```bash
CONFIG=config/smoke.toml \
CHECKPOINT_JOBS=1 CASE_JOBS=1 SPYDRPICK_JOBS=1 KOVAR_JOBS=1 \
bash scripts/run_full_pipeline.sh
```

The smoke population sizes, tolerance, and duration are not scientifically
interpretable.

## Outputs

Five neutral checkpoints are written below
`checkpoints_100kb_neutral10k_mu2e8/`. The 45 flat-generation case directories are
written as:

```text
runs_100kb_global_frequency_hgt_mu2e8/
  rep_####/
    cross_0|cross_0p002|cross_0p02/
      mode_0|mode_1|mode_2/
```

There are no `gen_0300`, `gen_0400`, or other generation subdirectories.

Important per-case outputs include:

- `equilibrium_status.tsv`: fixed sampling tick, first equilibrium tick, and final deviation;
- `focal_monitor.tsv`: every equilibrium check for p1-p4;
- `focal_endpoint.tsv`: target and observed A-B state frequencies at sampling;
- `selected_loci.tsv`: focal positions, mutation identifiers, mode, and cross-HGT;
- `simulation_tree.nwk`: rooted simulation genealogy at 90 kb;
- `spydrpick_all_pairs/mi_vs_distance.png` and MI summaries;
- `kovar_v083_simulation_tree/kovar_p_vs_distance.png` and distance-stratified QQ plot;
- KOVAR convergence and runtime diagnostics in `response_models.tsv` and
  `execution_metadata.tsv`.
- `lineage_pair_metrics.tsv.gz`: pair-level covariance decomposition and
  method-specific ranks for evaluation only.

Across-case focal tables and mode-by-HGT plots are written below
`results/spydrpick_all_pairs/` and `results/kovar_v083_simulation_tree/`.
The lineage-confounding category summary, category-count plot, and equal-budget
comparison plot are written below `results/lineage_confounding/`.
