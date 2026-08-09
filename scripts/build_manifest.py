#!/usr/bin/env python3
"""Expand a compact TOML design into checkpoint and continuation manifests."""

from __future__ import annotations

import argparse
from pathlib import Path

from simflow import (
    REPO_ROOT,
    deterministic_seed,
    load_config,
    probability_slug,
    repo_path,
    write_tsv_atomic,
)


def relative(path: Path) -> str:
    try:
        return path.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(path)


def build(config: dict) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    design = config["design"]
    genome = config["genome"]
    population = config["population"]
    loci = config["loci"]
    experiment = config["experiment"]
    postprocess = config["postprocess"]

    checkpoint_root = repo_path(config["paths"]["checkpoint_root"])
    run_root = repo_path(config["paths"]["run_root"])
    reference = relative(repo_path(config["paths"]["reference"]))
    master_seed = int(design["master_seed"])

    checkpoints: list[dict[str, object]] = []
    cases: list[dict[str, object]] = []

    for replicate in range(1, int(design["replicates"]) + 1):
        rep = f"rep_{replicate:04d}"
        checkpoint_dir = checkpoint_root / rep
        checkpoints.append(
            {
                "replicate": replicate,
                "checkpoint_id": rep,
                "seed": deterministic_seed(master_seed, "checkpoint", replicate),
                "out_dir": relative(checkpoint_dir),
                "reference": reference,
                "genome_length": genome["length"],
                "mutation_rate": genome["mutation_rate"],
                "tract_length": genome["mean_hgt_tract"],
                "within_hgt_probability": genome["within_hgt_probability"],
                "ancestral_size": population["ancestral_size"],
                "clade_size": population["clade_size"],
                "terminal_size": population["terminal_size"],
                "ancestral_generations": population["ancestral_generations"],
                "deep_generations": population["deep_clade_generations"],
                "terminal_generations": population["terminal_generations"],
                "global_frequency_min": loci["global_frequency_min"],
                "global_frequency_max": loci["global_frequency_max"],
                "lineage_frequency_min": loci["lineage_frequency_min"],
                "lineage_frequency_max": loci["lineage_frequency_max"],
                "minimum_distance": loci["minimum_distance"],
            }
        )

        for cross_probability in design["cross_hgt_probabilities"]:
            cross_probability = float(cross_probability)
            cross_slug = probability_slug(cross_probability)
            for mode in design["modes"]:
                mode = int(mode)
                case_id = f"{rep}__cross_{cross_slug}__mode_{mode}"
                cases.append(
                    {
                        "case_id": case_id,
                        "replicate": replicate,
                        "checkpoint_id": rep,
                        "checkpoint_dir": relative(checkpoint_dir),
                        "out_dir": relative(run_root / rep / f"cross_{cross_slug}" / f"mode_{mode}"),
                        "reference": reference,
                        "mode": mode,
                        "cross_hgt_probability": cross_probability,
                        "seed": deterministic_seed(
                            master_seed, "continuation", replicate, cross_probability, mode
                        ),
                        "genome_length": genome["length"],
                        "mutation_rate": genome["mutation_rate"],
                        "tract_length": genome["mean_hgt_tract"],
                        "within_hgt_probability": genome["within_hgt_probability"],
                        "terminal_size": population["terminal_size"],
                        "sample_per_terminal": population["sample_per_terminal"],
                        "experiment_generations": experiment["generations"],
                        "s_ab": experiment["s_ab"],
                        "s_cd": experiment["s_cd"],
                        "monitor_every": experiment["monitor_every"],
                        "ancestral_ne": postprocess["ancestral_ne"],
                        "oracle_tree_position": postprocess["oracle_tree_position"],
                    }
                )

    return checkpoints, cases


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    checkpoints, cases = build(config)
    manifest_root = repo_path(config["paths"]["manifest_root"])

    write_tsv_atomic(manifest_root / "checkpoints.tsv", checkpoints, list(checkpoints[0]))
    write_tsv_atomic(manifest_root / "cases.tsv", cases, list(cases[0]))
    print(f"[done] {len(checkpoints)} checkpoints -> {manifest_root / 'checkpoints.tsv'}")
    print(f"[done] {len(cases)} continuations -> {manifest_root / 'cases.tsv'}")


if __name__ == "__main__":
    main()
