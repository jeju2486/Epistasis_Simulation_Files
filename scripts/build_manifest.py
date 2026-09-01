#!/usr/bin/env python3
"""Expand the minimal two-mode design into checkpoint and case manifests."""

from __future__ import annotations

import argparse
from pathlib import Path

from simflow import REPO_ROOT, deterministic_seed, load_config, repo_path, write_tsv_atomic


MODE_LABEL = {
    1: "within_independent_pooled_dependent",
    2: "within_dependent_pooled_independent",
}


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
    frequency_dependence = config["frequency_dependence"]
    postprocess = config["postprocess"]

    genome_length = int(genome["length"])
    a_position = int(loci["a_position"])
    b_position = int(loci["b_position"])
    tree_position = int(postprocess["tree_position"])
    if not 0 <= a_position < b_position < genome_length:
        raise ValueError("require 0 <= A < B < genome length")
    if not 0 <= tree_position < genome_length:
        raise ValueError("tree position must lie within the genome")
    if tree_position in {a_position, b_position}:
        raise ValueError("tree position must not equal A or B")

    deep_split_tick = int(population["deep_split_tick"])
    terminal_split_tick = int(population["terminal_split_tick"])
    end_tick = int(population["end_tick"])
    if not 1 < deep_split_tick < terminal_split_tick < end_tick:
        raise ValueError("require 1 < deep split < terminal split < end tick")

    modes = [int(mode) for mode in design["modes"]]
    if not modes or any(mode not in MODE_LABEL for mode in modes):
        raise ValueError("modes must contain only 1 and 2")
    if len(set(modes)) != len(modes):
        raise ValueError("modes contains duplicates")

    strength = float(frequency_dependence["strength"])
    epsilon = float(frequency_dependence["epsilon"])
    if strength <= 0.0 or epsilon <= 0.0:
        raise ValueError("frequency-dependence strength and epsilon must be positive")

    checkpoint_root = repo_path(config["paths"]["checkpoint_root"])
    run_root = repo_path(config["paths"]["run_root"])
    reference = relative(repo_path(config["paths"]["reference"]))
    master_seed = int(design["master_seed"])
    checkpoints: list[dict[str, object]] = []
    cases: list[dict[str, object]] = []

    for replicate in range(1, int(design["replicates"]) + 1):
        rep = f"rep_{replicate:04d}"
        checkpoint_dir = checkpoint_root / rep
        checkpoint_id = rep
        checkpoints.append(
            {
                "replicate": replicate,
                "checkpoint_id": checkpoint_id,
                "seed": deterministic_seed(master_seed, "checkpoint", replicate),
                "out_dir": relative(checkpoint_dir),
                "reference": reference,
                "genome_length": genome_length,
                "mutation_rate": genome["mutation_rate"],
                "tract_length": genome["mean_hgt_tract"],
                "within_hgt_probability": genome["within_hgt_probability"],
                "ancestral_size": population["ancestral_size"],
                "clade_size": population["clade_size"],
                "terminal_size": population["terminal_size"],
                "deep_split_tick": deep_split_tick,
                "terminal_split_tick": terminal_split_tick,
                "a_position": a_position,
                "b_position": b_position,
            }
        )
        for mode in modes:
            case_id = f"{rep}__mode_{mode}"
            cases.append(
                {
                    "case_id": case_id,
                    "replicate": replicate,
                    "checkpoint_id": checkpoint_id,
                    "checkpoint_dir": relative(checkpoint_dir),
                    "out_dir": relative(run_root / rep / f"mode_{mode}"),
                    "reference": reference,
                    "mode": mode,
                    "regime": MODE_LABEL[mode],
                    "seed": deterministic_seed(master_seed, "continuation", replicate, mode),
                    "genome_length": genome_length,
                    "mutation_rate": genome["mutation_rate"],
                    "tract_length": genome["mean_hgt_tract"],
                    "within_hgt_probability": genome["within_hgt_probability"],
                    "terminal_size": population["terminal_size"],
                    "sample_per_terminal": population["sample_per_terminal"],
                    "end_tick": end_tick,
                    "fds_strength": strength,
                    "fds_epsilon": epsilon,
                    "ancestral_ne": postprocess["ancestral_ne"],
                    "tree_position": tree_position,
                    "a_position": a_position,
                    "b_position": b_position,
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
    print(f"[done] {len(cases)} cases -> {manifest_root / 'cases.tsv'}")


if __name__ == "__main__":
    main()
