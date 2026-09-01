#!/usr/bin/env python3
"""Expand the three-mode, three-cross-HGT benchmark manifests."""

from __future__ import annotations

import argparse
from pathlib import Path

from simflow import REPO_ROOT, deterministic_seed, load_config, repo_path, write_tsv_atomic


MODE_LABEL = {
    0: "balanced_independent_negative_control",
    1: "global_high_frequency_independent",
    2: "global_high_frequency_dependent",
}


def probability_label(value: float) -> str:
    """Return a stable path-safe decimal label such as 0, 0p002, or 0p02."""
    return format(value, ".15g").replace(".", "p")


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
    equilibrium = config["equilibrium"]
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
        raise ValueError("modes must contain only 0, 1, and 2")
    if len(set(modes)) != len(modes):
        raise ValueError("modes contains duplicates")

    cross_hgt_probabilities = [float(value) for value in design["cross_hgt_probabilities"]]
    if not cross_hgt_probabilities or any(
        value < 0.0 or value > 1.0 for value in cross_hgt_probabilities
    ):
        raise ValueError("cross-HGT probabilities must lie between zero and one")
    if len(set(cross_hgt_probabilities)) != len(cross_hgt_probabilities):
        raise ValueError("cross-HGT probabilities contain duplicates")
    within_hgt_probability = float(genome["within_hgt_probability"])
    if any(value + within_hgt_probability > 1.0 for value in cross_hgt_probabilities):
        raise ValueError("within- and cross-HGT probabilities must sum to at most one")

    strength = float(frequency_dependence["strength"])
    epsilon = float(frequency_dependence["epsilon"])
    if strength <= 0.0 or epsilon <= 0.0:
        raise ValueError("frequency-dependence strength and epsilon must be positive")

    monitor_every = int(equilibrium["monitor_every"])
    minimum_ticks = int(equilibrium["minimum_ticks"])
    stable_checks = int(equilibrium["stable_checks"])
    tolerance = float(equilibrium["tolerance"])
    if monitor_every < 1 or minimum_ticks < 0 or stable_checks < 1:
        raise ValueError("equilibrium timing values are invalid")
    if not 0.0 < tolerance < 1.0:
        raise ValueError("equilibrium tolerance must lie between zero and one")

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
                "within_hgt_probability": within_hgt_probability,
                "ancestral_size": population["ancestral_size"],
                "clade_size": population["clade_size"],
                "terminal_size": population["terminal_size"],
                "deep_split_tick": deep_split_tick,
                "terminal_split_tick": terminal_split_tick,
                "a_position": a_position,
                "b_position": b_position,
            }
        )
        for cross_hgt_probability in cross_hgt_probabilities:
            cross_label = probability_label(cross_hgt_probability)
            for mode in modes:
                case_id = f"{rep}__cross_{cross_label}__mode_{mode}"
                cases.append(
                    {
                        "case_id": case_id,
                        "replicate": replicate,
                        "checkpoint_id": checkpoint_id,
                        "checkpoint_dir": relative(checkpoint_dir),
                        "out_dir": relative(
                            run_root / rep / f"cross_{cross_label}" / f"mode_{mode}"
                        ),
                        "reference": reference,
                        "mode": mode,
                        "regime": MODE_LABEL[mode],
                        "cross_hgt_probability": cross_hgt_probability,
                        "cross_hgt_label": cross_label,
                        "seed": deterministic_seed(
                            master_seed,
                            "continuation",
                            replicate,
                            cross_hgt_probability,
                            mode,
                        ),
                        "genome_length": genome_length,
                        "mutation_rate": genome["mutation_rate"],
                        "tract_length": genome["mean_hgt_tract"],
                        "within_hgt_probability": within_hgt_probability,
                        "terminal_size": population["terminal_size"],
                        "sample_per_terminal": population["sample_per_terminal"],
                        "end_tick": end_tick,
                        "fds_strength": strength,
                        "fds_epsilon": epsilon,
                        "equilibrium_monitor_every": monitor_every,
                        "equilibrium_minimum_ticks": minimum_ticks,
                        "equilibrium_stable_checks": stable_checks,
                        "equilibrium_tolerance": tolerance,
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
