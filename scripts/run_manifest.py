#!/usr/bin/env python3
"""Run manifest rows with bounded local parallelism."""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import subprocess
import sys
from pathlib import Path

from simflow import REPO_ROOT, read_tsv


def checkpoint_command(row: dict[str, str]) -> list[str]:
    return [
        "bash", "scripts/run_checkpoint.sh",
        "--outdir", row["out_dir"],
        "--reference", row["reference"],
        "--seed", row["seed"],
        "--genome-length", row["genome_length"],
        "--mu", row["mutation_rate"],
        "--tract-length", row["tract_length"],
        "--within-hgt", row["within_hgt_probability"],
        "--ancestral-size", row["ancestral_size"],
        "--clade-size", row["clade_size"],
        "--terminal-size", row["terminal_size"],
        "--ancestral-generations", row["ancestral_generations"],
        "--deep-generations", row["deep_generations"],
        "--terminal-generations", row["terminal_generations"],
        "--global-freq-min", row["global_frequency_min"],
        "--global-freq-max", row["global_frequency_max"],
        "--lineage-freq-min", row["lineage_frequency_min"],
        "--lineage-freq-max", row["lineage_frequency_max"],
        "--minimum-distance", row["minimum_distance"],
    ]


def case_command(row: dict[str, str]) -> list[str]:
    return [
        "bash", "scripts/run_case.sh",
        "--outdir", row["out_dir"],
        "--checkpoint-dir", row["checkpoint_dir"],
        "--reference", row["reference"],
        "--seed", row["seed"],
        "--mode", row["mode"],
        "--cross-hgt", row["cross_hgt_probability"],
        "--genome-length", row["genome_length"],
        "--mu", row["mutation_rate"],
        "--tract-length", row["tract_length"],
        "--within-hgt", row["within_hgt_probability"],
        "--terminal-size", row["terminal_size"],
        "--sample-per-terminal", row["sample_per_terminal"],
        "--experiment-generations", row["experiment_generations"],
        "--s-ab", row["s_ab"],
        "--s-cd", row["s_cd"],
        "--monitor-every", row["monitor_every"],
        "--ancestral-ne", row["ancestral_ne"],
        "--oracle-tree-position", row["oracle_tree_position"],
    ]


def run_one(label: str, command: list[str], env: dict[str, str]) -> tuple[str, int, str]:
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    return label, completed.returncode, completed.stdout


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--stage", required=True, choices=("checkpoint", "case"))
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--slim-bin", default=os.environ.get("SLIM_BIN", "slim"))
    args = parser.parse_args()

    if args.jobs < 1:
        parser.error("--jobs must be at least one")

    rows = read_tsv(args.manifest)
    builder = checkpoint_command if args.stage == "checkpoint" else case_command
    label_key = "checkpoint_id" if args.stage == "checkpoint" else "case_id"
    env = os.environ.copy()
    env["SLIM_BIN"] = args.slim_bin

    failures = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(run_one, row[label_key], builder(row), env) for row in rows]
        for future in concurrent.futures.as_completed(futures):
            label, returncode, output = future.result()
            status = "done" if returncode == 0 else "FAILED"
            print(f"[{status}] {label}")
            if output.strip():
                print(output.rstrip())
            if returncode != 0:
                failures += 1

    if failures:
        raise SystemExit(f"{failures} {args.stage} job(s) failed")


if __name__ == "__main__":
    main()
