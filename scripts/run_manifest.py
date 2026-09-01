#!/usr/bin/env python3
"""Run checkpoint or continuation manifest rows with bounded parallelism."""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import subprocess

from simflow import REPO_ROOT, read_tsv


def checkpoint_command(row: dict[str, str]) -> list[str]:
    return [
        "bash", "scripts/run_checkpoint.sh", "--outdir", row["out_dir"],
        "--reference", row["reference"], "--seed", row["seed"],
        "--genome-length", row["genome_length"], "--mu", row["mutation_rate"],
        "--tract-length", row["tract_length"], "--within-hgt", row["within_hgt_probability"],
        "--ancestral-size", row["ancestral_size"], "--clade-size", row["clade_size"],
        "--terminal-size", row["terminal_size"], "--deep-split-tick", row["deep_split_tick"],
        "--terminal-split-tick", row["terminal_split_tick"],
        "--a-position", row["a_position"], "--b-position", row["b_position"],
    ]


def case_command(row: dict[str, str]) -> list[str]:
    return [
        "bash", "scripts/run_case.sh", "--outdir", row["out_dir"],
        "--checkpoint-dir", row["checkpoint_dir"], "--reference", row["reference"],
        "--seed", row["seed"], "--mode", row["mode"],
        "--genome-length", row["genome_length"], "--mu", row["mutation_rate"],
        "--tract-length", row["tract_length"], "--within-hgt", row["within_hgt_probability"],
        "--cross-hgt", row["cross_hgt_probability"],
        "--terminal-size", row["terminal_size"],
        "--sample-per-terminal", row["sample_per_terminal"],
        "--end-tick", row["end_tick"], "--fds-strength", row["fds_strength"],
        "--fds-epsilon", row["fds_epsilon"], "--ancestral-ne", row["ancestral_ne"],
        "--equilibrium-monitor-every", row["equilibrium_monitor_every"],
        "--equilibrium-minimum-ticks", row["equilibrium_minimum_ticks"],
        "--equilibrium-stable-checks", row["equilibrium_stable_checks"],
        "--equilibrium-tolerance", row["equilibrium_tolerance"],
        "--tree-position", row["tree_position"],
        "--a-position", row["a_position"], "--b-position", row["b_position"],
    ]


def run_one(label: str, command: list[str], env: dict[str, str]) -> tuple[str, int, str]:
    completed = subprocess.run(command, cwd=REPO_ROOT, env=env, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    return label, completed.returncode, completed.stdout


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--stage", required=True, choices=("checkpoint", "case"))
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--slim-bin", default=os.environ.get("SLIM_BIN", "slim"))
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
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
            print(f"[{'done' if returncode == 0 else 'FAILED'}] {label}")
            if output.strip():
                print(output.rstrip())
            failures += returncode != 0
    if failures:
        raise SystemExit(f"{failures} {args.stage} job(s) failed")


if __name__ == "__main__":
    main()
