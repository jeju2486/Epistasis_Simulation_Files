#!/usr/bin/env python3
"""Run KOVAR 0.8.3 across completed cases with bounded parallelism."""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import subprocess

from simflow import REPO_ROOT, read_tsv, repo_path


def run_one(row: dict[str, str], args: argparse.Namespace) -> tuple[str, int, str]:
    command = [
        os.environ.get("PYTHON_BIN", "python3"), "scripts/run_kovar_case.py",
        "--case-dir", row["out_dir"], "--threads", str(args.threads_per_case),
        "--min-maf", str(args.min_maf),
        "--min-cell-count", str(args.min_cell_count),
        "--spa-mode", args.spa_mode,
    ]
    if args.force:
        command.append("--force")
    completed = subprocess.run(
        command, cwd=REPO_ROOT, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False,
    )
    return row["case_id"], completed.returncode, completed.stdout


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="manifests/cases.tsv")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--threads-per-case", type=int, default=1)
    parser.add_argument("--min-maf", type=float, default=0.05)
    parser.add_argument("--min-cell-count", type=int, default=0)
    parser.add_argument("--spa-mode", choices=("off", "auto", "always"), default="auto")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if min(args.jobs, args.threads_per_case) < 1:
        parser.error("jobs and thread counts must be positive")

    rows = read_tsv(repo_path(args.manifest))
    failures = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(run_one, row, args) for row in rows]
        for future in concurrent.futures.as_completed(futures):
            label, returncode, output = future.result()
            print(f"[{'done' if returncode == 0 else 'FAILED'}] {label}")
            if output.strip():
                print(output.rstrip())
            failures += returncode != 0
    if failures:
        raise SystemExit(f"{failures} KOVAR case(s) failed")


if __name__ == "__main__":
    main()
