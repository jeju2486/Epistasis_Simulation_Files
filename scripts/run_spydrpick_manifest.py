#!/usr/bin/env python3
"""Run all-pair SpydrPick across completed cases with bounded parallelism."""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import subprocess
from pathlib import Path

from simflow import REPO_ROOT, read_tsv, repo_path


def run_one(row: dict[str, str], args: argparse.Namespace) -> tuple[str, int, str]:
    command = [
        os.environ.get("PYTHON_BIN", "python3"), "scripts/spydrpick_case.py",
        "--case-dir", row["out_dir"], "--spydrpick", args.spydrpick,
        "--threads", str(args.threads_per_case), "--min-maf", str(args.min_maf),
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
    parser.add_argument("--jobs", type=int, default=1, help="concurrent cases")
    parser.add_argument("--threads-per-case", type=int, default=1)
    parser.add_argument("--spydrpick", default=os.environ.get("SPYDRPICK_BIN", "SpydrPick"))
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--min-maf", type=float, default=0.05)
    args = parser.parse_args()
    if args.jobs < 1 or args.threads_per_case < 1:
        parser.error("job and thread counts must be at least one")

    rows = read_tsv(repo_path(args.manifest))
    incomplete = [row["case_id"] for row in rows if not (repo_path(row["out_dir"]) / "_SUCCESS").exists()]
    if incomplete:
        raise SystemExit(
            f"{len(incomplete)} case(s) are incomplete; first missing case: {incomplete[0]}"
        )
    cpus = os.cpu_count() or 1
    requested = args.jobs * args.threads_per_case
    if requested > cpus:
        print(f"[warning] requested {requested} compute threads on {cpus} visible CPUs")

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
        raise SystemExit(f"{failures} SpydrPick case(s) failed")


if __name__ == "__main__":
    main()
