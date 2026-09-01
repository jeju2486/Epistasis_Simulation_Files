#!/usr/bin/env python3
"""Run KOVAR 0.8.3 on the complete SpydrPick-eligible pair universe."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

from kovar_inputs import materialize_pairs, read_fasta


RESULT_NAME = "kovar_v083_simulation_tree"


def require(path: Path, nonempty: bool = True) -> None:
    if not path.exists() or (nonempty and path.is_file() and path.stat().st_size == 0):
        raise FileNotFoundError(f"missing required input: {path}")


def archive(path: Path, label: str) -> Path:
    destination = path.with_name(
        f"{path.name}.{label}.{time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())}"
    )
    path.replace(destination)
    return destination


def run_case(args: argparse.Namespace) -> None:
    case_dir = args.case_dir.resolve()
    spydrpick = case_dir / "spydrpick_all_pairs"
    binary = spydrpick / "all_snps.binary_ac.fa"
    compressed_pairs = spydrpick / "spydrpick.edges.gz"
    pairs = spydrpick / "kovar_pairs.tsv"
    tree = case_dir / "simulation_tree.nwk"
    for path in (
        case_dir / "_SUCCESS", spydrpick / "_SUCCESS", binary,
        compressed_pairs, tree,
    ):
        require(path, nonempty=path.name != "_SUCCESS")

    version = subprocess.run(
        [args.kovar, "--version"], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False,
    )
    if version.returncode or version.stdout.strip() != "KO-Variation 0.8.3":
        raise RuntimeError(f"expected KO-Variation 0.8.3, got: {version.stdout.strip()}")

    output = case_dir / RESULT_NAME
    if (output / "_SUCCESS").exists() and not args.force:
        print(f"[skip] {output}")
        return
    if output.exists() and args.force:
        print(f"[archive] {archive(output, 'replaced')}")

    pair_count = materialize_pairs(compressed_pairs, pairs)
    records = read_fasta(binary)
    samples = len(records)
    loci = len(records[0][1])
    expected_pairs = loci * (loci - 1) // 2
    if pair_count != expected_pairs:
        raise ValueError(f"found {pair_count} pairs; expected complete universe of {expected_pairs}")

    resume = output.exists() and (output / ".kovar_checkpoint").exists()
    if output.exists() and not resume:
        print(f"[archive] {archive(output, 'incomplete')}")

    command = [
        args.kovar, "--fasta", str(binary), "--pairs", str(pairs),
        "--tree", str(tree), "--out", str(output),
        "--min-maf", str(args.min_maf),
        "--min-cell-count", str(args.min_cell_count),
        "--spa-mode", args.spa_mode, "--threads", str(args.threads),
        "--no-progress",
    ]
    if resume:
        command.append("--resume")

    env = os.environ.copy()
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env.setdefault(variable, "1")
    log = case_dir / f"{RESULT_NAME}.log"
    with log.open("a" if resume else "w", encoding="utf-8") as handle:
        completed = subprocess.run(
            command, cwd=case_dir, env=env, stdout=handle,
            stderr=subprocess.STDOUT, text=True, check=False,
        )
    if completed.returncode:
        raise RuntimeError(f"KOVAR exited {completed.returncode}; see {log}")
    for filename in (
        "ko_variation.tsv", "response_models.tsv", "run_summary.txt",
        "execution_metadata.tsv",
    ):
        require(output / filename)

    metadata = {
        "version": version.stdout.strip(), "case_dir": str(case_dir),
        "samples": samples, "loci": loci, "pairs": pair_count,
        "candidate_universe": "all_MAF_eligible_unordered_pairs",
        "tree": str(tree), "resumed": resume, "kovar_command": command,
    }
    (output / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    shutil.copy2(log, output / "kovar.log")
    (output / "_SUCCESS").write_text("complete\n", encoding="utf-8")
    print(f"[done] {output}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", required=True, type=Path)
    parser.add_argument("--kovar", default=os.environ.get("KOVAR_BIN", "ko-variation"))
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--min-maf", type=float, default=0.05)
    parser.add_argument("--min-cell-count", type=int, default=5)
    parser.add_argument("--spa-mode", choices=("off", "auto", "always"), default="auto")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    if not 0.0 <= args.min_maf <= 0.5:
        parser.error("--min-maf must lie between zero and 0.5")
    if args.min_cell_count < 0:
        parser.error("--min-cell-count must be non-negative")
    run_case(args)


if __name__ == "__main__":
    main()
