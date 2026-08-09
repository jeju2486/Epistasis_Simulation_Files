#!/usr/bin/env python3
"""Run KOVAR 0.8.1 for one completed simulation case."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

from kovar_inputs import materialize_pairs, read_fasta


def require(path: Path, nonempty: bool = True) -> None:
    if not path.exists() or (nonempty and path.is_file() and path.stat().st_size == 0):
        raise FileNotFoundError(f"missing required input: {path}")


def run_checked(command: list[str], cwd: Path, log: Path, env: dict[str, str]) -> None:
    with log.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(
            command, cwd=cwd, env=env, stdout=handle, stderr=subprocess.STDOUT,
            text=True, check=False,
        )
    if completed.returncode:
        raise RuntimeError(f"command exited {completed.returncode}; see {log}")


def run_case(args: argparse.Namespace) -> None:
    case_dir = args.case_dir.resolve()
    spydrpick = case_dir / "spydrpick_all_pairs"
    require(case_dir / "_SUCCESS", nonempty=False)
    require(spydrpick / "_SUCCESS", nonempty=False)
    for path in (
        case_dir / "all_snps.fa", case_dir / "all_snps.positions.tsv",
        case_dir / "core_snps.fa", spydrpick / "spydrpick.edges.gz",
        spydrpick / "all_snps.binary_ac.fa",
    ):
        require(path)

    version = subprocess.run(
        [args.kovar, "--version"], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False,
    )
    if version.returncode or version.stdout.strip() != "KO-Variation 0.8.1":
        raise RuntimeError(f"expected KO-Variation 0.8.1, got: {version.stdout.strip()}")

    analysis_label = "kovar_v081" if args.tree_mode == "core" else f"kovar_v081_{args.tree_mode}"
    if args.max_pairs:
        analysis_label += f"_diagnostic_top_{args.max_pairs}"
    output = case_dir / analysis_label
    if (output / "_SUCCESS").exists() and not args.force:
        print(f"[skip] {output}")
        return
    if output.exists():
        stale = output.with_name(f"{output.name}.incomplete.{time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())}")
        output.replace(stale)
        print(f"[archive] {stale}")

    temporary = Path(tempfile.mkdtemp(prefix=output.name + ".tmp.", dir=case_dir))
    env = os.environ.copy()
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env.setdefault(variable, "1")
    try:
        binary = temporary / "all_snps.binary_ac.fa"
        pairs = temporary / "spydrpick_pairs.0based.tsv"
        shutil.copy2(spydrpick / "all_snps.binary_ac.fa", binary)
        binary_records = read_fasta(binary)
        samples = len(binary_records)
        loci = len(binary_records[0][1])
        pair_count = materialize_pairs(spydrpick / "spydrpick.edges.gz", pairs, args.max_pairs)

        tree: Path | None = None
        tree_command: list[str] | None = None
        if args.tree_mode == "core":
            prefix = temporary / "core_tree"
            tree_command = [
                args.iqtree, "-s", str(case_dir / "core_snps.fa"),
                "-pre", str(prefix), "-m", "GTR+ASC", "-seed", str(args.seed),
                "-nt", str(args.tree_threads), "-redo",
            ]
            run_checked(tree_command, temporary, temporary / "iqtree.log", env)
            tree = prefix.with_suffix(".treefile")
            require(tree)
        elif args.tree_mode == "oracle":
            tree = case_dir / "oracle_local_tree.nwk"
            require(tree)

        results = temporary / "results"
        command = [
            args.kovar, "--fasta", str(binary), "--pairs", str(pairs),
            "--out", str(results), "--direction-mode", "both",
            "--min-maf", str(args.min_maf), "--min-cell-count", str(args.min_cell_count),
            "--spa-mode", args.spa_mode, "--full-refit-p", str(args.full_refit_p),
            "--threads", str(args.threads), "--worker-chunk-size", str(args.worker_chunk_size),
            "--predictor-batch-size", str(args.predictor_batch_size), "--overwrite",
        ]
        if tree is not None:
            command.extend(["--tree", str(tree), "--tree-missing-length", "error"])
        run_checked(command, temporary, temporary / "kovar.log", env)
        for filename in ("ko_variation.tsv", "response_models.tsv", "run_summary.txt"):
            require(results / filename)

        # The exhaustive plain pair file can be reconstructed from the compressed
        # SpydrPick result and otherwise doubles pair-storage requirements.
        pairs.unlink()
        metadata = {
            "version": version.stdout.strip(), "case_dir": str(case_dir),
            "samples": samples, "loci": loci, "pairs": pair_count,
            "candidate_universe": "all_spydrpick_pairs" if args.max_pairs == 0 else "diagnostic_mi_prefix",
            "max_pairs": args.max_pairs, "tree_mode": args.tree_mode,
            "tree_command": tree_command, "kovar_command": command,
        }
        (temporary / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
        (temporary / "_SUCCESS").write_text("complete\n", encoding="utf-8")
        temporary.replace(output)
    except Exception:
        failed = output.with_name(f"{output.name}.failed.{time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())}")
        temporary.replace(failed)
        raise
    print(f"[done] {output}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", required=True, type=Path)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--kovar", default=os.environ.get("KOVAR_BIN", "ko-variation"))
    parser.add_argument("--iqtree", default=os.environ.get("IQTREE_BIN", "iqtree2"))
    parser.add_argument("--tree-mode", choices=("core", "oracle", "grm"), default="core")
    parser.add_argument("--tree-threads", type=int, default=1)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--min-maf", type=float, default=0.05)
    parser.add_argument("--min-cell-count", type=int, default=5)
    parser.add_argument("--spa-mode", choices=("off", "auto", "always"), default="off")
    parser.add_argument("--full-refit-p", type=float, default=0.0)
    parser.add_argument("--worker-chunk-size", type=int, default=1)
    parser.add_argument("--predictor-batch-size", type=int, default=256)
    parser.add_argument("--max-pairs", type=int, default=0)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if min(args.threads, args.tree_threads, args.worker_chunk_size, args.predictor_batch_size) < 1:
        parser.error("thread, chunk and batch values must be positive")
    if args.max_pairs < 0:
        parser.error("--max-pairs must be non-negative")
    run_case(args)


if __name__ == "__main__":
    main()
