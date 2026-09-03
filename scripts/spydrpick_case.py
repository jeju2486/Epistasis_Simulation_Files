#!/usr/bin/env python3
"""Run one threshold-free, ARACNE-free all-pair SpydrPick scan."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from kovar_inputs import write_maf_filtered_binary_fasta


OUTPUT_NAMES = {
    "default": "spydrpick_all_pairs",
    "none": "spydrpick_all_pairs_unweighted",
}


def output_name(sample_reweighting: str) -> str:
    try:
        return OUTPUT_NAMES[sample_reweighting]
    except KeyError as error:
        raise ValueError(f"unknown sample-reweighting mode: {sample_reweighting}") from error


def parse_edge(line: str) -> tuple[int, int, int, int, float]:
    fields = line.split()
    if len(fields) < 5:
        raise ValueError(f"invalid SpydrPick edge: {line.rstrip()}")
    return int(fields[0]), int(fields[1]), int(fields[2]), int(fields[3]), float(fields[4])


def pair_index(u: int, v: int, n_loci: int) -> int:
    u, v = sorted((u, v))
    if not (0 <= u < v < n_loci):
        raise ValueError(f"invalid unordered pair: {u}, {v}")
    return u * (2 * n_loci - u - 1) // 2 + (v - u - 1)


def normalize_edges(raw: Path, output: Path, positions: list[int]) -> tuple[int, int]:
    n_loci = len(positions)
    expected = n_loci * (n_loci - 1) // 2
    seen = bytearray(expected)
    count = 0
    previous_mi = float("inf")
    with raw.open(encoding="utf-8") as source, gzip.open(
        output, "wt", encoding="utf-8", newline=""
    ) as destination:
        for line in source:
            if not line.strip():
                continue
            raw_u, raw_v, _distance, aracne, mi = parse_edge(line)
            u, v = raw_u - 1, raw_v - 1
            if not (0 <= u < len(positions) and 0 <= v < len(positions) and u != v):
                raise ValueError(f"edge column outside alignment: {raw_u}, {raw_v}")
            if mi > previous_mi + 1e-12:
                raise ValueError("SpydrPick edges are not sorted by descending MI")
            if mi < -1e-12:
                raise ValueError(f"SpydrPick returned negative MI: {mi}")
            previous_mi = mi
            u, v = sorted((u, v))
            index = pair_index(u, v, n_loci)
            if seen[index]:
                raise ValueError(f"duplicate SpydrPick pair: {u}, {v}")
            seen[index] = 1
            distance = abs(positions[v] - positions[u])
            destination.write(f"{u} {v} {distance} {aracne} {mi:.17g}\n")
            count += 1
        # SpydrPick applies a strict MI > threshold comparison. With threshold
        # zero, exact-zero pairs are therefore absent rather than printed.
        for u in range(n_loci):
            for v in range(u + 1, n_loci):
                if not seen[pair_index(u, v, n_loci)]:
                    distance = abs(positions[v] - positions[u])
                    destination.write(f"{u} {v} {distance} 0 0\n")
    return expected, expected - count


def run_case(
    case_dir: Path,
    executable: str,
    threads: int,
    force: bool,
    min_maf: float,
    sample_reweighting: str,
) -> None:
    case_dir = case_dir.resolve()
    for path in (
        case_dir / "_SUCCESS", case_dir / "all_snps.fa",
        case_dir / "all_snps.positions.tsv",
    ):
        if not path.exists():
            raise FileNotFoundError(f"missing case input: {path}")
    output = case_dir / output_name(sample_reweighting)
    if (output / "_SUCCESS").exists() and not force:
        print(f"[skip] {output}")
        return
    if output.exists():
        shutil.rmtree(output)
    temporary = Path(tempfile.mkdtemp(prefix=output.name + ".tmp.", dir=case_dir))
    try:
        binary = temporary / "all_snps.binary_ac.fa"
        mapping = temporary / "eligible_loci.tsv"
        samples, input_loci, eligible_loci = write_maf_filtered_binary_fasta(
            case_dir / "all_snps.fa", case_dir / "all_snps.positions.tsv",
            binary, mapping, min_maf,
        )
        with mapping.open(encoding="utf-8", newline="") as handle:
            positions = [int(row["slim_position"]) for row in csv.DictReader(handle, delimiter="\t")]
        command = [
            executable, "--verbose", f"--threads={threads}", "--no-aracne",
            "--mi-threshold=0", "--no-filter-alignment",
        ]
        if sample_reweighting == "none":
            command.append("--no-sample-reweighting")
        command.append(str(binary))
        with (temporary / "spydrpick.log").open("w", encoding="utf-8") as log:
            completed = subprocess.run(command, cwd=temporary, stdout=log,
                                       stderr=subprocess.STDOUT, text=True, check=False)
        if completed.returncode:
            raise RuntimeError(f"SpydrPick exited {completed.returncode}")
        edges = list(temporary.glob("*.spydrpick_couplings.*-based.*edges"))
        if len(edges) != 1:
            raise RuntimeError(f"expected one SpydrPick edge file, found {len(edges)}")
        pair_count, restored_zero_pairs = normalize_edges(
            edges[0], temporary / "spydrpick.edges.gz", positions
        )
        edges[0].unlink()
        metadata = {
            "case_dir": str(case_dir), "command": command, "samples": samples,
            "input_loci": input_loci, "eligible_loci": eligible_loci,
            "pairs": pair_count, "min_maf": min_maf,
            "spydrpick_reported_pairs": pair_count - restored_zero_pairs,
            "zero_mi_pairs_restored": restored_zero_pairs,
            "sample_reweighting": sample_reweighting, "aracne": False,
        }
        (temporary / "run_metadata.json").write_text(
            json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
        )
        (temporary / "_SUCCESS").write_text("complete\n", encoding="utf-8")
        temporary.replace(output)
    except Exception:
        failed = output.with_name(output.name + ".failed")
        if failed.exists():
            shutil.rmtree(failed)
        temporary.replace(failed)
        raise
    print(f"[done] {output}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", required=True, type=Path)
    parser.add_argument("--spydrpick", default=os.environ.get("SPYDRPICK_BIN", "SpydrPick"))
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--min-maf", type=float, default=0.05)
    parser.add_argument(
        "--sample-reweighting", choices=tuple(OUTPUT_NAMES), default="default"
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    if not 0.0 <= args.min_maf <= 0.5:
        parser.error("--min-maf must lie between zero and 0.5")
    run_case(
        args.case_dir, args.spydrpick, args.threads, args.force,
        args.min_maf, args.sample_reweighting,
    )


if __name__ == "__main__":
    main()
