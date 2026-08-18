#!/usr/bin/env python3
"""Run threshold-free, ARACNE-free SpydrPick for one completed case."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from kovar_inputs import write_maf_filtered_binary_fasta


TRUTH_PAIRS = (("A", "B"), ("C", "D"))


def read_positions(path: Path) -> tuple[list[int], dict[int, int]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    columns = [int(row["alignment_column"]) for row in rows]
    if columns != list(range(len(rows))):
        raise ValueError("alignment columns must be contiguous and zero-based")
    positions = [int(row["vcf_position"]) - 1 for row in rows]
    if len(set(positions)) != len(positions):
        raise ValueError("all_snps.positions.tsv contains duplicate genomic positions")
    return positions, {position: column for column, position in enumerate(positions)}


def read_truth(path: Path) -> tuple[dict[str, int], set[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    truth = {row["label"]: int(row["position"]) for row in rows}
    if set(truth) != {"A", "B", "C", "D"}:
        raise ValueError("selected_loci.tsv must contain exactly A, B, C and D")
    seeded = {
        row["label"] for row in rows
        if row.get("seeded", "true").strip().lower() == "true"
    }
    return truth, seeded


def parse_edge(line: str) -> tuple[int, int, int, int, float]:
    fields = line.split()
    if len(fields) < 5:
        raise ValueError(f"invalid SpydrPick edge: {line.rstrip()}")
    return int(fields[0]), int(fields[1]), int(fields[2]), int(fields[3]), float(fields[4])


def compress_and_summarize(
    raw_edges: Path,
    compressed_edges: Path,
    truth_output: Path,
    positions: list[int],
    truth: dict[str, int],
    n_loci: int,
    all_positions: set[int] | None = None,
    seeded_labels: set[str] | None = None,
) -> dict[str, object]:
    if seeded_labels is None:
        seeded_labels = set(truth)
    position_to_column = {position: column for column, position in enumerate(positions)}
    targets = {
        tuple(sorted((truth[left], truth[right]))): left + right
        for left, right in TRUTH_PAIRS
        if left in seeded_labels and right in seeded_labels
        and truth[left] in position_to_column and truth[right] in position_to_column
    }
    found: dict[str, dict[str, object]] = {}
    edge_count = 0
    previous_mi = math.inf
    tie_start = 0
    tie_value: float | None = None

    with raw_edges.open(encoding="utf-8") as source, gzip.open(
        compressed_edges, "wt", encoding="utf-8", newline=""
    ) as destination:
        for line in source:
            if not line.strip():
                continue
            raw1, raw2, _raw_distance, aracne, mi = parse_edge(line)
            # Bioconda SpydrPick 1.2.0 lacks --mappings-list. Its default output
            # is one-based alignment-column indices, which we map back here.
            col1, col2 = raw1 - 1, raw2 - 1
            if not (0 <= col1 < n_loci and 0 <= col2 < n_loci):
                raise ValueError(f"edge column outside the alignment: {raw1}, {raw2}")
            pos1, pos2 = positions[col1], positions[col2]
            distance = abs(pos2 - pos1)
            if mi > previous_mi + 1e-12:
                raise ValueError("SpydrPick edges are not sorted by descending MI")
            edge_count += 1
            if tie_value is None or mi != tie_value:
                tie_value = mi
                tie_start = edge_count
            previous_mi = mi
            # Normalized schema remains SpydrPick's five fields, but positions
            # are zero-based alignment columns and distance is physical bp.
            destination.write(f"{col1} {col2} {distance} {aracne} {mi:.17g}\n")
            key = tuple(sorted((pos1, pos2)))
            if key in targets:
                label = targets[key]
                found[label] = {
                    "pair": label,
                    "u_column": col1,
                    "v_column": col2,
                    "u_position": pos1,
                    "v_position": pos2,
                    "physical_distance": abs(pos2 - pos1),
                    "mi": mi,
                    "rank_min": tie_start,
                }

    expected = n_loci * (n_loci - 1) // 2
    if edge_count != expected:
        raise ValueError(
            f"SpydrPick returned {edge_count} edges, but all {expected} pairs were expected; "
            "the run is not an all-pair result"
        )
    with truth_output.open("w", encoding="utf-8", newline="") as handle:
        fields = [
            "pair", "candidate_status", "u_column", "v_column", "u_position", "v_position",
            "physical_distance", "mi", "rank_min", "total_pairs", "rank_fraction",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for left, right in TRUTH_PAIRS:
            label = left + right
            if label in found:
                row = found[label]
                row["candidate_status"] = "eligible"
                row["total_pairs"] = edge_count
                row["rank_fraction"] = float(row["rank_min"]) / edge_count
            else:
                if left not in seeded_labels or right not in seeded_labels:
                    status = "not_seeded"
                else:
                    both_sampled = (
                        all_positions is not None
                        and truth[left] in all_positions
                        and truth[right] in all_positions
                    )
                    status = "maf_filtered" if both_sampled else "locus_absent"
                if (
                    left in seeded_labels and right in seeded_labels
                    and truth[left] in position_to_column and truth[right] in position_to_column
                ):
                    raise ValueError(f"eligible truth pair absent from all-pair output: {label}")
                row = {
                    "pair": label, "candidate_status": status,
                    "u_column": position_to_column.get(truth[left], ""),
                    "v_column": position_to_column.get(truth[right], ""),
                    "u_position": truth[left], "v_position": truth[right],
                    "physical_distance": abs(truth[right] - truth[left]),
                    "mi": "", "rank_min": "", "total_pairs": edge_count,
                    "rank_fraction": "",
                }
            writer.writerow(row)
    return {"n_loci": n_loci, "n_pairs": edge_count}


def run_case(case_dir: Path, executable: str, threads: int, force: bool, unweighted: bool, min_maf: float) -> None:
    case_dir = case_dir.resolve()
    required = [
        case_dir / "_SUCCESS", case_dir / "all_snps.fa",
        case_dir / "all_snps.positions.tsv", case_dir / "selected_loci.tsv",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("missing case input(s): " + ", ".join(missing))

    output = case_dir / ("spydrpick_all_pairs_unweighted" if unweighted else "spydrpick_all_pairs")
    if (output / "_SUCCESS").exists() and not force:
        print(f"[skip] {output}")
        return
    if output.exists():
        shutil.rmtree(output)

    all_positions, _original_position_to_column = read_positions(case_dir / "all_snps.positions.tsv")
    truth, seeded_labels = read_truth(case_dir / "selected_loci.tsv")
    absent = sorted({truth[label] for label in seeded_labels} - set(all_positions))
    if absent:
        print(f"[warning] truth positions absent from sampled SNP alignment: {absent}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=output.name + ".tmp.", dir=output.parent))
    try:
        binary_alignment = temporary / "all_snps.binary_ac.fa"
        eligible_loci = temporary / "eligible_loci.tsv"
        binary_samples, input_loci, binary_loci = write_maf_filtered_binary_fasta(
            case_dir / "all_snps.fa", case_dir / "all_snps.positions.tsv",
            binary_alignment, eligible_loci, min_maf,
        )
        with eligible_loci.open(encoding="utf-8", newline="") as handle:
            positions = [int(row["slim_position"]) for row in csv.DictReader(handle, delimiter="\t")]
        expected_pairs = binary_loci * (binary_loci - 1) // 2
        print(f"[input] {case_dir}: {input_loci} input loci, {binary_loci} pass MAF >= {min_maf}, {expected_pairs} pairs")
        command = [
            executable, "--verbose", f"--threads={threads}", "--no-aracne",
            "--mi-threshold=0", "--no-filter-alignment",
            str(binary_alignment),
        ]
        if unweighted:
            command.insert(-1, "--no-sample-reweighting")
        with (temporary / "spydrpick.log").open("w", encoding="utf-8") as log:
            completed = subprocess.run(
                command, cwd=temporary, stdout=log, stderr=subprocess.STDOUT,
                text=True, check=False,
            )
        if completed.returncode:
            raise RuntimeError(f"SpydrPick exited {completed.returncode}; see {temporary / 'spydrpick.log'}")
        edges = list(temporary.glob("*.spydrpick_couplings.*-based.*edges"))
        if len(edges) != 1:
            raise RuntimeError(f"expected one edges file, found {len(edges)}")
        metadata = compress_and_summarize(
            edges[0], temporary / "spydrpick.edges.gz", temporary / "truth_pairs.tsv",
            positions, truth, len(positions), set(all_positions), seeded_labels,
        )
        edges[0].unlink()
        metadata.update({
            "case_dir": str(case_dir), "command": command,
            "sample_reweighting": "off" if unweighted else "SpydrPick default",
            "aracne": False, "mi_threshold": 0,
            "min_maf": min_maf, "input_loci": input_loci,
            "eligible_loci": binary_loci,
            "locus_representation": "binary_reference_vs_any_nonreference_A_C",
            "binary_samples": binary_samples,
            "edge_schema": "zero-based alignment columns; physical distance in bp; ARACNE flag; MI",
        })
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
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--unweighted", action="store_true")
    parser.add_argument("--min-maf", type=float, default=0.05)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be at least one")
    if not 0.0 <= args.min_maf <= 0.5:
        parser.error("--min-maf must be between zero and 0.5")
    run_case(args.case_dir, args.spydrpick, args.threads, args.force, args.unweighted, args.min_maf)


if __name__ == "__main__":
    main()
