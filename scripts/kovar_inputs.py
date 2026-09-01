#!/usr/bin/env python3
"""Prepare compact binary FASTA and plain candidate pairs for KOVAR."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks).upper()))
                name = line[1:].split()[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError("FASTA sequence encountered before its header")
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    if not records:
        raise ValueError(f"no FASTA records in {path}")
    return records


def reference_states(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    columns = [int(row["alignment_column"]) for row in rows]
    if columns != list(range(len(rows))):
        raise ValueError("alignment columns must be contiguous and zero-based")
    states = [row["alleles"].split(",", 1)[0].upper() for row in rows]
    if any(len(state) != 1 or state not in "ACGT" for state in states):
        raise ValueError("position map contains a non-nucleotide reference state")
    return states


def write_binary_fasta(alignment: Path, positions: Path, output: Path) -> tuple[int, int]:
    records = read_fasta(alignment)
    references = reference_states(positions)
    expected_length = len(references)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        for name, sequence in records:
            if len(sequence) != expected_length:
                raise ValueError(f"alignment length mismatch for {name}")
            invalid = set(sequence) - set("ACGT")
            if invalid:
                raise ValueError(f"unsupported states for {name}: {sorted(invalid)}")
            binary = "".join("A" if observed == reference else "C" for observed, reference in zip(sequence, references, strict=True))
            handle.write(f">{name}\n")
            for start in range(0, len(binary), 80):
                handle.write(binary[start : start + 80] + "\n")
    return len(records), expected_length


def write_maf_filtered_binary_fasta(
    alignment: Path,
    positions: Path,
    output: Path,
    mapping_output: Path,
    min_maf: float = 0.05,
) -> tuple[int, int, int]:
    """Write a binary matrix after a predeclared marginal-frequency filter."""
    if not 0.0 <= min_maf <= 0.5:
        raise ValueError("min_maf must be between zero and 0.5")
    records = read_fasta(alignment)
    references = reference_states(positions)
    expected_length = len(references)
    binary_records: list[tuple[str, str]] = []
    presence_counts = [0] * expected_length
    for name, sequence in records:
        if len(sequence) != expected_length:
            raise ValueError(f"alignment length mismatch for {name}")
        invalid = set(sequence) - set("ACGT")
        if invalid:
            raise ValueError(f"unsupported states for {name}: {sorted(invalid)}")
        binary = "".join(
            "A" if observed == reference else "C"
            for observed, reference in zip(sequence, references, strict=True)
        )
        for column, state in enumerate(binary):
            presence_counts[column] += state == "C"
        binary_records.append((name, binary))

    n_samples = len(binary_records)
    prevalence = [count / n_samples for count in presence_counts]
    maf = [min(value, 1.0 - value) for value in prevalence]
    retained = [column for column, value in enumerate(maf) if value >= min_maf]
    if len(retained) < 2:
        raise ValueError(f"fewer than two loci pass min_maf={min_maf}")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        for name, binary in binary_records:
            filtered = "".join(binary[column] for column in retained)
            handle.write(f">{name}\n")
            for start in range(0, len(filtered), 80):
                handle.write(filtered[start : start + 80] + "\n")

    with positions.open(encoding="utf-8", newline="") as handle:
        position_rows = list(csv.DictReader(handle, delimiter="\t"))
    fields = ["filtered_column", "original_column", "slim_position", "prevalence", "maf"]
    with mapping_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for filtered_column, original_column in enumerate(retained):
            writer.writerow({
                "filtered_column": filtered_column,
                "original_column": original_column,
                "slim_position": int(position_rows[original_column]["vcf_position"]) - 1,
                "prevalence": prevalence[original_column],
                "maf": maf[original_column],
            })
    return n_samples, expected_length, len(retained)


def materialize_pairs(source: Path, output: Path) -> int:
    """Write all unique canonical pairs using only KOVAR's u/v contract."""
    pairs: set[tuple[int, int]] = set()
    previous_mi = float("inf")
    with gzip.open(source, "rt", encoding="utf-8") as reader:
        for line in reader:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) < 5:
                raise ValueError(f"invalid SpydrPick pair row: {line.rstrip()}")
            raw_u, raw_v, mi = int(fields[0]), int(fields[1]), float(fields[4])
            if raw_u < 0 or raw_v < 0 or raw_u == raw_v:
                raise ValueError(f"invalid KOVAR pair: {raw_u}, {raw_v}")
            if mi > previous_mi + 1e-12:
                raise ValueError("SpydrPick pairs are not sorted by descending MI")
            previous_mi = mi
            pair = tuple(sorted((raw_u, raw_v)))
            if pair in pairs:
                raise ValueError(f"duplicate unordered SpydrPick pair: {pair}")
            pairs.add(pair)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="\n") as writer:
        writer.write("u\tv\n")
        for u, v in sorted(pairs):
            writer.write(f"{u}\t{v}\n")
    if not pairs:
        raise ValueError(f"no pairs found in {source}")
    return len(pairs)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument("--positions", required=True, type=Path)
    parser.add_argument("--binary-output", required=True, type=Path)
    parser.add_argument("--pairs-gz", required=True, type=Path)
    parser.add_argument("--pairs-output", required=True, type=Path)
    args = parser.parse_args()
    samples, loci = write_binary_fasta(args.alignment, args.positions, args.binary_output)
    pairs = materialize_pairs(args.pairs_gz, args.pairs_output)
    print(f"[done] KOVAR inputs: samples={samples} loci={loci} pairs={pairs}")


if __name__ == "__main__":
    main()
