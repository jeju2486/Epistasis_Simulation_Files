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


def materialize_pairs(source: Path, output: Path, max_pairs: int = 0) -> int:
    if max_pairs < 0:
        raise ValueError("max_pairs must be non-negative")
    count = 0
    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(source, "rt", encoding="utf-8") as reader, output.open(
        "w", encoding="utf-8", newline="\n"
    ) as writer:
        for line in reader:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) < 5:
                raise ValueError(f"invalid SpydrPick pair row: {line.rstrip()}")
            u, v = int(fields[0]), int(fields[1])
            if u < 0 or v < 0 or u == v:
                raise ValueError(f"invalid KOVAR pair: {u}, {v}")
            writer.write(" ".join(fields[:5]) + "\n")
            count += 1
            if max_pairs and count >= max_pairs:
                break
    if count == 0:
        raise ValueError(f"no pairs found in {source}")
    return count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument("--positions", required=True, type=Path)
    parser.add_argument("--binary-output", required=True, type=Path)
    parser.add_argument("--pairs-gz", required=True, type=Path)
    parser.add_argument("--pairs-output", required=True, type=Path)
    parser.add_argument("--max-pairs", type=int, default=0)
    args = parser.parse_args()
    samples, loci = write_binary_fasta(args.alignment, args.positions, args.binary_output)
    pairs = materialize_pairs(args.pairs_gz, args.pairs_output, args.max_pairs)
    print(f"[done] KOVAR inputs: samples={samples} loci={loci} pairs={pairs}")


if __name__ == "__main__":
    main()
