#!/usr/bin/env python3
"""Validate the two uniformly seeded A-B focal-locus records."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def validate(
    path: Path,
    expected_positions: dict[str, int] | None = None,
    mode: int | None = None,
) -> None:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader((line for line in handle if line.strip()), delimiter="\t"))
    counts = Counter(row.get("label", "") for row in rows)
    if counts != Counter({"A": 1, "B": 1}):
        raise ValueError("selected_loci.tsv must contain exactly one A and one B record")
    if mode is not None and mode not in {1, 2}:
        raise ValueError("mode must be 1 or 2")
    mutation_ids: set[str] = set()
    for row in rows:
        label = row["label"]
        mutation_id = row.get("mutation_id", "").strip()
        if not mutation_id or mutation_id == "NA":
            raise ValueError(f"locus {label} lacks a mutation ID")
        mutation_ids.add(mutation_id)
        if mode is not None and int(row.get("mode", mode)) != mode:
            raise ValueError(f"locus {label} has mode inconsistent with mode {mode}")
        frequency = float(row.get("seeded_frequency", "nan"))
        if not 0.49 <= frequency <= 0.51:
            raise ValueError(f"locus {label} has unbalanced seeded_frequency={frequency}")
        if expected_positions is not None and int(row["position"]) != expected_positions[label]:
            raise ValueError(
                f"locus {label} is at {row['position']}, expected {expected_positions[label]}"
            )
    if len(mutation_ids) != 2:
        raise ValueError("A and B must have distinct mutation IDs")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=Path)
    parser.add_argument("--a-position", type=int)
    parser.add_argument("--b-position", type=int)
    parser.add_argument("--mode", type=int, choices=(1, 2))
    args = parser.parse_args()
    if (args.a_position is None) != (args.b_position is None):
        parser.error("supply both A and B positions or neither")
    expected = None
    if args.a_position is not None:
        expected = {"A": args.a_position, "B": args.b_position}
    try:
        validate(args.path, expected, args.mode)
    except (OSError, ValueError) as error:
        raise SystemExit(str(error)) from error


if __name__ == "__main__":
    main()
