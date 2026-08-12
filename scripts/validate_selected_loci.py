#!/usr/bin/env python3
"""Validate the four labelled standing variants selected by SLiM."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


EXPECTED_LABELS = {"A", "B", "C", "D"}


def validate(path: Path, expected_positions: dict[str, int] | None = None) -> None:
    with open(path, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader((line for line in handle if line.strip()), delimiter="\t"))

    counts = Counter(row.get("label", "") for row in rows)
    if set(counts) != EXPECTED_LABELS or any(counts[label] != 1 for label in EXPECTED_LABELS):
        raise ValueError("selected_loci.tsv must contain exactly one record for each of A, B, C, and D")

    mutation_ids: set[str] = set()
    for row in rows:
        if not row.get("mutation_id", "").strip() or not row.get("position", "").strip():
            raise ValueError(f"locus {row.get('label', '?')} lacks a mutation_id or position")
        mutation_ids.add(row["mutation_id"])
        if "seeding_design" in row and row["seeding_design"] != "balanced_16_haplotype_cycle":
            raise ValueError(f"locus {row['label']} has an unexpected seeding design")
        for column in (
            "global_frequency", "p1_frequency", "p2_frequency",
            "p3_frequency", "p4_frequency",
        ):
            if column in row and row[column]:
                frequency = float(row[column])
                if not 0.45 <= frequency <= 0.55:
                    raise ValueError(
                        f"locus {row['label']} has unbalanced initial {column}={frequency}"
                    )
        if expected_positions is not None:
            label = row["label"]
            if int(row["position"]) != expected_positions[label]:
                raise ValueError(
                    f"locus {label} is at {row['position']}, expected {expected_positions[label]}"
                )
    if len(mutation_ids) != 4:
        raise ValueError("A-D must have four distinct mutation IDs")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=Path)
    parser.add_argument("--a-position", type=int)
    parser.add_argument("--b-position", type=int)
    parser.add_argument("--c-position", type=int)
    parser.add_argument("--d-position", type=int)
    args = parser.parse_args()
    supplied = [args.a_position, args.b_position, args.c_position, args.d_position]
    if any(value is not None for value in supplied) and any(value is None for value in supplied):
        parser.error("supply all four expected positions or none")
    expected = None
    if all(value is not None for value in supplied):
        expected = dict(zip(("A", "B", "C", "D"), supplied, strict=True))
    try:
        validate(args.path, expected)
    except (OSError, ValueError) as error:
        raise SystemExit(str(error)) from error


if __name__ == "__main__":
    main()
