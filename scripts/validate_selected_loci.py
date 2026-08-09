#!/usr/bin/env python3
"""Validate the four labelled standing variants selected by SLiM."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


EXPECTED_LABELS = {"A", "B", "C", "D"}


def validate(path: Path) -> None:
    with open(path, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader((line for line in handle if line.strip()), delimiter="\t"))

    counts = Counter(row.get("label", "") for row in rows)
    if set(counts) != EXPECTED_LABELS or any(counts[label] != 1 for label in EXPECTED_LABELS):
        raise ValueError("selected_loci.tsv must contain exactly one record for each of A, B, C, and D")

    for row in rows:
        if not row.get("mutation_id", "").strip() or not row.get("position", "").strip():
            raise ValueError(f"locus {row.get('label', '?')} lacks a mutation_id or position")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=Path)
    args = parser.parse_args()
    try:
        validate(args.path)
    except (OSError, ValueError) as error:
        raise SystemExit(str(error)) from error


if __name__ == "__main__":
    main()
