#!/usr/bin/env python3
"""Validate mode-specific focal seeding records produced by SLiM."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


EXPECTED_LABELS = {"A", "B", "C", "D"}


MODE_SEEDED_LABELS = {0: {"A", "B", "C", "D"}, 1: {"A", "B"}, 2: {"C", "D"}}
MODE_ACTIVE_PAIR = {0: "AB_and_CD_neutral", 1: "AB", 2: "CD"}


def validate(
    path: Path,
    expected_positions: dict[str, int] | None = None,
    mode: int | None = None,
) -> None:
    with open(path, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader((line for line in handle if line.strip()), delimiter="\t"))

    counts = Counter(row.get("label", "") for row in rows)
    if set(counts) != EXPECTED_LABELS or any(counts[label] != 1 for label in EXPECTED_LABELS):
        raise ValueError("selected_loci.tsv must contain exactly one record for each of A, B, C, and D")

    if mode is not None and mode not in MODE_SEEDED_LABELS:
        raise ValueError("mode must be 0, 1 or 2")
    expected_seeded = MODE_SEEDED_LABELS.get(mode, EXPECTED_LABELS)
    mutation_ids: set[str] = set()
    for row in rows:
        label = row["label"]
        if not row.get("position", "").strip():
            raise ValueError(f"locus {label} lacks a position")
        seeded_text = row.get("seeded", "true").strip().lower()
        if seeded_text not in {"true", "false"}:
            raise ValueError(f"locus {label} has invalid seeded value {seeded_text!r}")
        seeded = seeded_text == "true"
        if mode is not None and seeded != (label in expected_seeded):
            raise ValueError(f"locus {label} has seeded={seeded_text}, inconsistent with mode {mode}")
        if mode is not None and row.get("active_pair", MODE_ACTIVE_PAIR[mode]) != MODE_ACTIVE_PAIR[mode]:
            raise ValueError(f"locus {label} has active_pair inconsistent with mode {mode}")

        mutation_id = row.get("mutation_id", "").strip()
        if seeded:
            if not mutation_id or mutation_id == "NA":
                raise ValueError(f"seeded locus {label} lacks a mutation_id")
            mutation_ids.add(mutation_id)
        elif mutation_id not in {"", "NA"}:
            raise ValueError(f"unseeded locus {label} unexpectedly has mutation_id={mutation_id}")

        design = row.get("seeding_design", "")
        expected_design = (
            "balanced_16_haplotype_cycle" if mode in (None, 0)
            else "balanced_four_haplotype_cycle"
        )
        if seeded and design and design != expected_design:
            raise ValueError(f"locus {label} has an unexpected seeding design")
        if not seeded and design and design != "not_seeded_for_mode":
            raise ValueError(f"locus {label} has an unexpected unseeded design")
        for column in (
            "global_frequency", "p1_frequency", "p2_frequency",
            "p3_frequency", "p4_frequency",
        ):
            if column in row and row[column]:
                frequency = float(row[column])
                if seeded and not 0.45 <= frequency <= 0.55:
                    raise ValueError(
                        f"locus {label} has unbalanced initial {column}={frequency}"
                    )
                if not seeded and frequency != 0.0:
                    raise ValueError(f"unseeded locus {label} has non-zero {column}={frequency}")
        if expected_positions is not None:
            if int(row["position"]) != expected_positions[label]:
                raise ValueError(
                    f"locus {label} is at {row['position']}, expected {expected_positions[label]}"
                )
    if len(mutation_ids) != len(expected_seeded):
        raise ValueError("seeded focal loci must have distinct mutation IDs")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=Path)
    parser.add_argument("--a-position", type=int)
    parser.add_argument("--b-position", type=int)
    parser.add_argument("--c-position", type=int)
    parser.add_argument("--d-position", type=int)
    parser.add_argument("--mode", type=int, choices=(0, 1, 2))
    args = parser.parse_args()
    supplied = [args.a_position, args.b_position, args.c_position, args.d_position]
    if any(value is not None for value in supplied) and any(value is None for value in supplied):
        parser.error("supply all four expected positions or none")
    expected = None
    if all(value is not None for value in supplied):
        expected = dict(zip(("A", "B", "C", "D"), supplied, strict=True))
    try:
        validate(args.path, expected, args.mode)
    except (OSError, ValueError) as error:
        raise SystemExit(str(error)) from error


if __name__ == "__main__":
    main()
