#!/usr/bin/env python3
"""Plot PAN-GWES-style physical distance versus MI for every case."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from itertools import combinations
from pathlib import Path

from simflow import read_tsv, repo_path
from spydrpick_case import parse_edge


FOCAL_LABELS = ("A", "B", "C", "D")
FOCAL_COLORS = {
    "AB": "#D55E00",
    "AC": "#E69F00",
    "AD": "#CC79A7",
    "BC": "#009E73",
    "BD": "#56B4E9",
    "CD": "#0072B2",
}


def eligible_positions(path: Path) -> list[int]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    columns = [int(row["filtered_column"]) for row in rows]
    if columns != list(range(len(rows))):
        raise ValueError("eligible-locus columns must be contiguous and zero-based")
    return [int(row["slim_position"]) for row in rows]


def focal_pair_columns(selected_loci: Path, positions: list[int]) -> dict[tuple[int, int], str]:
    focal = {row["label"]: int(row["position"]) for row in read_tsv(selected_loci)}
    if set(focal) != set(FOCAL_LABELS):
        raise ValueError("selected_loci.tsv must contain exactly A, B, C and D")
    position_to_column = {position: column for column, position in enumerate(positions)}
    pairs: dict[tuple[int, int], str] = {}
    for left, right in combinations(FOCAL_LABELS, 2):
        if focal[left] in position_to_column and focal[right] in position_to_column:
            u, v = sorted((position_to_column[focal[left]], position_to_column[focal[right]]))
            pairs[(u, v)] = left + right
    return pairs


def distance_mi_points(
    edges: Path,
    positions: list[int],
    focal_columns: dict[tuple[int, int], str],
    max_background_points: int,
) -> tuple[list[float], list[float], dict[str, tuple[float, float]]]:
    total_pairs = len(positions) * (len(positions) - 1) // 2
    stride = max(1, math.ceil(max(1, total_pairs - len(focal_columns)) / max_background_points))
    background_distance: list[float] = []
    background_mi: list[float] = []
    focal: dict[str, tuple[float, float]] = {}
    with gzip.open(edges, "rt", encoding="utf-8") as handle:
        for row_index, line in enumerate(handle):
            if not line.strip():
                continue
            col1, col2, distance, _aracne, mi = parse_edge(line)
            key = tuple(sorted((col1, col2)))
            distance_kb = distance / 1000.0
            if key in focal_columns:
                focal[focal_columns[key]] = (distance_kb, mi)
            elif row_index % stride == 0:
                background_distance.append(distance_kb)
                background_mi.append(mi)
    return background_distance, background_mi, focal


def plot_case(
    case_dir: Path,
    case_id: str,
    result_name: str,
    max_background_points: int,
) -> list[dict[str, str]]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    result = case_dir / result_name
    positions = eligible_positions(result / "eligible_loci.tsv")
    focal_columns = focal_pair_columns(case_dir / "selected_loci.tsv", positions)
    distances, mi_values, focal = distance_mi_points(
        result / "spydrpick.edges.gz", positions, focal_columns, max_background_points
    )

    fig, ax = plt.subplots(figsize=(7.6, 5.4), constrained_layout=True)
    ax.scatter(
        distances, mi_values, s=5, color="#777777", alpha=0.22,
        linewidths=0, rasterized=True, label="Other eligible pairs",
    )
    for pair in FOCAL_COLORS:
        if pair not in focal:
            continue
        distance, mi = focal[pair]
        ax.scatter(
            [distance], [mi], s=70, color=FOCAL_COLORS[pair],
            edgecolors="black", linewidths=0.6, label=pair, zorder=4,
        )
        ax.annotate(pair, (distance, mi), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set(
        xlabel="Physical distance between SNPs (kb)", ylabel="Mutual information (MI)",
        title=f"SpydrPick covariation by genomic distance: {case_id}",
    )
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.18, linewidth=0.6)
    ax.legend(frameon=False, ncols=2, fontsize=8)
    fig.savefig(result / "mi_vs_distance.png", dpi=220)
    plt.close(fig)
    return read_tsv(result / "truth_pairs.tsv")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="manifests/cases.tsv")
    parser.add_argument("--output-dir")
    parser.add_argument("--max-background-points", type=int, default=250_000)
    parser.add_argument("--unweighted", action="store_true")
    args = parser.parse_args()
    if args.max_background_points < 1:
        parser.error("--max-background-points must be positive")
    result_name = "spydrpick_all_pairs_unweighted" if args.unweighted else "spydrpick_all_pairs"
    aggregate: list[dict[str, str]] = []
    for case in read_tsv(repo_path(args.manifest)):
        case_dir = repo_path(case["out_dir"])
        if not (case_dir / result_name / "_SUCCESS").exists():
            raise SystemExit(f"missing completed SpydrPick result: {case_dir / result_name}")
        for row in plot_case(case_dir, case["case_id"], result_name, args.max_background_points):
            aggregate.append({**case, **row})

    output = repo_path(args.output_dir or f"results/{result_name}")
    output.mkdir(parents=True, exist_ok=True)
    fields = [
        "case_id", "replicate", "mode", "active_pair", "cross_hgt_probability",
        "experiment_generations", "pair", "candidate_status", "mi",
        "rank_min", "total_pairs", "rank_fraction", "u_column", "v_column",
        "u_position", "v_position", "physical_distance",
    ]
    with (output / "truth_pair_ranks.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(aggregate)
    print(f"[done] distance-MI plots and diagnostic rank table: {output}")


if __name__ == "__main__":
    main()
