#!/usr/bin/env python3
"""Plot physical distance versus directional KOVAR significance."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from plot_spydrpick import FOCAL_COLORS, eligible_positions, focal_pair_columns
from simflow import read_tsv, repo_path


DIRECTION_MARKERS = {"u_predicts_v": "^", "v_predicts_u": "v"}
DIRECTION_LABELS = {"u_predicts_v": "u predicts v", "v_predicts_u": "v predicts u"}
P_FLOOR = 1e-300


def kovar_points(
    result_table: Path,
    positions: list[int],
    focal_columns: dict[tuple[int, int], str],
    expected_pairs: int,
    max_background_points: int,
) -> tuple[dict[str, tuple[list[float], list[float]]], dict[str, list[tuple[str, float, float]]], int]:
    expected_rows = max(1, 2 * expected_pairs)
    stride = max(1, math.ceil(expected_rows / max_background_points))
    background = {direction: ([], []) for direction in DIRECTION_MARKERS}
    focal = {pair: [] for pair in FOCAL_COLORS}
    n_directional_tests = 0
    with result_table.open(encoding="utf-8", newline="") as handle:
        for row_index, row in enumerate(csv.DictReader(handle, delimiter="\t")):
            direction = row.get("direction", "")
            if direction not in DIRECTION_MARKERS:
                raise ValueError(f"unknown KOVAR direction: {direction!r}")
            try:
                p_value = float(row["p_primary"])
            except (KeyError, TypeError, ValueError):
                continue
            if not math.isfinite(p_value) or p_value < 0 or p_value > 1:
                continue
            if row.get("n_directional_tests") not in (None, "", "NA"):
                n_directional_tests = max(n_directional_tests, int(float(row["n_directional_tests"])))
            col1, col2 = int(row["u"]), int(row["v"])
            if not (0 <= col1 < len(positions) and 0 <= col2 < len(positions)):
                raise ValueError(f"KOVAR locus outside eligible alignment: {col1}, {col2}")
            key = tuple(sorted((col1, col2)))
            distance_kb = abs(positions[col2] - positions[col1]) / 1000.0
            score = -math.log10(max(p_value, P_FLOOR))
            if key in focal_columns:
                focal[focal_columns[key]].append((direction, distance_kb, score))
            elif row_index % stride == 0:
                background[direction][0].append(distance_kb)
                background[direction][1].append(score)
    return background, focal, n_directional_tests


def plot_case(case_dir: Path, case_id: str, result_name: str, max_background_points: int) -> None:
    result = case_dir / result_name
    spydrpick = case_dir / "spydrpick_all_pairs"
    positions = eligible_positions(spydrpick / "eligible_loci.tsv")
    focal_columns = focal_pair_columns(case_dir / "selected_loci.tsv", positions)
    metadata = json.loads((result / "run_metadata.json").read_text(encoding="utf-8"))
    background, focal, n_tests = kovar_points(
        result / "results" / "ko_variation.tsv", positions, focal_columns,
        int(metadata["pairs"]), max_background_points,
    )

    fig, ax = plt.subplots(figsize=(7.6, 5.4), constrained_layout=True)
    for direction, marker in DIRECTION_MARKERS.items():
        ax.scatter(
            background[direction][0], background[direction][1], s=6, marker=marker,
            color="#777777", alpha=0.18, linewidths=0, rasterized=True,
        )
    for pair, color in FOCAL_COLORS.items():
        for direction, distance, score in focal[pair]:
            ax.scatter(
                [distance], [score], s=72, marker=DIRECTION_MARKERS[direction],
                color=color, edgecolors="black", linewidths=0.6, zorder=4,
            )
            ax.annotate(pair, (distance, score), xytext=(4, 4), textcoords="offset points", fontsize=8)
    if n_tests > 0:
        threshold = -math.log10(0.05 / n_tests)
        ax.axhline(
            threshold, color="#333333", linestyle="--", linewidth=0.9,
            label="Bonferroni 0.05",
        )
    pair_handles = [
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=color,
               markeredgecolor="black", markersize=6, label=pair)
        for pair, color in FOCAL_COLORS.items() if focal[pair]
    ]
    direction_handles = [
        Line2D([0], [0], marker=marker, linestyle="", color="#777777",
               markersize=6, label=DIRECTION_LABELS[direction])
        for direction, marker in DIRECTION_MARKERS.items()
    ]
    threshold_handles, threshold_labels = ax.get_legend_handles_labels()
    handles = pair_handles + direction_handles + threshold_handles
    labels = [handle.get_label() for handle in pair_handles + direction_handles] + threshold_labels
    ax.legend(handles, labels, frameon=False, ncols=2, fontsize=8)
    ax.set(
        xlabel="Physical distance between SNPs (kb)",
        ylabel="−log10(KOVAR directional primary p-value)",
        title=f"KOVAR covariation by genomic distance: {case_id}",
    )
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.18, linewidth=0.6)
    fig.savefig(result / "kovar_p_vs_distance.png", dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="manifests/cases.tsv")
    parser.add_argument("--tree-mode", choices=("core", "oracle", "grm"), default="oracle")
    parser.add_argument("--result-name")
    parser.add_argument("--max-background-points", type=int, default=250_000)
    args = parser.parse_args()
    if args.max_background_points < 1:
        parser.error("--max-background-points must be positive")
    result_name = args.result_name or (
        "kovar_v081" if args.tree_mode == "core" else f"kovar_v081_{args.tree_mode}"
    )
    for case in read_tsv(repo_path(args.manifest)):
        case_dir = repo_path(case["out_dir"])
        result = case_dir / result_name
        if not (result / "_SUCCESS").exists():
            raise SystemExit(f"missing completed KOVAR result: {result}")
        plot_case(case_dir, case["case_id"], result_name, args.max_background_points)
    print(f"[done] KOVAR distance-significance plots: {result_name}")


if __name__ == "__main__":
    main()
