#!/usr/bin/env python3
"""Plot weighted all-pair SpydrPick MI and the focal A-B result."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path

from simflow import read_tsv, repo_path
from spydrpick_case import parse_edge


FOCAL_COLOR = "#D55E00"


def eligible_positions(path: Path) -> list[int]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    columns = [int(row["filtered_column"]) for row in rows]
    if columns != list(range(len(rows))):
        raise ValueError("eligible-locus columns must be contiguous and zero-based")
    return [int(row["slim_position"]) for row in rows]


def focal_pair_columns(selected_loci: Path, positions: list[int]) -> tuple[int, int]:
    focal = {row["label"]: int(row["position"]) for row in read_tsv(selected_loci)}
    if set(focal) != {"A", "B"}:
        raise ValueError("selected_loci.tsv must contain exactly A and B")
    position_to_column = {position: column for column, position in enumerate(positions)}
    if focal["A"] not in position_to_column or focal["B"] not in position_to_column:
        raise ValueError("A or B failed the shared MAF filter")
    return tuple(sorted((position_to_column[focal["A"]], position_to_column[focal["B"]])))


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return float("nan")
    index = probability * (len(ordered) - 1)
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] * (upper - index) + ordered[upper] * (index - lower)


def read_points(
    edges: Path, focal_columns: tuple[int, int], max_background_points: int,
) -> tuple[list[float], list[float], tuple[float, float, int], list[dict[str, float | int]]]:
    rows: list[tuple[float, float]] = []
    focal: tuple[float, float, int] | None = None
    bins: dict[int, list[float]] = {}
    with gzip.open(edges, "rt", encoding="utf-8") as handle:
        for rank, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            u, v, distance, _aracne, mi = parse_edge(line)
            distance_kb = distance / 1000.0
            if tuple(sorted((u, v))) == focal_columns:
                focal = (distance_kb, mi, rank)
            else:
                rows.append((distance_kb, mi))
            bins.setdefault(int(distance_kb), []).append(mi)
    if focal is None:
        raise ValueError("A-B pair is absent from the complete SpydrPick output")
    stride = max(1, math.ceil(len(rows) / max_background_points))
    sampled = rows[::stride]
    summaries = [
        {
            "bin_start_kb": start, "bin_end_kb": start + 1, "n_pairs": len(values),
            "median_mi": quantile(values, 0.5), "q95_mi": quantile(values, 0.95),
            "q99_mi": quantile(values, 0.99),
        }
        for start, values in sorted(bins.items())
    ]
    return [x for x, _ in sampled], [y for _, y in sampled], focal, summaries


def draw_plot(result: Path, case_id: str, distances: list[float], values: list[float],
              focal: tuple[float, float, int], maximum_kb: float | None, filename: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.6, 5.4), constrained_layout=True)
    ax.scatter(distances, values, s=5, color="#777777", alpha=0.22,
               linewidths=0, rasterized=True, label="Other eligible pairs")
    if maximum_kb is None or focal[0] <= maximum_kb:
        ax.scatter([focal[0]], [focal[1]], s=75, color=FOCAL_COLOR,
                   edgecolors="black", linewidths=0.6, label="AB", zorder=4)
        ax.annotate("AB", (focal[0], focal[1]), xytext=(4, 4),
                    textcoords="offset points", fontsize=9)
    ax.axvline(1.0, color="#333333", linestyle="--", linewidth=0.8, label="1 kb")
    ax.set(xlabel="Physical distance between SNPs (kb)", ylabel="Mutual information (MI)",
           title=f"SpydrPick covariation by genomic distance: {case_id}")
    ax.set_xlim(0, maximum_kb if maximum_kb is not None else None)
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.18, linewidth=0.6)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(result / filename, dpi=220)
    plt.close(fig)


def plot_case(case_dir: Path, case_id: str, max_background_points: int) -> dict[str, str | int | float]:
    result = case_dir / "spydrpick_all_pairs"
    positions = eligible_positions(result / "eligible_loci.tsv")
    focal_columns = focal_pair_columns(case_dir / "selected_loci.tsv", positions)
    distances, values, focal, summaries = read_points(
        result / "spydrpick.edges.gz", focal_columns, max_background_points
    )
    draw_plot(result, case_id, distances, values, focal, None, "mi_vs_distance.png")
    draw_plot(result, case_id, distances, values, focal, 5.0, "mi_vs_distance_0_5kb.png")
    with (result / "mi_distance_quantiles.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["bin_start_kb", "bin_end_kb", "n_pairs", "median_mi", "q95_mi", "q99_mi"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summaries)
    return {
        "u_column": focal_columns[0], "v_column": focal_columns[1],
        "physical_distance_kb": focal[0], "mi": focal[1], "mi_order_rank": focal[2],
        "eligible_loci": len(positions), "total_pairs": len(positions) * (len(positions) - 1) // 2,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="manifests/cases.tsv")
    parser.add_argument("--output-dir", default="results/spydrpick_all_pairs")
    parser.add_argument("--max-background-points", type=int, default=250_000)
    args = parser.parse_args()
    if args.max_background_points < 1:
        parser.error("--max-background-points must be positive")
    aggregate: list[dict[str, str | int | float]] = []
    for case in read_tsv(repo_path(args.manifest)):
        case_dir = repo_path(case["out_dir"])
        if not (case_dir / "spydrpick_all_pairs" / "_SUCCESS").exists():
            raise SystemExit(f"missing completed SpydrPick result: {case_dir}")
        aggregate.append({**case, **plot_case(case_dir, case["case_id"], args.max_background_points)})
    output = repo_path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    fields = list(aggregate[0]) if aggregate else []
    with (output / "focal_ab.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(aggregate)
    if aggregate:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(5.6, 4.6), constrained_layout=True)
        for mode in (1, 2):
            subset = [row for row in aggregate if int(row["mode"]) == mode]
            ax.scatter([mode] * len(subset), [float(row["mi"]) for row in subset],
                       s=45, color=FOCAL_COLOR, edgecolors="black", linewidths=0.5)
        ax.set(xticks=[1, 2], xticklabels=["Mode 1", "Mode 2"],
               ylabel="A-B mutual information",
               title="SpydrPick focal result across replicates")
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", alpha=0.2)
        fig.savefig(output / "focal_ab_by_mode.png", dpi=220)
        plt.close(fig)
    print(f"[done] SpydrPick plots and focal summary: {output}")


if __name__ == "__main__":
    main()
