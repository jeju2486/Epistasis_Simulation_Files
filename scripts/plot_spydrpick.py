#!/usr/bin/env python3
"""Create per-case all-pair MI heatmaps and an aggregate truth-rank plot."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from simflow import read_tsv, repo_path
from spydrpick_case import parse_edge


def plot_case(case_dir: Path, genome_length: int, bins: int, result_name: str) -> list[dict[str, str]]:
    result = case_dir / result_name
    truth = read_tsv(result / "truth_pairs.tsv")
    total_pairs = int(truth[0]["total_pairs"])
    rank_stride = max(1, math.ceil(total_pairs / 100_000))
    sampled_ranks: list[int] = []
    sampled_mi: list[float] = []
    maximum = np.full((bins, bins), np.nan)
    with gzip.open(result / "spydrpick.edges.gz", "rt", encoding="utf-8") as handle:
        for rank, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            p1, p2, _distance, _aracne, mi = parse_edge(line)
            if rank == 1 or rank == total_pairs or rank % rank_stride == 0:
                sampled_ranks.append(rank)
                sampled_mi.append(mi)
            i = min(bins - 1, p1 * bins // genome_length)
            j = min(bins - 1, p2 * bins // genome_length)
            if math.isnan(maximum[i, j]) or mi > maximum[i, j]:
                maximum[i, j] = mi
                maximum[j, i] = mi

    fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    image = ax.imshow(
        maximum, origin="lower", extent=(0, genome_length, 0, genome_length),
        cmap="viridis", aspect="equal", interpolation="nearest",
    )
    for row, color in zip(truth, ("#e45756", "#f2cf5b"), strict=True):
        if row["candidate_status"] != "eligible":
            continue
        x, y = int(row["u_position"]), int(row["v_position"])
        ax.scatter([x], [y], s=80, facecolors="none", edgecolors=color, linewidths=2, label=row["pair"])
    ax.set(xlabel="Locus 1 position (bp)", ylabel="Locus 2 position (bp)", title=f"All-pair SpydrPick MI: {case_dir.name}")
    ax.legend(frameon=False)
    fig.colorbar(image, ax=ax, label="Maximum MI within genomic bin pair")
    fig.savefig(result / "all_pair_mi_heatmap.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
    ax.plot(np.asarray(sampled_ranks) / total_pairs, sampled_mi, color="#4c78a8", linewidth=1)
    for row, color in zip(truth, ("#e45756", "#f2cf5b"), strict=True):
        if row["candidate_status"] != "eligible":
            continue
        ax.scatter(float(row["rank_fraction"]), float(row["mi"]), s=65, color=color, label=row["pair"], zorder=3)
    ax.set_xscale("log")
    ax.set(xlabel="MI rank / all pairs (lower is better)", ylabel="Mutual information", title=f"SpydrPick all-pair rank curve: {case_dir.name}")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)
    fig.savefig(result / "mi_rank_curve.png", dpi=180)
    plt.close(fig)
    return truth


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="manifests/cases.tsv")
    parser.add_argument("--output-dir")
    parser.add_argument("--bins", type=int, default=150)
    parser.add_argument("--unweighted", action="store_true")
    args = parser.parse_args()
    if args.bins < 10:
        parser.error("--bins must be at least 10")
    result_name = "spydrpick_all_pairs_unweighted" if args.unweighted else "spydrpick_all_pairs"
    aggregate: list[dict[str, str]] = []
    for case in read_tsv(repo_path(args.manifest)):
        case_dir = repo_path(case["out_dir"])
        if not (case_dir / result_name / "_SUCCESS").exists():
            raise SystemExit(f"missing completed SpydrPick result: {case_dir / result_name}")
        for row in plot_case(case_dir, int(case["genome_length"]), args.bins, result_name):
            aggregate.append({**case, **row})

    output = repo_path(args.output_dir or f"results/{result_name}")
    output.mkdir(parents=True, exist_ok=True)
    fields = [
        "case_id", "replicate", "mode", "cross_hgt_probability", "pair", "candidate_status", "mi",
        "rank_min", "total_pairs", "rank_fraction", "u_column", "v_column",
        "u_position", "v_position", "physical_distance",
    ]
    with (output / "truth_pair_ranks.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(aggregate)

    groups: dict[tuple[str, str, str], list[float]] = {}
    for row in aggregate:
        if row["candidate_status"] != "eligible":
            continue
        key = (row["mode"], row["cross_hgt_probability"], row["pair"])
        groups.setdefault(key, []).append(float(row["rank_fraction"]))
    keys = sorted(groups, key=lambda key: (int(key[0]), float(key[1]), key[2]))
    if not groups:
        raise SystemExit("no eligible AB or CD truth pairs were available to plot")
    data = [groups[key] for key in keys]
    labels = [f"m{mode}\nhgt={hgt}\n{pair}" for mode, hgt, pair in keys]
    fig, ax = plt.subplots(figsize=(max(8, len(keys) * 0.75), 5.2), constrained_layout=True)
    ax.boxplot(data, tick_labels=labels, showmeans=True)
    ax.set_yscale("log")
    ax.set(ylabel="MI rank / all pairs (lower is better)", title="SpydrPick rank of truth pairs")
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(output / "truth_pair_rank_summary.png", dpi=180)
    plt.close(fig)
    print(f"[done] plots and summary: {output}")


if __name__ == "__main__":
    main()
