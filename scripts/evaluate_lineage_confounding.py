#!/usr/bin/env python3
"""Evaluate lineage-driven covariance without exposing labels to either method."""

from __future__ import annotations

import argparse
from array import array
import csv
import gzip
import heapq
import math
from pathlib import Path

from kovar_inputs import read_fasta
from plot_spydrpick import eligible_positions, focal_pair_columns
from run_kovar_case import RESULT_NAME
from simflow import read_tsv, repo_path
from spydrpick_case import parse_edge


P_FLOOR = 1e-300
COVARIANCE_MINIMUM = 0.01
COMPONENT_FRACTION = 0.80
FOCAL_WINDOW_BP = 5_000
TOP_FRACTION = 0.01


def triangular_index(u: int, v: int, n_loci: int) -> int:
    if u == v:
        raise ValueError("self-pairs have no triangular index")
    u, v = sorted((u, v))
    if not (0 <= u < v < n_loci):
        raise ValueError(f"pair outside locus universe: {u}, {v}")
    return u * (2 * n_loci - u - 1) // 2 + (v - u - 1)


def decompose_covariance(
    n11: int,
    n10: int,
    n01: int,
    n00: int,
    population_u: list[float],
    population_v: list[float],
    population_weights: list[float],
) -> tuple[float, float, float, float]:
    total_n = n11 + n10 + n01 + n00
    if total_n <= 0:
        raise ValueError("joint table is empty")
    if not (len(population_u) == len(population_v) == len(population_weights)):
        raise ValueError("population-frequency vectors have unequal lengths")
    p_u = (n11 + n10) / total_n
    p_v = (n11 + n01) / total_n
    total = n11 / total_n - p_u * p_v
    between = sum(
        weight * (freq_u - p_u) * (freq_v - p_v)
        for freq_u, freq_v, weight in zip(
            population_u, population_v, population_weights, strict=True
        )
    )
    within = total - between
    fraction = abs(between) / (abs(between) + abs(within) + 1e-15)
    return total, within, between, fraction


def classify_pair(
    *,
    is_focal: bool,
    focal_proximal: bool,
    distance_bp: int,
    total_covariance: float,
    within_covariance: float,
    lineage_fraction: float,
) -> str:
    if is_focal:
        return "focal_AB"
    if focal_proximal:
        return "focal_proximal"
    if distance_bp <= FOCAL_WINDOW_BP:
        return "short_distance"
    if abs(total_covariance) >= COVARIANCE_MINIMUM and lineage_fraction >= COMPONENT_FRACTION:
        return "lineage_driven"
    within_fraction = abs(within_covariance) / (
        abs(within_covariance) + abs(total_covariance - within_covariance) + 1e-15
    )
    if abs(total_covariance) >= COVARIANCE_MINIMUM and within_fraction >= COMPONENT_FRACTION:
        return "within_population"
    return "other_distant"


def population_frequencies(
    binary_fasta: Path, sample_map: Path
) -> tuple[list[str], dict[str, list[float]], list[float], int]:
    mapping = {row["analysis_label"]: row["population"] for row in read_tsv(sample_map)}
    records = read_fasta(binary_fasta)
    n_loci = len(records[0][1])
    populations = sorted(set(mapping.values()))
    counts = {population: [0] * n_loci for population in populations}
    sizes = {population: 0 for population in populations}
    for sample, sequence in records:
        if sample not in mapping:
            raise ValueError(f"sample {sample!r} is absent from sample_names.tsv")
        population = mapping[sample]
        sizes[population] += 1
        target = counts[population]
        for column, state in enumerate(sequence):
            target[column] += state == "C"
    if any(size == 0 for size in sizes.values()):
        raise ValueError("one or more populations have no samples")
    total = sum(sizes.values())
    frequencies = {
        population: [count / sizes[population] for count in counts[population]]
        for population in populations
    }
    weights = [sizes[population] / total for population in populations]
    return populations, frequencies, weights, n_loci


def load_spydrpick(edges: Path, n_loci: int) -> tuple[array, array]:
    expected = n_loci * (n_loci - 1) // 2
    mi = array("d", [math.nan]) * expected
    ranks = array("I", [0]) * expected
    with gzip.open(edges, "rt", encoding="utf-8") as handle:
        for rank, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            u, v, _distance, _aracne, value = parse_edge(line)
            index = triangular_index(u, v, n_loci)
            if ranks[index]:
                raise ValueError(f"duplicate SpydrPick pair: {u}, {v}")
            mi[index] = value
            ranks[index] = rank
    if any(rank == 0 for rank in ranks):
        raise ValueError("SpydrPick output does not contain the complete pair universe")
    return mi, ranks


def kovar_top_threshold(path: Path) -> tuple[float, int]:
    def values():
        with path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                try:
                    value = float(row["p_primary"])
                except (KeyError, TypeError, ValueError):
                    continue
                if math.isfinite(value) and 0 <= value <= 1:
                    yield value

    finite_count = 0
    for _value in values():
        finite_count += 1
    if finite_count == 0:
        return -1.0, 0
    keep = max(1, math.ceil(TOP_FRACTION * finite_count))
    smallest = heapq.nsmallest(keep, values())
    return smallest[-1], finite_count


def evaluate_case(case: dict[str, str], output_name: str) -> list[dict[str, object]]:
    case_dir = repo_path(case["out_dir"])
    spydrpick = case_dir / "spydrpick_all_pairs"
    kovar = case_dir / RESULT_NAME
    for marker in (case_dir / "_SUCCESS", spydrpick / "_SUCCESS", kovar / "_SUCCESS"):
        if not marker.exists():
            raise FileNotFoundError(f"missing completed stage: {marker}")

    positions = eligible_positions(spydrpick / "eligible_loci.tsv")
    focal_columns = focal_pair_columns(case_dir / "selected_loci.tsv", positions)
    focal_positions = [positions[column] for column in focal_columns]
    populations, frequencies, weights, n_loci = population_frequencies(
        spydrpick / "all_snps.binary_ac.fa", case_dir / "sample_names.tsv"
    )
    if n_loci != len(positions):
        raise ValueError("eligible position map and binary FASTA have different lengths")
    mi, ranks = load_spydrpick(spydrpick / "spydrpick.edges.gz", n_loci)
    expected_pairs = n_loci * (n_loci - 1) // 2
    spyd_top_limit = max(1, math.ceil(TOP_FRACTION * expected_pairs))
    result_table = kovar / "ko_variation.tsv"
    kovar_threshold, finite_kovar = kovar_top_threshold(result_table)

    fields = [
        "case_id", "u", "v", "u_position", "v_position", "distance_bp", "category",
        "total_covariance", "within_population_covariance", "between_population_covariance",
        "lineage_fraction", "spydrpick_mi", "spydrpick_rank", "spydrpick_top_1pct",
        "kovar_p_primary", "kovar_top_1pct", "kovar_bonferroni",
    ]
    metrics_path = case_dir / output_name
    summaries: dict[str, dict[str, float | int]] = {}
    with gzip.open(metrics_path, "wt", encoding="utf-8", newline="") as output, result_table.open(
        encoding="utf-8", newline=""
    ) as source:
        writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in csv.DictReader(source, delimiter="\t"):
            u, v = int(row["u"]), int(row["v"])
            index = triangular_index(u, v, n_loci)
            try:
                counts = tuple(int(float(row[name])) for name in ("n11", "n10", "n01", "n00"))
            except (KeyError, TypeError, ValueError):
                continue
            population_u = [frequencies[population][u] for population in populations]
            population_v = [frequencies[population][v] for population in populations]
            total, within, between, lineage_fraction = decompose_covariance(
                *counts, population_u, population_v, weights
            )
            distance = abs(positions[v] - positions[u])
            is_focal = tuple(sorted((u, v))) == focal_columns
            focal_proximal = any(
                abs(position - focal) <= FOCAL_WINDOW_BP
                for position in (positions[u], positions[v])
                for focal in focal_positions
            )
            category = classify_pair(
                is_focal=is_focal, focal_proximal=focal_proximal, distance_bp=distance,
                total_covariance=total, within_covariance=within,
                lineage_fraction=lineage_fraction,
            )
            try:
                p_value = float(row["p_primary"])
            except (KeyError, TypeError, ValueError):
                p_value = math.nan
            spyd_top = ranks[index] <= spyd_top_limit
            kovar_top = (
                kovar_threshold >= 0 and math.isfinite(p_value) and p_value <= kovar_threshold
            )
            try:
                bonferroni = int(float(row.get("bonferroni_significant", "0") or 0))
            except (TypeError, ValueError):
                bonferroni = 0
            writer.writerow({
                "case_id": case["case_id"], "u": u, "v": v,
                "u_position": positions[u], "v_position": positions[v],
                "distance_bp": distance, "category": category,
                "total_covariance": total,
                "within_population_covariance": within,
                "between_population_covariance": between,
                "lineage_fraction": lineage_fraction, "spydrpick_mi": mi[index],
                "spydrpick_rank": ranks[index], "spydrpick_top_1pct": int(spyd_top),
                "kovar_p_primary": p_value if math.isfinite(p_value) else "NA",
                "kovar_top_1pct": int(kovar_top), "kovar_bonferroni": bonferroni,
            })
            summary = summaries.setdefault(category, {
                "pairs": 0, "mi_sum": 0.0, "spydrpick_top_1pct": 0,
                "finite_kovar": 0, "kovar_neglog10_sum": 0.0,
                "kovar_top_1pct": 0, "kovar_bonferroni": 0,
            })
            summary["pairs"] += 1
            summary["mi_sum"] += mi[index]
            summary["spydrpick_top_1pct"] += int(spyd_top)
            if math.isfinite(p_value):
                summary["finite_kovar"] += 1
                summary["kovar_neglog10_sum"] += -math.log10(max(p_value, P_FLOOR))
                summary["kovar_top_1pct"] += int(kovar_top)
                summary["kovar_bonferroni"] += bonferroni

    output_rows: list[dict[str, object]] = []
    for category, summary in sorted(summaries.items()):
        pairs = int(summary["pairs"])
        finite = int(summary["finite_kovar"])
        output_rows.append({
            **case, "category": category, "pairs": pairs,
            "mean_spydrpick_mi": float(summary["mi_sum"]) / pairs,
            "spydrpick_top_1pct": int(summary["spydrpick_top_1pct"]),
            "spydrpick_top_1pct_rate": int(summary["spydrpick_top_1pct"]) / pairs,
            "finite_kovar": finite,
            "mean_kovar_neglog10_p": (
                float(summary["kovar_neglog10_sum"]) / finite if finite else "NA"
            ),
            "kovar_top_1pct": int(summary["kovar_top_1pct"]),
            "kovar_top_1pct_rate": int(summary["kovar_top_1pct"]) / finite if finite else "NA",
            "kovar_bonferroni": int(summary["kovar_bonferroni"]),
            "kovar_bonferroni_rate": int(summary["kovar_bonferroni"]) / finite if finite else "NA",
            "finite_kovar_universe": finite_kovar,
        })
    return output_rows


def plot_summary(rows: list[dict[str, object]], output: Path) -> None:
    lineage = [row for row in rows if row["category"] == "lineage_driven"]
    if not lineage:
        return
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    modes = sorted({int(row["mode"]) for row in lineage})
    fig, axes = plt.subplots(1, len(modes), figsize=(4.3 * len(modes), 4.2),
                             sharey=True, constrained_layout=True)
    if len(modes) == 1:
        axes = [axes]
    for ax, mode in zip(axes, modes, strict=True):
        subset = [row for row in lineage if int(row["mode"]) == mode]
        cross_values = sorted({float(row["cross_hgt_probability"]) for row in subset})
        for method, field, color, marker in (
            ("SpydrPick", "spydrpick_top_1pct_rate", "#0072B2", "o"),
            ("KOVAR", "kovar_top_1pct_rate", "#D55E00", "s"),
        ):
            for cross_index, cross in enumerate(cross_values):
                values = [
                    float(row[field]) for row in subset
                    if math.isclose(float(row["cross_hgt_probability"]), cross)
                    and row[field] != "NA"
                ]
                ax.scatter([cross_index] * len(values), values, color=color, marker=marker,
                           s=35, alpha=0.8, label=method if cross_index == 0 else None)
        ax.set(xticks=range(len(cross_values)),
               xticklabels=[f"{value:g}" for value in cross_values],
               xlabel="Cross-HGT probability", title=f"Mode {mode}")
        ax.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("Fraction of lineage-driven pairs in method top 1%")
    axes[0].legend(frameon=False)
    fig.savefig(output / "lineage_driven_top1_by_hgt.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, len(modes), figsize=(4.3 * len(modes), 4.2),
                             sharey=True, constrained_layout=True)
    if len(modes) == 1:
        axes = [axes]
    for ax, mode in zip(axes, modes, strict=True):
        subset = [row for row in lineage if int(row["mode"]) == mode]
        cross_values = sorted({float(row["cross_hgt_probability"]) for row in subset})
        for cross_index, cross in enumerate(cross_values):
            values = [
                int(row["pairs"]) for row in subset
                if math.isclose(float(row["cross_hgt_probability"]), cross)
            ]
            ax.scatter([cross_index] * len(values), values, color="#5F5F5F", s=35, alpha=0.8)
        ax.set(xticks=range(len(cross_values)),
               xticklabels=[f"{value:g}" for value in cross_values],
               xlabel="Cross-HGT probability", title=f"Mode {mode}")
        ax.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("Number of lineage-driven distant pairs")
    fig.savefig(output / "lineage_driven_pair_count_by_hgt.png", dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="manifests/cases.tsv")
    parser.add_argument("--output-dir", default="results/lineage_confounding")
    parser.add_argument("--pair-metrics-name", default="lineage_pair_metrics.tsv.gz")
    args = parser.parse_args()
    rows: list[dict[str, object]] = []
    for case in read_tsv(repo_path(args.manifest)):
        print(f"[evaluate] {case['case_id']}")
        rows.extend(evaluate_case(case, args.pair_metrics_name))
    output = repo_path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with (output / "pair_category_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    plot_summary(rows, output)
    print(f"[done] lineage-confounding evaluation: {output}")


if __name__ == "__main__":
    main()
