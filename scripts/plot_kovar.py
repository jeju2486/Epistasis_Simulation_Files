#!/usr/bin/env python3
"""Plot unordered KOVAR 0.8.3 results, QQ curves, and the focal A-B test."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from plot_spydrpick import FOCAL_COLOR, eligible_positions, focal_pair_columns
from run_kovar_case import RESULT_NAME
from simflow import read_tsv, repo_path


P_FLOOR = 1e-300


def distance_bin(distance_kb: float) -> str:
    if distance_kb <= 1:
        return "0-1 kb"
    if distance_kb <= 5:
        return "1-5 kb"
    if distance_kb <= 20:
        return "5-20 kb"
    return ">20 kb"


def read_results(table: Path, positions: list[int], focal_columns: tuple[int, int]) -> tuple[list[dict], dict]:
    points: list[dict] = []
    focal: dict | None = None
    with table.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                p_value = float(row["p_primary"])
            except (KeyError, TypeError, ValueError):
                continue
            if not math.isfinite(p_value) or not 0 <= p_value <= 1:
                continue
            u, v = int(row["u"]), int(row["v"])
            if not (0 <= u < len(positions) and 0 <= v < len(positions)):
                raise ValueError(f"KOVAR locus outside eligible alignment: {u}, {v}")
            distance = abs(positions[v] - positions[u]) / 1000.0
            point = {
                "u": u, "v": v, "distance_kb": distance, "bin": distance_bin(distance),
                "p": p_value, "score": -math.log10(max(p_value, P_FLOOR)),
                "bonferroni": int(float(row.get("bonferroni_significant", "0") or 0)),
                "n_tests": int(float(row.get("n_tests", "0") or 0)),
            }
            points.append(point)
            if tuple(sorted((u, v))) == focal_columns:
                focal = point
    if focal is None:
        raise ValueError("A-B has no finite KOVAR primary p-value")
    return points, focal


def plot_distance(result: Path, case_id: str, points: list[dict], focal: dict,
                  max_background_points: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    stride = max(1, math.ceil(len(points) / max_background_points))
    background = [point for point in points if point is not focal][::stride]
    fig, ax = plt.subplots(figsize=(7.6, 5.4), constrained_layout=True)
    ax.scatter([p["distance_kb"] for p in background], [p["score"] for p in background],
               s=6, color="#777777", alpha=0.2, linewidths=0, rasterized=True,
               label="Other tested pairs")
    ax.scatter([focal["distance_kb"]], [focal["score"]], s=75, color=FOCAL_COLOR,
               edgecolors="black", linewidths=0.6, zorder=4, label="AB")
    ax.annotate("AB", (focal["distance_kb"], focal["score"]), xytext=(4, 4),
                textcoords="offset points", fontsize=9)
    n_tests = max((p["n_tests"] for p in points), default=0)
    if n_tests:
        ax.axhline(-math.log10(0.05 / n_tests), color="#333333", linestyle="--",
                   linewidth=0.9, label="Bonferroni 0.05")
    ax.set(xlabel="Physical distance between SNPs (kb)",
           ylabel="−log10(KOVAR primary p-value)",
           title=f"KOVAR covariation by genomic distance: {case_id}")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.18, linewidth=0.6)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(result / "kovar_p_vs_distance.png", dpi=220)
    plt.close(fig)


def plot_qq(result: Path, case_id: str, points: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.2, 5.5), constrained_layout=True)
    maximum = 0.0
    for label in ("0-1 kb", "1-5 kb", "5-20 kb", ">20 kb"):
        observed = sorted(point["p"] for point in points if point["bin"] == label)
        if not observed:
            continue
        expected = [(index + 0.5) / len(observed) for index in range(len(observed))]
        x = [-math.log10(value) for value in expected]
        y = [-math.log10(max(value, P_FLOOR)) for value in observed]
        maximum = max(maximum, max(x), max(y))
        ax.plot(x, y, marker=".", markersize=2.5, linewidth=0.7, label=f"{label} (n={len(observed)})")
    ax.plot([0, maximum], [0, maximum], color="#333333", linestyle="--", linewidth=0.8)
    ax.set(xlabel="Expected −log10(p)", ylabel="Observed −log10(p)",
           title=f"KOVAR QQ by physical distance: {case_id}")
    ax.grid(alpha=0.18, linewidth=0.6)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(result / "kovar_qq_by_distance.png", dpi=220)
    plt.close(fig)


def plot_case(case_dir: Path, case_id: str, max_background_points: int) -> tuple[dict, list[dict]]:
    result = case_dir / RESULT_NAME
    spydrpick = case_dir / "spydrpick_all_pairs"
    positions = eligible_positions(spydrpick / "eligible_loci.tsv")
    focal_columns = focal_pair_columns(case_dir / "selected_loci.tsv", positions)
    points, focal = read_results(result / "ko_variation.tsv", positions, focal_columns)
    plot_distance(result, case_id, points, focal, max_background_points)
    plot_qq(result, case_id, points)
    counts = []
    for label in ("0-1 kb", "1-5 kb", "5-20 kb", ">20 kb"):
        subset = [point for point in points if point["bin"] == label]
        counts.append({"distance_bin": label, "tested": len(subset),
                       "bonferroni_significant": sum(point["bonferroni"] for point in subset)})
    return {
        "u_column": focal["u"], "v_column": focal["v"],
        "physical_distance_kb": focal["distance_kb"], "p_primary": focal["p"],
        "neglog10_p": focal["score"], "bonferroni_significant": focal["bonferroni"],
        "n_tests": focal["n_tests"],
    }, counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="manifests/cases.tsv")
    parser.add_argument("--output-dir", default="results/kovar_v083_simulation_tree")
    parser.add_argument("--max-background-points", type=int, default=250_000)
    args = parser.parse_args()
    if args.max_background_points < 1:
        parser.error("--max-background-points must be positive")
    focal_rows: list[dict] = []
    count_rows: list[dict] = []
    for case in read_tsv(repo_path(args.manifest)):
        case_dir = repo_path(case["out_dir"])
        if not (case_dir / RESULT_NAME / "_SUCCESS").exists():
            raise SystemExit(f"missing completed KOVAR result: {case_dir / RESULT_NAME}")
        focal, counts = plot_case(case_dir, case["case_id"], args.max_background_points)
        focal_rows.append({**case, **focal})
        count_rows.extend({**case, **row} for row in counts)
    output = repo_path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    for filename, rows in (("focal_ab.tsv", focal_rows), ("significant_by_distance.tsv", count_rows)):
        fields = list(rows[0]) if rows else []
        with (output / filename).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    if focal_rows:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(5.6, 4.6), constrained_layout=True)
        for mode in (1, 2):
            subset = [row for row in focal_rows if int(row["mode"]) == mode]
            ax.scatter([mode] * len(subset), [float(row["neglog10_p"]) for row in subset],
                       s=45, color=FOCAL_COLOR, edgecolors="black", linewidths=0.5)
        ax.set(xticks=[1, 2], xticklabels=["Mode 1", "Mode 2"],
               ylabel="A-B −log10(KOVAR primary p-value)",
               title="KOVAR focal result across replicates")
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", alpha=0.2)
        fig.savefig(output / "focal_ab_by_mode.png", dpi=220)
        plt.close(fig)
    print(f"[done] KOVAR plots and summaries: {output}")


if __name__ == "__main__":
    main()
