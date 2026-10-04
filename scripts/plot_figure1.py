#!/usr/bin/env python3
"""Export manuscript Figure 1 panels from completed simulation analyses.

MI is never treated as a P value. Panel A requires an explicit null-pair list
and separately calibrated marginal P values. No simulation or model is rerun.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
from pathlib import Path

import numpy as np

from evaluate_lineage_confounding import (
    classify_pair, decompose_covariance, population_frequencies, triangular_index,
)
from plot_spydrpick import eligible_positions, focal_pair_columns
from run_kovar_case import result_name, sha256_file
from simflow import read_tsv, repo_path
from spydrpick_case import output_name, parse_edge

COLORS = {"MI": "#0072B2", "KOVAR": "#D55E00"}
CATEGORIES = {"other_distant": 0, "focal_AB": 1, "lineage_driven": 2,
              "focal_proximal": 3, "short_distance": 4, "within_population": 5,
              "unclassified": 6}
P_FLOOR = 1e-300


def average_ranks(values: np.ndarray, *, descending: bool = False) -> np.ndarray:
    """One-based average ranks; nonfinite values remain unavailable, not worst."""
    ranks = np.full(len(values), np.nan)
    finite = np.flatnonzero(np.isfinite(values))
    ordered = finite[np.argsort(-values[finite] if descending else values[finite])]
    if not len(ordered):
        return ranks
    starts = np.r_[0, np.flatnonzero(np.diff(values[ordered]) != 0) + 1]
    ends = np.r_[starts[1:], len(ordered)]
    ranks[ordered] = np.repeat((starts + ends + 1) / 2, ends - starts)
    return ranks


def valid_p(value: str) -> float:
    try:
        p = float(value)
    except (ValueError, TypeError):
        return math.nan
    return p if math.isfinite(p) and 0 <= p <= 1 else math.nan


def write_rows(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]),
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_case(case: dict, spa_mode: str, weighting: str) -> dict:
    """Join by unordered zero-based pair identity, validate the complete universe."""
    case_dir = repo_path(case["out_dir"])
    marginal = case_dir / output_name(weighting)
    kovar = case_dir / result_name(spa_mode)
    for stage in (case_dir, marginal, kovar):
        if not (stage / "_SUCCESS").exists():
            raise FileNotFoundError(f"missing completed stage: {stage}")
    positions = np.asarray(eligible_positions(marginal / "eligible_loci.tsv"))
    default_positions = eligible_positions(
        case_dir / "spydrpick_all_pairs" / "eligible_loci.tsv")
    if positions.tolist() != default_positions:
        raise ValueError(f"{case['case_id']}: MI and KOVAR position maps differ")
    if len(set(positions)) != len(positions):
        raise ValueError("eligible locus positions are not unique")
    n = len(positions)
    u, v = np.triu_indices(n, 1)
    total = len(u)
    if total == 0:
        raise ValueError(f"{case['case_id']}: no eligible pairs")
    mi, p = np.full(total, np.nan), np.full(total, np.nan)
    seen = np.zeros(total, dtype=bool)
    with gzip.open(marginal / "spydrpick.edges.gz", "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            a, b, distance, _aracne, value = parse_edge(line)
            index = triangular_index(a, b, n)
            if seen[index] or not math.isfinite(value) or value < 0:
                raise ValueError(f"invalid or duplicate marginal pair: {a}, {b}")
            if distance != abs(int(positions[a]) - int(positions[b])):
                raise ValueError(f"marginal distance disagrees with position map: {a}, {b}")
            seen[index], mi[index] = True, value
    if not seen.all():
        raise ValueError("marginal results are not the complete eligible pair universe")
    selected = {row["label"]: int(row["position"])
                for row in read_tsv(case_dir / "selected_loci.tsv")}
    focal = focal_pair_columns(case_dir / "selected_loci.tsv", positions.tolist())
    focal_index = triangular_index(*focal, n) if focal is not None else None
    pops, freq, weights, n_binary = population_frequencies(
        marginal / "all_snps.binary_ac.fa", case_dir / "sample_names.tsv")
    if n_binary != n:
        raise ValueError("binary alignment disagrees with eligible position map")
    category = np.full(total, CATEGORIES["unclassified"], dtype=np.uint8)
    seen[:] = False
    with (kovar / "ko_variation.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            a, b = int(row["u"]), int(row["v"])
            index = triangular_index(a, b, n)
            if seen[index]:
                raise ValueError(f"duplicate KOVAR pair: {a}, {b}")
            seen[index], p[index] = True, valid_p(row.get("p_primary"))
            try:
                raw_counts = [float(row[name]) for name in ("n11", "n10", "n01", "n00")]
                if any(not math.isfinite(x) or x < 0 or not x.is_integer() for x in raw_counts):
                    raise ValueError("invalid joint counts")
                counts = [int(x) for x in raw_counts]
                total_cov, within, _between, fraction = decompose_covariance(
                    *counts, [freq[pop][a] for pop in pops],
                    [freq[pop][b] for pop in pops], weights)
            except (KeyError, TypeError, ValueError):
                continue  # Preserve unavailable tests; do not manufacture lineage labels.
            label = classify_pair(
                is_focal=index == focal_index,
                focal_proximal=any(abs(int(pos) - anchor) <= 5_000
                                   for pos in (positions[a], positions[b])
                                   for anchor in selected.values()),
                distance_bp=abs(int(positions[b]) - int(positions[a])),
                total_covariance=total_cov, within_covariance=within,
                lineage_fraction=fraction)
            category[index] = CATEGORIES[label]
    if not seen.all():
        raise ValueError("KOVAR results do not contain the identical eligible pair universe; "
                         "failed tests must have rows with unavailable P values")
    if focal_index is not None:
        category[focal_index] = CATEGORIES["focal_AB"]
    finite = np.isfinite(p)
    mi_rank = average_ranks(mi, descending=True)
    p_rank = average_ranks(p)
    shared_mi_rank = np.full(total, np.nan)
    shared_mi_rank[finite] = average_ranks(mi[finite], descending=True)
    inputs = [marginal / "eligible_loci.tsv", marginal / "spydrpick.edges.gz",
              marginal / "all_snps.binary_ac.fa", kovar / "ko_variation.tsv",
              case_dir / "selected_loci.tsv", case_dir / "sample_names.tsv"]
    return {"case": case, "u": u, "v": v, "positions": positions,
            "distance": np.abs(positions[v] - positions[u]), "mi": mi, "p": p,
            "mi_rank": mi_rank, "p_rank": p_rank, "shared_mi_rank": shared_mi_rank,
            "category": category, "focal": focal_index, "n_pairs": total,
            "n_finite": int(finite.sum()),
            "input_sha256": {str(path): sha256_file(path) for path in inputs}}


def save(fig, output: Path, stem: str, formats: list[str]) -> list[str]:
    import matplotlib.pyplot as plt
    files = []
    for fmt in formats:
        path = output / f"{stem}.{fmt}"
        # SVG must retain editable points as vectors, not embedded raster images.
        rasterized = [artist for artist in fig.findobj() if artist.get_rasterized()]
        if fmt == "svg":
            for artist in rasterized:
                artist.set_rasterized(False)
        fig.savefig(path, dpi=300)
        if fmt == "svg":
            for artist in rasterized:
                artist.set_rasterized(True)
        files.append(path.name)
    plt.close(fig)
    return files


def sample_indices(indices: np.ndarray, maximum: int) -> np.ndarray:
    """Deterministic display-only thinning; statistics always use every pair."""
    if len(indices) <= maximum:
        return indices
    return indices[np.linspace(0, len(indices) - 1, maximum, dtype=int)]


def panel_b(data: dict, output: Path, args) -> list[str]:
    import matplotlib.pyplot as plt
    files = []
    for method, key, ylabel in (("MI", "mi", "Mutual information (MI)"),
                                 ("KOVAR", "p", "−log10(KOVAR primary P)")):
        y = data[key] if key == "mi" else -np.log10(np.maximum(data[key], P_FLOOR))
        indices = sample_indices(np.flatnonzero(np.isfinite(y)), args.max_points)
        fig, ax = plt.subplots(figsize=(6.2, 4.5), layout="constrained")
        ax.scatter(data["distance"][indices] / 1000, y[indices], s=5,
                   color="#777777", alpha=.25, linewidths=0, rasterized=True)
        focal = data["focal"]
        if focal is not None and math.isfinite(y[focal]):
            ax.scatter(data["distance"][focal] / 1000, y[focal], marker="*", s=110,
                       color="#CC79A7", edgecolor="black", linewidth=.5, label="A–B", zorder=4)
        ax.axvspan(0, args.distal_bp / 1000, color="#EEEEEE", zorder=0)
        ax.axvline(args.distal_bp / 1000, color="black", linestyle="--", linewidth=.8,
                   label=f"≤{args.distal_bp / 1000:g} kb excluded from distal network")
        ax.set(xlabel="Genomic separation (kb)", ylabel=ylabel,
               title=f"B  {data['case']['case_id']} · {method}", xlim=(0, None), ylim=(0, None))
        ax.legend(frameon=False, fontsize=8)
        files += save(fig, output, f"B_{data['case']['case_id']}_{method}", args.formats)
    return files


def panel_e(data: dict, output: Path, args, panel: str = "E") -> list[str]:
    import matplotlib.pyplot as plt
    finite = np.isfinite(data["p"])
    fig, ax = plt.subplots(figsize=(5.5, 5.1), layout="constrained")
    highlighted = (data["category"] == CATEGORIES["lineage_driven"]) & (data["distance"] > args.distal_bp)
    for mask, color, label, marker, size in (
        (finite & ~highlighted, "#888888", "Other pairs", ".", 8),
        (finite & highlighted, COLORS["MI"], "Distal lineage-associated", "o", 9)):
        indices = sample_indices(np.flatnonzero(mask), args.max_points)
        ax.scatter(data["shared_mi_rank"][indices], data["p_rank"][indices],
                   s=size, color=color, alpha=.35, marker=marker, linewidths=0,
                   rasterized=True, label=label)
    focal = data["focal"]
    if focal is not None and finite[focal]:
        ax.scatter(data["shared_mi_rank"][focal], data["p_rank"][focal],
                   marker="*", s=130, color="#CC79A7", edgecolor="black", linewidth=.5,
                   label="A–B", zorder=5)
    maximum = max(2, data["n_finite"])
    ax.plot([1, maximum], [1, maximum], "--", color="black", linewidth=.8)
    ax.set(xscale="log", yscale="log", xlim=(.8, maximum * 1.1), ylim=(.8, maximum * 1.1),
           xlabel="Marginal MI rank (1 = highest)", ylabel="KOVAR P rank (1 = lowest P)",
           title=f"{panel}  {data['case']['case_id']}")
    ax.text(.02, .02, f"Same {data['n_finite']:,} finite-test pairs ranked by both methods\n"
            f"{data['n_pairs'] - data['n_finite']:,} unavailable KOVAR tests reported separately",
            transform=ax.transAxes, fontsize=8, va="bottom")
    ax.legend(frameon=False, fontsize=8)
    return save(fig, output, f"{panel}_{data['case']['case_id']}_identical_pair_ranks", args.formats)


def top_mask(values: np.ndarray, fraction: float, *, descending: bool,
             budget: int | None = None) -> np.ndarray:
    finite = np.isfinite(values)
    mask = np.zeros(len(values), dtype=bool)
    count = int(finite.sum())
    if count:
        budget = min(count, max(1, math.ceil(fraction * count) if budget is None else budget))
        ordered = np.sort(values[finite])
        cutoff = ordered[-budget] if descending else ordered[budget - 1]
        mask = finite & (values >= cutoff if descending else values <= cutoff)
    return mask


def summarize(data: dict, args) -> tuple[dict, dict]:
    case = data["case"]
    finite = np.isfinite(data["p"])
    # Focal recovery uses the full candidate budget; failed tests cannot recover.
    mi_top = top_mask(data["mi"], args.top_fraction, descending=True)
    p_top = top_mask(data["p"], args.top_fraction, descending=False,
                     budget=math.ceil(args.top_fraction * data["n_pairs"]))
    focal = data["focal"]
    status = ("not_maf_eligible" if focal is None else
              "available" if finite[focal] else "kovar_p_unavailable")
    common = {"case_id": case["case_id"], "replicate": case["replicate"],
              "mode": case["mode"], "cross_hgt_probability": case["cross_hgt_probability"],
              "n_candidate_pairs": data["n_pairs"], "n_finite_kovar": data["n_finite"],
              "focal_status": status}
    focal_row = {**common, "mi": data["mi"][focal] if focal is not None else "NA",
                 "p_primary": data["p"][focal] if status == "available" else "NA",
                 "full_mi_rank": data["mi_rank"][focal] if focal is not None else "NA",
                 "shared_mi_rank": data["shared_mi_rank"][focal] if status == "available" else "NA",
                 "kovar_rank": data["p_rank"][focal] if status == "available" else "NA",
                 "mi_top_recovered": int(mi_top[focal]) if focal is not None else 0,
                 "kovar_top_recovered": int(p_top[focal]) if focal is not None else 0,
                 "kovar_bonferroni_recovered": int(data["p"][focal] <= args.alpha / data["n_pairs"])
                    if status == "available" else 0}
    lineage = (data["category"] == CATEGORIES["lineage_driven"]) & (data["distance"] > args.distal_bp)
    universe = finite & (data["distance"] > args.distal_bp) & (data["category"] != CATEGORIES["unclassified"])
    # Panel D applies its distal filter before imposing the discovery budget.
    marginal_top = top_mask(np.where(universe, data["mi"], np.nan), args.top_fraction, descending=True)
    adjusted_top = top_mask(np.where(universe, data["p"], np.nan), args.top_fraction, descending=False)
    base_rate = (lineage & universe).sum() / universe.sum() if universe.any() else math.nan
    lineage_row = {**common, "distal_pairs": int(universe.sum()),
                   "distal_lineage_pairs": int((lineage & universe).sum()),
                   "unclassified_pairs": int((data["category"] == CATEGORIES["unclassified"]).sum())}
    for method, mask in (("mi", marginal_top), ("kovar", adjusted_top)):
        n_top = int(mask.sum())
        n_lineage = int((mask & lineage).sum())
        enrichment = n_lineage / n_top / base_rate if n_top and base_rate > 0 else math.nan
        lineage_row.update({f"{method}_top_pairs": n_top, f"{method}_top_lineage": n_lineage,
                            f"{method}_lineage_enrichment": enrichment})
    return focal_row, lineage_row


def condition(row: dict) -> tuple[int, float]:
    return int(row["mode"]), float(row["cross_hgt_probability"])


def panel_c(rows: list[dict], output: Path, args) -> list[str]:
    import matplotlib.pyplot as plt
    groups = sorted({condition(row) for row in rows})
    labels = [f"M{mode}\nHGT={hgt:g}" for mode, hgt in groups]
    files, recovery = [], []
    fig, ax = plt.subplots(figsize=(max(6, len(groups) * .85), 4.5), layout="constrained")
    for j, (method, field) in enumerate((("MI", "shared_mi_rank"), ("KOVAR", "kovar_rank"))):
        for x, group in enumerate(groups):
            selected = [row for row in rows if condition(row) == group]
            values = [float(row[field]) / int(row["n_finite_kovar"]) * 100
                      for row in selected if row[field] != "NA" and int(row["n_finite_kovar"]) > 0]
            jitter = np.linspace(-.045, .045, len(values)) if len(values) > 1 else np.zeros(len(values))
            ax.scatter(x + (j - .5) * .22 + jitter, values, s=32,
                       color=COLORS[method], alpha=.8, label=method if x == 0 else None)
            recovered = sum(int(row[f"{method.lower()}_top_recovered"]) for row in selected)
            recovery.append({"mode": group[0], "cross_hgt_probability": group[1],
                             "method": method, "n_cases": len(selected),
                             "n_focal_available": sum(row["focal_status"] == "available" for row in selected),
                             "n_top_recovered": recovered, "recovery_fraction": recovered / len(selected)})
    ax.set(xticks=range(len(groups)), xticklabels=labels, ylabel="A–B rank / shared finite pairs (%)",
           title="C  Focal A–B rank across simulation conditions", ylim=(0, 105),
           xlim=(-.5, len(groups) - .5))
    ax.grid(axis="y", alpha=.2)
    ax.legend(frameon=False)
    files += save(fig, output, "C_focal_rank", args.formats)
    fig, ax = plt.subplots(figsize=(max(6, len(groups) * .85), 4.5), layout="constrained")
    for method in COLORS:
        subset = [row for row in recovery if row["method"] == method]
        ax.plot(range(len(groups)), [row["recovery_fraction"] for row in subset],
                marker="o", color=COLORS[method], linewidth=1, label=method)
    ax.set(xticks=range(len(groups)), xticklabels=labels, ylim=(-.03, 1.05),
           ylabel=f"Fraction of cases with A–B in top {100 * args.top_fraction:g}%",
           title="C  Focal A–B recovery", xlim=(-.5, len(groups) - .5))
    ax.text(.02, .5, "Denominator: all selected cases", transform=ax.transAxes, fontsize=8)
    ax.legend(frameon=False, loc="upper left")
    ax.grid(axis="y", alpha=.2)
    files += save(fig, output, "C_focal_recovery", args.formats)
    write_rows(output / "C_recovery.tsv", recovery)
    return files + ["C_recovery.tsv"]


def panel_d(rows: list[dict], output: Path, args) -> list[str]:
    import matplotlib.pyplot as plt
    groups = sorted({condition(row) for row in rows})
    fig, ax = plt.subplots(figsize=(max(6, len(groups) * .85), 4.5), layout="constrained")
    for j, method in enumerate(COLORS):
        for x, group in enumerate(groups):
            subset = [row for row in rows if condition(row) == group]
            values = [float(row[f"{method.lower()}_lineage_enrichment"]) for row in subset]
            values = [value for value in values if math.isfinite(value)]
            jitter = np.linspace(-.045, .045, len(values)) if len(values) > 1 else np.zeros(len(values))
            ax.scatter(x + (j - .5) * .22 + jitter, values, color=COLORS[method],
                       s=32, alpha=.8, label=method if x == 0 else None)
    ax.axhline(1, color="black", linestyle="--", linewidth=.8)
    ax.set(xticks=range(len(groups)),
           xticklabels=[f"M{m}\nHGT={h:g}" for m, h in groups],
           ylabel="Lineage-associated enrichment in top-ranked distal pairs",
           title=f"D  Nonfocal lineage-associated pairs (>{args.distal_bp / 1000:g} kb)",
           ylim=(0, None), xlim=(-.5, len(groups) - .5))
    if not any(math.isfinite(float(row[f"{method.lower()}_lineage_enrichment"]))
               for row in rows for method in COLORS):
        ax.text(.5, .5, "No finite enrichment estimates\nSee D_lineage_cases.tsv for pair counts",
                transform=ax.transAxes, ha="center", fontsize=9)
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=.2)
    return save(fig, output, "D_lineage_enrichment", args.formats)


def keyed_table(path: Path, pvalues: bool) -> dict:
    result = {}
    required = {"case_id", "u", "v"} | ({"p_marginal"} if pvalues else set())
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: expected columns {sorted(required)}")
        for row in reader:
            a, b = sorted((int(row["u"]), int(row["v"])))
            key = row["case_id"], a, b
            if a < 0 or a == b or key in result:
                raise ValueError(f"{path}: invalid or duplicate pair {key}")
            value = valid_p(row["p_marginal"]) if pvalues else True
            if pvalues and not math.isfinite(value):
                raise ValueError(f"{path}: invalid calibrated marginal P for {key}")
            result[key] = value
    if not result:
        raise ValueError(f"{path}: empty input")
    return result


def null_values(data: dict, null_pairs: dict, marginal_p: dict) -> tuple[list[float], list[float]]:
    values, adjusted = [], []
    case_id = data["case"]["case_id"]
    n = len(data["positions"])
    for key in null_pairs:
        if key[0] != case_id:
            continue
        index = triangular_index(key[1], key[2], n)
        if int(data["case"]["mode"]) == 2 and index == data["focal"]:
            raise ValueError("the implanted mode-2 A–B pair cannot be a calibration null")
        if key not in marginal_p:
            raise ValueError(f"missing independently calibrated marginal P: {key}")
        if math.isfinite(data["p"][index]):
            values.append(marginal_p[key])
            adjusted.append(data["p"][index])
    return values, adjusted


def panel_a(groups: dict, output: Path, args) -> list[str]:
    import matplotlib.pyplot as plt
    files, rows = [], []
    for (mode, hgt), (marginal, adjusted) in sorted(groups.items()):
        if not marginal:
            continue
        fig, ax = plt.subplots(figsize=(5.2, 4.8), layout="constrained")
        maximum = 1.
        for method, values in (("MI", marginal), ("KOVAR", adjusted)):
            ordered = np.sort(values)
            expected = (np.arange(len(ordered)) + .5) / len(ordered)
            x, y = -np.log10(expected), -np.log10(np.maximum(ordered, P_FLOOR))
            ax.plot(x, y, color=COLORS[method], linewidth=1.1, marker=".", markersize=3,
                    label=f"{'Calibrated marginal' if method == 'MI' else 'KOVAR'} (n={len(values):,})")
            maximum = max(maximum, float(x.max()), float(y.max()))
            rows.extend({"mode": mode, "cross_hgt_probability": hgt, "method": method,
                         "expected_p": e, "observed_p": p} for e, p in zip(expected, ordered))
        ax.plot([0, maximum], [0, maximum], "--", color="black", linewidth=.8)
        ax.set(xlabel="Expected −log10(P)", ylabel="Observed −log10(P)",
               title=f"A  Predeclared null pairs · mode {mode} · cross-HGT {hgt:g}")
        ax.legend(frameon=False, fontsize=8)
        slug = f"mode_{mode}_hgt_{hgt:g}".replace(".", "p")
        files += save(fig, output, f"A_qq_{slug}", args.formats)
    if not rows:
        raise ValueError("no null pairs have finite paired marginal/KOVAR P values")
    write_rows(output / "A_qq.tsv", rows)
    return files + ["A_qq.tsv"]


def export_pairs(data: dict, output: Path) -> str:
    path = output / f"pairs_{data['case']['case_id']}.tsv.gz"
    fields = ["u", "v", "u_position", "v_position", "distance_bp", "category", "mi",
              "p_primary", "mi_rank_full", "mi_rank_shared_finite", "kovar_rank", "kovar_status"]
    labels = {value: key for key, value in CATEGORIES.items()}
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for i, (u, v) in enumerate(zip(data["u"], data["v"])):
            writer.writerow([u, v, data["positions"][u], data["positions"][v],
                             data["distance"][i], labels[data["category"][i]], data["mi"][i],
                             *[data[key][i] if math.isfinite(data[key][i]) else "NA"
                               for key in ("p", "mi_rank", "shared_mi_rank", "p_rank")],
                             "finite" if math.isfinite(data["p"][i]) else "unavailable"])
    return path.name


def recovery_by_budget(data: dict, args) -> list[dict]:
    """Focal recovery at equal distal-candidate budgets, including cutoff ties."""
    distal = data["distance"] > args.distal_bp
    n_candidates = int(distal.sum())
    focal = data["focal"]
    status = ("not_maf_eligible" if focal is None else "not_distal" if not distal[focal]
              else "available" if math.isfinite(data["p"][focal]) else "kovar_p_unavailable")
    rows = []
    for percent in args.budget_percent:
        budget = max(1, math.ceil(percent / 100 * n_candidates)) if n_candidates else 0
        for method, values, descending in (("MI", data["mi"], True), ("KOVAR", data["p"], False)):
            selected = top_mask(np.where(distal, values, np.nan), percent / 100,
                                descending=descending, budget=budget)
            rows.append({"case_id": data["case"]["case_id"],
                         "replicate": data["case"]["replicate"], "mode": data["case"]["mode"],
                         "cross_hgt_probability": data["case"]["cross_hgt_probability"],
                         "method": method, "budget_percent": percent, "nominal_budget": budget,
                         "n_distal_candidates": n_candidates,
                         "n_finite_distal_kovar": int((distal & np.isfinite(data["p"])).sum()),
                         "n_selected_including_ties": int(selected.sum()), "focal_status": status,
                         "recovered": int(selected[focal]) if focal is not None else 0})
    return rows


def panel_recovery(rows: list[dict], output: Path, args) -> list[str]:
    import matplotlib.pyplot as plt
    hgt_values = sorted({float(row["cross_hgt_probability"]) for row in rows})
    fig, axes = plt.subplots(1, len(hgt_values), figsize=(3.1 * len(hgt_values), 3.7),
                             sharey=True, squeeze=False, layout="constrained")
    summaries = []
    for ax, hgt in zip(axes[0], hgt_values):
        subset = [row for row in rows if float(row["cross_hgt_probability"]) == hgt]
        for method in COLORS:
            fractions = []
            for budget in args.budget_percent:
                group = [row for row in subset if row["method"] == method and row["budget_percent"] == budget]
                recovered = sum(row["recovered"] for row in group)
                fraction = recovered / len(group)
                fractions.append(fraction)
                summaries.append({"mode": args.main_mode, "cross_hgt_probability": hgt,
                                  "method": method, "budget_percent": budget, "n_replicates": len(group),
                                  "n_recovered": recovered, "recovery_fraction": fraction})
            ax.plot(args.budget_percent, fractions, color=COLORS[method],
                    marker="o" if method == "MI" else "s", markersize=4, linewidth=1.2, label=method)
        n_replicates = len({row["replicate"] for row in subset})
        ax.set(xscale="log", xlabel="Distal-pair discovery budget (%)", ylim=(-.04, 1.04),
               title=f"Cross-HGT = {hgt:g} · n = {n_replicates}",
               xticks=args.budget_percent, xticklabels=[f"{x:g}" for x in args.budget_percent])
        ax.tick_params(axis="x", labelrotation=45)
        ax.grid(axis="y", alpha=.2)
    axes[0][0].set_ylabel("Fraction of replicates recovering A–B")
    axes[0][0].legend(frameon=False)
    fig.suptitle(f"D  Focal A–B recovery · mode {args.main_mode}")
    write_rows(output / "D_recovery_by_budget.tsv", summaries)
    write_rows(output / "D_recovery_per_replicate.tsv", rows)
    return save(fig, output, "D_focal_recovery_by_budget", args.formats) + [
        "D_recovery_by_budget.tsv", "D_recovery_per_replicate.tsv"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default="manifests/cases.tsv")
    parser.add_argument("--case-id", action="append", help="exact case ID; repeat to select cases")
    parser.add_argument("--spa-mode", choices=("off", "auto"), required=True)
    parser.add_argument("--sample-reweighting", choices=("default", "none"), default="default")
    parser.add_argument("--panels", nargs="+", choices=list("ABCD"), default=list("ABCD"))
    parser.add_argument("--main-mode", type=int, choices=(0, 1, 2), default=2)
    parser.add_argument("--example-replicate", type=int, default=1)
    parser.add_argument("--example-hgt", type=float, default=.002)
    parser.add_argument("--supplementary", action="store_true", help="also export per-case plots and control summaries for all manifest cases")
    parser.add_argument("--budget-percent", nargs="+", type=float, default=[.01, .05, .1, .5, 1, 2, 5])
    parser.add_argument("--output-dir")
    parser.add_argument("--formats", nargs="+", choices=("pdf", "svg", "png"), default=["png", "svg"])
    parser.add_argument("--distal-bp", type=int, default=10_000)
    parser.add_argument("--max-points", type=int, default=100_000)
    parser.add_argument("--top-fraction", type=float, default=.01)
    parser.add_argument("--alpha", type=float, default=.05)
    parser.add_argument("--null-pairs", help="TSV with case_id,u,v; independently defined null membership")
    parser.add_argument("--marginal-pvalues", help="TSV with case_id,u,v,p_marginal; already calibrated P values")
    parser.add_argument("--calibration-note", help="null definition and marginal calibration provenance")
    args = parser.parse_args()
    if any(not math.isfinite(x) or not 0 < x <= 100 for x in args.budget_percent):
        parser.error("budget-percent values must be finite and in (0, 100]")
    args.budget_percent = sorted(set(args.budget_percent))
    if args.distal_bp < 0 or args.max_points < 1 or not 0 < args.top_fraction <= 1 or not 0 < args.alpha < 1:
        parser.error("require distal-bp >= 0, max-points > 0, 0 < top-fraction <= 1 and 0 < alpha < 1")
    calibration_ready = bool(args.null_pairs and args.marginal_pvalues and args.calibration_note)
    if any((args.null_pairs, args.marginal_pvalues, args.calibration_note)) and not calibration_ready:
        parser.error("supply null-pairs, marginal-pvalues and calibration-note together")
    panels = set(args.panels)
    if panels == {"A"} and not calibration_ready:
        parser.error("panel A requires null-pairs, marginal-pvalues and calibration-note")
    cases = read_tsv(repo_path(args.manifest))
    ids = [row["case_id"] for row in cases]
    if len(set(ids)) != len(ids) or any(not re.fullmatch(r"[A-Za-z0-9_.-]+", x) for x in ids):
        parser.error("manifest case IDs must be unique safe filenames")
    if args.case_id:
        missing = set(args.case_id) - set(ids)
        if missing:
            parser.error(f"unknown case IDs: {sorted(missing)}")
        cases = [row for row in cases if row["case_id"] in args.case_id]
    if not cases:
        parser.error("no cases selected")
    main_cases = [row for row in cases if int(row["mode"]) == args.main_mode]
    if not main_cases:
        parser.error("no cases for main-mode in the selected manifest")
    if len({(int(row["replicate"]), float(row["cross_hgt_probability"])) for row in main_cases}) != len(main_cases):
        parser.error("duplicate replicate/HGT conditions in main-mode")
    examples = [row for row in main_cases if int(row["replicate"]) == args.example_replicate
                and math.isclose(float(row["cross_hgt_probability"]), args.example_hgt, rel_tol=0, abs_tol=1e-12)]
    if panels & {"A", "B", "C"} and len(examples) != 1:
        parser.error("the predefined example is absent or ambiguous; set example-replicate/example-hgt explicitly")
    example_id = examples[0]["case_id"] if examples else None
    processing_cases = cases if args.supplementary else main_cases
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans",
                               "pdf.fonttype": 42, "svg.fonttype": "none"})
    output = repo_path(args.output_dir or f"results/figure1_redesigned_spa_{args.spa_mode}_{args.sample_reweighting}")
    output.mkdir(parents=True, exist_ok=True)
    null_pairs = keyed_table(repo_path(args.null_pairs), False) if calibration_ready else {}
    marginal_p = keyed_table(repo_path(args.marginal_pvalues), True) if calibration_ready else {}
    manifest_ids = set(ids)
    if any(key[0] not in manifest_ids for key in null_pairs.keys() | marginal_p.keys()):
        parser.error("calibration tables contain a case not present in the manifest")
    focal_rows, lineage_rows, recovery_rows, files, input_hashes = [], [], [], [], {}
    skipped = {}
    for case in processing_cases:
        print(f"[figure1] {case['case_id']}", flush=True)
        data = read_case(case, args.spa_mode, args.sample_reweighting)
        input_hashes.update(data["input_sha256"])
        focal, lineage = summarize(data, args)
        focal_rows.append(focal)
        lineage_rows.append(lineage)
        if int(case["mode"]) == args.main_mode and "D" in panels:
            recovery_rows.extend(recovery_by_budget(data, args))
        if case["case_id"] == example_id:
            files.append(export_pairs(data, output))
            if "B" in panels:
                files += panel_b(data, output, args)
            if "C" in panels:
                files += panel_e(data, output, args, panel="C")
            if "A" in panels:
                if not calibration_ready:
                    skipped["A"] = "Needs independent null membership and calibrated marginal P values"
                else:
                    marginal, adjusted = null_values(data, null_pairs, marginal_p)
                    if not marginal:
                        raise ValueError("the predefined example has no finite paired null P values")
                    files += panel_a({condition(case): (marginal, adjusted)}, output, args)
        if args.supplementary:
            supplementary = output / "supplementary" / case["case_id"]
            supplementary.mkdir(parents=True, exist_ok=True)
            generated = [export_pairs(data, supplementary)]
            generated += panel_b(data, supplementary, args)
            generated += panel_e(data, supplementary, args)
            if calibration_ready:
                marginal, adjusted = null_values(data, null_pairs, marginal_p)
                if marginal:
                    generated += panel_a({condition(case): (marginal, adjusted)}, supplementary, args)
                else:
                    skipped[f"supplementary/{case['case_id']}/A"] = "No paired finite null P values supplied"
            else:
                skipped[f"supplementary/{case['case_id']}/A"] = "No independent null/calibration inputs supplied"
            files.extend(str((supplementary / name).relative_to(output)) for name in generated)
    write_rows(output / "focal_cases.tsv", focal_rows)
    write_rows(output / "lineage_cases.tsv", lineage_rows)
    files += ["focal_cases.tsv", "lineage_cases.tsv"]
    if recovery_rows:
        files += panel_recovery(recovery_rows, output, args)
    if args.supplementary:
        supplementary = output / "supplementary"
        generated = panel_c(focal_rows, supplementary, args) + panel_d(lineage_rows, supplementary, args)
        files.extend(str((supplementary / name).relative_to(output)) for name in generated)
    for path in [repo_path(args.manifest)] + ([repo_path(args.null_pairs), repo_path(args.marginal_pvalues)] if calibration_ready else []):
        input_hashes[str(path)] = sha256_file(path)
    report = {"settings": vars(args), "cases": [row["case_id"] for row in processing_cases],
              "example_case_id": example_id, "main_cases": [row["case_id"] for row in main_cases],
              "input_sha256": input_hashes,
              "files": files, "skipped_panels": skipped,
              "rank_rule": "one-based average ties; unavailable tests have no rank; C reranks both methods on identical finite-test pairs",
              "lineage_rule": "existing covariance-component evaluator; >distal-bp filter before supplementary enrichment selection",
              "display_only_thinning": True,
              "p_plotting_floor": P_FLOOR,
              "inference": "D summarizes focal recovery across replicates, not conventional precision; no assumed outcome"}
    (output / "figure1_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for panel, reason in skipped.items():
        print(f"[pending] panel {panel}: {reason}")
    print(f"[done] Figure 1 panels: {output}")


if __name__ == "__main__":
    main()
