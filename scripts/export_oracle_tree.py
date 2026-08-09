#!/usr/bin/env python3
"""Export one labelled local true-ancestry tree from a SLiM tree sequence.

This is deliberately an oracle diagnostic, not the realistic primary tree.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pyslim
import tskit


def sample_labels(path: Path) -> dict[int, str]:
    with open(path, encoding="utf-8", newline="") as handle:
        return {
            int(row["pedigree_id"]): row["analysis_label"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def export(
    trees_path: Path,
    sample_map_path: Path,
    output_path: Path,
    position: float,
    ancestral_ne: float,
    seed: int,
) -> int:
    labels_by_pedigree = sample_labels(sample_map_path)
    ts = tskit.load(trees_path)
    ts = pyslim.recapitate(ts, ancestral_Ne=ancestral_ne, random_seed=seed)

    old_nodes: list[int] = []
    old_labels: dict[int, str] = {}
    for individual in ts.individuals():
        metadata = individual.metadata or {}
        pedigree_id = metadata.get("pedigree_id")
        if pedigree_id not in labels_by_pedigree:
            continue
        non_null_nodes = [node for node in individual.nodes if node != tskit.NULL]
        if len(non_null_nodes) != 1:
            raise ValueError(
                f"Expected one haploid node for pedigree {pedigree_id}, found {len(non_null_nodes)}"
            )
        node = non_null_nodes[0]
        old_nodes.append(node)
        old_labels[node] = labels_by_pedigree[pedigree_id]

    if len(old_nodes) != len(labels_by_pedigree):
        raise ValueError(
            f"Mapped {len(old_nodes)} tree individuals but sample map contains {len(labels_by_pedigree)}"
        )

    simplified, node_map = ts.simplify(
        samples=old_nodes,
        keep_input_roots=True,
        map_nodes=True,
    )
    node_labels = {int(node_map[old]): label for old, label in old_labels.items()}
    focal_position = min(max(0.0, position), simplified.sequence_length - 1e-9)
    tree = simplified.at(focal_position)
    if tree.num_roots != 1:
        raise ValueError(f"Oracle tree has {tree.num_roots} roots after recapitation")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(tree.as_newick(node_labels=node_labels) + "\n", encoding="utf-8")
    return len(node_labels)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trees", required=True, type=Path)
    parser.add_argument("--sample-map", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--position", required=True, type=float)
    parser.add_argument("--ancestral-ne", required=True, type=float)
    parser.add_argument("--seed", required=True, type=int)
    args = parser.parse_args()
    count = export(
        args.trees,
        args.sample_map,
        args.output,
        args.position,
        args.ancestral_ne,
        args.seed,
    )
    print(f"[done] wrote oracle local tree with {count} tips to {args.output}")


if __name__ == "__main__":
    main()
