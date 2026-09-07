#!/usr/bin/env python3
"""Infer and midpoint-root one case tree with IQ-TREE 2."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import tempfile
import time
from pathlib import Path

from kovar_inputs import read_fasta
from simflow import read_tsv


RESULT_NAME = "iqtree_phylogeny"
ROOTED_TREE_NAME = "iqtree_rooted.nwk"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def archive(path: Path, label: str) -> Path:
    destination = path.with_name(
        f"{path.name}.{label}.{time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())}"
    )
    path.replace(destination)
    return destination


def write_variable_alignment(
    alignment: Path,
    positions: Path,
    selected_loci: Path,
    output: Path,
    site_map: Path,
) -> tuple[int, int, int]:
    """Remove focal and invariant columns before applying IQ-TREE's +ASC model."""
    records = read_fasta(alignment)
    input_sites = len(records[0][1])
    if any(len(sequence) != input_sites for _name, sequence in records):
        raise ValueError("IQ-TREE input sequences have unequal lengths")
    invalid = sorted(set().union(*(set(sequence) for _name, sequence in records)) - set("ACGT"))
    if invalid:
        raise ValueError(f"IQ-TREE input contains unsupported states: {invalid}")

    with positions.open(encoding="utf-8", newline="") as handle:
        position_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(position_rows) != input_sites:
        raise ValueError("alignment and position map have different site counts")
    columns = [int(row["alignment_column"]) for row in position_rows]
    if columns != list(range(input_sites)):
        raise ValueError("alignment columns must be contiguous and zero-based")

    focal_positions = {int(row["position"]) for row in read_tsv(selected_loci)}
    retained: list[int] = []
    for column, row in enumerate(position_rows):
        slim_position = int(row["vcf_position"]) - 1
        states = {sequence[column] for _name, sequence in records}
        if slim_position not in focal_positions and len(states) > 1:
            retained.append(column)
    if not retained:
        raise ValueError("no nonfocal variable sites remain for IQ-TREE")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        for name, sequence in records:
            filtered = "".join(sequence[column] for column in retained)
            handle.write(f">{name}\n{filtered}\n")

    fields = ["iqtree_column", "original_column", "slim_position"]
    with site_map.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for iqtree_column, original_column in enumerate(retained):
            writer.writerow({
                "iqtree_column": iqtree_column,
                "original_column": original_column,
                "slim_position": int(position_rows[original_column]["vcf_position"]) - 1,
            })
    return len(records), input_sites, len(retained)


def midpoint_root(treefile: Path, output: Path, expected_names: list[str]) -> None:
    # Import lazily so alignment-preparation tests do not require Biopython.
    from Bio import Phylo

    tree = Phylo.read(treefile, "newick")
    tree.root_at_midpoint()
    observed = [tip.name for tip in tree.get_terminals()]
    if len(observed) != len(set(observed)):
        raise ValueError("IQ-TREE output contains duplicate tip labels")
    if set(observed) != set(expected_names):
        missing = sorted(set(expected_names) - set(observed))[:5]
        extra = sorted(set(observed) - set(expected_names))[:5]
        raise ValueError(f"IQ-TREE tip mismatch; missing={missing}, extra={extra}")
    Phylo.write(tree, output, "newick")


def run_case(
    case_dir: Path,
    executable: str,
    threads: int,
    seed: int,
    model: str,
    force: bool,
) -> None:
    case_dir = case_dir.resolve()
    alignment = case_dir / "all_snps.fa"
    positions = case_dir / "all_snps.positions.tsv"
    selected = case_dir / "selected_loci.tsv"
    for path in (case_dir / "_SUCCESS", alignment, positions, selected):
        if not path.exists() or (path.is_file() and path.name != "_SUCCESS" and path.stat().st_size == 0):
            raise FileNotFoundError(f"missing IQ-TREE input: {path}")

    output = case_dir / RESULT_NAME
    requested = {
        "model": model,
        "seed": seed,
        "input_sha256": {
            path.name: sha256_file(path) for path in (alignment, positions, selected)
        },
    }
    if (output / "_SUCCESS").exists() and not force:
        metadata_path = output / "run_metadata.json"
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            metadata = None
        observed = (
            {key: metadata.get(key) for key in requested}
            if isinstance(metadata, dict) else None
        )
        if observed == requested:
            print(f"[skip] {output}")
            return
        print(f"[archive] {archive(output, 'settings_mismatch')}")
    elif output.exists():
        print(f"[archive] {archive(output, 'replaced' if force else 'incomplete')}")

    temporary = Path(tempfile.mkdtemp(prefix=RESULT_NAME + ".tmp.", dir=case_dir))
    try:
        variable_alignment = temporary / "nonfocal_variable_snps.fa"
        site_map = temporary / "nonfocal_variable_sites.tsv"
        samples, input_sites, retained_sites = write_variable_alignment(
            alignment, positions, selected, variable_alignment, site_map
        )
        sample_names = [name for name, _sequence in read_fasta(variable_alignment)]
        prefix = temporary / "iqtree"
        command = [
            executable, "-s", str(variable_alignment), "-st", "DNA", "-m", model,
            "-nt", str(threads), "-seed", str(seed), "-keep-ident",
            "-pre", str(prefix), "-quiet",
        ]
        with (temporary / "iqtree.stdout.log").open("w", encoding="utf-8") as log:
            completed = subprocess.run(
                command, cwd=temporary, stdout=log, stderr=subprocess.STDOUT,
                text=True, check=False,
            )
        if completed.returncode:
            raise RuntimeError(f"IQ-TREE exited {completed.returncode}")
        treefile = prefix.with_suffix(".treefile")
        if not treefile.exists() or treefile.stat().st_size == 0:
            raise FileNotFoundError(f"IQ-TREE did not produce {treefile}")
        rooted = temporary / ROOTED_TREE_NAME
        midpoint_root(treefile, rooted, sample_names)
        metadata = {
            "case_dir": str(case_dir), "command": command, "model": model,
            "rooting": "midpoint", "samples": samples, "input_sites": input_sites,
            "nonfocal_variable_sites": retained_sites, "seed": seed,
            "input_sha256": requested["input_sha256"],
        }
        (temporary / "run_metadata.json").write_text(
            json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
        )
        (temporary / "_SUCCESS").write_text("complete\n", encoding="utf-8")
        temporary.replace(output)
    except Exception:
        failed = output.with_name(
            f"{output.name}.failed.{time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())}"
        )
        temporary.replace(failed)
        raise
    print(f"[done] IQ-TREE phylogeny: {output}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", required=True, type=Path)
    parser.add_argument("--iqtree", default=os.environ.get("IQTREE_BIN", "iqtree"))
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--model", default="GTR+ASC")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    if args.seed < 1:
        parser.error("--seed must be positive")
    run_case(args.case_dir, args.iqtree, args.threads, args.seed, args.model, args.force)


if __name__ == "__main__":
    main()
