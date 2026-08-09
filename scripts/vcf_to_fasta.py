#!/usr/bin/env python3
"""Convert haploid SLiM VCF SNPs to a compact aligned FASTA."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_sample_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    with open(path, encoding="utf-8", newline="") as handle:
        return {
            row["vcf_label"]: row["analysis_label"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def load_excluded_positions(path: Path | None) -> set[int]:
    if path is None:
        return set()
    with open(path, encoding="utf-8", newline="") as handle:
        # SLiM positions are zero-based; VCF positions are one-based.
        return {int(row["position"]) + 1 for row in csv.DictReader(handle, delimiter="\t")}


def genotype_base(sample_field: str, alleles: list[str]) -> str:
    genotype = sample_field.split(":", 1)[0]
    if genotype in {".", "./.", ".|."}:
        return "N"
    token = genotype.replace("|", "/").split("/", 1)[0]
    try:
        allele = alleles[int(token)]
    except (ValueError, IndexError):
        return "N"
    return allele.upper() if len(allele) == 1 and allele.upper() in "ACGT" else "N"


def convert(
    vcf: Path,
    output: Path,
    positions_output: Path,
    sample_map: Path | None = None,
    excluded_loci: Path | None = None,
) -> tuple[int, int]:
    mapping = load_sample_map(sample_map)
    excluded = load_excluded_positions(excluded_loci)
    sample_names: list[str] | None = None
    analysis_names: list[str] | None = None
    sequences: list[list[str]] = []
    positions: list[tuple[str, int, str, str]] = []

    with open(vcf, encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("##"):
                continue
            if raw.startswith("#CHROM"):
                fields = raw.rstrip("\n").split("\t")
                sample_names = fields[9:]
                analysis_names = [mapping.get(name, name.replace(":", "_")) for name in sample_names]
                if len(set(analysis_names)) != len(analysis_names):
                    raise ValueError("Analysis sample labels are not unique")
                sequences = [[] for _ in sample_names]
                continue
            if raw.startswith("#") or not raw.strip():
                continue
            if sample_names is None:
                raise ValueError("VCF has data before the #CHROM header")

            fields = raw.rstrip("\n").split("\t")
            chrom, pos_text, variant_id, ref, alt = fields[:5]
            pos = int(pos_text)
            if pos in excluded:
                continue
            alleles = [ref, *alt.split(",")]
            if any(len(allele) != 1 or allele.upper() not in "ACGT" for allele in alleles):
                continue

            for index, sample_field in enumerate(fields[9:]):
                sequences[index].append(genotype_base(sample_field, alleles))
            positions.append((chrom, pos, variant_id, ",".join(alleles)))

    if sample_names is None or analysis_names is None:
        raise ValueError("VCF has no #CHROM header")
    if not positions:
        raise ValueError("No biallelic or multiallelic SNP positions remain after filtering")

    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8", newline="\n") as handle:
        for name, sequence in zip(analysis_names, sequences, strict=True):
            joined = "".join(sequence)
            handle.write(f">{name}\n")
            for offset in range(0, len(joined), 80):
                handle.write(joined[offset : offset + 80] + "\n")

    with open(positions_output, "w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["alignment_column", "chrom", "vcf_position", "variant_id", "alleles"])
        for column, position in enumerate(positions):
            writer.writerow([column, *position])

    return len(sample_names), len(positions)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--positions-output", required=True, type=Path)
    parser.add_argument("--sample-map", type=Path)
    parser.add_argument("--exclude-loci", type=Path)
    args = parser.parse_args()

    samples, sites = convert(
        args.vcf,
        args.output,
        args.positions_output,
        args.sample_map,
        args.exclude_loci,
    )
    print(f"[done] wrote {samples} samples and {sites} SNP columns to {args.output}")


if __name__ == "__main__":
    main()
