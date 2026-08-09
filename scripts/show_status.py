#!/usr/bin/env python3
"""Summarize manifest completion from atomic _SUCCESS markers."""

from __future__ import annotations

import argparse
from collections import Counter

from simflow import read_tsv, repo_path


def status(out_dir: str) -> str:
    directory = repo_path(out_dir)
    if (directory / "_SUCCESS").exists():
        return "complete"
    if directory.exists():
        return "incomplete"
    return "pending"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifests", nargs="+")
    args = parser.parse_args()

    for manifest in args.manifests:
        rows = read_tsv(manifest)
        counts = Counter(status(row["out_dir"]) for row in rows)
        print(
            f"{manifest}: total={len(rows)} complete={counts['complete']} "
            f"incomplete={counts['incomplete']} pending={counts['pending']}"
        )


if __name__ == "__main__":
    main()
