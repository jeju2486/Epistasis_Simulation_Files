#!/usr/bin/env python3
"""Create the deterministic ancestral reference defined by a TOML config."""

from __future__ import annotations

import argparse
import random
from pathlib import Path

from simflow import load_config, repo_path


def write_reference(path: Path, length: int, seed: int, force: bool = False) -> None:
    if path.exists() and not force:
        sequence = "".join(
            line.strip() for line in path.read_text(encoding="utf-8").splitlines() if not line.startswith(">")
        )
        if len(sequence) != length:
            raise SystemExit(
                f"Existing reference has length {len(sequence)}, expected {length}; use --force to replace it."
            )
        print(f"[skip] valid reference already exists: {path}")
        return

    rng = random.Random(seed)
    sequence = "".join(rng.choices("ACGT", k=length))
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with open(temporary, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(">simulated_reference\n")
        for offset in range(0, length, 80):
            handle.write(sequence[offset : offset + 80] + "\n")
    temporary.replace(path)
    print(f"[done] wrote {length}-bp reference: {path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    write_reference(
        repo_path(config["paths"]["reference"]),
        int(config["genome"]["length"]),
        int(config["genome"]["reference_seed"]),
        args.force,
    )


if __name__ == "__main__":
    main()
