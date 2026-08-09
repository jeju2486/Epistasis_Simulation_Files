"""Shared helpers for the manifest-driven simulation workflow."""

from __future__ import annotations

import csv
import hashlib
import os
import tempfile
import tomllib
from pathlib import Path
from typing import Iterable, Mapping


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_config(path: str | os.PathLike[str]) -> dict:
    with open(path, "rb") as handle:
        return tomllib.load(handle)


def repo_path(value: str | os.PathLike[str]) -> Path:
    path = Path(value)
    return path if path.is_absolute() else REPO_ROOT / path


def deterministic_seed(master_seed: int, *parts: object) -> int:
    material = "|".join([str(master_seed), *(str(part) for part in parts)])
    digest = hashlib.sha256(material.encode("utf-8")).digest()
    seed = int.from_bytes(digest[:8], "big") & ((1 << 62) - 1)
    return seed or 1


def probability_slug(value: float) -> str:
    if value == 0:
        return "0"
    return format(value, ".12g").replace("-", "m").replace(".", "p").replace("+", "")


def write_tsv_atomic(path: Path, rows: Iterable[Mapping[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    temporary.replace(path)


def read_tsv(path: str | os.PathLike[str]) -> list[dict[str, str]]:
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))
