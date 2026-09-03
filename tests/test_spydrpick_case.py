from __future__ import annotations

import gzip
import sys
import tempfile
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from spydrpick_case import normalize_edges, output_name  # noqa: E402


class SpydrPickSummaryTests(unittest.TestCase):
    def test_weighting_modes_have_distinct_output_names(self) -> None:
        self.assertEqual(output_name("default"), "spydrpick_all_pairs")
        self.assertEqual(output_name("none"), "spydrpick_all_pairs_unweighted")
        with self.assertRaisesRegex(ValueError, "unknown sample-reweighting mode"):
            output_name("invalid")

    def test_all_pairs_are_zero_based_and_retain_physical_distance(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "raw.edges"
            raw.write_text("1 2 1 0 0.9\n1 3 2 0 0.8\n2 3 1 0 0.7\n", encoding="utf-8")
            output = root / "edges.gz"
            self.assertEqual(normalize_edges(raw, output, [10, 20, 40]), (3, 0))
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                rows = handle.readlines()
            self.assertEqual(rows[0].split()[:3], ["0", "1", "10"])
            self.assertEqual(rows[1].split()[:3], ["0", "2", "30"])

    def test_missing_pairs_are_restored_as_exact_zero_mi(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "raw.edges"
            raw.write_text("1 2 1 0 0.9\n", encoding="utf-8")
            output = root / "edges.gz"
            self.assertEqual(normalize_edges(raw, output, [10, 20, 30]), (3, 2))
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                rows = [line.split() for line in handle]
            self.assertEqual(rows, [
                ["0", "1", "10", "0", "0.90000000000000002"],
                ["0", "2", "20", "0", "0"],
                ["1", "2", "10", "0", "0"],
            ])

    def test_duplicate_pair_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "raw.edges"
            raw.write_text("1 2 1 0 0.9\n2 1 1 0 0.8\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate SpydrPick pair"):
                normalize_edges(raw, root / "edges.gz", [10, 20])


if __name__ == "__main__":
    unittest.main()
