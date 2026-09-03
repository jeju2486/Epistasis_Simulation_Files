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
            self.assertEqual(normalize_edges(raw, output, [10, 20, 40]), 3)
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                rows = handle.readlines()
            self.assertEqual(rows[0].split()[:3], ["0", "1", "10"])
            self.assertEqual(rows[1].split()[:3], ["0", "2", "30"])

    def test_missing_pair_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "raw.edges"
            raw.write_text("1 2 1 0 0.9\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "expected all 3 pairs"):
                normalize_edges(raw, root / "edges.gz", [10, 20, 30])


if __name__ == "__main__":
    unittest.main()
