from __future__ import annotations

import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from spydrpick_case import compress_and_summarize  # noqa: E402


class SpydrPickSummaryTests(unittest.TestCase):
    def test_all_pairs_and_truth_ranks(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "raw.edges"
            raw.write_text(
                "10 20 10 1 0.9\n10 30 20 1 0.8\n20 30 10 1 0.7\n"
                "10 40 30 1 0.6\n20 40 20 1 0.5\n30 40 10 1 0.4\n",
                encoding="utf-8",
            )
            metadata = compress_and_summarize(
                raw, root / "edges.gz", root / "truth.tsv",
                {10: 0, 20: 1, 30: 2, 40: 3},
                {"A": 10, "B": 20, "C": 30, "D": 40}, 4,
            )
            self.assertEqual(metadata["n_pairs"], 6)
            with (root / "truth.tsv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["pair"], "AB")
            self.assertEqual(rows[0]["candidate_status"], "eligible")
            self.assertEqual(rows[0]["rank_min"], "1")
            self.assertEqual(rows[1]["pair"], "CD")
            self.assertEqual(rows[1]["rank_min"], "6")
            with gzip.open(root / "edges.gz", "rt", encoding="utf-8") as handle:
                self.assertEqual(sum(1 for _ in handle), 6)

    def test_missing_pair_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "raw.edges"
            raw.write_text("10 20 10 1 0.9\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "all 3 pairs"):
                compress_and_summarize(
                    raw, root / "edges.gz", root / "truth.tsv",
                    {10: 0, 20: 1, 30: 2},
                    {"A": 10, "B": 20, "C": 20, "D": 30}, 3,
                )

    def test_absent_truth_locus_is_labelled(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "raw.edges"
            raw.write_text("10 20 10 1 0.9\n", encoding="utf-8")
            compress_and_summarize(
                raw, root / "edges.gz", root / "truth.tsv",
                {10: 0, 20: 1}, {"A": 10, "B": 99, "C": 10, "D": 20}, 2,
            )
            with (root / "truth.tsv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["candidate_status"], "locus_absent")
            self.assertEqual(rows[1]["candidate_status"], "eligible")


if __name__ == "__main__":
    unittest.main()
