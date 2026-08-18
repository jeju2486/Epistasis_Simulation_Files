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
                "1 2 1 1 0.9\n1 3 2 1 0.8\n2 3 1 1 0.7\n"
                "1 4 3 1 0.6\n2 4 2 1 0.5\n3 4 1 1 0.4\n",
                encoding="utf-8",
            )
            metadata = compress_and_summarize(
                raw, root / "edges.gz", root / "truth.tsv",
                [10, 20, 30, 40],
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
                normalized = handle.readlines()
            self.assertEqual(len(normalized), 6)
            self.assertEqual(normalized[0].split()[:3], ["0", "1", "10"])

    def test_missing_pair_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "raw.edges"
            raw.write_text("1 2 1 1 0.9\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "all 3 pairs"):
                compress_and_summarize(
                    raw, root / "edges.gz", root / "truth.tsv",
                    [10, 20, 30],
                    {"A": 10, "B": 20, "C": 20, "D": 30}, 3,
                )

    def test_absent_truth_locus_is_labelled(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "raw.edges"
            raw.write_text("1 2 1 1 0.9\n", encoding="utf-8")
            compress_and_summarize(
                raw, root / "edges.gz", root / "truth.tsv",
                [10, 20], {"A": 10, "B": 99, "C": 10, "D": 20}, 2,
            )
            with (root / "truth.tsv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["candidate_status"], "locus_absent")
            self.assertEqual(rows[1]["candidate_status"], "eligible")

    def test_maf_filtered_truth_locus_is_distinguished(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "raw.edges"
            raw.write_text("1 2 1 1 0.9\n", encoding="utf-8")
            compress_and_summarize(
                raw, root / "edges.gz", root / "truth.tsv", [10, 20],
                {"A": 10, "B": 30, "C": 10, "D": 20}, 2, {10, 20, 30},
            )
            with (root / "truth.tsv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["candidate_status"], "maf_filtered")

    def test_inactive_mode_pair_is_labelled_not_seeded(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "raw.edges"
            raw.write_text("1 2 1 1 0.9\n", encoding="utf-8")
            compress_and_summarize(
                raw, root / "edges.gz", root / "truth.tsv", [10, 20],
                {"A": 10, "B": 20, "C": 30, "D": 40}, 2,
                {10, 20}, {"A", "B"},
            )
            with (root / "truth.tsv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["candidate_status"], "eligible")
            self.assertEqual(rows[1]["candidate_status"], "not_seeded")


if __name__ == "__main__":
    unittest.main()
