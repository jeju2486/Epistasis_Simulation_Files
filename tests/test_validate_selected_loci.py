from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from validate_selected_loci import validate  # noqa: E402


HEADER = "label\tmutation_id\tposition\n"


class SelectedLociValidationTests(unittest.TestCase):
    def test_accepts_four_labels_with_blank_records(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "selected_loci.tsv"
            path.write_text(
                HEADER + "\nA\t1\t10\n\nB\t2\t20\n\nC\t3\t30\n\nD\t4\t40\n",
                encoding="utf-8",
            )
            validate(path)

    def test_rejects_duplicate_and_missing_labels(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "selected_loci.tsv"
            path.write_text(HEADER + "A\t1\t10\nB\t2\t20\nC\t3\t30\nC\t4\t40\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "exactly one record"):
                validate(path)

    def test_validates_predeclared_positions(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "selected_loci.tsv"
            path.write_text(
                HEADER + "A\t1\t10\nB\t2\t20\nC\t3\t60\nD\t4\t85\n",
                encoding="utf-8",
            )
            expected = {"A": 10, "B": 20, "C": 60, "D": 85}
            validate(path, expected)
            expected["D"] = 86
            with self.assertRaisesRegex(ValueError, "expected 86"):
                validate(path, expected)

    def test_rejects_unbalanced_seed_frequency(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "selected_loci.tsv"
            path.write_text(
                "label\tmutation_id\tposition\tseeding_design\tglobal_frequency\n"
                "A\t1\t10\tbalanced_16_haplotype_cycle\t0.8\n"
                "B\t2\t20\tbalanced_16_haplotype_cycle\t0.5\n"
                "C\t3\t60\tbalanced_16_haplotype_cycle\t0.5\n"
                "D\t4\t85\tbalanced_16_haplotype_cycle\t0.5\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "unbalanced initial"):
                validate(path)


if __name__ == "__main__":
    unittest.main()
