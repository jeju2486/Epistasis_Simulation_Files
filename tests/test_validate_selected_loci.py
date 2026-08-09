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


if __name__ == "__main__":
    unittest.main()
