from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from validate_selected_loci import validate  # noqa: E402


class SelectedLociValidationTests(unittest.TestCase):
    def write(self, root: Path, a_frequency: float = 0.5, labels=("A", "B")) -> Path:
        path = root / "selected_loci.tsv"
        rows = ["label\tmutation_id\tposition\tseeded_frequency\tmode"]
        values = {"A": (10, 10, a_frequency), "B": (11, 50, 0.5), "C": (12, 70, 0.5)}
        rows.extend(f"{label}\t{values[label][0]}\t{values[label][1]}\t{values[label][2]}\t1" for label in labels)
        path.write_text("\n".join(rows) + "\n", encoding="utf-8")
        return path

    def test_accepts_two_balanced_loci(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            validate(self.write(Path(directory)), {"A": 10, "B": 50}, mode=1)

    def test_accepts_mode_zero(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = self.write(Path(directory))
            path.write_text(path.read_text().replace("\t1\n", "\t0\n"), encoding="utf-8")
            validate(path, {"A": 10, "B": 50}, mode=0)

    def test_rejects_extra_or_missing_label(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "exactly one A and one B"):
                validate(self.write(Path(directory), labels=("A", "B", "C")))

    def test_rejects_unbalanced_seed_frequency(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "unbalanced"):
                validate(self.write(Path(directory), a_frequency=0.6))


if __name__ == "__main__":
    unittest.main()
