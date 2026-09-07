from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from infer_iqtree_case import write_variable_alignment  # noqa: E402


class IqtreeInputTests(unittest.TestCase):
    def test_alignment_excludes_focal_and_invariant_sites(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            alignment = root / "all_snps.fa"
            positions = root / "all_snps.positions.tsv"
            selected = root / "selected_loci.tsv"
            output = root / "nonfocal.fa"
            site_map = root / "site_map.tsv"
            alignment.write_text(
                ">s1\nAACG\n>s2\nACCG\n>s3\nAGTG\n", encoding="utf-8"
            )
            positions.write_text(
                "alignment_column\tchrom\tvcf_position\tvariant_id\talleles\n"
                "0\t1\t11\t.\tA,C\n1\t1\t21\t.\tA,C,G\n"
                "2\t1\t31\t.\tC,T\n3\t1\t41\t.\tG,A\n",
                encoding="utf-8",
            )
            selected.write_text(
                "label\tposition\nA\t20\nB\t40\n", encoding="utf-8"
            )

            self.assertEqual(
                write_variable_alignment(alignment, positions, selected, output, site_map),
                (3, 4, 1),
            )
            self.assertEqual(output.read_text(encoding="utf-8"), ">s1\nC\n>s2\nC\n>s3\nT\n")
            with site_map.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows, [{
                "iqtree_column": "0", "original_column": "2", "slim_position": "30"
            }])

    def test_alignment_rejects_no_remaining_nonfocal_variation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "all.fa").write_text(">s1\nA\n>s2\nC\n", encoding="utf-8")
            (root / "positions.tsv").write_text(
                "alignment_column\tvcf_position\n0\t21\n", encoding="utf-8"
            )
            (root / "selected.tsv").write_text(
                "label\tposition\nA\t20\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "no nonfocal variable sites"):
                write_variable_alignment(
                    root / "all.fa", root / "positions.tsv", root / "selected.tsv",
                    root / "out.fa", root / "map.tsv",
                )


if __name__ == "__main__":
    unittest.main()
