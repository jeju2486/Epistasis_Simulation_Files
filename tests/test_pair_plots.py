from __future__ import annotations

import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from plot_kovar import read_results  # noqa: E402
from plot_spydrpick import focal_pair_columns, read_points  # noqa: E402


class PairPlotTests(unittest.TestCase):
    def test_ab_is_identified_and_never_thinned(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "selected_loci.tsv"
            selected.write_text("label\tposition\nA\t10\nB\t50\n", encoding="utf-8")
            positions = [10, 20, 50]
            focal_columns = focal_pair_columns(selected, positions)
            self.assertEqual(focal_columns, (0, 2))
            edges = root / "edges.gz"
            with gzip.open(edges, "wt", encoding="utf-8") as handle:
                handle.write("0 2 40 0 0.9\n0 1 10 0 0.8\n1 2 30 0 0.7\n")
            _distance, _mi, focal, bins = read_points(edges, focal_columns, 1)
            self.assertEqual(focal, (0.04, 0.9, 1))
            self.assertTrue(bins)

    def test_kovar_reads_one_unordered_focal_result(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            table = Path(temporary) / "ko_variation.tsv"
            fields = ["u", "v", "p_primary", "bonferroni_significant", "n_tests"]
            with table.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerow({"u": 0, "v": 1, "p_primary": 0.001,
                                 "bonferroni_significant": 1, "n_tests": 2})
                writer.writerow({"u": 0, "v": 2, "p_primary": 0.5,
                                 "bonferroni_significant": 0, "n_tests": 2})
            points, focal = read_results(table, [10, 50, 90], (0, 1))
            self.assertEqual(len(points), 2)
            self.assertEqual(focal["u"], 0)
            self.assertAlmostEqual(focal["score"], 3.0)
            self.assertEqual(focal["bonferroni"], 1)


if __name__ == "__main__":
    unittest.main()
