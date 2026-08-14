from __future__ import annotations

import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from plot_kovar import kovar_points  # noqa: E402
from plot_spydrpick import distance_mi_points, focal_pair_columns  # noqa: E402


class PairPlotTests(unittest.TestCase):
    def test_all_six_focal_pairs_are_identified_and_never_thinned(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "selected_loci.tsv"
            selected.write_text(
                "label\tposition\nA\t10\nB\t20\nC\t60\nD\t85\n",
                encoding="utf-8",
            )
            positions = [10, 20, 60, 85, 90]
            focal_columns = focal_pair_columns(selected, positions)
            self.assertEqual(set(focal_columns.values()), {"AB", "AC", "AD", "BC", "BD", "CD"})

            edges = root / "edges.gz"
            with gzip.open(edges, "wt", encoding="utf-8") as handle:
                for u, v, mi in (
                    (0, 1, 0.9), (0, 2, 0.8), (0, 3, 0.7), (1, 2, 0.6),
                    (1, 3, 0.5), (2, 3, 0.4), (0, 4, 0.3), (1, 4, 0.2),
                    (2, 4, 0.1), (3, 4, 0.05),
                ):
                    handle.write(f"{u} {v} {abs(positions[v] - positions[u])} 0 {mi}\n")
            _distance, _mi, focal = distance_mi_points(
                edges, positions, focal_columns, max_background_points=1
            )
            self.assertEqual(set(focal), set(focal_columns.values()))
            self.assertEqual(focal["AB"], (0.01, 0.9))
            self.assertEqual(focal["CD"], (0.025, 0.4))

    def test_kovar_keeps_both_focal_directions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            table = root / "ko_variation.tsv"
            fields = ["u", "v", "direction", "p_primary", "n_directional_tests"]
            with table.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerow({"u": 0, "v": 1, "direction": "u_predicts_v", "p_primary": 0.01, "n_directional_tests": 4})
                writer.writerow({"u": 0, "v": 1, "direction": "v_predicts_u", "p_primary": 0.001, "n_directional_tests": 4})
                writer.writerow({"u": 0, "v": 2, "direction": "u_predicts_v", "p_primary": 0.5, "n_directional_tests": 4})
                writer.writerow({"u": 0, "v": 2, "direction": "v_predicts_u", "p_primary": 0.5, "n_directional_tests": 4})
            background, focal, n_tests = kovar_points(
                table, [10, 20, 60], {(0, 1): "AB"}, expected_pairs=2,
                max_background_points=10,
            )
            self.assertEqual(n_tests, 4)
            self.assertEqual([row[0] for row in focal["AB"]], ["u_predicts_v", "v_predicts_u"])
            self.assertAlmostEqual(focal["AB"][0][2], 2.0)
            self.assertAlmostEqual(focal["AB"][1][2], 3.0)
            self.assertEqual(len(background["u_predicts_v"][0]), 1)
            self.assertEqual(len(background["v_predicts_u"][0]), 1)


if __name__ == "__main__":
    unittest.main()
