from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from evaluate_lineage_confounding import (  # noqa: E402
    classify_pair,
    decompose_covariance,
    evaluate_case,
    kovar_top_threshold,
    triangular_index,
)


class LineageEvaluationTests(unittest.TestCase):
    def test_triangular_index_covers_each_unordered_pair_once(self) -> None:
        indices = {
            triangular_index(u, v, 5)
            for u in range(5)
            for v in range(u + 1, 5)
        }
        self.assertEqual(indices, set(range(10)))
        self.assertEqual(triangular_index(4, 1, 5), triangular_index(1, 4, 5))

    def test_covariance_decomposition_identifies_lineage_component(self) -> None:
        total, within, between, fraction = decompose_covariance(
            2, 0, 0, 2, [0.0, 1.0], [0.0, 1.0], [0.5, 0.5]
        )
        self.assertAlmostEqual(total, 0.25)
        self.assertAlmostEqual(within, 0.0)
        self.assertAlmostEqual(between, 0.25)
        self.assertAlmostEqual(fraction, 1.0)

    def test_covariance_decomposition_identifies_within_population_component(self) -> None:
        total, within, between, fraction = decompose_covariance(
            2, 0, 0, 2, [0.5, 0.5], [0.5, 0.5], [0.5, 0.5]
        )
        self.assertAlmostEqual(total, 0.25)
        self.assertAlmostEqual(within, 0.25)
        self.assertAlmostEqual(between, 0.0)
        self.assertAlmostEqual(fraction, 0.0)

    def test_pair_classification_separates_focal_linkage_and_lineage(self) -> None:
        common = dict(total_covariance=0.1, within_covariance=0.0, lineage_fraction=1.0)
        self.assertEqual(classify_pair(
            is_focal=True, focal_proximal=True, distance_bp=40_000, **common
        ), "focal_AB")
        self.assertEqual(classify_pair(
            is_focal=False, focal_proximal=True, distance_bp=40_000, **common
        ), "focal_proximal")
        self.assertEqual(classify_pair(
            is_focal=False, focal_proximal=False, distance_bp=1_000, **common
        ), "short_distance")
        self.assertEqual(classify_pair(
            is_focal=False, focal_proximal=False, distance_bp=40_000, **common
        ), "lineage_driven")

    def test_kovar_top_threshold_ignores_filtered_rows(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ko_variation.tsv"
            path.write_text(
                "u\tv\tp_primary\n0\t1\tNA\n0\t2\t0.2\n1\t2\t0.01\n",
                encoding="utf-8",
            )
            threshold, finite = kovar_top_threshold(path)
            self.assertEqual(finite, 2)
            self.assertEqual(threshold, 0.01)

    def test_evaluation_continues_when_focal_loci_are_not_maf_eligible(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            case_dir = Path(directory)
            spydrpick = case_dir / "spydrpick_all_pairs"
            kovar = case_dir / "kovar_v083_iqtree_spa_off"
            spydrpick.mkdir()
            kovar.mkdir()
            for marker in (case_dir / "_SUCCESS", spydrpick / "_SUCCESS", kovar / "_SUCCESS"):
                marker.touch()
            (case_dir / "selected_loci.tsv").write_text(
                "label\tposition\nA\t10000\nB\t50000\n", encoding="utf-8"
            )
            (case_dir / "sample_names.tsv").write_text(
                "analysis_label\tpopulation\n"
                "s1\tp1\ns2\tp1\ns3\tp2\ns4\tp2\n",
                encoding="utf-8",
            )
            (spydrpick / "eligible_loci.tsv").write_text(
                "filtered_column\toriginal_column\tslim_position\tprevalence\tmaf\n"
                "0\t0\t20000\t0.5\t0.5\n1\t1\t30000\t0.5\t0.5\n",
                encoding="utf-8",
            )
            (spydrpick / "all_snps.binary_ac.fa").write_text(
                ">s1\nAA\n>s2\nCC\n>s3\nAC\n>s4\nCA\n", encoding="utf-8"
            )
            import gzip

            with gzip.open(spydrpick / "spydrpick.edges.gz", "wt", encoding="utf-8") as handle:
                handle.write("0 1 10000 0 0.0\n")
            (kovar / "ko_variation.tsv").write_text(
                "u\tv\tn11\tn10\tn01\tn00\tp_primary\tbonferroni_significant\n"
                "0\t1\t1\t1\t1\t1\t0.5\t0\n",
                encoding="utf-8",
            )

            rows = evaluate_case(
                {"case_id": "missing_focal", "out_dir": str(case_dir)},
                "metrics.tsv.gz", "off",
            )
            self.assertEqual(len(rows), 1)
            self.assertNotEqual(rows[0]["category"], "focal_AB")
            self.assertTrue((case_dir / "metrics.tsv.gz").exists())


if __name__ == "__main__":
    unittest.main()
