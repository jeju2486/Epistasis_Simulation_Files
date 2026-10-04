from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from plot_figure1 import average_ranks, null_values, read_case, summarize


def fixture(root: Path, *, filtered_focal=False) -> tuple[dict, Path]:
    case = root / "case"
    marginal = case / "spydrpick_all_pairs"
    kovar = case / "kovar_v083_iqtree_spa_off"
    for stage in (case, marginal, kovar):
        stage.mkdir(parents=True, exist_ok=True)
        (stage / "_SUCCESS").touch()
    positions = [10000, 20000, 50000, 80000]
    (marginal / "eligible_loci.tsv").write_text(
        "filtered_column\tslim_position\n" + "".join(f"{i}\t{x}\n" for i, x in enumerate(positions)), encoding="utf-8")
    (case / "selected_loci.tsv").write_text(
        f"label\tposition\nA\t{10001 if filtered_focal else 10000}\nB\t50000\n", encoding="utf-8")
    (case / "sample_names.tsv").write_text(
        "analysis_label\tpopulation\nx1\tp1\nx2\tp1\nx3\tp2\nx4\tp2\n", encoding="utf-8")
    (marginal / "all_snps.binary_ac.fa").write_text(
        ">x1\nAAAA\n>x2\nAACC\n>x3\nCCAA\n>x4\nCCCC\n", encoding="utf-8")
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    with gzip.open(marginal / "spydrpick.edges.gz", "wt", encoding="utf-8") as handle:
        for i, (u, v) in enumerate(pairs):
            handle.write(f"{u} {v} {positions[v]-positions[u]} 0 {[.6,.5,.3,.3,0,0][i]}\n")
    path = kovar / "ko_variation.tsv"
    # Deliberately reverse row order/orientation: joins must use pair identities.
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["u", "v", "p_primary", "n11", "n10", "n01", "n00"], delimiter="\t")
        writer.writeheader()
        for i in reversed(range(6)):
            u, v = pairs[i]
            writer.writerow({"u": v, "v": u, "p_primary": [.5,.001,.3,.3,"NA",0][i],
                             "n11": 1, "n10": 1, "n01": 1, "n00": 1})
    row = {"case_id": "test_case", "out_dir": str(case), "replicate": "1", "mode": "2", "cross_hgt_probability": "0"}
    return row, path


class Figure1Tests(unittest.TestCase):
    def test_ties_and_unavailable_pvalues(self):
        values = np.array([0., 0., .5, np.nan, 1.])
        np.testing.assert_allclose(average_ranks(values), [1.5, 1.5, 3, np.nan, 4], equal_nan=True)
        np.testing.assert_allclose(average_ranks(values, descending=True), [3.5, 3.5, 2, np.nan, 1], equal_nan=True)

    def test_identity_join_zero_mi_and_failed_tests_are_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            row, _ = fixture(Path(tmp))
            data = read_case(row, "off", "default")
            self.assertEqual(data["n_pairs"], 6)
            self.assertEqual(data["n_finite"], 5)
            self.assertEqual(data["focal"], 1)
            self.assertEqual(data["p"][1], .001)
            self.assertTrue(np.isnan(data["shared_mi_rank"][4]))
            self.assertEqual(data["mi"][5], 0)
            self.assertEqual(data["p_rank"][2], data["p_rank"][3])
            self.assertEqual(data["shared_mi_rank"][1], 2)

    def test_mismatched_universe_fails_instead_of_silent_intersection(self):
        with tempfile.TemporaryDirectory() as tmp:
            row, path = fixture(Path(tmp))
            lines = path.read_text().splitlines()
            path.write_text("\n".join(lines[:-1]) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "identical eligible pair universe"):
                read_case(row, "off", "default")

    def test_filtered_focal_is_unrecovered_not_dropped(self):
        with tempfile.TemporaryDirectory() as tmp:
            row, _ = fixture(Path(tmp), filtered_focal=True)
            data = read_case(row, "off", "default")
            focal, _ = summarize(data, SimpleNamespace(top_fraction=.01, alpha=.05, distal_bp=10000))
            self.assertEqual(focal["focal_status"], "not_maf_eligible")
            self.assertEqual(focal["mi_top_recovered"], 0)
            self.assertEqual(focal["kovar_top_recovered"], 0)
            self.assertEqual(focal["kovar_rank"], "NA")

    def test_null_pairs_must_have_independent_marginal_calibration(self):
        with tempfile.TemporaryDirectory() as tmp:
            row, _ = fixture(Path(tmp))
            data = read_case(row, "off", "default")
            with self.assertRaisesRegex(ValueError, "missing independently calibrated"):
                null_values(data, {("test_case", 0, 1): True}, {})

    def test_lineage_enrichment_uses_same_classified_distal_universe(self):
        data = {"case": {"case_id": "x", "replicate": "1", "mode": "2", "cross_hgt_probability": "0"},
                "n_pairs": 6, "n_finite": 5, "focal": None,
                "mi": np.array([.99, .9, .8, .7, .6, .5]),
                "p": np.array([.001, .5, .1, .2, np.nan, .01]),
                "mi_rank": np.arange(1, 7), "shared_mi_rank": np.array([1, 2, 3, 4, np.nan, 5]),
                "p_rank": np.array([1, 5, 3, 4, np.nan, 2]),
                "distance": np.array([10000, 20000, 30000, 40000, 50000, 60000]),
                "category": np.array([2, 2, 0, 0, 2, 6])}
        _, row = summarize(data, SimpleNamespace(top_fraction=.2, alpha=.05, distal_bp=10000))
        # Excludes exact 10 kb, unavailable P, and unclassified joint counts.
        self.assertEqual(row["distal_pairs"], 3)
        self.assertEqual(row["distal_lineage_pairs"], 1)
        self.assertEqual(row["mi_lineage_enrichment"], 3)
        self.assertEqual(row["kovar_lineage_enrichment"], 0)

    def test_cli_renders_all_panels_with_explicit_null_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            row, _ = fixture(root)
            manifest = root / "cases.tsv"
            with manifest.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
                writer.writeheader()
                writer.writerow(row)
            (root / "null.tsv").write_text("case_id\tu\tv\ntest_case\t0\t1\n", encoding="utf-8")
            (root / "p.tsv").write_text("case_id\tu\tv\tp_marginal\ntest_case\t0\t1\t0.5\n", encoding="utf-8")
            output = root / "plots"
            command = [sys.executable, str(SCRIPTS / "plot_figure1.py"), "--manifest", str(manifest),
                       "--spa-mode", "off", "--output-dir", str(output), "--max-points", "1",
                       "--null-pairs", str(root / "null.tsv"), "--marginal-pvalues", str(root / "p.tsv"),
                       "--calibration-note", "synthetic wiring fixture, not biological null evidence"]
            subprocess.run(command, check=True, capture_output=True, text=True)
            report = json.loads((output / "figure1_report.json").read_text())
            self.assertEqual(report["skipped_panels"], {})
            self.assertEqual(len(list(output.glob("*.png"))), 7)
            self.assertEqual(len(list(output.glob("*.svg"))), 7)
            self.assertEqual(len(list(output.glob("*.pdf"))), 7)
            for name in report["files"]:
                self.assertGreater((output / name).stat().st_size, 0)
            # Read exported data: unavailable KOVAR pair is retained, not omitted.
            with gzip.open(output / "pairs_test_case.tsv.gz", "rt", encoding="utf-8") as handle:
                pairs = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(pairs), 6)
            self.assertEqual(sum(row["kovar_status"] == "unavailable" for row in pairs), 1)
            # No calibration input: default generates B–E and explicitly marks A pending.
            pending = root / "pending"
            subprocess.run(command[:command.index("--null-pairs")] +
                           ["--output-dir", str(pending), "--formats", "png"],
                           check=True, capture_output=True, text=True)
            report = json.loads((pending / "figure1_report.json").read_text())
            self.assertIn("A", report["skipped_panels"])
            self.assertFalse(list(pending.glob("A_*.png")))


if __name__ == "__main__":
    unittest.main()
