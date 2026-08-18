from __future__ import annotations

import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from build_manifest import build  # noqa: E402
from simflow import deterministic_seed, probability_slug, to_msprime_seed  # noqa: E402


class ManifestTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = {
            "paths": {
                "reference": "inputs/reference.fa",
                "checkpoint_root": "checkpoints",
                "run_root": "runs",
            },
            "design": {
                "replicates": 2,
                "modes": [0, 1, 2],
                "cross_hgt_probabilities": [0.0, 0.002, 0.02],
                "master_seed": 123,
            },
            "genome": {
                "length": 100,
                "mutation_rate": 1e-7,
                "mean_hgt_tract": 5,
                "within_hgt_probability": 0.02,
            },
            "population": {
                "ancestral_size": 100,
                "clade_size": 50,
                "terminal_size": 25,
                "sample_per_terminal": 5,
                "ancestral_generations": 10,
                "deep_clade_generations": 4,
                "terminal_generations": 4,
            },
            "loci": {
                "a_position": 10,
                "b_position": 20,
                "c_position": 60,
                "d_position": 85,
            },
            "experiment": {
                "generations": [10, 20],
                "s_ab": 0.01,
                "s_cd": 0.015,
                "monitor_every": 2,
            },
            "postprocess": {"ancestral_ne": 100, "oracle_tree_position": 50},
        }

    def test_expansion_is_paired(self) -> None:
        checkpoints, cases = build(self.config)
        self.assertEqual(len(checkpoints), 6)
        self.assertEqual(len(cases), 36)
        self.assertEqual(
            {row["checkpoint_id"] for row in cases},
            {
                "rep_0001__cross_0", "rep_0001__cross_0p002", "rep_0001__cross_0p02",
                "rep_0002__cross_0", "rep_0002__cross_0p002", "rep_0002__cross_0p02",
            },
        )
        self.assertEqual({row["mode"] for row in cases}, {0, 1, 2})
        self.assertEqual({row["active_pair"] for row in cases}, {"AB_and_CD_neutral", "AB", "CD"})
        self.assertEqual({row["cross_hgt_probability"] for row in cases}, {0.0, 0.002, 0.02})
        self.assertEqual({row["experiment_generations"] for row in cases}, {10, 20})
        self.assertEqual(
            {row["cross_hgt_probability"] for row in checkpoints}, {0.0, 0.002, 0.02}
        )
        self.assertEqual(checkpoints[0]["a_position"], 10)
        self.assertEqual(cases[0]["d_position"], 85)

    def test_equal_focal_pair_distances_are_rejected(self) -> None:
        self.config["loci"]["d_position"] = 70
        with self.assertRaisesRegex(ValueError, "distances must differ"):
            build(self.config)

    def test_seeds_are_stable_and_unique(self) -> None:
        _, first = build(self.config)
        _, second = build(self.config)
        self.assertEqual([row["seed"] for row in first], [row["seed"] for row in second])
        self.assertEqual(len({row["seed"] for row in first}), len(first))
        self.assertEqual(deterministic_seed(1, "a"), deterministic_seed(1, "a"))

    def test_probability_slug(self) -> None:
        self.assertEqual(probability_slug(0.0), "0")
        self.assertEqual(probability_slug(0.002), "0p002")
        self.assertEqual(probability_slug(0.02), "0p02")

    def test_msprime_seed_is_deterministic_and_in_range(self) -> None:
        self.assertEqual(to_msprime_seed(1), 1)
        self.assertEqual(to_msprime_seed((1 << 32) - 1), 1)
        mapped = to_msprime_seed(3163853192524712447)
        self.assertGreater(mapped, 0)
        self.assertLess(mapped, 1 << 32)
        self.assertEqual(mapped, to_msprime_seed(3163853192524712447))


if __name__ == "__main__":
    unittest.main()
