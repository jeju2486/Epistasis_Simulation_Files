from __future__ import annotations

import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from build_manifest import build  # noqa: E402
from simflow import deterministic_seed, to_msprime_seed  # noqa: E402


class ManifestTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = {
            "paths": {"reference": "inputs/reference.fa", "checkpoint_root": "checkpoints",
                      "run_root": "runs", "manifest_root": "manifests"},
            "design": {"replicates": 2, "modes": [0, 1, 2],
                       "cross_hgt_probabilities": [0.0, 0.002, 0.02], "master_seed": 123},
            "genome": {"length": 100, "mutation_rate": 1e-7, "mean_hgt_tract": 5,
                       "within_hgt_probability": 0.02},
            "population": {"ancestral_size": 100, "clade_size": 50, "terminal_size": 25,
                           "sample_per_terminal": 5, "deep_split_tick": 10,
                           "terminal_split_tick": 20, "end_tick": 40},
            "loci": {"a_position": 10, "b_position": 50},
            "frequency_dependence": {"strength": 0.25, "epsilon": 0.001},
            "equilibrium": {"monitor_every": 5, "minimum_ticks": 10,
                            "stable_checks": 2, "tolerance": 0.1},
            "postprocess": {"ancestral_ne": 100, "tree_position": 90},
        }

    def test_expansion_has_one_checkpoint_and_nine_cases_per_replicate(self) -> None:
        checkpoints, cases = build(self.config)
        self.assertEqual(len(checkpoints), 2)
        self.assertEqual(len(cases), 18)
        self.assertEqual({row["checkpoint_id"] for row in cases}, {"rep_0001", "rep_0002"})
        self.assertEqual({row["mode"] for row in cases}, {0, 1, 2})
        self.assertEqual({row["cross_hgt_probability"] for row in cases}, {0.0, 0.002, 0.02})
        self.assertEqual(
            {row["regime"] for row in cases},
            {"balanced_independent_negative_control", "global_high_frequency_independent",
             "global_high_frequency_dependent"},
        )
        self.assertEqual(cases[0]["out_dir"], "runs/rep_0001/cross_0/mode_0")
        self.assertEqual(cases[0]["case_id"], "rep_0001__cross_0__mode_0")
        self.assertEqual(cases[4]["cross_hgt_label"], "0p002")

    def test_invalid_timeline_is_rejected(self) -> None:
        self.config["population"]["end_tick"] = 20
        with self.assertRaisesRegex(ValueError, "deep split"):
            build(self.config)

    def test_invalid_mode_is_rejected(self) -> None:
        self.config["design"]["modes"] = [0, 3]
        with self.assertRaisesRegex(ValueError, "only 0, 1, and 2"):
            build(self.config)

    def test_invalid_hgt_sum_is_rejected(self) -> None:
        self.config["design"]["cross_hgt_probabilities"] = [0.99]
        with self.assertRaisesRegex(ValueError, "sum to at most one"):
            build(self.config)

    def test_seeds_are_stable_and_unique(self) -> None:
        _, first = build(self.config)
        _, second = build(self.config)
        self.assertEqual([row["seed"] for row in first], [row["seed"] for row in second])
        self.assertEqual(len({row["seed"] for row in first}), len(first))
        self.assertEqual(deterministic_seed(1, "a"), deterministic_seed(1, "a"))

    def test_msprime_seed_is_deterministic_and_in_range(self) -> None:
        mapped = to_msprime_seed(3163853192524712447)
        self.assertGreater(mapped, 0)
        self.assertLess(mapped, 1 << 32)
        self.assertEqual(mapped, to_msprime_seed(3163853192524712447))


if __name__ == "__main__":
    unittest.main()
