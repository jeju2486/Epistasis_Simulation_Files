import math
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]


class SlimScriptTests(unittest.TestCase):
    def test_target_tables_have_the_intended_dependence_structure(self) -> None:
        mode0 = (0.25, 0.25, 0.25, 0.25)
        mode1 = (0.25, 0.25, 0.25, 0.25)
        mode2 = (0.45, 0.05, 0.05, 0.45)
        for table in (mode0, mode1):
            expected_11 = (table[2] + table[3]) * (table[1] + table[3])
            self.assertAlmostEqual(table[3], expected_11)
        self.assertEqual(mode1[2] + mode1[3], mode2[2] + mode2[3])
        self.assertEqual(mode1[1] + mode1[3], mode2[1] + mode2[3])
        self.assertAlmostEqual((mode2[0] * mode2[3]) / (mode2[1] * mode2[2]), 81.0)
        expected_mi_bits = sum(cell * math.log2(cell / 0.25) for cell in mode2)
        self.assertAlmostEqual(expected_mi_bits, 0.5310044064107189)

    def test_seeding_uses_four_haplotypes_and_slim_52_mutation_api(self) -> None:
        script = (REPO_ROOT / "slim" / "continue_experiment.slim").read_text()
        self.assertIn("haplotypeOrder = c(0, 3, 1, 2);", script)
        self.assertIn("shuffled.tag = haplotypeOrder[integerMod(seqAlong(shuffled), 4)];", script)
        self.assertIn("originSubpop=NULL, nucleotide=derived", script)
        self.assertIn("targets = carriers.haplosomesNonNull;", script)
        self.assertNotIn("C_POSITION", script)

    def test_modes_encode_the_predeclared_joint_targets(self) -> None:
        script = (REPO_ROOT / "slim" / "continue_experiment.slim").read_text()
        self.assertIn("c(0.25, 0.25, 0.25, 0.25)", script)
        self.assertIn("c(0.45, 0.05, 0.05, 0.45)", script)
        self.assertIn("FDS_STRENGTH * log((target + FDS_EPSILON) / (freqs + FDS_EPSILON))", script)

    def test_mode_zero_disables_focal_frequency_correction(self) -> None:
        script = (REPO_ROOT / "slim" / "continue_experiment.slim").read_text()
        self.assertIn("if (MODE == 0) {", script)
        self.assertIn('status = "neutral_no_focal_control";', script)
        neutral_block = script.split("if (MODE == 0) {", 1)[1].split("}", 1)[0]
        self.assertIn("juveniles.fitnessScaling", neutral_block)
        self.assertIn("next;", neutral_block)

    def test_cross_hgt_and_equilibrium_monitor_are_explicit(self) -> None:
        checkpoint = (REPO_ROOT / "slim" / "build_checkpoint.slim").read_text()
        self.assertNotIn("HGT_P_CROSS", checkpoint)
        continuation = (REPO_ROOT / "slim" / "continue_experiment.slim").read_text()
        self.assertIn("HGT_P_CROSS", continuation)
        self.assertIn("otherPops", continuation)
        self.assertIn("EQUILIBRIUM_MONITOR_EVERY", continuation)
        self.assertIn('status = "not_reached_by_sampling";', continuation)
        self.assertIn('status = "equilibrium_reached";', continuation)
        self.assertIn('status = "neutral_no_focal_control";', continuation)
        self.assertIn("sim.setValue(\"first_equilibrium_tick\", community.tick);", continuation)
        self.assertNotIn("stop(\"Equilibrium was not reached", continuation)

    def test_sampling_occurs_only_at_fixed_end_tick(self) -> None:
        script = (REPO_ROOT / "slim" / "continue_experiment.slim").read_text()
        self.assertIn("if (community.tick >= END_TICK) {", script)
        self.assertIn("finishRun(maximumDeviation);", script)
        monitor_block = script.split("if (shouldMonitor) {", 1)[1].split(
            "if (community.tick >= END_TICK)", 1
        )[0]
        self.assertNotIn("finishRun", monitor_block)

    def test_tree_sequence_time_unit_matches_msprime(self) -> None:
        for filename in ("build_checkpoint.slim", "continue_experiment.slim"):
            script = (REPO_ROOT / "slim" / filename).read_text()
            self.assertIn('initializeTreeSeq(timeUnit="generations");', script)


if __name__ == "__main__":
    unittest.main()
