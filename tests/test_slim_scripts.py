from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]


class SlimScriptTests(unittest.TestCase):
    def test_target_tables_have_the_intended_dependence_structure(self) -> None:
        mode1_high = (0.0225, 0.1275, 0.1275, 0.7225)
        mode1_low = tuple(reversed(mode1_high))
        for table in (mode1_high, mode1_low):
            expected_11 = (table[2] + table[3]) * (table[1] + table[3])
            self.assertAlmostEqual(table[3], expected_11)
        pooled_mode1 = tuple((left + right) / 2 for left, right in zip(mode1_high, mode1_low))
        self.assertNotAlmostEqual(pooled_mode1[3], 0.25)

        mode2_left = (0.25, 0.05, 0.45, 0.25)
        mode2_right = (0.25, 0.45, 0.05, 0.25)
        for table in (mode2_left, mode2_right):
            self.assertGreater(table[0] * table[3], table[1] * table[2])
        pooled_mode2 = tuple((left + right) / 2 for left, right in zip(mode2_left, mode2_right))
        self.assertEqual(pooled_mode2, (0.25, 0.25, 0.25, 0.25))

    def test_seeding_uses_four_haplotypes_and_slim_52_mutation_api(self) -> None:
        script = (REPO_ROOT / "slim" / "continue_experiment.slim").read_text()
        self.assertIn("haplotypeOrder = c(0, 3, 1, 2);", script)
        self.assertIn("shuffled.tag = haplotypeOrder[integerMod(seqAlong(shuffled), 4)];", script)
        self.assertIn("originSubpop=NULL, nucleotide=derived", script)
        self.assertIn("targets = carriers.haplosomesNonNull;", script)
        self.assertNotIn("C_POSITION", script)

    def test_modes_encode_the_predeclared_joint_targets(self) -> None:
        script = (REPO_ROOT / "slim" / "continue_experiment.slim").read_text()
        self.assertIn("c(0.0225, 0.1275, 0.1275, 0.7225)", script)
        self.assertIn("c(0.7225, 0.1275, 0.1275, 0.0225)", script)
        self.assertIn("c(0.25, 0.05, 0.45, 0.25)", script)
        self.assertIn("c(0.25, 0.45, 0.05, 0.25)", script)
        self.assertIn("FDS_STRENGTH * log((target + FDS_EPSILON) / (freqs + FDS_EPSILON))", script)

    def test_no_cross_hgt_migration_or_periodic_monitor(self) -> None:
        for filename in ("build_checkpoint.slim", "continue_experiment.slim"):
            script = (REPO_ROOT / "slim" / filename).read_text()
            self.assertNotIn("HGT_P_CROSS", script)
            self.assertNotIn("MONITOR_EVERY", script)
        continuation = (REPO_ROOT / "slim" / "continue_experiment.slim").read_text()
        self.assertNotIn("otherPops", continuation)

    def test_tree_sequence_time_unit_matches_msprime(self) -> None:
        for filename in ("build_checkpoint.slim", "continue_experiment.slim"):
            script = (REPO_ROOT / "slim" / filename).read_text()
            self.assertIn('initializeTreeSeq(timeUnit="generations");', script)


if __name__ == "__main__":
    unittest.main()
