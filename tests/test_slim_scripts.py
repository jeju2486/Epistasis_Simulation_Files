from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]


class SlimScriptTests(unittest.TestCase):
    def test_experimental_seeding_uses_slim_52_mutation_api(self) -> None:
        script = (REPO_ROOT / "slim" / "continue_experiment.slim").read_text()

        self.assertIn("originSubpop=NULL, nucleotide=derived", script)
        self.assertNotIn("NULL, NULL, derived", script)
        self.assertIn("balanced_four_haplotype_cycle", script)
        self.assertIn("targets = carriers.haplosomesNonNull;", script)

    def test_haplotype_cycle_pairs_complementary_states(self) -> None:
        script = (REPO_ROOT / "slim" / "continue_experiment.slim").read_text()

        self.assertIn(
            "haplotypeOrder = c(0, 15, 1, 14, 2, 13, 3, 12,", script
        )
        self.assertIn("haplotypeOrder = c(0, 3, 1, 2);", script)
        self.assertIn(
            "shuffled.tag = haplotypeOrder[integerMod(seqAlong(shuffled), cycleSize)];",
            script,
        )
        for order, loci in (
            ([0, 15, 1, 14, 2, 13, 3, 12, 4, 11, 5, 10, 6, 9, 7, 8], 4),
            ([0, 3, 1, 2], 2),
        ):
            for population_size in range(1, 101):
                states = [order[index % len(order)] for index in range(population_size)]
                for bit in range(loci):
                    carriers = sum((state >> bit) & 1 for state in states)
                    self.assertLessEqual(abs(population_size - 2 * carriers), 1)

    def test_cross_hgt_is_active_in_neutral_checkpoint(self) -> None:
        script = (REPO_ROOT / "slim" / "build_checkpoint.slim").read_text()
        self.assertIn("HGT_P_CROSS", script)
        self.assertIn("otherPops = sim.subpopulations[sim.subpopulations != subpop];", script)
        self.assertIn('defineConstant("DEEP_SPLIT_TICK", 1 + ANCESTRAL_GENERATIONS);', script)

    def test_tree_sequence_time_unit_matches_msprime(self) -> None:
        for filename in ("build_checkpoint.slim", "continue_experiment.slim"):
            script = (REPO_ROOT / "slim" / filename).read_text()
            self.assertIn('initializeTreeSeq(timeUnit="generations");', script)
