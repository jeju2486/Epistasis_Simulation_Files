from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]


class SlimScriptTests(unittest.TestCase):
    def test_checkpoint_seeding_uses_slim_52_mutation_api(self) -> None:
        script = (REPO_ROOT / "slim" / "build_checkpoint.slim").read_text()

        self.assertIn("originSubpop=NULL, nucleotide=derived", script)
        self.assertNotIn("NULL, NULL, derived", script)
        self.assertIn("targets = carriers.haplosomesNonNull;", script)

    def test_tree_sequence_time_unit_matches_msprime(self) -> None:
        for filename in ("build_checkpoint.slim", "continue_experiment.slim"):
            script = (REPO_ROOT / "slim" / filename).read_text()
            self.assertIn('initializeTreeSeq(timeUnit="generations");', script)
