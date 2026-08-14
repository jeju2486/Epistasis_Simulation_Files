from __future__ import annotations

import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

try:
    import tskit

    from export_oracle_tree import generations_time_units
except ModuleNotFoundError:
    tskit = None
    generations_time_units = None


@unittest.skipIf(tskit is None, "tree-sequence dependencies are unavailable")
class OracleTreeTimeUnitTests(unittest.TestCase):
    @staticmethod
    def tree_sequence(time_units: str):
        tables = tskit.TableCollection(sequence_length=1)
        tables.time_units = time_units
        return tables.tree_sequence()

    def test_relabels_legacy_ticks_without_rescaling_times(self) -> None:
        ts = self.tree_sequence("ticks")
        converted = generations_time_units(ts)
        self.assertEqual(converted.time_units, "generations")
        self.assertEqual(converted.tables.nodes, ts.tables.nodes)

    def test_preserves_generation_label(self) -> None:
        ts = self.tree_sequence("generations")
        self.assertIs(generations_time_units(ts), ts)

    def test_rejects_unknown_units(self) -> None:
        with self.assertRaisesRegex(ValueError, "Expected tree time units"):
            generations_time_units(self.tree_sequence("years"))
