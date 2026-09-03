from __future__ import annotations

import gzip
import argparse
import sys
import tempfile
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from kovar_inputs import materialize_pairs, write_binary_fasta, write_maf_filtered_binary_fasta  # noqa: E402
from run_kovar_case import requested_settings, require  # noqa: E402


class KovarInputTests(unittest.TestCase):
    def test_minimum_cell_count_zero_is_part_of_restart_signature(self) -> None:
        settings = requested_settings(argparse.Namespace(
            min_maf=0.05, min_cell_count=0, spa_mode="auto"
        ))
        self.assertEqual(settings["min_cell_count"], 0)
        self.assertEqual(settings["candidate_source"],
                         "spydrpick_all_pairs_default_weighting")

    def test_empty_success_marker_is_valid(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            marker = Path(directory) / "_SUCCESS"
            marker.touch()
            require(marker, nonempty=False)

    def test_binary_encoding_uses_vcf_reference_state(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta, positions, output = root / "all.fa", root / "positions.tsv", root / "binary.fa"
            fasta.write_text(">s1\nATG\n>s2\nGCA\n", encoding="utf-8")
            positions.write_text("alignment_column\tchrom\tvcf_position\tvariant_id\talleles\n"
                                 "0\t1\t10\t.\tA,G\n1\t1\t20\t.\tC,T\n2\t1\t30\t.\tG,A\n", encoding="utf-8")
            self.assertEqual(write_binary_fasta(fasta, positions, output), (2, 3))
            self.assertEqual(output.read_text(encoding="utf-8"), ">s1\nACA\n>s2\nCAC\n")

    def test_pair_materialization_is_canonical_uv_only(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root / "pairs.gz", root / "pairs.tsv"
            with gzip.open(source, "wt", encoding="utf-8") as handle:
                handle.write("2 0 100 0 0.8\n1 3 200 0 0.7\n")
            self.assertEqual(materialize_pairs(source, output), 2)
            self.assertEqual(output.read_text(encoding="utf-8"), "u\tv\n0\t2\n1\t3\n")

    def test_maf_filter_writes_shared_locus_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta, positions = root / "all.fa", root / "positions.tsv"
            output, mapping = root / "filtered.fa", root / "eligible.tsv"
            fasta.write_text(">s1\nAAA\n>s2\nCAA\n>s3\nCCA\n>s4\nCCA\n", encoding="utf-8")
            positions.write_text("alignment_column\tchrom\tvcf_position\tvariant_id\talleles\n"
                                 "0\t1\t10\t.\tA,C\n1\t1\t20\t.\tA,C\n2\t1\t30\t.\tA,C\n", encoding="utf-8")
            self.assertEqual(write_maf_filtered_binary_fasta(
                fasta, positions, output, mapping, min_maf=0.25), (4, 3, 2))
            self.assertIn("0\t0\t9\t0.75\t0.25", mapping.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
