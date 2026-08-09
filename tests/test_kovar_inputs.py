from __future__ import annotations

import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from kovar_inputs import materialize_pairs, write_binary_fasta  # noqa: E402
from run_kovar_case import require  # noqa: E402


class KovarInputTests(unittest.TestCase):
    def test_empty_success_marker_is_valid(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            marker = Path(directory) / "_SUCCESS"
            marker.touch()
            require(marker, nonempty=False)

    def test_binary_encoding_uses_vcf_reference_state(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "all.fa"
            positions = root / "positions.tsv"
            output = root / "binary.fa"
            fasta.write_text(">s1\nATG\n>s2\nGCA\n", encoding="utf-8")
            positions.write_text(
                "alignment_column\tchrom\tvcf_position\tvariant_id\talleles\n"
                "0\t1\t10\t.\tA,G\n1\t1\t20\t.\tC,T\n2\t1\t30\t.\tG,A\n",
                encoding="utf-8",
            )
            self.assertEqual(write_binary_fasta(fasta, positions, output), (2, 3))
            self.assertEqual(output.read_text(encoding="utf-8"), ">s1\nACA\n>s2\nCAC\n")

    def test_pair_materialization_and_diagnostic_limit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "pairs.gz"
            with gzip.open(source, "wt", encoding="utf-8") as handle:
                handle.write("0 2 100 1 0.8\n1 3 200 1 0.7\n")
            output = root / "pairs.tsv"
            self.assertEqual(materialize_pairs(source, output, max_pairs=1), 1)
            self.assertEqual(output.read_text(encoding="utf-8"), "0 2 100 1 0.8\n")


if __name__ == "__main__":
    unittest.main()
