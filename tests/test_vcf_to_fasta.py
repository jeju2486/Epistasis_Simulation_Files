from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from vcf_to_fasta import convert  # noqa: E402


VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tp3:i0\tp4:i1
1\t2\t1\tA\tC\t.\tPASS\t.\tGT\t0\t1
1\t5\t2\tG\tT\t.\tPASS\t.\tGT\t1\t0
1\t9\t3\tA\tAT\t.\tPASS\t.\tGT\t1\t0
"""


class VcfToFastaTests(unittest.TestCase):
    def test_conversion_and_exclusion(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            vcf = root / "test.vcf"
            sample_map = root / "samples.tsv"
            loci = root / "loci.tsv"
            output = root / "out.fa"
            positions = root / "positions.tsv"
            vcf.write_text(VCF, encoding="utf-8")
            sample_map.write_text(
                "vcf_label\tanalysis_label\tpopulation\tpedigree_id\tindividual_index\n"
                "p3:i0\tp3_i0\tp1\t10\t0\n"
                "p4:i1\tp4_i1\tp2\t11\t1\n",
                encoding="utf-8",
            )
            # Position 1 in SLiM is position 2 in VCF.
            loci.write_text(
                "label\tmutation_id\tposition\tglobal_frequency\tp1_frequency\tp2_frequency\tp3_frequency\tp4_frequency\n"
                "A\t1\t1\t0.5\t0.5\t0.5\t0.5\t0.5\n",
                encoding="utf-8",
            )

            samples, sites = convert(vcf, output, positions, sample_map, loci)
            self.assertEqual((samples, sites), (2, 1))
            self.assertEqual(output.read_text(encoding="utf-8"), ">p3_i0\nT\n>p4_i1\nG\n")
            with open(positions, encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["vcf_position"], "5")


if __name__ == "__main__":
    unittest.main()
