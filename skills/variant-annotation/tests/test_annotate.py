#!/usr/bin/env python3
"""Unit tests for annotate_vcf.py.

Focuses on the coordinate-conversion logic (the most error-prone part)
and VCF parsing edge cases.

Run:
    cd skills/variant-annotation
    python3 -m pytest tests/test_annotate.py -v
    # or without pytest:
    python3 tests/test_annotate.py
"""

import os
import sys
import tempfile
import unittest

# Ensure the parent directory is on the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from annotate_vcf import format_vep_region, parse_vcf


class TestFormatVepRegion(unittest.TestCase):
    """Test VCF → VEP region coordinate conversion.

    These are the most accuracy-critical transformations in the pipeline.
    Off-by-one errors here silently annotate the wrong genomic position.
    """

    # ---- SNVs ----

    def test_snv_basic(self):
        """SNV: single base change, coordinates unchanged."""
        result = format_vep_region("1", 11796321, "G", "A")
        self.assertEqual(result, "1 11796321 11796321 G/A 1")

    def test_snv_chr_prefix(self):
        """SNV: 'chr' prefix is stripped."""
        result = format_vep_region("chr10", 94942290, "C", "T")
        self.assertEqual(result, "10 94942290 94942290 C/T 1")

    def test_snv_chrX(self):
        """SNV: chrX → X."""
        result = format_vep_region("chrX", 100, "A", "G")
        self.assertEqual(result, "X 100 100 A/G 1")

    # ---- Deletions ----

    def test_deletion_single_base(self):
        """Deletion: VCF REF=AG, ALT=A at pos 100 → deleted G at pos 101.

        VCF:  pos=100, REF=AG, ALT=A  (padding base A at 100, G deleted at 101)
        VEP:  chrom 101 101 - 1
        """
        result = format_vep_region("1", 100, "AG", "A")
        self.assertEqual(result, "1 101 101 G/- 1")

    def test_deletion_multi_base(self):
        """Deletion: VCF REF=ACGT, ALT=A → 3 bases deleted at pos 101-103.

        VCF:  pos=100, REF=ACGT, ALT=A  (padding A at 100, CGT deleted 101-103)
        VEP:  chrom 101 103 - 1
        """
        result = format_vep_region("2", 100, "ACGT", "A")
        self.assertEqual(result, "2 101 103 CGT/- 1")

    # ---- Insertions ----

    def test_insertion_single_base(self):
        """Insertion: VCF REF=A, ALT=AG at pos 100 → G inserted after 100.

        VCF:  pos=100, REF=A, ALT=AG  (anchor A at 100, G inserted after it)
        VEP:  chrom 101 100 G 1  (start > end signals insertion)
        """
        result = format_vep_region("3", 100, "A", "AG")
        self.assertEqual(result, "3 101 100 -/G 1")

    def test_insertion_multi_base(self):
        """Insertion: VCF REF=T, ALT=TCCC at pos 200 → CCC inserted after 200.

        VEP:  chrom 201 200 CCC 1
        """
        result = format_vep_region("5", 200, "T", "TCCC")
        self.assertEqual(result, "5 201 200 -/CCC 1")

    # ---- Complex / MNV ----

    def test_complex_substitution(self):
        """Complex: REF=AC, ALT=GT at pos 50 → replace AC with GT.

        No shared prefix, so start=50, end=51.
        VEP:  chrom 50 51 GT 1
        """
        result = format_vep_region("7", 50, "AC", "GT")
        self.assertEqual(result, "7 50 51 AC/GT 1")

    def test_complex_with_shared_prefix(self):
        """Complex: REF=ACG, ALT=ATT at pos 50 → shared prefix A.

        After stripping shared A: REF=CG, ALT=TT, start=51, end=52.
        VEP:  chrom 51 52 TT 1
        """
        result = format_vep_region("7", 50, "ACG", "ATT")
        self.assertEqual(result, "7 51 52 CG/TT 1")

    # ---- Edge cases ----

    def test_empty_ref_raises(self):
        """Empty REF should raise ValueError."""
        with self.assertRaises(ValueError):
            format_vep_region("1", 100, "", "A")

    def test_empty_alt_raises(self):
        """Empty ALT should raise ValueError."""
        with self.assertRaises(ValueError):
            format_vep_region("1", 100, "A", "")


class TestParseVcf(unittest.TestCase):
    """Test VCF parsing."""

    def _write_vcf(self, lines: list[str]) -> str:
        """Write lines to a temp file and return the path."""
        fd, path = tempfile.mkstemp(suffix=".vcf")
        with os.fdopen(fd, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        return path

    def test_basic_parsing(self):
        path = self._write_vcf([
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "1\t100\t.\tA\tG\t50\tPASS\t.",
        ])
        variants = parse_vcf(path)
        os.unlink(path)
        self.assertEqual(len(variants), 1)
        self.assertEqual(variants[0], ("1", 100, "A", "G"))

    def test_multi_allelic_split(self):
        """Multi-allelic ALTs (comma-separated) must be split."""
        path = self._write_vcf([
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "1\t200\t.\tA\tC,T\t50\tPASS\t.",
        ])
        variants = parse_vcf(path)
        os.unlink(path)
        self.assertEqual(len(variants), 2)
        self.assertEqual(variants[0], ("1", 200, "A", "C"))
        self.assertEqual(variants[1], ("1", 200, "A", "T"))

    def test_skips_star_alt(self):
        """Spanning deletion marker (*) should be skipped."""
        path = self._write_vcf([
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "1\t300\t.\tA\t*\t50\tPASS\t.",
        ])
        variants = parse_vcf(path)
        os.unlink(path)
        self.assertEqual(len(variants), 0)

    def test_skips_dot_alt(self):
        """Missing ALT (.) should be skipped."""
        path = self._write_vcf([
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "1\t400\t.\tA\t.\t50\tPASS\t.",
        ])
        variants = parse_vcf(path)
        os.unlink(path)
        self.assertEqual(len(variants), 0)

    def test_header_lines_skipped(self):
        """All lines starting with # are skipped."""
        path = self._write_vcf([
            "##fileformat=VCFv4.2",
            "##INFO=<ID=DP,Number=1,Type=Integer>",
            "##custom_header=something",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "1\t500\t.\tC\tT\t50\tPASS\t.",
        ])
        variants = parse_vcf(path)
        os.unlink(path)
        self.assertEqual(len(variants), 1)

    def test_uppercase_conversion(self):
        """REF and ALT should be uppercased."""
        path = self._write_vcf([
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "1\t600\t.\ta\tg\t50\tPASS\t.",
        ])
        variants = parse_vcf(path)
        os.unlink(path)
        self.assertEqual(variants[0], ("1", 600, "A", "G"))


if __name__ == "__main__":
    unittest.main()
