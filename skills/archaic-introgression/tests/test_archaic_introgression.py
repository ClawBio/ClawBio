"""Tests for archaic_introgression.py -- 24 tests across 10 test classes."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import archaic_introgression as ai

SKILL_DIR = Path(__file__).resolve().parent.parent
DEMO_MODERN = SKILL_DIR / "examples" / "demo_modern.vcf"
DEMO_ARCHAIC = SKILL_DIR / "examples" / "demo_archaic.vcf"


# -----------------------------------------------------------------------
# TestToolDiscovery (4 tests)
# -----------------------------------------------------------------------


class TestToolDiscovery:
    """Test find_tool() and list_available_methods()."""

    def test_find_tool_nonexistent(self):
        """A tool that does not exist returns None."""
        assert ai.find_tool("definitely_not_a_real_tool_xyz") is None

    def test_find_tool_returns_string_or_none(self):
        """find_tool always returns str or None."""
        result = ai.find_tool("python3")
        assert result is None or isinstance(result, str)

    def test_find_tool_sprime_jar_missing(self):
        """When SPRIME_JAR points to a nonexistent file and sprime not on PATH."""
        with patch.object(ai, "SPRIME_JAR", "/nonexistent/sprime.jar"):
            with patch("shutil.which", return_value=None):
                assert ai.find_tool("sprime") is None

    def test_list_available_methods_always_has_fallback(self):
        """list_available_methods includes ibdmix_fallback when no binary exists."""
        with patch("shutil.which", return_value=None):
            with patch.object(ai, "SPRIME_JAR", "/nonexistent/sprime.jar"):
                methods = ai.list_available_methods()
                assert "ibdmix_fallback" in methods


# -----------------------------------------------------------------------
# TestVCFParsing (5 tests)
# -----------------------------------------------------------------------


class TestVCFParsing:
    """Test VCF parsing functions."""

    def test_parse_vcf_samples_modern(self):
        """Modern VCF has 3 samples."""
        samples = ai.parse_vcf_samples(DEMO_MODERN)
        assert samples == ["SAMPLE01", "SAMPLE02", "SAMPLE03"]

    def test_parse_vcf_samples_archaic(self):
        """Archaic VCF has 1 sample."""
        samples = ai.parse_vcf_samples(DEMO_ARCHAIC)
        assert samples == ["NEANDERTAL"]

    def test_parse_vcf_variants_count(self):
        """Demo VCF has 10 variants."""
        variants = ai.parse_vcf_variants(DEMO_MODERN)
        assert len(variants) == 10

    def test_parse_vcf_variants_fields(self):
        """Variant records contain expected fields."""
        variants = ai.parse_vcf_variants(DEMO_MODERN)
        first = variants[0]
        assert first["chrom"] == "chr22"
        assert first["pos"] == "16050075"
        assert first["ref"] == "A"
        assert first["alt"] == "G"
        assert "AA=A" in first["info"]

    def test_parse_vcf_genotypes_shape(self):
        """Genotype array has shape (3 samples, 10 variants)."""
        genos = ai.parse_vcf_genotypes(DEMO_MODERN)
        assert genos.shape == (3, 10)


# -----------------------------------------------------------------------
# TestSharedVariants (2 tests)
# -----------------------------------------------------------------------


class TestSharedVariants:
    """Test shared position detection."""

    def test_shared_positions_count(self):
        """All 10 positions are shared between modern and archaic demo VCFs."""
        shared = ai.get_shared_positions(DEMO_MODERN, DEMO_ARCHAIC)
        assert len(shared) == 10

    def test_shared_positions_sorted(self):
        """Shared positions are sorted numerically."""
        shared = ai.get_shared_positions(DEMO_MODERN, DEMO_ARCHAIC)
        assert shared == sorted(shared)


# -----------------------------------------------------------------------
# TestIBDmixRunner (2 tests)
# -----------------------------------------------------------------------


class TestIBDmixRunner:
    """Test IBDmix runner (pure-Python fallback)."""

    def test_run_ibdmix_returns_segments(self, tmp_path):
        """run_ibdmix returns a list of IntrogressionSegment objects."""
        segments = ai.run_ibdmix(DEMO_MODERN, DEMO_ARCHAIC, tmp_path, lod_threshold=1.0)
        assert isinstance(segments, list)
        for seg in segments:
            assert isinstance(seg, ai.IntrogressionSegment)

    def test_run_ibdmix_with_sample_filter(self, tmp_path):
        """Filtering to one sample returns segments only for that sample."""
        segments = ai.run_ibdmix(
            DEMO_MODERN, DEMO_ARCHAIC, tmp_path,
            lod_threshold=1.0, samples=["SAMPLE01"]
        )
        for seg in segments:
            assert seg.sample == "SAMPLE01"


# -----------------------------------------------------------------------
# TestSprimeRunner (1 test)
# -----------------------------------------------------------------------


class TestSprimeRunner:
    """Test Sprime command builder."""

    def test_build_sprime_no_binary(self):
        """build_sprime_command returns empty list when tool not found."""
        with patch.object(ai, "SPRIME_JAR", "/nonexistent/sprime.jar"):
            with patch("shutil.which", return_value=None):
                cmd = ai.build_sprime_command("modern.vcf", "outgroup.vcf", "out")
                assert cmd == []


# -----------------------------------------------------------------------
# TestHmmixRunner (1 test)
# -----------------------------------------------------------------------


class TestHmmixRunner:
    """Test hmmix command builder."""

    def test_build_hmmix_no_binary(self):
        """build_hmmix_command returns empty list when tool not found."""
        with patch("shutil.which", return_value=None):
            cmd = ai.build_hmmix_command("modern.vcf", "/tmp/output")
            assert cmd == []


# -----------------------------------------------------------------------
# TestIntrogressionSegment (3 tests)
# -----------------------------------------------------------------------


class TestIntrogressionSegment:
    """Test IntrogressionSegment dataclass."""

    def test_length_property(self):
        """Length is end - start."""
        seg = ai.IntrogressionSegment(
            sample="S1", chrom="chr1", start=1000, end=2000,
            score=3.5, num_variants=5,
        )
        assert seg.length == 1000

    def test_to_dict_includes_length(self):
        """to_dict() includes the computed length field."""
        seg = ai.IntrogressionSegment(
            sample="S1", chrom="chr1", start=1000, end=2000,
            score=3.5, num_variants=5,
        )
        d = seg.to_dict()
        assert d["length"] == 1000
        assert d["sample"] == "S1"
        assert d["chrom"] == "chr1"

    def test_to_bed_format(self):
        """to_bed() returns a valid BED line with 0-based start."""
        seg = ai.IntrogressionSegment(
            sample="S1", chrom="chr1", start=1000, end=2000,
            score=3.5, num_variants=5,
        )
        bed = seg.to_bed()
        fields = bed.split("\t")
        assert fields[0] == "chr1"
        assert fields[1] == "999"  # 0-based
        assert fields[2] == "2000"
        assert fields[3] == "S1"
        assert fields[4] == "3.50"


# -----------------------------------------------------------------------
# TestSummaryStats (1 test)
# -----------------------------------------------------------------------


class TestSummaryStats:
    """Test compute_summary()."""

    def test_summary_from_segments(self):
        """Summary computes correct per-sample statistics."""
        segments = [
            ai.IntrogressionSegment(
                sample="S1", chrom="chr1", start=1000, end=2000,
                score=4.0, num_variants=5,
            ),
            ai.IntrogressionSegment(
                sample="S1", chrom="chr1", start=3000, end=5000,
                score=6.0, num_variants=8,
            ),
            ai.IntrogressionSegment(
                sample="S2", chrom="chr1", start=1000, end=1500,
                score=3.5, num_variants=3,
            ),
        ]
        summary = ai.compute_summary(segments)
        assert summary["S1"]["total_segments"] == 2
        assert summary["S1"]["total_length_bp"] == 3000
        assert summary["S1"]["mean_segment_length"] == 1500.0
        assert summary["S2"]["total_segments"] == 1
        assert summary["S2"]["total_length_bp"] == 500


# -----------------------------------------------------------------------
# TestCLI (3 tests)
# -----------------------------------------------------------------------


class TestCLI:
    """Test CLI entry point."""

    def test_demo_mode(self, tmp_path):
        """--demo runs successfully and writes output files."""
        outdir = tmp_path / "demo_out"
        rc = ai.main(["--demo", "--output", str(outdir)])
        assert rc == 0
        assert (outdir / "introgression_results.json").exists()
        assert (outdir / "segments.bed").exists()

    def test_demo_json_structure(self, tmp_path):
        """Demo JSON output has expected structure."""
        outdir = tmp_path / "demo_json"
        ai.main(["--demo", "--output", str(outdir)])
        with open(outdir / "introgression_results.json") as fh:
            data = json.load(fh)
        assert "method" in data
        assert "segments" in data
        assert "summary" in data
        assert "lod_threshold" in data
        assert "num_samples" in data

    def test_missing_input_returns_error(self):
        """Missing --input without --demo causes an error."""
        with pytest.raises(SystemExit):
            ai.main(["--archaic", "foo.vcf"])


# -----------------------------------------------------------------------
# TestEdgeCases (2 tests)
# -----------------------------------------------------------------------


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_vcf(self, tmp_path):
        """A VCF with only a header returns no samples."""
        empty_vcf = tmp_path / "empty.vcf"
        empty_vcf.write_text(
            "##fileformat=VCFv4.1\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"
        )
        samples = ai.parse_vcf_samples(empty_vcf)
        assert samples == []

    def test_genotypes_missing_values(self, tmp_path):
        """Missing genotypes (./.) are encoded as -1."""
        vcf = tmp_path / "missing.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.1\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr1\t100\t.\tA\tG\t100\tPASS\t.\tGT\t./.\n"
        )
        genos = ai.parse_vcf_genotypes(vcf)
        assert genos.shape == (1, 1)
        assert genos[0, 0] == -1
