"""Tests for pipeline.py -- 15 tests across 7 test classes."""

from __future__ import annotations

import sys
from math import ceil
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import pipeline as pl

SKILL_DIR = Path(__file__).resolve().parent.parent
DEMO_ARCHAIC = SKILL_DIR / "examples" / "demo_archaic.vcf"


# -----------------------------------------------------------------------
# Fixtures
# -----------------------------------------------------------------------


@pytest.fixture
def eigenstrat_dir(tmp_path):
    """Create a minimal EIGENSTRAT dataset with 5 individuals, 20 SNPs.

    10 SNPs on chr1 (positions 100..1000, step 100)
    10 SNPs on chr22 (positions 16050000..16050900, step 100)
    Text format .geno.
    """
    # .ind file
    ind_path = tmp_path / "test.ind"
    ind_lines = [
        "IND01 M PopA",
        "IND02 F PopA",
        "IND03 M PopB",
        "IND04 F PopB",
        "IND05 M PopC",
    ]
    ind_path.write_text("\n".join(ind_lines) + "\n")

    # .snp file
    snp_path = tmp_path / "test.snp"
    snp_lines = []
    for i in range(10):
        pos = (i + 1) * 100
        snp_lines.append(f"rs{i+1:04d} 1 0.0 {pos} A G")
    for i in range(10):
        pos = 16050000 + i * 100
        snp_lines.append(f"rs{i+11:04d} 22 0.0 {pos} C T")
    snp_path.write_text("\n".join(snp_lines) + "\n")

    # .geno file (text format: each line is one SNP, digits = samples)
    # 5 individuals, 20 SNPs
    np.random.seed(42)
    geno_path = tmp_path / "test.geno"
    lines = []
    for _ in range(20):
        row = "".join(str(np.random.choice([0, 1, 2])) for _ in range(5))
        lines.append(row)
    geno_path.write_text("\n".join(lines) + "\n")

    return {
        "ind": ind_path,
        "snp": snp_path,
        "geno": geno_path,
        "n_samples": 5,
        "n_snps": 20,
    }


@pytest.fixture
def archaic_vcf(tmp_path):
    """Create a minimal archaic VCF matching chr22 positions from eigenstrat_dir."""
    vcf_path = tmp_path / "archaic.vcf"
    header = (
        "##fileformat=VCFv4.1\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNEANDERTAL\n"
    )
    records = []
    for i in range(10):
        pos = 16050000 + i * 100
        gt = np.random.choice(["0/0", "0/1", "1/1"])
        records.append(f"22\t{pos}\trs{i+11:04d}\tC\tT\t100\tPASS\t.\tGT\t{gt}")
    vcf_path.write_text(header + "\n".join(records) + "\n")
    return vcf_path


# -----------------------------------------------------------------------
# TestReadEigenstratInd (2 tests)
# -----------------------------------------------------------------------


class TestReadEigenstratInd:
    """Test read_eigenstrat_ind()."""

    def test_ind_count(self, eigenstrat_dir):
        """Reads correct number of individuals."""
        inds = pl.read_eigenstrat_ind(eigenstrat_dir["ind"])
        assert len(inds) == 5

    def test_ind_fields(self, eigenstrat_dir):
        """Each individual has name, sex, population."""
        inds = pl.read_eigenstrat_ind(eigenstrat_dir["ind"])
        assert inds[0]["name"] == "IND01"
        assert inds[0]["sex"] == "M"
        assert inds[0]["population"] == "PopA"
        assert inds[3]["sex"] == "F"
        assert inds[3]["population"] == "PopB"


# -----------------------------------------------------------------------
# TestReadEigenstratSnp (3 tests)
# -----------------------------------------------------------------------


class TestReadEigenstratSnp:
    """Test read_eigenstrat_snp()."""

    def test_snp_count_all(self, eigenstrat_dir):
        """Reads all 20 SNPs when no chromosome filter."""
        snps = pl.read_eigenstrat_snp(eigenstrat_dir["snp"])
        assert len(snps) == 20

    def test_snp_count_filtered(self, eigenstrat_dir):
        """Reads 10 SNPs when filtered to chr1."""
        snps = pl.read_eigenstrat_snp(eigenstrat_dir["snp"], chrom="1")
        assert len(snps) == 10

    def test_list_chromosomes(self, eigenstrat_dir):
        """list_chromosomes returns chroms in order of first appearance."""
        chroms = pl.list_chromosomes(eigenstrat_dir["snp"])
        assert chroms == ["1", "22"]


# -----------------------------------------------------------------------
# TestReadEigenstratGeno (3 tests)
# -----------------------------------------------------------------------


class TestReadEigenstratGeno:
    """Test read_eigenstrat_geno() for text and binary formats."""

    def test_text_geno_shape(self, eigenstrat_dir):
        """Text .geno returns correct shape."""
        genos = pl.read_eigenstrat_geno(
            eigenstrat_dir["geno"],
            eigenstrat_dir["n_samples"],
            eigenstrat_dir["n_snps"],
        )
        assert genos.shape == (5, 20)

    def test_text_geno_values(self, eigenstrat_dir):
        """Text .geno values are in {0, 1, 2, -1}."""
        genos = pl.read_eigenstrat_geno(
            eigenstrat_dir["geno"],
            eigenstrat_dir["n_samples"],
            eigenstrat_dir["n_snps"],
        )
        unique = set(np.unique(genos))
        assert unique.issubset({-1, 0, 1, 2})

    def test_binary_geno_roundtrip(self, tmp_path):
        """Binary packed .geno reads back correctly."""
        n_samples = 7
        n_snps = 5
        # Create known genotype data
        np.random.seed(99)
        expected = np.random.choice([0, 1, 2], size=(n_samples, n_snps)).astype(np.int8)

        # Write binary format
        bytes_per_row = ceil(n_samples / 4)
        bin_path = tmp_path / "test_binary.geno"
        with open(bin_path, "wb") as fh:
            # Magic header
            fh.write(b"GENO")

            # Write each SNP row
            for snp_idx in range(n_snps):
                row_bytes = bytearray(bytes_per_row)
                for s_idx in range(n_samples):
                    byte_idx = s_idx // 4
                    bit_pos = s_idx % 4
                    geno = int(expected[s_idx, snp_idx])
                    row_bytes[byte_idx] |= (geno & 0x03) << (bit_pos * 2)
                fh.write(bytes(row_bytes))

        result = pl.read_eigenstrat_geno(bin_path, n_samples, n_snps)
        np.testing.assert_array_equal(result, expected)


# -----------------------------------------------------------------------
# TestExtractArchaicGenotypes (2 tests)
# -----------------------------------------------------------------------


class TestExtractArchaicGenotypes:
    """Test extract_archaic_genotypes()."""

    def test_extract_from_demo_vcf(self):
        """Extracts genotypes from demo archaic VCF."""
        genos = pl.extract_archaic_genotypes(DEMO_ARCHAIC)
        assert len(genos) == 10
        # All values should be valid genotypes
        for pos, gt in genos.items():
            assert gt in {0, 1, 2, -1}

    def test_extract_with_position_filter(self):
        """Filtering to specific positions works."""
        genos = pl.extract_archaic_genotypes(
            DEMO_ARCHAIC, positions=[16050075, 16050115]
        )
        assert len(genos) == 2
        assert 16050075 in genos
        assert 16050115 in genos


# -----------------------------------------------------------------------
# TestBuildGenotypeTable (2 tests)
# -----------------------------------------------------------------------


class TestBuildGenotypeTable:
    """Test build_genotype_table()."""

    def test_table_shape(self, eigenstrat_dir, archaic_vcf):
        """Table has n_samples + 1 rows (archaic + modern)."""
        table, positions, names = pl.build_genotype_table(
            eigenstrat_dir["ind"],
            eigenstrat_dir["snp"],
            eigenstrat_dir["geno"],
            archaic_vcf,
            chrom="22",
        )
        assert table.shape[0] == 6  # 1 archaic + 5 modern
        assert len(names) == 5

    def test_table_positions_shared(self, eigenstrat_dir, archaic_vcf):
        """Returned positions are the intersection of modern and archaic."""
        table, positions, names = pl.build_genotype_table(
            eigenstrat_dir["ind"],
            eigenstrat_dir["snp"],
            eigenstrat_dir["geno"],
            archaic_vcf,
            chrom="22",
        )
        # positions should be a subset of both modern and archaic positions
        assert len(positions) <= 10  # At most 10 shared


# -----------------------------------------------------------------------
# TestRunIBDmixChromosome (1 test)
# -----------------------------------------------------------------------


class TestRunIBDmixChromosome:
    """Test run_ibdmix_on_table()."""

    def test_returns_segments(self, eigenstrat_dir, archaic_vcf):
        """run_ibdmix_on_table returns IntrogressionSegment list."""
        table, positions, names = pl.build_genotype_table(
            eigenstrat_dir["ind"],
            eigenstrat_dir["snp"],
            eigenstrat_dir["geno"],
            archaic_vcf,
            chrom="22",
        )
        if len(positions) > 0:
            segments = pl.run_ibdmix_on_table(
                table, positions, names, chrom="22", lod_threshold=0.5
            )
            assert isinstance(segments, list)
            for seg in segments:
                assert hasattr(seg, "sample")
                assert hasattr(seg, "chrom")
                assert seg.chrom == "22"


# -----------------------------------------------------------------------
# TestPipelineChromosome (2 tests)
# -----------------------------------------------------------------------


class TestPipelineChromosome:
    """Test run_chromosome() end-to-end."""

    def test_run_chromosome_returns_list(self, eigenstrat_dir, archaic_vcf):
        """run_chromosome returns a list."""
        segments = pl.run_chromosome(
            eigenstrat_dir["ind"],
            eigenstrat_dir["snp"],
            eigenstrat_dir["geno"],
            archaic_vcf,
            chrom="22",
            lod_threshold=0.5,
        )
        assert isinstance(segments, list)

    def test_run_chromosome_no_overlap(self, eigenstrat_dir, archaic_vcf):
        """Querying a chromosome with no data returns empty list."""
        segments = pl.run_chromosome(
            eigenstrat_dir["ind"],
            eigenstrat_dir["snp"],
            eigenstrat_dir["geno"],
            archaic_vcf,
            chrom="99",  # nonexistent chromosome
            lod_threshold=0.5,
        )
        assert segments == []
