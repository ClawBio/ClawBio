from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import multiomics_locus_integrator as mli  # noqa: E402


def test_demo_cli_writes_the_documented_output_contract(tmp_path: Path) -> None:
    output_dir = tmp_path / "out"
    rc = mli.main(["--demo", "--output", str(output_dir)])
    assert rc == 0
    assert (output_dir / "report.md").exists()
    assert (output_dir / "result.json").exists()
    assert (output_dir / "tables" / "evidence_table.tsv").exists()
    assert (output_dir / "reproducibility" / "commands.sh").exists()
    assert (output_dir / "reproducibility" / "environment.yml").exists()
    assert (output_dir / "reproducibility" / "checksums.sha256").exists()


def test_demo_covers_every_tier(tmp_path: Path) -> None:
    output_dir = tmp_path / "out"
    mli.main(["--demo", "--output", str(output_dir)])
    result = json.loads((output_dir / "result.json").read_text())
    tiers = {r["tier"] for r in result["loci"]}
    assert tiers == {"high", "partial", "unsupported", "not_assayed"}


def test_not_assayed_is_distinct_from_unsupported() -> None:
    loci = [("chr1", 100, 101, "L1")]
    accessibility = [("chr1", 1000, 1100, "peak1")]  # far from the locus, defines a narrow span
    genes = [("chr1", 50, 200, "GENE_X")]
    result = mli.analyze(loci, accessibility, [], [], genes)
    row = result["loci"][0]
    assert row["tier"] == "not_assayed"
    assert "accessibility" in row["layers_not_assayed"]
    assert row["layers_checked"] == []


def test_true_absence_is_unsupported_not_not_assayed() -> None:
    loci = [("chr1", 100, 101, "L1")]
    # Two peaks on either side of the locus: their combined span covers the
    # locus (so it is "assayed"), but neither overlaps it directly.
    accessibility = [("chr1", 0, 50, "peak1"), ("chr1", 200, 300, "peak2")]
    genes = [("chr1", 0, 50, "GENE_X")]
    result = mli.analyze(loci, accessibility, [], [], genes)
    row = result["loci"][0]
    assert row["tier"] == "unsupported"
    assert row["layers_checked"] == ["accessibility"]
    assert row["layers_supported"] == []


def test_contact_gene_disagreement_is_flagged_not_resolved() -> None:
    loci = [("chr1", 100, 110, "L1")]
    genes = [("chr1", 0, 50, "GENE_NEAR"), ("chr1", 1000, 1050, "GENE_FAR")]
    loops = [("chr1", 90, 120, "chr1", 1010, 1040, "LOOP_1")]
    result = mli.analyze(loci, [], [], loops, genes)
    row = result["loci"][0]
    assert row["target_gene_nearest"] == "GENE_NEAR"
    assert row["target_gene_contact"] == "GENE_FAR"
    assert row["nearest_vs_contact_disagreement"] is True


def test_contact_support_is_independent_of_a_known_gene_at_the_far_anchor() -> None:
    # PR #678 review (manuelcorpas): a locus genuinely sitting inside a real
    # loop anchor must count as contact-supported even when --gene-targets is
    # too sparse to name the gene at the loop's far end. "Is there a loop
    # here?" and "do we know the gene?" are two separate questions; a thin
    # gene list must never hide real chromatin-contact evidence.
    loci = [("chr1", 100, 110, "L1")]
    genes: list[tuple[str, int, int, str]] = []  # sparse: nothing here at all
    loops = [("chr1", 90, 120, "chr1", 1000, 1050, "LOOP_1")]
    result = mli.analyze(loci, [], [], loops, genes)
    row = result["loci"][0]
    assert "contact" in row["layers_supported"]
    assert row["tier"] != "unsupported"
    assert row["target_gene_contact"] is None  # gene genuinely unknown, reported separately


def test_gene_targets_tss_window_plus_strand(tmp_path: Path) -> None:
    # '+' strand TSS = BED start; window = [TSS-1000, TSS+250).
    bed = tmp_path / "genes.bed"
    bed.write_text("chr1\t5000\t8000\tGENE_PLUS\t+\n")
    genes = mli.parse_gene_targets(bed)
    assert genes == [("chr1", 4000, 5250, "GENE_PLUS")]


def test_gene_targets_tss_window_minus_strand(tmp_path: Path) -> None:
    # '-' strand TSS = BED end; window = [TSS-250, TSS+1000).
    bed = tmp_path / "genes.bed"
    bed.write_text("chr1\t5000\t8000\tGENE_MINUS\t-\n")
    genes = mli.parse_gene_targets(bed)
    assert genes == [("chr1", 7750, 9000, "GENE_MINUS")]


def test_gene_targets_requires_strand(tmp_path: Path) -> None:
    bed = tmp_path / "genes.bed"
    bed.write_text("chr1\t5000\t8000\tGENE_NOSTRAND\n")  # only 4 fields, no strand
    with pytest.raises(ValueError, match="strand"):
        mli.parse_gene_targets(bed)


def test_gene_targets_rejects_invalid_strand(tmp_path: Path) -> None:
    bed = tmp_path / "genes.bed"
    bed.write_text("chr1\t5000\t8000\tGENE_BAD\t.\n")
    with pytest.raises(ValueError, match="strand"):
        mli.parse_gene_targets(bed)


def test_gene_targets_window_never_goes_negative(tmp_path: Path) -> None:
    bed = tmp_path / "genes.bed"
    bed.write_text("chr1\t500\t2000\tGENE_NEAR_START\t+\n")  # TSS=500, 500-1000 would be negative
    genes = mli.parse_gene_targets(bed)
    assert genes[0][1] == 0  # clamped, not negative


def test_high_tier_requires_two_distinct_supported_layers() -> None:
    loci = [("chr1", 100, 110, "L1")]
    genes = [("chr1", 0, 50, "GENE_X")]
    loops = [("chr1", 90, 120, "chr1", 10, 40, "LOOP_1")]
    result = mli.analyze(loci, [], [], loops, genes)
    row = result["loci"][0]
    # contact alone (weight 2 of 2 checked = score 1.0) must NOT reach "high"
    # because only one distinct layer type is supported.
    assert row["score"] == 1.0
    assert row["tier"] == "partial"
