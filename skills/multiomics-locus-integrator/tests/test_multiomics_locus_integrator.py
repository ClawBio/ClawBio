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
    genes = [("chr1", 50, 200, "GENE_X", 50)]
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
    genes = [("chr1", 0, 50, "GENE_X", 0)]
    result = mli.analyze(loci, accessibility, [], [], genes)
    row = result["loci"][0]
    assert row["tier"] == "unsupported"
    assert row["layers_checked"] == ["accessibility"]
    assert row["layers_supported"] == []


def test_contact_gene_disagreement_is_flagged_not_resolved() -> None:
    loci = [("chr1", 100, 110, "L1")]
    genes = [("chr1", 0, 50, "GENE_NEAR", 0), ("chr1", 1000, 1050, "GENE_FAR", 1000)]
    loops = [("chr1", 90, 120, "chr1", 1010, 1040, "LOOP_1")]
    result = mli.analyze(loci, [], [], loops, genes)
    row = result["loci"][0]
    assert row["target_gene_nearest"] == "GENE_NEAR"
    assert row["target_gene_contact"] == "GENE_FAR"
    assert row["nearest_vs_contact_disagreement"] is True


def test_contact_support_is_independent_of_a_known_gene_at_the_far_anchor() -> None:
    # a locus genuinely sitting inside a real
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
    assert genes == [("chr1", 4000, 5250, "GENE_PLUS", 5000)]  # tss=5000 (BED start, + strand)


def test_gene_targets_tss_window_minus_strand(tmp_path: Path) -> None:
    # '-' strand TSS = BED end - 1 (BED end is exclusive, so the actual
    # last base / TSS is one bp lower); window = [TSS-250+1, TSS+1000+1)
    # -- the +1 on both bounds keeps the TSS base itself classified as
    # the first base of the downstream segment, matching the plus-strand
    # convention (see the window-construction comment in
    # parse_gene_targets for why a direct mirror would be 1bp off).
    bed = tmp_path / "genes.bed"
    bed.write_text("chr1\t5000\t8000\tGENE_MINUS\t-\n")
    genes = mli.parse_gene_targets(bed)
    assert genes == [("chr1", 7750, 9000, "GENE_MINUS", 7999)]  # tss=7999 (BED end-1, - strand)


def test_gene_targets_requires_strand(tmp_path: Path) -> None:
    # Unrelated fixture directory name: this
    # test's own name contains "strand", and tmp_path embeds the test
    # function name, so match="strand" could coincidentally match the
    # path instead of the real message -- same class of issue already
    # fixed for the space-delimited tests in round 4. Verified this one
    # already passed for the right reason (the message genuinely says
    # "...gene_symbol, strand) needs..."), but removing the ambiguity
    # entirely is cheap and matches the established convention.
    bed_dir = tmp_path / "unrelated_fixture_dir"
    bed_dir.mkdir()
    bed = bed_dir / "unrelated_fixture_name.bed"
    bed.write_text("chr1\t5000\t8000\tGENE_NOSTRAND\n")  # only 4 fields, no strand
    with pytest.raises(ValueError, match="strand"):
        mli.parse_gene_targets(bed)


def test_gene_targets_rejects_invalid_strand(tmp_path: Path) -> None:
    # Same unrelated-fixture-name reasoning as above.
    bed_dir = tmp_path / "unrelated_fixture_dir"
    bed_dir.mkdir()
    bed = bed_dir / "unrelated_fixture_name.bed"
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
    genes = [("chr1", 0, 50, "GENE_X", 0)]
    loops = [("chr1", 90, 120, "chr1", 10, 40, "LOOP_1")]
    result = mli.analyze(loci, [], [], loops, genes)
    row = result["loci"][0]
    # contact alone (weight 2 of 2 checked = score 1.0) must NOT reach "high"
    # because only one distinct layer type is supported.
    assert row["score"] == 1.0
    assert row["tier"] == "partial"


def test_all_contacted_genes_reported_not_just_the_first() -> None:
    # gene_contact = contacts[0][0]
    # silently dropped every contacted gene past the first, so
    # nearest_vs_contact_disagreement could read False even when a loop
    # genuinely contacts a different gene than the nearest one.
    loci = [('chr1', 100, 110, 'L1')]
    genes = [('chr1', 50, 100, 'GENE_NEAR', 50), ('chr1', 2000, 2050, 'GENE_FAR', 2000)]
    loops = [
        ('chr1', 90, 120, 'chr1', 50, 100, 'LOOP_TO_NEAR'),
        ('chr1', 90, 120, 'chr1', 2000, 2050, 'LOOP_TO_FAR'),
    ]
    result = mli.analyze(loci, [], [], loops, genes)
    row = result['loci'][0]
    assert 'GENE_NEAR' in row['target_gene_contact']
    assert 'GENE_FAR' in row['target_gene_contact']
    assert row['nearest_vs_contact_disagreement'] is True


def test_parse_bed_space_delimited_gets_a_specific_error(tmp_path: Path) -> None:
    # a space-delimited BED previously raised the
    # generic ">=3 tab-separated fields" count error even though the real
    # problem (wrong delimiter) is detectable and worth naming directly.
    # NOTE: pytest's tmp_path fixture embeds the TEST FUNCTION NAME in its
    # generated directory path, so an earlier version of this test
    # (naively matching the substring "space") passed for the wrong
    # reason -- it matched "..._space_delimited_..." in the tmp_path
    # itself, never the actual exception text. Using an unrelated,
    # unique directory name here closes that loophole.
    bed_dir = tmp_path / "unrelated_fixture_dir"
    bed_dir.mkdir()
    bed = bed_dir / "unrelated_fixture_name.bed"
    bed.write_text("chr1 100 200 PEAK1\n")
    with pytest.raises(ValueError, match="space-separated"):
        mli.parse_bed(bed)


def test_parse_gene_targets_space_delimited_gets_a_specific_error(tmp_path: Path) -> None:
    bed_dir = tmp_path / "unrelated_fixture_dir"
    bed_dir.mkdir()
    bed = bed_dir / "unrelated_fixture_name.bed"
    bed.write_text("chr1 100 200 GENE1 +\n")
    with pytest.raises(ValueError, match="space-separated"):
        mli.parse_gene_targets(bed)


def test_evidence_table_has_no_leading_comment_line(tmp_path: Path) -> None:
    # the round-4 comment line is gone now that
    # the distance column is genuine TSS-point distance (no caveat left
    # to document), and a plain csv.reader -- unlike pandas
    # read_csv(comment="#") -- does not skip "#" lines, a real risk for
    # a chaining partner reading this table.
    output_dir = tmp_path / "out"
    mli.main(["--demo", "--output", str(output_dir)])
    lines = (output_dir / "tables" / "evidence_table.tsv").read_text().splitlines()
    assert lines[0].split("\t")[0] == "locus_id"  # header is the real first line


def test_parse_bedpe_space_delimited_gets_a_specific_error(tmp_path: Path) -> None:
    # parse_bedpe still raised the plain field-count
    # error, missing the space-separated hint parse_bed/parse_gene_targets
    # already give.
    bed_dir = tmp_path / "unrelated_fixture_dir"
    bed_dir.mkdir()
    bedpe = bed_dir / "unrelated_fixture_name.bedpe"
    bedpe.write_text("chr1 90 120 chr1 50 100 LOOP_1\n")
    with pytest.raises(ValueError, match="space-separated"):
        mli.parse_bedpe(bedpe)


def test_nearest_gene_ranks_by_true_tss_not_window_edge(tmp_path: Path) -> None:
    # the window-edge method could rank the WRONG
    # gene as nearest. GENE_A's window edge is closer to this locus than
    # GENE_B's, but GENE_B's true TSS is genuinely closer than GENE_A's.
    #   GENE_A: + strand, TSS=2000, window=[1000, 2250) -> edge dist to 900 is 100
    #   GENE_B: + strand, TSS=500,  window=[0, 750)     -> edge dist to 900 is 150
    #   -> old window-edge ranking would have picked GENE_A (100 < 150), wrong.
    #   -> true TSS distance: GENE_A |900-2000|=1100, GENE_B |900-500|=400
    #   -> correct ranking picks GENE_B.
    bed = tmp_path / "genes.bed"
    bed.write_text(
        "chr1\t2000\t3000\tGENE_A\t+\n"
        "chr1\t500\t1000\tGENE_B\t+\n"
    )
    genes = mli.parse_gene_targets(bed)
    name, dist = mli.nearest_gene("chr1", 900, 901, genes)
    assert name == "GENE_B"
    assert dist == 400


def test_gene_targets_accepts_standard_bed6(tmp_path: Path) -> None:
    # a standard BED6 file (chrom, start, end,
    # name, score, strand) should be readable directly -- strand lives in
    # column 6 there, with a score in column 5 this skill doesn't use.
    bed = tmp_path / "genes.bed"
    bed.write_text("chr1\t5000\t8000\tGENE_BED6\t960\t+\n")
    genes = mli.parse_gene_targets(bed)
    assert genes == [("chr1", 4000, 5250, "GENE_BED6", 5000)]


def test_gene_targets_bed6_minus_strand(tmp_path: Path) -> None:
    bed = tmp_path / "genes.bed"
    bed.write_text("chr1\t5000\t8000\tGENE_BED6\t960\t-\n")
    genes = mli.parse_gene_targets(bed)
    assert genes == [("chr1", 7750, 9000, "GENE_BED6", 7999)]


def test_nearest_gene_flags_when_distant(tmp_path: Path) -> None:
    # a sparse gene-targets file can leave the
    # only candidate gene megabases away; the distance column shows the
    # raw number, but nothing flagged that it should be read with
    # suspicion. target_gene_nearest_is_distant makes that explicit.
    loci = [("chr1", 100, 110, "L1")]
    genes = [("chr1", 2_000_000, 2_000_050, "GENE_FAR", 2_000_025)]  # ~2Mb away
    loops = []
    result = mli.analyze(loci, [], [], loops, genes)
    row = result["loci"][0]
    assert row["target_gene_nearest_distance_bp"] > mli.NEAREST_GENE_FAR_THRESHOLD_BP
    assert row["target_gene_nearest_is_distant"] is True


def test_nearest_gene_not_flagged_when_close() -> None:
    loci = [("chr1", 100, 110, "L1")]
    genes = [("chr1", 0, 50, "GENE_NEAR", 25)]
    result = mli.analyze(loci, [], [], [], genes)
    row = result["loci"][0]
    assert row["target_gene_nearest_is_distant"] is False


def test_layer_coverage_reported_in_summary() -> None:
    # one stray far-off interval can silently
    # stretch a layer's observed span and mask a real coverage gap
    # elsewhere; surfacing the actual span per layer lets a reader spot
    # that directly instead of having to infer it from which loci read
    # not_assayed.
    loci = [("chr1", 100, 101, "L1")]
    accessibility = [("chr1", 50, 150, "peak1")]
    histone = [("chr1", 200, 300, "peak2")]
    loops = [("chr1", 90, 120, "chr1", 500, 600, "LOOP_1")]
    genes = [("chr1", 0, 50, "GENE_X", 0)]
    result = mli.analyze(loci, accessibility, histone, loops, genes)
    coverage = result["summary"]["layer_coverage"]
    assert coverage["accessibility"] == "chr1:50-150"
    assert coverage["histone"] == "chr1:200-300"
    assert coverage["contact"] == "chr1:90-600"


def test_contact_not_assayed_uses_per_loop_span_not_chromosome_wide_merge() -> None:
    # Coverage used one merged span per chromosome
    # across every loop's anchors, so a single distant loop widened the
    # "assayed" region for a locus nowhere near any actual loop. A locus
    # sitting between two unrelated loops on the same chromosome -- inside
    # their combined chromosome-wide bounding box, but outside either
    # loop's own span -- must read not_assayed, not checked-and-unsupported.
    loci = [("chr1", 5000, 5001, "L1")]
    loops = [
        ("chr1", 90, 120, "chr1", 150, 180, "LOOP_NEAR_START"),
        ("chr1", 9000, 9050, "chr1", 9100, 9150, "LOOP_NEAR_END"),
    ]
    genes: list[tuple[str, int, int, str, int]] = []
    result = mli.analyze(loci, [], [], loops, genes)
    row = result["loci"][0]
    assert "contact" in row["layers_not_assayed"]
    assert "contact" not in row["layers_checked"]


def test_contact_still_assayed_within_a_single_loops_own_span() -> None:
    # Sanity check: a locus between the two anchors of ONE loop (not
    # overlapping either anchor directly) is still correctly assayed and
    # unsupported, not not_assayed.
    loci = [("chr1", 500, 501, "L1")]
    loops = [("chr1", 90, 120, "chr1", 900, 950, "LOOP_1")]
    genes: list[tuple[str, int, int, str, int]] = []
    result = mli.analyze(loci, [], [], loops, genes)
    row = result["loci"][0]
    assert "contact" in row["layers_checked"]
    assert "contact" not in row["layers_supported"]


def test_coverage_bed_gives_a_tighter_not_assayed_than_the_track_span() -> None:
    # Default behaviour (no coverage file) treats the whole span between a
    # layer's first and last peak as "assayed" -- see
    # test_true_absence_is_unsupported_not_not_assayed. An explicit
    # --accessibility-coverage-bed states the REAL assayed footprint (e.g.
    # a capture panel's target regions), so gaps between its entries
    # genuinely do mean "never measured" and should read not_assayed, not
    # unsupported, even though the peak track's own span would cover it.
    loci = [("chr1", 5000, 5001, "L1")]
    accessibility = [("chr1", 0, 50, "peak1"), ("chr1", 9000, 9050, "peak2")]
    accessibility_coverage = [("chr1", 0, 100, "panel_a"), ("chr1", 9000, 9100, "panel_b")]
    genes: list[tuple[str, int, int, str, int]] = []
    result = mli.analyze(
        loci, accessibility, [], [], genes,
        accessibility_coverage=accessibility_coverage,
    )
    row = result["loci"][0]
    assert "accessibility" in row["layers_not_assayed"]
    assert "accessibility" not in row["layers_checked"]


def test_coverage_bed_still_assayed_within_its_own_block() -> None:
    # Sanity check: a locus inside one of the coverage file's own entries
    # is assayed, whether or not a peak was actually called there.
    loci = [("chr1", 50, 51, "L1")]
    accessibility = [("chr1", 9000, 9050, "peak_far_away")]
    accessibility_coverage = [("chr1", 0, 100, "panel_a")]
    genes: list[tuple[str, int, int, str, int]] = []
    result = mli.analyze(
        loci, accessibility, [], [], genes,
        accessibility_coverage=accessibility_coverage,
    )
    row = result["loci"][0]
    assert "accessibility" in row["layers_checked"]
    assert "accessibility" not in row["layers_supported"]


def test_coverage_bed_omitted_keeps_default_span_behaviour() -> None:
    # Omitting the coverage argument entirely must reproduce the existing,
    # deliberate default (chromosome-wide track span) exactly --
    # the new parameter must never change behaviour for every existing
    # caller that does not pass it.
    loci = [("chr1", 100, 101, "L1")]
    accessibility = [("chr1", 0, 50, "peak1"), ("chr1", 200, 300, "peak2")]
    genes = [("chr1", 0, 50, "GENE_X", 0)]
    result = mli.analyze(loci, accessibility, [], [], genes)
    row = result["loci"][0]
    assert row["tier"] == "unsupported"


def test_coverage_bed_applies_per_layer_independently() -> None:
    # accessibility_coverage must not affect histone's not_assayed
    # decision, and vice versa -- each layer's coverage override is
    # independent.
    loci = [("chr1", 5000, 5001, "L1")]
    accessibility = [("chr1", 0, 50, "peak1"), ("chr1", 9000, 9050, "peak2")]
    histone = [("chr1", 0, 50, "peak1"), ("chr1", 9000, 9050, "peak2")]
    accessibility_coverage = [("chr1", 0, 100, "panel_a"), ("chr1", 9000, 9100, "panel_b")]
    genes: list[tuple[str, int, int, str, int]] = []
    result = mli.analyze(
        loci, accessibility, histone, [], genes,
        accessibility_coverage=accessibility_coverage,
    )
    row = result["loci"][0]
    assert "accessibility" in row["layers_not_assayed"]
    assert "histone" in row["layers_checked"]  # unaffected, still the default span


def test_chrom_naming_normalizes_chr_prefix_for_overlap() -> None:
    # "5" and "chr5" are the same chromosome under a UCSC-vs-plain naming
    # difference. A locus on "5" against a peak track on "chr5" used to
    # compare unequal and silently read as no overlap at all.
    loci = [("5", 100, 200, "L1")]
    accessibility = [("chr5", 100, 200, "peak1")]
    genes: list[tuple[str, int, int, str, int]] = []
    result = mli.analyze(loci, accessibility, [], [], genes)
    row = result["loci"][0]
    assert "accessibility" in row["layers_supported"]


def test_chrom_naming_normalizes_for_not_assayed_span() -> None:
    # Same naming difference, but for the chromosome-wide span check
    # (not just the direct overlap check): a locus on "5" must be
    # "checked" against an accessibility track whose whole span is on
    # "chr5", not read as not_assayed just because the strings differ.
    loci = [("5", 500, 501, "L1")]
    accessibility = [("chr5", 0, 1000, "peak1")]
    genes: list[tuple[str, int, int, str, int]] = []
    result = mli.analyze(loci, accessibility, [], [], genes)
    row = result["loci"][0]
    assert "accessibility" in row["layers_checked"]
    assert "accessibility" not in row["layers_not_assayed"]


def test_chrom_naming_normalizes_for_contact() -> None:
    loci = [("5", 100, 110, "L1")]
    loops = [("chr5", 90, 120, "chr5", 1000, 1050, "LOOP_1")]
    genes = [("chr5", 1000, 1050, "GENE_X", 1000)]
    result = mli.analyze(loci, [], [], loops, genes)
    row = result["loci"][0]
    assert "contact" in row["layers_supported"]
    assert row["target_gene_contact"] == "GENE_X"


def test_chrom_naming_normalizes_for_nearest_gene() -> None:
    loci = [("5", 100, 101, "L1")]
    genes = [("chr5", 0, 50, "GENE_X", 0)]
    result = mli.analyze(loci, [], [], [], genes)
    row = result["loci"][0]
    assert row["target_gene_nearest"] == "GENE_X"


def test_chrom_naming_normalizes_for_coverage_bed_override() -> None:
    loci = [("5", 50, 51, "L1")]
    accessibility = [("chr5", 9000, 9050, "peak_far_away")]
    accessibility_coverage = [("chr5", 0, 100, "panel_a")]
    genes: list[tuple[str, int, int, str, int]] = []
    result = mli.analyze(
        loci, accessibility, [], [], genes,
        accessibility_coverage=accessibility_coverage,
    )
    row = result["loci"][0]
    assert "accessibility" in row["layers_checked"]


def test_contact_default_span_flag_set_without_coverage_override() -> None:
    loci = [("chr1", 500, 501, "L1")]
    loops = [("chr1", 90, 120, "chr1", 900, 950, "LOOP_1")]
    genes: list[tuple[str, int, int, str, int]] = []
    result = mli.analyze(loci, [], [], loops, genes)
    assert result["summary"]["contact_uses_default_span"] is True


def test_contact_default_span_flag_false_with_coverage_override() -> None:
    loci = [("chr1", 500, 501, "L1")]
    loops = [("chr1", 90, 120, "chr1", 900, 950, "LOOP_1")]
    contact_coverage = [("chr1", 0, 1000, "capture_region")]
    genes: list[tuple[str, int, int, str, int]] = []
    result = mli.analyze(loci, [], [], loops, genes, contact_coverage=contact_coverage)
    assert result["summary"]["contact_uses_default_span"] is False


def test_contact_default_span_flag_false_without_any_loops() -> None:
    loci = [("chr1", 500, 501, "L1")]
    genes: list[tuple[str, int, int, str, int]] = []
    result = mli.analyze(loci, [], [], [], genes)
    assert result["summary"]["contact_uses_default_span"] is False


def test_report_legend_explains_contact_default_span_only_when_default(tmp_path: Path) -> None:
    loci_bed = tmp_path / "loci.bed"
    loci_bed.write_text("chr1\t500\t501\tL1\n")
    loops_bedpe = tmp_path / "loops.bedpe"
    loops_bedpe.write_text("chr1\t90\t120\tchr1\t900\t950\tLOOP_1\n")
    genes_bed = tmp_path / "genes.bed"
    genes_bed.write_text("chr1\t0\t50\tGENE_X\t+\n")
    output_dir = tmp_path / "out"

    mli.main([
        "--loci", str(loci_bed),
        "--contact-evidence", str(loops_bedpe),
        "--gene-targets", str(genes_bed),
        "--output", str(output_dir),
    ])
    report = (output_dir / "report.md").read_text()
    assert "WHOLE region between a loop's two anchors" in report


def test_contact_coverage_bed_overrides_the_per_loop_span() -> None:
    # Default (per-loop-span) behaviour: a locus between two unrelated
    # loops reads not_assayed -- see
    # test_contact_not_assayed_uses_per_loop_span_not_chromosome_wide_merge.
    # An explicit --contact-coverage-bed states the real capture/assayed
    # footprint of the 3D-contact assay (e.g. a Capture-C target region),
    # so a gap inside it is a true negative (checked, unsupported), not
    # missing data, even though no loop was called there specifically.
    loci = [("chr1", 5000, 5001, "L1")]
    loops = [
        ("chr1", 90, 120, "chr1", 150, 180, "LOOP_NEAR_START"),
        ("chr1", 9000, 9050, "chr1", 9100, 9150, "LOOP_NEAR_END"),
    ]
    contact_coverage = [("chr1", 0, 10000, "capture_region")]
    genes: list[tuple[str, int, int, str, int]] = []
    result = mli.analyze(loci, [], [], loops, genes, contact_coverage=contact_coverage)
    row = result["loci"][0]
    assert "contact" in row["layers_checked"]
    assert "contact" not in row["layers_not_assayed"]
    assert "contact" not in row["layers_supported"]


def test_contact_coverage_bed_omitted_keeps_default_per_loop_behaviour() -> None:
    loci = [("chr1", 5000, 5001, "L1")]
    loops = [
        ("chr1", 90, 120, "chr1", 150, 180, "LOOP_NEAR_START"),
        ("chr1", 9000, 9050, "chr1", 9100, 9150, "LOOP_NEAR_END"),
    ]
    genes: list[tuple[str, int, int, str, int]] = []
    result = mli.analyze(loci, [], [], loops, genes)
    row = result["loci"][0]
    assert "contact" in row["layers_not_assayed"]


def test_contact_coverage_bed_cli_flag_produces_the_tighter_check(tmp_path: Path) -> None:
    loci_bed = tmp_path / "loci.bed"
    loci_bed.write_text("chr1\t5000\t5001\tL1\n")
    loops_bedpe = tmp_path / "loops.bedpe"
    loops_bedpe.write_text(
        "chr1\t90\t120\tchr1\t150\t180\tLOOP_NEAR_START\n"
        "chr1\t9000\t9050\tchr1\t9100\t9150\tLOOP_NEAR_END\n"
    )
    coverage_bed = tmp_path / "capture.bed"
    coverage_bed.write_text("chr1\t0\t10000\tcapture_region\n")
    genes_bed = tmp_path / "genes.bed"
    genes_bed.write_text("chr1\t0\t50\tGENE_X\t+\n")
    output_dir = tmp_path / "out"

    rc = mli.main([
        "--loci", str(loci_bed),
        "--contact-evidence", str(loops_bedpe),
        "--contact-coverage-bed", str(coverage_bed),
        "--gene-targets", str(genes_bed),
        "--output", str(output_dir),
    ])
    assert rc == 0
    result = json.loads((output_dir / "result.json").read_text())
    row = result["loci"][0]
    assert "contact" in row["layers_checked"]
    assert "contact" not in row["layers_not_assayed"]


def test_chrom_naming_normalizes_mitochondrial_mt_vs_chrm() -> None:
    # "MT" (Ensembl/RefSeq-style, no "chr" prefix) and "chrM" (UCSC-style)
    # name the same chromosome, but stripping a leading "chr" alone
    # leaves "chrM" -> "M" while "MT" stays "MT" -- still unequal. Both
    # must normalize to the same key.
    loci = [("MT", 100, 200, "L1")]
    accessibility = [("chrM", 100, 200, "peak1")]
    genes: list[tuple[str, int, int, str, int]] = []
    result = mli.analyze(loci, accessibility, [], [], genes)
    row = result["loci"][0]
    assert "accessibility" in row["layers_supported"]


def test_tss_window_treats_the_tss_base_consistently_across_strands(tmp_path: Path) -> None:
    # The TSS base itself must be classified as the first base of the
    # DOWNSTREAM segment on both strands -- not downstream on one strand
    # and upstream on the other, which a direct mirror of the plus-strand
    # formula would do. Locus = exactly the TSS base on each gene.
    bed = tmp_path / "genes.bed"
    bed.write_text(
        "chr1\t5000\t8000\tGENE_PLUS\t+\n"
        "chr1\t5000\t8000\tGENE_MINUS\t-\n"
    )
    genes = mli.parse_gene_targets(bed)
    plus_window = next(g for g in genes if g[3] == "GENE_PLUS")
    minus_window = next(g for g in genes if g[3] == "GENE_MINUS")
    plus_tss = plus_window[4]
    minus_tss = minus_window[4]
    # Both windows' total width must still be UP+DOWN=1250, unaffected
    # by which base is classified which way.
    assert plus_window[2] - plus_window[1] == 1250
    assert minus_window[2] - minus_window[1] == 1250
    # On the plus strand, the TSS base starts the downstream segment:
    # window_start is exactly UP bases before tss.
    assert plus_tss - plus_window[1] == 1000
    # On the minus strand, after the fix, the TSS base also starts its
    # strand's downstream segment: window_start is exactly DOWN-1 bases
    # before tss (not DOWN, which would put tss one base into upstream).
    assert minus_tss - minus_window[1] == 249


def test_multiomics_locus_cli_registry_allows_all_coverage_bed_flags() -> None:
    # clawbio.cli's "multiomics-locus" alias must pass --accessibility-
    # coverage-bed / --histone-coverage-bed / --contact-coverage-bed
    # through "clawbio run" the same way it already does for the tracks
    # they override -- missing from either set here means the runner
    # either rejects the flag outright or passes it through as a bare
    # string instead of resolving it as a path.
    import importlib
    cli = importlib.import_module("clawbio.cli")
    entry = cli.SKILLS["multiomics-locus"]
    for flag in ("--accessibility-coverage-bed", "--histone-coverage-bed", "--contact-coverage-bed"):
        assert flag in entry["allowed_extra_flags"], flag
        assert flag in entry["extra_path_flags"], flag


def test_coverage_bed_cli_flags_produce_the_tighter_not_assayed(tmp_path: Path) -> None:
    loci_bed = tmp_path / "loci.bed"
    loci_bed.write_text("chr1\t5000\t5001\tL1\n")
    accessibility_bed = tmp_path / "accessibility.bed"
    accessibility_bed.write_text("chr1\t0\t50\tpeak1\nchr1\t9000\t9050\tpeak2\n")
    coverage_bed = tmp_path / "coverage.bed"
    coverage_bed.write_text("chr1\t0\t100\tpanel_a\nchr1\t9000\t9100\tpanel_b\n")
    genes_bed = tmp_path / "genes.bed"
    genes_bed.write_text("chr1\t0\t50\tGENE_X\t+\n")
    output_dir = tmp_path / "out"

    rc = mli.main([
        "--loci", str(loci_bed),
        "--accessibility-track", str(accessibility_bed),
        "--accessibility-coverage-bed", str(coverage_bed),
        "--gene-targets", str(genes_bed),
        "--output", str(output_dir),
    ])
    assert rc == 0
    result = json.loads((output_dir / "result.json").read_text())
    row = result["loci"][0]
    assert "accessibility" in row["layers_not_assayed"]


def test_gene_targets_bed6_unstranded_names_the_right_column(tmp_path: Path) -> None:
    # A BED6 line with strand "." in column 6 fell
    # back to column 5 (the score column) as if IT were strand, so the
    # error named a nonsense value ("960") with no indication it came
    # from the wrong column. Column 6 is strand whenever the line has
    # >=6 fields (BED6-shaped); falling back to column 5 is only correct
    # for this skill's OWN 5-field format, never as a second guess after
    # a 6+-field line's real strand column fails to parse.
    bed_dir = tmp_path / "unrelated_fixture_dir"
    bed_dir.mkdir()
    bed = bed_dir / "unrelated_fixture_name.bed"
    bed.write_text("chr1\t5000\t8000\tGENE\t960\t.\n")
    with pytest.raises(ValueError, match="column 6"):
        mli.parse_gene_targets(bed)
