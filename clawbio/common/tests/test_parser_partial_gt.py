"""Partial VCF GT values must not become callable haploid genotypes."""
from __future__ import annotations

from pathlib import Path

from clawbio.common.parsers import detect_format, genotypes_to_simple, parse_genetic_file


def test_public_vcf_parser_drops_partial_missing_calls_and_keeps_complete_ploidies(
    tmp_path: Path,
) -> None:
    path = tmp_path / "calls.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "1\t1\trs_partial_ref\tC\tT\t.\t.\t.\tGT\t0/.\n"
        "1\t2\trs_partial_alt\tC\tT\t.\t.\t.\tGT\t./1\n"
        "1\t3\trs_missing\tC\tT\t.\t.\t.\tGT\t./.\n"
        "1\t4\trs_haploid\tC\tT\t.\t.\t.\tGT\t1\n"
        "1\t5\trs_polyploid\tC\tT\t.\t.\t.\tGT\t0/1/1\n",
        encoding="utf-8",
    )
    assert detect_format(path) == "vcf"
    calls = genotypes_to_simple(parse_genetic_file(path))
    assert "rs_partial_ref" not in calls
    assert "rs_partial_alt" not in calls
    assert "rs_missing" not in calls
    assert calls["rs_haploid"] == "T"
    assert calls["rs_polyploid"] == "CTT"
