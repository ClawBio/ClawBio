"""BOM-prefixed vendor exports must parse through the public genetic-file API."""
from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from clawbio.common.parsers import detect_format, parse_genetic_file


@pytest.mark.parametrize(
    ("filename", "content", "expected_format", "rsid", "genotype"),
    [
        (
            "ancestry.tsv",
            "rsid\tchromosome\tposition\tallele1\tallele2\nrs100\t1\t101\tA\tG\n",
            "ancestry",
            "rs100",
            "AG",
        ),
        (
            "myheritage.csv",
            "RSID,CHROMOSOME,POSITION,RESULT\nrs200,2,202,CT\n",
            "myheritage",
            "rs200",
            "CT",
        ),
    ],
)
def test_public_parser_accepts_bom_prefixed_vendor_exports(
    tmp_path: Path,
    filename: str,
    content: str,
    expected_format: str,
    rsid: str,
    genotype: str,
) -> None:
    path = tmp_path / filename
    path.write_text(content, encoding="utf-8-sig")
    assert detect_format(path) == expected_format
    assert parse_genetic_file(path)[rsid].genotype == genotype


def test_public_parser_accepts_bom_prefixed_gzip_ancestry_export(tmp_path: Path) -> None:
    path = tmp_path / "ancestry.tsv.gz"
    with gzip.open(path, "wt", encoding="utf-8-sig") as handle:
        handle.write("rsid\tchromosome\tposition\tallele1\tallele2\nrs300\t3\t303\tC\tT\n")
    assert detect_format(path) == "ancestry"
    assert parse_genetic_file(path)["rs300"].genotype == "CT"
