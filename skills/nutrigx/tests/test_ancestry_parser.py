"""
test_ancestry_parser.py — AncestryDNA parsing for NutriGx Advisor.

Every test in test_nutrigx.py uses fmt="23andme", so the AncestryDNA path had no
coverage. Kept in its own file so it does not collide with other changes to
test_nutrigx.py.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

def test_parse_ancestry_skips_comments_and_joins_alleles(tmp_path):
    from parse_input import parse_ancestry

    f = tmp_path / "ancestry.txt"
    f.write_text(
        "#AncestryDNA raw data download\n"
        "#another comment\n"
        "rsid\tchromosome\tposition\tallele1\tallele2\n"
        "rs1801133\t1\t11856378\tC\tT\n"
        "rs4988235\t2\t136608646\tA\tG\n"
        "notanrsid\t3\t300\tA\tA\n",
        encoding="utf-8",
    )
    assert parse_ancestry(str(f)) == {"rs1801133": "CT", "rs4988235": "AG"}


def test_parse_ancestry_tolerates_truncated_rows(tmp_path):
    from parse_input import parse_ancestry

    f = tmp_path / "truncated.txt"
    f.write_text(
        "rsid\tchromosome\tposition\tallele1\tallele2\n"
        "rs1801133\t1\t11856378\tC\tT\n"
        "rs4988235\t2\t136608646\tA\n",          # allele2 missing
        encoding="utf-8",
    )
    table = parse_ancestry(str(f))
    assert table["rs1801133"] == "CT"
    assert table["rs4988235"] == "A"
