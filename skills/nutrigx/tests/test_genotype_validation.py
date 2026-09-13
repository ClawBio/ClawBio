"""
test_genotype_validation.py — genotype calls are validated on every entry path.

Covers the whitelist itself, no_call status for rejected calls, the CLI (which
parses through clawbio.common.parsers), the API (caller-built dicts) and the
report guard.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from extract_genotypes import extract_snp_genotypes

SKILL_DIR = Path(__file__).parent.parent
PANEL = SKILL_DIR / "data" / "snp_panel.json"


def load_panel():
    with open(PANEL) as f:
        return json.load(f)


def test_clean_genotype_accepts_nucleotides_and_rejects_everything_else():
    from parse_input import clean_genotype

    assert clean_genotype("ct") == "CT"
    assert clean_genotype(" AG ") == "AG"
    assert clean_genotype("DI") == "DI"
    for bad in ("--", "00", "", None, "C`T", "<b>", "A|G", "CT; rm -rf"):
        assert clean_genotype(bad) is None, bad


def test_invalid_calls_are_reported_as_no_call_not_untested():
    from parse_input import clean_genotype_table

    panel = load_panel()
    table = clean_genotype_table({"rs1801133": "C`T</code>", "rs4988235": "--"})
    calls = extract_snp_genotypes(table, panel)
    assert calls["rs1801133"]["status"] == "no_call"
    assert calls["rs4988235"]["status"] == "no_call"


def test_cli_path_applies_the_whitelist():
    """nutrigx.py parses through clawbio.common.parsers, which does not validate."""
    src = (SKILL_DIR / "nutrigx.py").read_text()
    assert "clean_genotype_table(genotypes_to_simple(records))" in src


def test_api_path_applies_the_whitelist():
    import api

    result = api.run({"rs1801133": "<script>", "rs4988235": "AG", "rs9939609": "TT"})
    assert result["snp_calls"]["rs1801133"]["status"] == "no_call"


def test_safe_display_genotype_strips_table_breaking_characters():
    from generate_report import safe_display_genotype

    out = safe_display_genotype("A|G`<x>")
    for ch in ("|", "`", "<", ">"):
        assert ch not in out
    assert safe_display_genotype("") == "--"


def test_skill_parsers_validate_through_parse_genetic_file(tmp_path):
    from parse_input import parse_genetic_file

    f = tmp_path / "hostile.txt"
    f.write_text(
        "rsid\tchromosome\tposition\tgenotype\n"
        "rs1801133\t1\t1\tC`T</code>\n"
        "rs4988235\t2\t2\tAG\n",
        encoding="utf-8",
    )
    table = parse_genetic_file(str(f), fmt="23andme")
    assert table == {"rs1801133": "", "rs4988235": "AG"}
