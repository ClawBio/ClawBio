"""Regression tests for phenotype labels in the PharmGx terminal preview."""

from contextlib import redirect_stdout
from io import StringIO
import re

import pytest
from hypothesis import example, given, strategies as st

from clawbio import cli


def _gene_line(phenotype: str, diplotype: str = "*1/*1") -> str:
    report = (
        "## Gene Profiles\n\n"
        "| Gene | Full Name | Diplotype | Phenotype |\n"
        "|------|-----------|-----------|-----------|\n"
        f"| TESTGENE | Synthetic test gene | {diplotype} | {phenotype} |\n\n"
        "---\n"
    )
    output = StringIO()
    with redirect_stdout(output):
        cli.format_pharmgx_preview(report, "report.md")
    return next(
        line for line in output.getvalue().splitlines() if "TESTGENE" in line
    )


def _plain(text: str) -> str:
    return re.sub(r"\x1b\[[0-9;]*m", "", text)


@pytest.mark.parametrize(
    "phenotype",
    [
        "Ultrarapid Metabolizer",
        "Rapid Metabolizer",
        "Intermediate Metabolizer",
        "Poor Metabolizer",
        "Normal Metabolizer",
        "Normal Function",
        "Intermediate Warfarin Sensitivity",
        "High Warfarin Sensitivity",
        "Poor Function",
        "Decreased Function",
        "CYP3A5 Expressor",
        "Intermediate Expressor",
        "CYP3A5 Non-expressor",
        "Reduced MTHFR enzyme activity (677CT)",
        "Strongly reduced MTHFR enzyme activity (677TT)",
        "Indeterminate (incomplete coverage: 3/4 SNPs tested)",
        "Indeterminate (rs8175347 not tested)",
        "Indeterminate (not genotyped)",
        "Indeterminate (not in panel)",
        "Indeterminate (phase ambiguity: *2 + *17)",
        "Unknown (unmapped diplotype)",
        "Unrecognized phenotype",
    ],
)
def test_preview_preserves_complete_phenotype(phenotype):
    assert _plain(_gene_line(phenotype)).split("*1/*1", 1)[1].strip() == phenotype


@given(
    category=st.sampled_from(["Indeterminate", "Unknown", "Unmapped"]),
    reason=st.text(alphabet="abcdefghijklmnopqrstuvwxyz 0123456789:/", max_size=40),
)
@example(category="Indeterminate", reason="not genotyped")
@example(category="Indeterminate", reason="incomplete coverage: 3/4")
def test_uncertain_phenotype_is_never_replaced_with_normal(category, reason):
    phenotype = f"{category} ({reason})"
    displayed = _plain(_gene_line(phenotype)).split("*1/*1", 1)[1].strip()
    assert displayed == phenotype
    assert displayed != "Normal"


@pytest.mark.parametrize(
    ("phenotype", "color"),
    [
        ("Normal Metabolizer", "GREEN"),
        ("Normal Function", "GREEN"),
        ("Normal Warfarin Sensitivity", "GREEN"),
        ("Ultrarapid Metabolizer", "YELLOW"),
        ("Rapid Metabolizer", "YELLOW"),
        ("CYP3A5 Expressor", "YELLOW"),
        ("Reduced MTHFR enzyme activity (677CT)", "YELLOW"),
        ("Indeterminate (Normal/Normal; incomplete coverage)", "YELLOW"),
        ("Unknown (High sensitivity unconfirmed)", "YELLOW"),
        ("Unmapped diplotype", "YELLOW"),
        ("Unrecognized phenotype", "YELLOW"),
        ("Poor Function", "RED"),
        ("High Warfarin Sensitivity", "RED"),
        ("CYP3A5 Non-expressor", "DIM"),
    ],
)
def test_preview_uses_explicit_phenotype_colors(monkeypatch, phenotype, color):
    colors = {"GREEN": "\x1b[32m", "YELLOW": "\x1b[33m",
              "RED": "\x1b[31m", "DIM": "\x1b[2m"}
    for name, value in colors.items():
        monkeypatch.setattr(cli, name, value)
    monkeypatch.setattr(cli, "RESET", "\x1b[0m")
    line = _gene_line(phenotype)
    assert f"{colors[color]}{phenotype}\x1b[0m" in line


def test_empty_phenotype_is_unknown_not_normal():
    assert _plain(_gene_line("")).split("*1/*1", 1)[1].strip() == "Unknown"


def test_preview_retains_incomplete_coverage_warning_and_diplotype_note():
    line = _plain(_gene_line(
        "Indeterminate (incomplete coverage: 3/4 SNPs tested)",
        "Normal/Normal (3/4 SNPs tested)",
    ))
    assert "Normal/Normal" in line
    assert "Indeterminate (incomplete coverage: 3/4 SNPs tested)" in line
    assert "(3/4 SNPs tested)" in line
