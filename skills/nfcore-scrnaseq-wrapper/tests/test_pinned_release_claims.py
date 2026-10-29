"""The checks this wrapper leaves to nf-core/scrnaseq must be in the pinned release.

samplesheet_builder.py and preflight.py name what scrnaseq 4.1.0 refuses on its own
(nf-schema over assets/schema_input.json and nextflow_schema.json, and
PIPELINE_INITIALISATION). These tests read that release, vendored in
tests/fixtures/nfcore_pinned or the tree under CLAWBIO_NFCORE_UPSTREAM_DIR, so a
claim that stops being true fails here instead of letting bad input reach a task.
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path

import pytest

from clawbio.common.tests.nfcore_test_imports import load_tested_module

_SKILL_DIR = Path(__file__).resolve().parent.parent
schemas = load_tested_module(_SKILL_DIR, "schemas")
preflight = load_tested_module(_SKILL_DIR, "preflight")

_RELEASE = f"scrnaseq-{schemas.DEFAULT_PIPELINE_VERSION}"
_UPSTREAM = os.environ.get("CLAWBIO_NFCORE_UPSTREAM_DIR")
_TREE = Path(_UPSTREAM) / _RELEASE if _UPSTREAM else _SKILL_DIR.parents[1] / "tests" / "fixtures" / "nfcore_pinned" / _RELEASE
_SAMPLESHEET = json.loads((_TREE / "assets" / "schema_input.json").read_text(encoding="utf-8"))["items"]
_COLUMNS = _SAMPLESHEET["properties"]
_PIPELINE_UTILS = (_TREE / "subworkflows" / "local" / "utils_nfcore_scrnaseq_pipeline" / "main.nf").read_text(encoding="utf-8")


def _parameters() -> dict[str, dict]:
    schema = json.loads((_TREE / "nextflow_schema.json").read_text(encoding="utf-8"))
    found: dict[str, dict] = {}
    for group in (schema.get("$defs") or schema.get("definitions") or {}).values():
        found.update(group.get("properties", {}))
    return found


def test_schema_requires_sample_and_both_fastqs():
    assert {"sample", "fastq_1", "fastq_2"} <= set(_SAMPLESHEET["required"])


@pytest.mark.parametrize("column", ["fastq_1", "fastq_2"])
def test_schema_checks_each_fastq_exists_and_is_named_right(column):
    assert _COLUMNS[column]["format"] == "file-path"
    assert _COLUMNS[column]["exists"] is True
    assert _COLUMNS[column]["pattern"]


def test_schema_types_expected_cells_and_enumerates_the_type_columns():
    assert _COLUMNS["expected_cells"]["type"] == "integer"
    assert _COLUMNS["sample_type"]["enum"] == ["atac", "gex"]
    assert _COLUMNS["feature_type"]["enum"] == ["gex", "vdj", "ab", "crispr", "cmo"]


@pytest.mark.parametrize("rule,pattern", [
    ("a cellrangerarc ATAC row without fastq_barcode",
     r"if \(!fastq_2 \|\| \(meta\.sample_type == \"atac\" && !meta\.fastq_barcode\)\) \{\s*"
     r"error\(\"Please check input samplesheet -> cellrangerarc requires both paired-end reads and barcode fastq files"),
    ("a cellrangerarc sample_type other than gex or atac",
     r"def valid_sample_types = \[\"gex\", \"atac\"\].{0,200}?"
     r"error\(\"Please check input samplesheet -> The property 'sample_type' is required and can only be 'gex' or 'atac'"),
    ("runs of one sample with different endedness",
     r"def validateInputSamplesheet\(input\) \{.{0,300}?meta\.single_end \}\.unique\(\)\.size == 1.{0,80}?"
     r"error\(\"Please check input samplesheet -> Multiple runs of a sample must be of the same datatype"),
])
def test_pipeline_initialisation_refuses_what_the_wrapper_left_to_it(rule, pattern):
    assert re.search(pattern, _PIPELINE_UTILS, re.S), rule


def test_reference_paths_left_to_nf_schema_are_existence_checked_there():
    parameters = _parameters()
    left = [name for name in schemas.ALL_REFERENCE_PATH_FIELDS if name not in preflight._PATHS_THE_SCHEMA_DOES_NOT_CHECK]
    assert left and [name for name in left if parameters.get(name, {}).get("exists") is not True] == []


def test_reference_paths_the_wrapper_still_checks_are_not_existence_checked_there():
    parameters = _parameters()
    assert [name for name in preflight._PATHS_THE_SCHEMA_DOES_NOT_CHECK if parameters.get(name, {}).get("exists") is True] == []
