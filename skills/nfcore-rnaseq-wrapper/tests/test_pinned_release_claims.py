"""The checks this wrapper leaves to nf-core/rnaseq must be in the pinned release.

samplesheet_builder.py and preflight.py name what rnaseq 3.26.0 refuses on its own
(nf-schema over assets/schema_input.json and nextflow_schema.json, and
checkSamplesAfterGrouping). These tests read that release, vendored in
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

_RELEASE = f"rnaseq-{schemas.DEFAULT_PIPELINE_VERSION}"
_UPSTREAM = os.environ.get("CLAWBIO_NFCORE_UPSTREAM_DIR")
_TREE = Path(_UPSTREAM) / _RELEASE if _UPSTREAM else _SKILL_DIR.parents[1] / "tests" / "fixtures" / "nfcore_pinned" / _RELEASE
_SAMPLESHEET = json.loads((_TREE / "assets" / "schema_input.json").read_text(encoding="utf-8"))["items"]
_COLUMNS = _SAMPLESHEET["properties"]
_PIPELINE_UTILS = (_TREE / "subworkflows" / "local" / "utils_nfcore_rnaseq_pipeline" / "main.nf").read_text(encoding="utf-8")


def _parameters() -> dict[str, dict]:
    schema = json.loads((_TREE / "nextflow_schema.json").read_text(encoding="utf-8"))
    found: dict[str, dict] = {}
    for group in (schema.get("$defs") or schema.get("definitions") or {}).values():
        found.update(group.get("properties", {}))
    return found


def test_schema_requires_sample_fastq_1_and_strandedness():
    assert {"sample", "fastq_1", "strandedness"} <= set(_SAMPLESHEET["required"])
    assert _COLUMNS["sample"]["pattern"] == r"^\S+$"


def test_schema_enumerates_strandedness():
    assert _COLUMNS["strandedness"]["enum"] == ["forward", "reverse", "unstranded", "auto"]


@pytest.mark.parametrize("column", ["fastq_1", "fastq_2", "genome_bam", "transcriptome_bam"])
def test_schema_checks_every_file_column_exists_and_is_named_right(column):
    assert _COLUMNS[column]["format"] == "file-path"
    assert _COLUMNS[column]["exists"] is True
    assert _COLUMNS[column]["pattern"]


def test_schema_requires_the_input_to_be_an_existing_csv():
    assert _parameters()["input"]["pattern"] == r"^\S+\.csv$"
    assert _parameters()["input"]["exists"] is True


def test_grouping_refuses_one_sample_with_two_strandedness_values():
    assert re.search(
        r"def checkSamplesAfterGrouping\(input\) \{.{0,600}?"
        r"meta\.strandedness \}\.unique\(\)\.size == 1.{0,80}?"
        r"error\(\"Please check input samplesheet -> Multiple runs of a sample must have the same strandedness",
        _PIPELINE_UTILS, re.S,
    )


def test_reference_paths_left_to_nf_schema_are_existence_checked_there():
    parameters = _parameters()
    left = preflight._EXPLICIT_REFERENCE_FIELDS
    assert left and [name for name in left if parameters.get(name, {}).get("exists") is not True] == []


def test_reference_paths_the_wrapper_still_checks_are_not_existence_checked_there():
    parameters = _parameters()
    assert [name for name in preflight._CONTAMINANT_PATH_FIELDS if parameters.get(name, {}).get("exists") is True] == []
