"""The checks this wrapper leaves to nf-core/sarek must be in the pinned release.

samplesheet_builder.py and preflight.py name what Sarek 3.8.1 refuses on its own
(nf-schema over assets/schema_input.json and nextflow_schema.json, and
subworkflows/local/samplesheet_to_channel). These tests read that release, vendored
in tests/fixtures/nfcore_pinned or the tree under CLAWBIO_NFCORE_UPSTREAM_DIR, so a
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
samplesheet_builder = load_tested_module(_SKILL_DIR, "samplesheet_builder")

_RELEASE = f"sarek-{schemas.DEFAULT_PIPELINE_VERSION}"
_UPSTREAM = os.environ.get("CLAWBIO_NFCORE_UPSTREAM_DIR")
_TREE = Path(_UPSTREAM) / _RELEASE if _UPSTREAM else _SKILL_DIR.parents[1] / "tests" / "fixtures" / "nfcore_pinned" / _RELEASE
_SAMPLESHEET = json.loads((_TREE / "assets" / "schema_input.json").read_text(encoding="utf-8"))["items"]
_COLUMNS = _SAMPLESHEET["properties"]
_TO_CHANNEL = (_TREE / "subworkflows" / "local" / "samplesheet_to_channel" / "main.nf").read_text(encoding="utf-8")


def _parameters() -> dict[str, dict]:
    schema = json.loads((_TREE / "nextflow_schema.json").read_text(encoding="utf-8"))
    found: dict[str, dict] = {}
    for group in (schema.get("$defs") or schema.get("definitions") or {}).values():
        found.update(group.get("properties", {}))
    return found


def test_schema_requires_patient_and_sample_with_no_blanks():
    assert {"patient", "sample"} <= set(_SAMPLESHEET["required"])
    for column in ("patient", "sample", "lane"):
        assert _COLUMNS[column]["pattern"] == r"^\S+$", column


def test_schema_enumerates_sex_and_status():
    assert _COLUMNS["sex"]["enum"] == ["XX", "XY", "NA"]
    assert _COLUMNS["status"]["type"] == "integer" and _COLUMNS["status"]["enum"] == [0, 1]


def test_schema_types_contamination_as_a_number():
    """Only the type: the wrapper still refuses values outside [0, 1]."""
    assert _COLUMNS["contamination"]["type"] == "number"


@pytest.mark.parametrize("column", samplesheet_builder._PATH_COLUMNS)
def test_schema_checks_every_file_column_exists_and_is_named_right(column):
    assert _COLUMNS[column]["format"] == "file-path"
    assert _COLUMNS[column]["exists"] is True
    assert _COLUMNS[column]["pattern"]


def test_schema_pairs_read_columns_and_refuses_repeated_rows():
    assert _SAMPLESHEET["dependentRequired"] == {"fastq_2": ["fastq_1"], "spring_2": ["spring_1"]}
    assert set(_SAMPLESHEET["uniqueEntries"]) == {"patient", "sample", "lane"}


@pytest.mark.parametrize("rule,pattern", [
    ("sex for every sample with ascat, controlfreec or varlociraptor",
     r"contains\('ascat'\) \|\| tools\.split\(','\)\.contains\('controlfreec'\) \|\| "
     r"tools\.split\(','\)\.contains\('varlociraptor'\).{0,200}?if \(it\[0\]\.sex == 'NA'\) \{\s*"
     r"error\(\"Please specify sex information for each sample"),
    ("contamination for every tumour with varlociraptor",
     r"contains\('varlociraptor'\)\).{0,200}?if \(it\[0\]\.status == 1 && !it\[0\]\.containsKey\('contamination'\)\) \{\s*"
     r"error\(\"Please specify contamination information for each tumor sample"),
    ("two normals with exactly one tumour",
     r"if \(status1_count == 1 && status0_count > 1\) \{.{0,300}?"
     r"error\(\"Execution halted due to sample status inconsistency"),
    ("duplicate patient-sample-status-lane rows",
     r"if \(combination_list\.size\(\) > 1\) \{.{0,800}?"
     r"error\(\"Execution halted due to duplicate patient-sample-status-lane combination"),
    ("a sample under two patients",
     r"if \(unique_patients\.size\(\) > 1\) \{.{0,400}?"
     r"error\(\"Execution halted due to sample status inconsistency"),
    ("FASTQ rows that do not fit --step", r"error\(\"Samplesheet contains fastq files but step is `\$\{step\}`"),
    ("Spring rows that do not fit --step", r"error\(\"Samplesheet contains spring files .{0,80}? but step is `\$\{step\}`"),
    ("BAM rows that do not fit --step", r"error\(\"Samplesheet contains bam files but step is `\$\{step\}`"),
    ("CRAM rows that do not fit --step", r"error\(\"Samplesheet contains cram files but step is `\$\{step\}`"),
    ("VCF rows that do not fit --step", r"error\(\"Samplesheet contains vcf files but step is `\$\{step\}`"),
    ("tumour-only tools on a normal-only sheet",
     r"error\('The sample-sheet only contains normal-samples, but the following tools"),
    ("normal-requiring tools on a tumour-only sheet",
     r"error\('The sample-sheet only contains tumor-samples, but the following tools"),
])
def test_samplesheet_to_channel_refuses_what_the_wrapper_left_to_it(rule, pattern):
    assert re.search(pattern, _TO_CHANNEL, re.S), rule


def test_reference_paths_left_to_nf_schema_are_existence_checked_there():
    parameters = _parameters()
    left = [name for name in preflight.REFERENCE_PATH_PARAMS if name not in preflight._PATHS_THE_SCHEMA_DOES_NOT_CHECK]
    assert left and [name for name in left if parameters.get(name, {}).get("exists") is not True] == []


def test_reference_paths_the_wrapper_still_checks_are_not_existence_checked_there():
    parameters = _parameters()
    assert [name for name in preflight._PATHS_THE_SCHEMA_DOES_NOT_CHECK if parameters.get(name, {}).get("exists") is True] == []
