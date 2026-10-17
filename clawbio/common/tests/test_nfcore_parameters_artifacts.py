"""Exact parameter bytes and portable artifact-writing contracts."""

import csv
import importlib
import json
from pathlib import Path

import pytest
import yaml
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st


def test_parameter_header_and_key_order_are_stable(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.parameters")
    expected = (
        "# 'input' and 'outdir' are relative to the Nextflow launch directory (output_dir).\n"
        "# Use reproducibility/commands.sh for the recommended replay path.\n"
        "# Manual replay: cd <output_dir> && nextflow run nf-core/example -params-file reproducibility/params.yaml\n"
        "second: false\nfirst: 7\nempty: null\n"
    )
    values = {"second": False, "first": 7, "empty": None}
    assert (
        api.serialize_params_yaml(values, pipeline_name="nf-core/example") == expected
    )
    path = api.write_params_yaml(
        values, output_dir=tmp_path, pipeline_name="nf-core/example"
    )
    assert path.read_bytes() == expected.encode()



def test_long_reference_paths_stay_on_one_line():
    """PyYAML wraps long scalars at 80 columns by default; a reference path with a
    space past that column then spans two lines, which the line-based, stdlib-only
    remap_paths.py cannot read whole. Every value must stay on its own line."""
    api = importlib.import_module("clawbio.common.nfcore.parameters")
    path = "/Users/dana/Library/Mobile_Documents/com~apple~CloudDocs/Project_Data/refs/GRCh38 release v2/genome.fa"
    text = api.serialize_params_yaml({"fasta": path, "outdir": "upstream/results"}, pipeline_name="nf-core/x")
    assert f"\nfasta: {path}\n" in text
    assert yaml.safe_load(text)["fasta"] == path

@settings(deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(st.lists(st.text(alphabet='abc é,"', max_size=20), max_size=5))
def test_csv_rows_round_trip_with_lf(tmp_path, values):
    api = importlib.import_module("clawbio.common.nfcore.artifacts")
    path = tmp_path / "rows.csv"
    rows = [{"sample": value, "read": "read.fastq.gz"} for value in values]
    api.write_csv_lf(path, ["sample", "read"], rows)
    assert b"\r" not in path.read_bytes()
    with path.open(newline="") as handle:
        assert list(csv.DictReader(handle)) == rows


def test_json_writer_respects_newline_and_unicode(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.artifacts")
    path = tmp_path / "record.json"
    api.write_json_lf(path, {"name": "é"}, trailing_newline=False)
    assert path.read_bytes() == b'{\n  "name": "\\u00e9"\n}'
    api.write_json_lf(path, {"name": "é"}, trailing_newline=True)
    assert path.read_bytes().endswith(b"\n")
    assert json.loads(path.read_text()) == {"name": "é"}


def test_safe_config_basename_has_stable_fallback():
    api = importlib.import_module("clawbio.common.nfcore.artifacts")
    assert api.safe_config_basename(" a é.config ") == "a_.config"
    assert api.safe_config_basename("...") == "nextflow.config"


def test_json_snapshots_preserve_names_and_newline_policy(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.artifacts")
    api.write_json_snapshots(tmp_path, {"one.json": {"z": 1}, "two.json": {"a": False}})
    assert {p.name for p in tmp_path.iterdir()} == {"one.json", "two.json"}
    assert (tmp_path / "one.json").read_bytes() == b'{\n  "z": 1\n}'
    assert (tmp_path / "two.json").read_bytes() == b'{\n  "a": false\n}'


def test_sorted_json_writer_supports_legacy_path_conversion(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.artifacts")
    path = tmp_path / "record.json"
    api.write_json_lf(
        path,
        {"z": 1, "a": Path("file")},
        trailing_newline=True,
        sort_keys=True,
        default=str,
    )
    assert path.read_bytes() == b'{\n  "a": "file",\n  "z": 1\n}\n'


_IGNORING_CONTRACT = {
    "ignored_parameters": ["genomes", "freebayes_filter", "vep_cache_version"],
    "parameters": {
        "freebayes_filter": {"type": "string"},
        "vep_cache_version": {"type": "string"},
        "aligner": {"type": "string", "enum": ["bwa-mem", "dragmap"]},
    },
}


def test_upstream_ignored_parameters_are_not_validated():
    """nf-schema skips validation.defaultIgnoreParams. nf-core/sarek 3.8.1 lists
    freebayes_filter and vep_cache_version there because its schema types them as
    strings while its own defaults are integers (30 and 113)."""
    from clawbio.common.nfcore.parameters import validate_against_contract

    validate_against_contract({"freebayes_filter": 30, "vep_cache_version": 113}, _IGNORING_CONTRACT)
    validate_against_contract({"freebayes_filter": 20}, _IGNORING_CONTRACT, patterns=False)
    with pytest.raises(ValueError):
        validate_against_contract({"aligner": "star"}, _IGNORING_CONTRACT)
    with pytest.raises(ValueError):
        validate_against_contract({"freebayes_filter": 30}, {**_IGNORING_CONTRACT, "ignored_parameters": []})


@given(
    st.one_of(
        st.integers(), st.floats(allow_nan=True, allow_infinity=True), st.text(max_size=8),
        st.booleans(), st.lists(st.integers(), max_size=3),
    )
)
def test_any_value_of_an_ignored_parameter_passes_and_others_still_fail(value):
    from clawbio.common.nfcore.parameters import validate_against_contract

    validate_against_contract({"freebayes_filter": value}, _IGNORING_CONTRACT)
    with pytest.raises(ValueError):
        validate_against_contract({"freebayes_filter": value, "aligner": "star"}, _IGNORING_CONTRACT)


@given(st.integers(1, 512), st.integers(1, 2048), st.one_of(st.none(), st.fixed_dictionaries(
    {"cpus": st.integers(1, 64), "memory_gb": st.integers(1, 256)})))
def test_profile_limits_only_ever_lower_capacity(cpus, memory_gb, limits):
    from clawbio.common.nfcore.resources import cap_to_profile_limits

    bounded = cap_to_profile_limits(cpus, memory_gb, limits)
    assert bounded[0] <= cpus and bounded[1] <= memory_gb
    if limits is None:
        assert bounded == (cpus, memory_gb)
    else:
        assert bounded == (min(cpus, limits["cpus"]), min(memory_gb, limits["memory_gb"]))
