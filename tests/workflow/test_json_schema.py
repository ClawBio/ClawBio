"""Validate pipeline.schema.json against the real demo skills.

This pins the schema to actual usage: any change that breaks variant_annotation
or gwas_lookup is caught here. Cross-step semantics (cycles, unresolved
placeholders, etc.) live in the python validator and are tested separately in
test_pipeline_load.py — JSON Schema only covers shape.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

jsonschema = pytest.importorskip("jsonschema")


REPO_ROOT = Path(__file__).resolve().parents[2]
VARIANT_ANNOTATION = Path(__file__).resolve().parent / "fixtures" / "variant_annotation"
GWAS_LOOKUP = REPO_ROOT / "skills" / "gwas-lookup"
ANALYZE_FASTA = REPO_ROOT / "skills" / "analyze-fasta"
RNASEQ_DE = REPO_ROOT / "skills" / "rnaseq-de"
CLAW_ANCESTRY_PCA = REPO_ROOT / "skills" / "claw-ancestry-pca"
BIGQUERY_PUBLIC = REPO_ROOT / "skills" / "bigquery-public"
CLAW_METAGENOMICS = REPO_ROOT / "skills" / "claw-metagenomics"
CLINPGX = REPO_ROOT / "skills" / "clinpgx"
CRISPR_SCREEN_TRIAGE = REPO_ROOT / "skills" / "crispr-screen-triage"
DIFF_VISUALIZER = REPO_ROOT / "skills" / "diff-visualizer"
EQUITY_SCORER = REPO_ROOT / "skills" / "equity-scorer"
GENOME_COMPARE = REPO_ROOT / "skills" / "genome-compare"
GWAS_PRS = REPO_ROOT / "skills" / "gwas-prs"
MARKER_DOMINANCE_MAPPER = REPO_ROOT / "skills" / "marker-dominance-mapper"
NUTRIGX = REPO_ROOT / "skills" / "nutrigx"
PHARMGX_REPORTER = REPO_ROOT / "skills" / "pharmgx-reporter"
PHYLOGENETICS_BUILDER = REPO_ROOT / "skills" / "phylogenetics-builder"
PROFILE_REPORT = REPO_ROOT / "skills" / "profile-report"
SAMPLE_QC_TRIAGE = REPO_ROOT / "skills" / "sample-qc-triage"
UKB_NAVIGATOR = REPO_ROOT / "skills" / "ukb-navigator"
SCRNA_EMBEDDING = REPO_ROOT / "skills" / "scrna-embedding"
SCRNA_ORCHESTRATOR = REPO_ROOT / "skills" / "scrna-orchestrator"
BIOCONDUCTOR_BRIDGE = REPO_ROOT / "skills" / "bioconductor-bridge"
GALAXY_BRIDGE = REPO_ROOT / "skills" / "galaxy-bridge"
NFCORE_RNASEQ_WRAPPER = REPO_ROOT / "skills" / "nfcore-rnaseq-wrapper"
NFCORE_SAREK_WRAPPER = REPO_ROOT / "skills" / "nfcore-sarek-wrapper"
NFCORE_SCRNASEQ_WRAPPER = REPO_ROOT / "skills" / "nfcore-scrnaseq-wrapper"

SCHEMA_PATH = REPO_ROOT / "docs" / "clawbio2" / "pipeline.schema.json"


@pytest.fixture(scope="module")
def validator():
    schema = json.loads(SCHEMA_PATH.read_text())
    jsonschema.Draft7Validator.check_schema(schema)
    return jsonschema.Draft7Validator(schema)


@pytest.mark.parametrize(
    "skill_dir",
    [
        VARIANT_ANNOTATION,
        GWAS_LOOKUP,
        ANALYZE_FASTA,
        RNASEQ_DE,
        CLAW_ANCESTRY_PCA,
        BIGQUERY_PUBLIC,
        CLAW_METAGENOMICS,
        CLINPGX,
        CRISPR_SCREEN_TRIAGE,
        DIFF_VISUALIZER,
        EQUITY_SCORER,
        GENOME_COMPARE,
        GWAS_PRS,
        MARKER_DOMINANCE_MAPPER,
        NUTRIGX,
        PHARMGX_REPORTER,
        PHYLOGENETICS_BUILDER,
        PROFILE_REPORT,
        SAMPLE_QC_TRIAGE,
        UKB_NAVIGATOR,
        SCRNA_EMBEDDING,
        SCRNA_ORCHESTRATOR,
        BIOCONDUCTOR_BRIDGE,
        GALAXY_BRIDGE,
        NFCORE_RNASEQ_WRAPPER,
        NFCORE_SAREK_WRAPPER,
        NFCORE_SCRNASEQ_WRAPPER,
    ],
)
def test_real_skills_pass_json_schema(validator, skill_dir: Path):
    data = yaml.safe_load((skill_dir / "pipeline.yaml").read_text())
    errors = list(validator.iter_errors(data))
    assert not errors, "\n".join(f"{list(e.absolute_path)}: {e.message}" for e in errors)


def test_schema_rejects_missing_required_top_level(validator):
    bad = {"skill": "x", "version": 1}  # no steps
    assert list(validator.iter_errors(bad))


def test_schema_rejects_wrong_version(validator):
    bad = {"skill": "x", "version": 2, "steps": [{"id": "a", "summary": "s", "command_template": "echo"}]}
    errors = list(validator.iter_errors(bad))
    assert any("version" in str(e.absolute_path) or "const" in e.message for e in errors)


def test_schema_rejects_empty_steps(validator):
    bad = {"skill": "x", "version": 1, "steps": []}
    assert list(validator.iter_errors(bad))


def test_schema_rejects_step_missing_required(validator):
    bad = {"skill": "x", "version": 1, "steps": [{"id": "a", "summary": "s"}]}  # no command_template
    assert list(validator.iter_errors(bad))


def test_schema_rejects_unknown_step_field(validator):
    bad = {
        "skill": "x", "version": 1,
        "steps": [{"id": "a", "summary": "s", "command_template": "echo", "when": "always"}],
    }
    assert list(validator.iter_errors(bad))


def test_schema_rejects_bad_id_format(validator):
    bad = {
        "skill": "x", "version": 1,
        "steps": [{"id": "1-bad", "summary": "s", "command_template": "echo"}],
    }
    assert list(validator.iter_errors(bad))
