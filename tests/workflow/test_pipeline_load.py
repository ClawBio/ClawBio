"""Pipeline loader + validator tests.

Loads the real demo skills (no synthetic minimal yaml when an existing skill
already exercises the codepath) and stresses the validator with negative cases
written as actual yaml files on disk.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from clawbio.workflow.pipeline import PipelineError, load, render_command

REPO_ROOT = Path(__file__).resolve().parents[2]
VARIANT_ANNOTATION = Path(__file__).resolve().parent / "fixtures" / "variant_annotation"
GWAS_LOOKUP = REPO_ROOT / "skills" / "gwas-lookup"


# ---------------------------------------------------------------------------
# Loading real skills
# ---------------------------------------------------------------------------

def test_load_variant_annotation():
    p = load(VARIANT_ANNOTATION)
    assert p.skill == "variant_annotation"
    assert p.version == 1
    assert p.step_ids == ["normalize", "annotate"]
    assert {i.name for i in p.inputs} == {"vcf_path", "reference"}
    assert p.step("annotate").prereqs == ["normalize"]
    assert p.step("annotate").requires == {"normalized_vcf": "present"}


def test_load_gwas_lookup():
    """The 11-step skill — proves the schema accepts a real DAG with fan-out
    and a heavy `merge` step that gates on all upstream emits."""
    p = load(GWAS_LOOKUP)
    assert p.skill == "gwas-lookup"
    assert len(p.steps) == 11

    merge = p.step("merge")
    # 8 parallel fetches + the resolver must all be done before merge
    assert "resolve_variant" in merge.prereqs
    assert len(merge.prereqs) == 9
    # merge gates on every upstream emit
    assert set(merge.requires.keys()) == {
        "variant_json",
        "gwas_catalog_json",
        "open_targets_json",
        "open_targets_credsets_json",
        "pheweb_ukb_json",
        "finngen_json",
        "pheweb_bbj_json",
        "gtex_json",
        "eqtl_catalogue_json",
    }
    assert merge.emits == ["merged_json"]


def test_load_missing_yaml(tmp_path: Path):
    with pytest.raises(PipelineError, match="no pipeline.yaml"):
        load(tmp_path)


# ---------------------------------------------------------------------------
# Validator negative cases — each writes a real yaml file
# ---------------------------------------------------------------------------

def _write(tmp_path: Path, body: str) -> Path:
    (tmp_path / "pipeline.yaml").write_text(body)
    return tmp_path


def test_reject_unknown_schema_version(tmp_path: Path):
    skill = _write(tmp_path, "skill: s\nversion: 2\nsteps:\n  - id: a\n    summary: x\n    command_template: 'echo a'\n")
    with pytest.raises(PipelineError, match="unsupported pipeline schema version"):
        load(skill)


def test_reject_duplicate_step_ids(tmp_path: Path):
    skill = _write(tmp_path, """
skill: s
version: 1
steps:
  - {id: a, summary: x, command_template: 'echo a'}
  - {id: a, summary: x, command_template: 'echo a2'}
""")
    with pytest.raises(PipelineError, match="duplicate step ids"):
        load(skill)


def test_reject_unknown_prereq(tmp_path: Path):
    skill = _write(tmp_path, """
skill: s
version: 1
steps:
  - {id: a, summary: x, command_template: 'echo a', prereqs: [ghost]}
""")
    with pytest.raises(PipelineError, match="unknown prereq 'ghost'"):
        load(skill)


def test_reject_forward_prereq(tmp_path: Path):
    """A step cannot reference a later step as prereq (would form a cycle)."""
    skill = _write(tmp_path, """
skill: s
version: 1
steps:
  - {id: a, summary: x, command_template: 'echo a', prereqs: [b]}
  - {id: b, summary: x, command_template: 'echo b'}
""")
    with pytest.raises(PipelineError, match="not declared earlier"):
        load(skill)


def test_reject_unknown_placeholder(tmp_path: Path):
    skill = _write(tmp_path, """
skill: s
version: 1
steps:
  - {id: a, summary: x, command_template: 'echo {nope}'}
""")
    with pytest.raises(PipelineError, match=r"unknown placeholder \{nope\}"):
        load(skill)


def test_reject_undeclared_requires_emit(tmp_path: Path):
    skill = _write(tmp_path, """
skill: s
version: 1
steps:
  - {id: a, summary: x, command_template: 'echo a'}
  - id: b
    summary: x
    command_template: 'echo b'
    requires: {phantom_emit: present}
""")
    with pytest.raises(PipelineError, match="requires 'phantom_emit'"):
        load(skill)


def test_reject_requires_inputs_undeclared(tmp_path: Path):
    skill = _write(tmp_path, """
skill: s
version: 1
inputs:
  - {name: real}
steps:
  - {id: a, summary: x, command_template: 'echo {real}', requires_inputs: [missing]}
""")
    with pytest.raises(PipelineError, match="undeclared input 'missing'"):
        load(skill)


def test_empty_steps_rejected(tmp_path: Path):
    skill = _write(tmp_path, "skill: s\nversion: 1\nsteps: []\n")
    with pytest.raises(PipelineError, match="non-empty list"):
        load(skill)


# ---------------------------------------------------------------------------
# render_command
# ---------------------------------------------------------------------------

def test_render_command_substitutes_inputs_and_emits():
    p = load(VARIANT_ANNOTATION)
    cmd = render_command(
        p.step("annotate"),
        inputs={"vcf_path": "/in.vcf", "reference": "/ref.fa"},
        emitted={"normalized_vcf": "/run/out/normalized.vcf"},
        out_dir="/run/out",
        skill_dir=str(VARIANT_ANNOTATION),
    )
    assert "/run/out/normalized.vcf" in cmd
    assert "/run/out/annotated.vcf" in cmd
    assert "{" not in cmd  # nothing unresolved


def test_render_command_unresolved_raises():
    p = load(VARIANT_ANNOTATION)
    with pytest.raises(PipelineError, match="unresolved placeholder"):
        render_command(
            p.step("annotate"),
            inputs={},
            emitted={},
            out_dir="/run/out",
            skill_dir=str(VARIANT_ANNOTATION),
        )
