"""Unit tests for proposer.py. LLM backend is mocked."""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from skills.clawpathy_autoresearch.proposer import (
    ProposalResult,
    _build_proposer_prompt,
    _extract_proposal,
    propose_edit,
)


def test_extract_proposal_from_fenced_json():
    text = '```json\n{"label": "add MHC collapse", "new_skill": "# Skill\\nStep 1"}\n```'
    result = _extract_proposal(text)
    assert result["label"] == "add MHC collapse"
    assert "Step 1" in result["new_skill"]


def test_extract_proposal_handles_preamble():
    text = 'Here is the edit:\n{"label": "x", "new_skill": "y"}\ntrailing'
    result = _extract_proposal(text)
    assert result["label"] == "x"


def test_extract_proposal_bad_json_returns_none():
    assert _extract_proposal("not json at all") is None


def test_build_proposer_prompt_includes_history_and_skill():
    prompt = _build_proposer_prompt(
        workspace_name="GWAS Repro",
        current_skill="# Pipeline\n1. clump at 142kb",
        history_labels=["baseline", "window 500→400kb"],
        last_score=0.23,
    )
    assert "GWAS Repro" in prompt
    assert "142kb" in prompt
    assert "window 500" in prompt
    assert "0.23" in prompt


def test_propose_edit_writes_new_skill(tmp_path: Path):
    skill_path = tmp_path / "SKILL.md"
    skill_path.write_text("# Original\n1. Do thing")

    fake_response = (
        '{"label": "add MHC collapse step", '
        '"new_skill": "# Original\\n1. Do thing\\n2. Collapse MHC"}'
    )
    with patch(
        "skills.clawpathy_autoresearch.proposer._call_llm",
        return_value=fake_response,
    ):
        result = propose_edit(
            workspace_name="test",
            skill_path=skill_path,
            history_labels=["baseline"],
            last_score=0.5,
        )

    assert isinstance(result, ProposalResult)
    assert result.label == "add MHC collapse step"
    assert "Collapse MHC" in skill_path.read_text()


def test_propose_edit_handles_llm_failure(tmp_path: Path):
    skill_path = tmp_path / "SKILL.md"
    skill_path.write_text("# Original")

    with patch(
        "skills.clawpathy_autoresearch.proposer._call_llm",
        side_effect=RuntimeError("backend down"),
    ):
        result = propose_edit(
            workspace_name="t",
            skill_path=skill_path,
            history_labels=[],
            last_score=0.5,
        )

    assert result.error is not None
    assert result.label == "no-op (proposer failed)"
    assert skill_path.read_text() == "# Original"  # unchanged
