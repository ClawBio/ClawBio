"""Whitelist invariant for examples/nebius_agent.py."""

from __future__ import annotations

import importlib.util
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "examples" / "nebius_agent.py"


def _load():
    spec = importlib.util.spec_from_file_location("nebius_agent", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_every_skill_is_demo_pinned():
    mod = _load()
    for name, spec in mod.SKILLS.items():
        args = spec["args"]
        assert "--demo" in args, name
        assert "--input" not in args, name


def test_whitelist_rejects_input_path(monkeypatch):
    mod = _load()
    poisoned = {
        "bad": {
            "script": "skills/vcf-annotator/vcf_annotator.py",
            "args": ["--input", "patient.vcf"],
            "report": "report.md",
            "description": "should not load",
        }
    }
    monkeypatch.setitem(mod.__dict__, "SKILLS", poisoned)
    try:
        mod._require_demo_only_whitelist(poisoned)
    except SystemExit as exc:
        assert "Offending: bad" in str(exc)
    else:
        raise AssertionError("expected SystemExit")
