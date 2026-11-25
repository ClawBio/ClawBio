"""
Tests for clawbio.common.portable_commands
Run with: pytest tests/test_portable_commands.py -v
"""

import importlib.util
import tempfile
from pathlib import Path


# Load portable_commands directly to avoid clawbio/__init__.py → runner.py chain
# (runner.py loads clawbio.py which uses `str | None`, Python 3.10+ only)
_mod_path = Path(__file__).parent.parent / "common" / "portable_commands.py"
_spec = importlib.util.spec_from_file_location("clawbio.common.portable_commands", _mod_path)
_mod = importlib.util.module_from_spec(_spec)  # type: ignore[arg-type]
_spec.loader.exec_module(_mod)  # type: ignore[union-attr]
build_portable_commands_sh = _mod.build_portable_commands_sh
write_portable_commands_sh = _mod.write_portable_commands_sh


# ── build_portable_commands_sh tests ──────────────────────────────────────────

class TestBuildPortableCommandsSh:

    def test_contains_shebang(self):
        result = build_portable_commands_sh(
            skill_name="lit-synthesizer",
            script_name="lit_synthesizer.py",
            args={"--query": "CRISPR"},
        )
        assert result.startswith("#!/usr/bin/env bash")

    def test_boolean_flag_no_value(self):
        result = build_portable_commands_sh(
            skill_name="lit-synthesizer",
            script_name="lit_synthesizer.py",
            args={"--demo": None},
        )
        assert "--demo" in result
        # Should not have "--demo None"
        assert "--demo None" not in result

    def test_uses_portable_python_interpreter(self):
        """The replay command must invoke a portable interpreter (`${PYTHON:-python3}`),
        never a bare `python`. Modern macOS and many Linux distributions ship only
        `python3` (PEP 394), so a bare `python "$SKILL_SCRIPT"` replay fails with
        `python: command not found`. Enforced at the source template so every skill
        that builds a bundle inherits the fix (no per-skill post-generation patch)."""
        result = build_portable_commands_sh(
            skill_name="lit-synthesizer",
            script_name="lit_synthesizer.py",
            args={"--query": "CRISPR"},
        )
        assert '"${PYTHON:-python3}" "$SKILL_SCRIPT"' in result
        assert "\npython \"$SKILL_SCRIPT\"" not in result

    def test_list_value_expands_to_repeated_flags(self):
        """A list value must emit the flag once per element, not a Python repr."""
        result = build_portable_commands_sh(
            skill_name="nfcore-rnaseq-wrapper",
            script_name="nfcore_rnaseq_wrapper.py",
            args={"--nextflow-config": ["/path/a.config", "/path/b.config"]},
        )
        assert "--nextflow-config /path/a.config" in result
        assert "--nextflow-config /path/b.config" in result
        assert "'/path/a.config'" not in result  # Python list repr must not appear


# ── write_portable_commands_sh tests ──────────────────────────────────────────

class TestWritePortableCommandsSh:

    def test_creates_commands_sh(self):
        with tempfile.TemporaryDirectory() as tmp:
            repro = Path(tmp) / "reproducibility"
            write_portable_commands_sh(
                repro_dir=repro,
                skill_name="lit-synthesizer",
                script_name="lit_synthesizer.py",
                args={"--query": "CRISPR", "--output": "./report"},
            )
            assert (repro / "commands.sh").exists()
