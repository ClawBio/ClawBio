"""Own this wrapper's bare modules without replacing the public package.

The real clawbio package imports the launcher without heavy runtime SDKs. Keep
it intact so shared, wrapper and public API tests can be collected together.
"""
from __future__ import annotations

import sys
from pathlib import Path

# ── Cross-skill bare-module isolation (shared verbatim with nfcore-sarek /
# nfcore-scrnaseq conftests) ─────────────────────────────────────────────────
# Several nf-core wrapper skills ship modules under bare top-level names
# (``errors``, ``schemas``, …). A single pytest session that collects more than
# one wrapper suite can shadow them across skills. We force THIS skill's dir to
# the front of sys.path and ensure each bare name resolves to this skill's
# *canonical* module object — restoring the same object (not reimporting) when a
# foreign copy is cached, so ``pytest.raises(SkillError)`` class identity holds.
_ISO_SKILL_DIR = Path(__file__).resolve().parent.parent
_ISO_TESTS_DIR = Path(__file__).resolve().parent
_ISO_BARE_MODULES = (
    "errors",
    "schemas",
    "preflight",
    "params_builder",
    "command_builder",
    "samplesheet_builder",
    "outputs_parser",
    "provenance",
    "reporting",
    "executor",
    "pipeline_source",
    "remap_paths",
    "repro_commands",
    "nfcore_4_1_0_contract",
    "_isolated_imports",
)
_ISO_CANONICAL: dict[str, object] = {}


def _claim_local_modules() -> None:
    for path in (_ISO_TESTS_DIR, _ISO_SKILL_DIR):
        if str(path) in sys.path:
            sys.path.remove(str(path))
        sys.path.insert(0, str(path))
    for name in _ISO_BARE_MODULES:
        module = sys.modules.get(name)
        if module is not None:
            module_file = Path(getattr(module, "__file__", "") or "")
            local = (
                module_file == _ISO_SKILL_DIR / f"{name}.py"
                or _ISO_SKILL_DIR in module_file.parents
            )
            if local:
                _ISO_CANONICAL[name] = module
                continue
        if name in _ISO_CANONICAL:
            sys.modules[name] = _ISO_CANONICAL[name]
        elif module is not None:
            sys.modules.pop(name, None)


_claim_local_modules()


def pytest_collectstart(collector):  # noqa: ARG001 - pytest hook signature
    _claim_local_modules()


def pytest_runtest_setup(item):  # noqa: ARG001 - pytest hook signature
    _claim_local_modules()
