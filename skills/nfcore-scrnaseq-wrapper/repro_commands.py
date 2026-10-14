"""Portable, nextflow-direct ``commands.sh`` generator for nfcore-scrnaseq-wrapper.

Mirrors the nf-core/sarek reproducibility model: the replay command runs
``nextflow run <pipeline> -r <version> -params-file reproducibility/params.yaml``
directly, rather than re-invoking the ClawBio python wrapper. This decouples
replay from wrapper-code stability and removes every machine-specific absolute
path from the script — the only machine-dependent inputs (data + reference
paths) live in params.yaml / the samplesheet and are rewritten by remap_paths.py.

The script self-anchors: it resolves the output directory from its own location
(``reproducibility/commands.sh`` → output_dir is the parent of the bundle dir),
cds there, and uses paths relative to the output directory throughout.
"""

from __future__ import annotations

import shlex
import sys
from pathlib import Path
from typing import Any

from clawbio.common.nfcore import build_runtime_resource_guard
from clawbio.common.textio import write_text_lf

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))
sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("nfcore_4_1_0_contract")

from nfcore_4_1_0_contract import MACOS_DEMO_RESOURCE_LIMITS, STAR_ALIGN_BASE_EXT_ARGS

_HEADER = """\
#!/usr/bin/env bash
# ClawBio reproducibility bundle — portable replay command (nextflow-direct)
# Generated: {generated_at}
# Skill: nfcore-scrnaseq-wrapper
#
# How to replay (from anywhere):
#   bash reproducibility/commands.sh
#
# This runs the pinned nf-core/scrnaseq pipeline directly against the bundled
# params.yaml. No ClawBio installation is required to replay.

set -euo pipefail

# ── Locate the run directory ──────────────────────────────────────────────────
# This script lives in <output_dir>/reproducibility/. Resolve <output_dir> from
# the script's own location so replay works regardless of where the bundle was
# copied, then run Nextflow from there (params.yaml uses output-relative paths).
BUNDLE_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
OUTPUT_DIR="$(cd "$BUNDLE_DIR/.." && pwd)"
cd "$OUTPUT_DIR"

if [[ ! -f "reproducibility/params.yaml" ]]; then
  echo "ERROR: reproducibility/params.yaml not found next to this script." >&2
  exit 1
fi
"""

_LOCAL_WARNING = """\

# ── WARNING: non-portable pipeline source ─────────────────────────────────────
# This run used a LOCAL pipeline checkout, which is NOT portable to other
# machines. The path below only exists on the original machine. To replay
# elsewhere, edit the 'nextflow run' line to use the remote form instead:
#     nextflow run nf-core/scrnaseq -r {resolved_version} \\
# Local checkout path was: {source_ref}
"""

# Emitted (only when macos_docker_config=True) before the run command. The macOS
# Docker workarounds must apply based on the REPLAY host, not the host that
# generated the bundle: a macOS-generated bundle replayed on Linux must NOT load
# macOS config (it would force --platform linux/amd64 — wrong on arm64 Linux —
# and alter STAR/resource behaviour), and a Linux-generated bundle replayed on
# macOS must still load it (or the run fails with VirtioFS EDEADLK / FIFO errors).
# A plain, always-assigned string variable is used (not a bash array) because
# macOS ships bash 3.2, where `"${arr[@]}"` on an empty array errors under set -u.
_MACOS_GUARD = """\
# Apply macOS-only Docker workarounds (VirtioFS EDEADLK, Apple Silicon
# --platform, STAR FIFOs, 4 h test cap) ONLY when replaying on macOS, so the
# same bundle runs identically on Linux and on macOS.
EXTRA_CONFIG=""
if [[ "$(uname -s)" == "Darwin" ]]; then
  EXTRA_CONFIG="-c reproducibility/macos_docker.config"
fi

"""

# Idempotent replay: reuse the previous run's cache when replaying against an
# already-completed output dir, but start fresh on a first run. Nextflow keeps its
# session history in the launch dir's .nextflow/, so its presence means a prior run
# happened here. -resume is a no-op (with only a warning) on a fresh dir, but the
# guard keeps first-run output clean and mirrors the nfcore-rnaseq bundle. A plain
# assigned string (not an array) keeps macOS's bash 3.2 happy under `set -u`.
_RESUME_GUARD = """\
RESUME=""
if [[ -d ".nextflow" ]]; then
  RESUME="-resume"
fi

"""

# Preserve the recorded ceilings while bounding them to the replay Docker VM.
# The checksummed source config stays unchanged on every platform.
_RESOURCE_GUARD = build_runtime_resource_guard()

# Nextflow config with the macOS + Apple-Silicon Docker workarounds. Always
# shipped inside docker bundles; applied at replay only under the uname guard
# above. See nfcore_scrnaseq_wrapper for the full rationale of each fix.
#
# Two concerns are kept strictly separate (audit H-1):
#   * Workarounds that are ALWAYS safe (VirtioFS stageInMode, Apple-Silicon
#     --platform, STAR FIFO --outTmpDir) — emitted for every macOS+docker run.
#   * resourceLimits copied from the nf-core test profile — emitted ONLY for
#     ``--demo`` runs. Applying a 4-CPU / 15 GB ceiling to real datasets would
#     starve STAR (a human index needs far more than 15 GB) and silently break
#     production runs on the wrapper's primary platform.
_MACOS_CONFIG_HEADER = """\
// macOS + Docker workaround for VirtioFS EDEADLK (errno 35), ARM64 hosts,
// and STAR FIFOs.
process {
    stageInMode = "copy"
"""

# Emitted only for demo runs. See MACOS_DEMO_RESOURCE_LIMITS in the 4.1.0 contract.
_MACOS_DEMO_RESOURCE_BLOCK = """\
    // DEMO/TEST ONLY: the nf-core test profile caps tasks at 1 h, too short for
    // STAR genome generation under emulation, so the cap is raised to 4 h. These
    // ceilings match the small test machine and MUST NOT be applied to real data.
    resourceLimits = [
        cpus: {cpus},
        memory: '{memory}',
        time: '{time}'
    ]
"""

# STAR FIFO workaround. The base flags are the pinned STAR_ALIGN_BASE_EXT_ARGS
# (verbatim from conf/modules.config); only --outTmpDir is added here so the
# override never silently drops an upstream flag. Do NOT use a closure that
# references task.ext.args — that causes a StackOverflowError.
_MACOS_CONFIG_STAR_AND_DOCKER = """\
    // STAR creates FIFOs in _STARtmp; VirtioFS does not support FIFOs.
    // --outTmpDir /tmp/star_tmp routes _STARtmp to the container's /tmp
    // (Linux tmpfs), which does support FIFOs.
    withName: '.*STAR_ALIGN.*' {{
        ext.args = {{ "{star_args} --outTmpDir /tmp/star_tmp" }}
    }}
}}
docker {{
    runOptions = "--platform linux/amd64"
}}
"""


def build_macos_docker_config(*, demo: bool) -> str:
    """Return the macOS+Docker Nextflow config text.

    The always-safe workarounds (stageInMode, --platform, STAR --outTmpDir) are
    always included. The test-profile ``resourceLimits`` ceilings are included
    only when ``demo`` is True (audit H-1).
    """
    config = _MACOS_CONFIG_HEADER
    if demo:
        config += _MACOS_DEMO_RESOURCE_BLOCK.format(**MACOS_DEMO_RESOURCE_LIMITS)
    config += _MACOS_CONFIG_STAR_AND_DOCKER.format(star_args=STAR_ALIGN_BASE_EXT_ARGS)
    return config


def write_macos_docker_config(output_dir: Path, *, demo: bool = False) -> Path:
    """Write reproducibility/macos_docker.config (LF, byte-stable). Returns its path.

    Shared by the live run (applied only on macOS+docker) and the bundle (always
    shipped for docker so commands.sh can apply it on macOS replay hosts). ``demo``
    gates the test-profile resource ceilings (see build_macos_docker_config).
    """
    config_path = Path(output_dir) / "reproducibility" / "macos_docker.config"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    write_text_lf(config_path, build_macos_docker_config(demo=demo))
    return config_path


def _pipeline_invocation(pipeline_source: dict[str, Any]) -> str:
    source_kind = str(pipeline_source["source_kind"])
    source_ref = str(pipeline_source["source_ref"])
    resolved_version = str(pipeline_source["resolved_version"])
    if source_kind == "local_checkout":
        # Emit the absolute checkout path (only valid on the original machine);
        # a loud warning block is appended separately by the caller.
        return shlex.quote(source_ref)
    return f"{shlex.quote(source_ref)} -r {shlex.quote(resolved_version)}"


def build_nextflow_commands_sh(
    *,
    pipeline_source: dict[str, Any],
    profile: str,
    demo: bool,
    macos_docker_config: bool,
    nextflow_version: str | None,
    generated_at: str,
    extra_configs: list[str] | None = None,
    work_dir: str = "upstream/work",
    resource_limits_config: bool = False,
) -> str:
    """Return the full contents of a portable, nextflow-direct ``commands.sh``.

    Parameters
    ----------
    pipeline_source:
        The resolved pipeline source dict (``source_kind``, ``source_ref``,
        ``resolved_version``). Remote sources are fully portable; local
        checkouts emit a warning block.
    profile:
        Backend profile (e.g. ``docker``). For demo runs the ``test`` profile
        is prepended (``test,docker``) to match the live invocation.
    demo:
        Prepend the ``test`` profile when True.
    macos_docker_config:
        When True, emit a ``uname``-gated block that applies
        ``reproducibility/macos_docker.config`` only when the *replay* host is
        macOS. The config file itself is always shipped in docker bundles, so
        the same bundle runs identically on Linux and macOS.
    nextflow_version:
        When given, pin the Nextflow engine via ``export NXF_VER=<version>`` so
        replay uses the exact engine that produced the run (the pipeline itself
        is already pinned with ``-r``). The engine version can change execution,
        so this is required for byte-identical reproduction.
    generated_at:
        Timestamp string for the header comment.
    """
    effective_profile = f"test,{profile}" if demo else profile
    invocation = _pipeline_invocation(pipeline_source)

    run_lines = [
        f"nextflow run {invocation} \\",
        f"  -profile {effective_profile} \\",
        "  -params-file reproducibility/params.yaml \\",
        f"  -work-dir {shlex.quote(work_dir)}",
    ]
    # Append optional flags as continuation lines on the last fixed line.
    # $EXTRA_CONFIG / $RESUME are unquoted so they word-split to a flag or to
    # nothing; both are always assigned, so they are safe under `set -u`.
    tail: list[str] = []
    tail.append("  $RESUME")
    if macos_docker_config:
        tail.append("  $EXTRA_CONFIG")
    if resource_limits_config:
        tail.append("  $RESOURCE_CONFIG")
    for config_path in extra_configs or []:
        tail.append(f"  -c {shlex.quote(config_path)}")
    if tail:
        run_lines[-1] = run_lines[-1] + " \\"
        for i, line in enumerate(tail):
            run_lines.append(line + (" \\" if i < len(tail) - 1 else ""))

    parts = [
        _HEADER.format(generated_at=generated_at),
        "\n# ── Replay command ────────────────────────────────────────────────────────────\n",
        _RESUME_GUARD,
    ]
    if resource_limits_config:
        parts.append(_RESOURCE_GUARD)
    if nextflow_version:
        # Pin the Nextflow engine version (the pipeline is pinned via -r). Nextflow
        # auto-fetches this exact version, so the engine cannot drift on replay.
        parts.append(
            "# Pin the Nextflow engine to the version that produced this run.\n"
            f'export NXF_VER="{nextflow_version}"\n\n'
        )
    if macos_docker_config:
        parts.append(_MACOS_GUARD)
    parts.append("\n".join(run_lines))
    parts.append("\n")
    if str(pipeline_source["source_kind"]) == "local_checkout":
        parts.append(
            _LOCAL_WARNING.format(
                resolved_version=str(pipeline_source["resolved_version"]),
                source_ref=str(pipeline_source["source_ref"]),
            )
        )
    return "".join(parts)
