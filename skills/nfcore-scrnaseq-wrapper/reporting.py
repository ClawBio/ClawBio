from __future__ import annotations

from clawbio.common.nfcore.artifacts import safe_config_basename, stage_config_files

import json
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from clawbio.common.report import (
    generate_report_footer,
    generate_report_header,
    write_result_json,
)
from clawbio.common.textio import write_text_lf

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))
sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("repro_commands", "schemas")

from repro_commands import build_nextflow_commands_sh, write_macos_docker_config
from schemas import (
    ALL_REFERENCE_PATH_FIELDS,
    SKILL_ALIAS,
    SKILL_DIR,
    SKILL_NAME,
    SKILL_VERSION,
    profile_includes,
)

# Repository-relative (not an absolute path resolved on the generating machine): an
# absolute path baked into the report would not exist on any other machine. Matches
# the `clawbio.py run ...` convention documented across the repo. Invoked with
# `python3` for portability to python3-only systems (PEP 394).
_CLAWBIO_SCRIPT = "clawbio.py"


def write_report(
    output_dir: Path,
    *,
    args,
    pipeline_source: dict[str, object],
    preflight_result: dict[str, Any],
    parsed_outputs: dict[str, Any],
    command_str: str,
) -> Path:
    lines = build_report_lines(
        output_dir,
        args=args,
        pipeline_source=pipeline_source,
        preflight_result=preflight_result,
        parsed_outputs=parsed_outputs,
        command_str=command_str,
    )
    report_path = output_dir / "report.md"
    write_text_lf(report_path, "\n".join(lines))
    return report_path


def build_report_lines(
    output_dir: Path,
    *,
    args,
    pipeline_source: dict[str, object],
    preflight_result: dict[str, Any],
    parsed_outputs: dict[str, Any],
    command_str: str,
) -> list[str]:
    header = generate_report_header(
        "nf-core/scrnaseq Wrapper Report",
        SKILL_NAME,
        SKILL_VERSION,
        extra_metadata={
            "Preset": args.preset,
            "Profile": args.profile,
            "Pipeline source": str(pipeline_source["source_kind"]),
            "Pipeline ref": str(pipeline_source["resolved_version"]),
        },
    )
    preferred_h5ad = str(parsed_outputs.get("preferred_h5ad", ""))
    sample_count = _reported_sample_count(preflight_result, parsed_outputs)
    return [
        header,
        "## Summary",
        "",
        f"- Preset: `{args.preset}`",
        f"- Effective aligner: `{parsed_outputs.get('aligner_effective', '')}`",
        f"- Pipeline source: `{pipeline_source['source_kind']}`",
        f"- Pipeline ref: `{pipeline_source['resolved_version']}`",
        f"- Output root: `{output_dir}`",
        "",
        "## Preflight",
        "",
        f"- Java: `{preflight_result['java']['version']}`",
        f"- Nextflow: `{preflight_result['nextflow']['version']}`",
        f"- Backend profile: `{args.profile}`",
        f"- Samples: `{sample_count}`",
        "",
        "## Outputs",
        "",
        f"- Preferred h5ad: `{preferred_h5ad or 'not available'}`",
        f"- MultiQC report: `{parsed_outputs.get('multiqc_report', '') or 'not found'}`",
        f"- Pipeline info: `{parsed_outputs.get('pipeline_info_dir', '') or 'not found'}`",
        f"- CellBender detected: `{parsed_outputs.get('cellbender_used', False)}`",
        "",
        "## Reproducibility",
        "",
        f"- Command: `{command_str}`",
        f"- Repro bundle: `{output_dir / 'reproducibility'}`",
        "",
        *_build_handoff_lines(preferred_h5ad),
        generate_report_footer(),
    ]


def _build_handoff_lines(preferred_h5ad: str) -> list[str]:
    if preferred_h5ad:
        return [
            "## Next Steps",
            "",
            f"- `python3 {_CLAWBIO_SCRIPT} run scrna --input {preferred_h5ad} --output <dir>`",
            f"- `python3 {_CLAWBIO_SCRIPT} run scrna-embedding --input {preferred_h5ad} --output <dir>`",
            "",
        ]
    return [
        "## Next Steps",
        "",
        "- No canonical `.h5ad` was selected automatically. Inspect `result.json` and the upstream outputs before chaining downstream skills.",
        "",
    ]


def _reported_sample_count(
    preflight_result: dict[str, Any], parsed_outputs: dict[str, Any]
) -> int:
    samplesheet = preflight_result.get("samplesheet", {})
    if isinstance(samplesheet, dict):
        sample_count = samplesheet.get("sample_count")
        if isinstance(sample_count, str) and sample_count.isdigit():
            sample_count = int(sample_count)
        # Only trust a *positive* local samplesheet count. Under --demo the
        # samplesheet is supplied remotely by the upstream test profile, so the
        # local count is 0 — fall back to the samples detected in the outputs
        # rather than reporting 0 for a run that produced real samples.
        if isinstance(sample_count, int) and sample_count > 0:
            return sample_count
    samples_detected = parsed_outputs.get("samples_detected", [])
    return len(samples_detected) if isinstance(samples_detected, list) else 0


_PORTABILITY_NOTICE = """\

# ── Portability notice ────────────────────────────────────────────────────────
# Replaying on a different machine? The pipeline itself is fetched from nf-core,
# but your data and reference paths are machine-specific. Before replaying:
#
#   1. Remap FASTQ paths in the samplesheet:
#        python3 reproducibility/remap_paths.py --old /original/prefix --new /new/prefix
#
#   2. Remap reference/index paths in params.yaml (--fasta, --gtf, indexes, ...):
#        python3 reproducibility/remap_paths.py --refs-old /original/refs --refs-new /new/refs
#
#   3. Verify everything resolves on this machine:
#        python3 reproducibility/remap_paths.py --verify
#
# The output directory is auto-detected from this script's location — no
# --output edit is needed when you move the whole bundle.
"""

_REMAP_SCRIPT_SRC = SKILL_DIR / "remap_paths.py"


def write_repro_commands(
    output_dir: Path,
    *,
    args,
    pipeline_source: dict[str, object],
    nextflow_version: str | None = None,
) -> None:
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)
    # Docker bundles always ship macos_docker.config; commands.sh applies it only
    # on macOS replay hosts (uname-gated), so the same bundle is OS-independent.
    macos_docker_config = profile_includes(args.profile, "docker")
    if macos_docker_config:
        # Match the live run: test-profile resource ceilings only for demo bundles.
        write_macos_docker_config(output_dir, demo=bool(getattr(args, "demo", False)))
    copied_extra_configs = _copy_user_nextflow_configs(
        repro_dir, getattr(args, "extra_config", []) or []
    )
    # The host-scaled resourceLimits cap is written into the bundle by the live run
    # (only for non-macOS docker real runs). Ship+apply it on replay when present so a
    # from-scratch reproduction on the generating host does not re-abort on memory.
    resource_limits_config = (repro_dir / "resource_limits.config").is_file()
    script = build_nextflow_commands_sh(
        pipeline_source=pipeline_source,
        profile=args.profile,
        demo=bool(getattr(args, "demo", False)),
        macos_docker_config=macos_docker_config,
        nextflow_version=nextflow_version,
        generated_at=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        extra_configs=copied_extra_configs,
        work_dir=_replay_work_dir(args),
        resource_limits_config=resource_limits_config,
    )
    if not getattr(args, "demo", False):
        script += _PORTABILITY_NOTICE
    commands_sh = repro_dir / "commands.sh"
    write_text_lf(commands_sh, script)
    _write_remap_script(repro_dir)


def _copy_user_nextflow_configs(repro_dir: Path, config_paths: list[str]) -> list[str]:
    if not config_paths:
        return []
    config_dir = repro_dir / "nextflow_configs"
    config_dir.mkdir(parents=True, exist_ok=True)
    copied: list[str] = []
    for index, raw_path in enumerate(config_paths, start=1):
        source = Path(raw_path).expanduser().resolve()
        destination = (
            config_dir / f"config_{index:02d}_{_safe_config_basename(source.name)}"
        )
        stage_config_files([source], destination_dir=destination.parent, filenames=[destination.name])
        copied.append(f"reproducibility/nextflow_configs/{destination.name}")
    return copied


def _replay_work_dir(args) -> str:
    raw_work_dir = getattr(args, "work_dir", None)
    if not raw_work_dir:
        return "upstream/work"
    work_dir = str(raw_work_dir).strip()
    if "://" in work_dir:
        return work_dir
    return Path(work_dir).expanduser().resolve().as_posix()


def _safe_config_basename(name: str) -> str:
    return safe_config_basename(name)


_REMAP_AUTOGEN_START = "# >>> AUTO-GENERATED reference keys (do not edit by hand) >>>"
_REMAP_AUTOGEN_END = "# <<< AUTO-GENERATED reference keys <<<"
_REMAP_AUTOGEN_RE = re.compile(
    re.escape(_REMAP_AUTOGEN_START) + r".*?" + re.escape(_REMAP_AUTOGEN_END),
    re.DOTALL,
)


def _render_reference_keys_block() -> str:
    """Render the _PARAMS_REFERENCE_KEYS tuple from the single canonical source."""
    keys = "".join(f'    "{field}",\n' for field in ALL_REFERENCE_PATH_FIELDS)
    return (
        f"{_REMAP_AUTOGEN_START}\n"
        f"_PARAMS_REFERENCE_KEYS = (\n{keys})\n"
        f"{_REMAP_AUTOGEN_END}"
    )


def _write_remap_script(repro_dir: Path) -> None:
    """Ship remap_paths.py into the bundle, regenerating its reference-key block
    from schemas.ALL_REFERENCE_PATH_FIELDS so the standalone bundle copy can never
    drift from the canonical list (the bundle has no access to schemas at replay)."""
    repro_dir.mkdir(parents=True, exist_ok=True)
    source = _REMAP_SCRIPT_SRC.read_text(encoding="utf-8")
    regenerated, n = _REMAP_AUTOGEN_RE.subn(
        lambda _m: _render_reference_keys_block(), source, count=1
    )
    if n != 1:
        # Sentinels missing/duplicated — fail loud rather than ship a stale copy.
        raise RuntimeError(
            f"Expected exactly one AUTO-GENERATED reference-keys block in {_REMAP_SCRIPT_SRC}, found {n}."
        )
    write_text_lf(repro_dir / "remap_paths.py", regenerated)


def write_result(
    output_dir: Path,
    *,
    args,
    pipeline_source: dict[str, object],
    parsed_outputs: dict[str, Any],
    command_str: str,
) -> Path:
    output_artifacts = build_output_artifacts(parsed_outputs)
    summary = build_result_summary(
        args=args,
        pipeline_source=pipeline_source,
        parsed_outputs=parsed_outputs,
        output_artifacts=output_artifacts,
    )
    data = build_result_data(
        parsed_outputs=parsed_outputs,
        output_artifacts=output_artifacts,
        command_str=command_str,
    )
    return write_result_json(
        output_dir,
        skill=SKILL_ALIAS,
        version=SKILL_VERSION,
        # Minimal shared cross-wrapper contract: a successful run is ok/"ok".
        status="ok",
        ok=True,
        summary=summary,
        data=data,
    )


def build_output_artifacts(parsed_outputs: dict[str, Any]) -> dict[str, object]:
    return {
        "preferred_h5ad": parsed_outputs.get("preferred_h5ad", ""),
        "multiqc_report": parsed_outputs.get("multiqc_report", ""),
        "pipeline_info_dir": parsed_outputs.get("pipeline_info_dir", ""),
        "official_outputs": parsed_outputs.get("official_outputs", {}),
        "h5ad_candidates": parsed_outputs.get("h5ad_candidates", []),
        "rds_candidates": parsed_outputs.get("rds_candidates", []),
    }


def build_result_summary(
    *,
    args,
    pipeline_source: dict[str, object],
    parsed_outputs: dict[str, Any],
    output_artifacts: dict[str, object],
) -> dict[str, object]:
    return {
        "preset": args.preset,
        "aligner_effective": parsed_outputs.get("aligner_effective", ""),
        "pipeline_source_kind": pipeline_source["source_kind"],
        "pipeline_version_or_commit": pipeline_source["resolved_version"],
        "profile": args.profile,
        "resume_used": bool(args.resume),
        "multiqc_report": parsed_outputs.get("multiqc_report", ""),
        "pipeline_info_dir": parsed_outputs.get("pipeline_info_dir", ""),
        "preferred_h5ad": parsed_outputs.get("preferred_h5ad", ""),
        "handoff_available": parsed_outputs.get("handoff_available", False),
        "samples_detected": len(parsed_outputs.get("samples_detected", [])),
        "cellbender_used": parsed_outputs.get("cellbender_used", False),
        "output_artifacts": output_artifacts,
    }


def build_result_data(
    *,
    parsed_outputs: dict[str, Any],
    output_artifacts: dict[str, object],
    command_str: str,
) -> dict[str, object]:
    return {
        "canonical_skill_name": SKILL_NAME,
        "cli_alias": SKILL_ALIAS,
        "command": command_str,
        "output_artifacts": output_artifacts,
        "outputs": parsed_outputs,
    }


def write_check_result(output_dir: Path, payload: dict[str, object]) -> Path:
    path = output_dir / "check_result.json"
    write_text_lf(path, json.dumps(payload, indent=2))
    return path
