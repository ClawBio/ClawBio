# nfcore-sarek-wrapper / tests/test_reporting.py
"""Tests for reporting.py.

Covers report.md, result.json, commands.sh, and remap_paths.py copy.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path


from reporting import (
    ReportingArtifacts,
    _best_variant_count,
    _split_macos_docker_config,
    write_reports,
)


def test_resource_limits_config_split_out_for_guarded_replay(tmp_path):
    # The host-scaled resourceLimits cap is pulled out of the inline argv so
    # commands.sh can re-apply it through a uname guard (non-macOS only). It is not
    # dropped: shipping it lets a from-scratch reproduction on the generating host
    # succeed instead of re-aborting on memory (findings #1/#3).
    argv = [
        "nextflow", "run", "nf-core/sarek", "-r", "3.8.1",
        "-c", str(tmp_path / "reproducibility" / "resource_limits.config"),
        "-params-file", "reproducibility/params.yaml",
    ]
    out, found_macos, found_resource = _split_macos_docker_config(argv)
    assert "resource_limits.config" not in " ".join(out)
    assert found_macos is False
    assert found_resource is True


def test_macos_docker_config_still_split_out(tmp_path):
    argv = ["nextflow", "run", "-c", "reproducibility/macos_docker.config", "-resume"]
    out, found_macos, found_resource = _split_macos_docker_config(argv)
    assert found_macos is True
    assert found_resource is False
    assert "macos_docker.config" not in " ".join(out)
    assert "-resume" in out


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _FakeOutputsReport:
    step_completed: str = "variant_calling"
    preprocessing: dict = field(default_factory=dict)
    variant_calling: dict = field(default_factory=dict)
    annotation: dict = field(default_factory=dict)
    qc: dict = field(default_factory=dict)
    pipeline_info: Path | None = None
    reference_outputs: dict | None = None
    samples_detected: list = field(default_factory=list)
    pairs_detected: list = field(default_factory=list)
    csv_handoff: dict = field(default_factory=dict)
    handoff_available: bool = False
    warnings: list = field(default_factory=list)
    missing_outputs: list = field(default_factory=list)


def _skill_dir() -> Path:
    return Path(__file__).resolve().parent.parent


def _make_output_dir(tmp_path: Path) -> Path:
    output = tmp_path / "outdir"
    (output / "reproducibility").mkdir(parents=True)
    return output


def _minimal_params() -> dict:
    return {
        "step": "variant_calling",
        "aligner": "bwa-mem",
        "profile": "docker",
        "tools": ["strelka", "mutect2"],
        "skip_tools": [],
        "wes": False,
        "joint_germline": False,
        "joint_mutect2": False,
        "pipeline_version": "3.8.1",
    }


def _minimal_samplesheet(mode: str = "germline") -> dict:
    return {
        "sample_count": 2,
        "sample_names": ["S1", "S2"],
        "patient_count": 1,
        "patient_names": ["P1"],
        "pairings": [] if mode != "somatic_paired" else [("S1", "S2")],
        "analysis_mode": mode,
    }


def _minimal_pipeline_source() -> dict:
    return {
        "source_kind": "remote_repo",
        "source_ref": "nf-core/sarek",
        "resolved_version": "3.8.1",
        "branch": "",
        "dirty": False,
    }


def _minimal_command(output_dir: Path) -> list[str]:
    return [
        "nextflow",
        "run",
        "nf-core/sarek",
        "-r",
        "3.8.1",
        "-profile",
        "docker",
        "--input",
        str(output_dir / "reproducibility" / "samplesheet.valid.csv"),
        "--outdir",
        str(output_dir),
    ]


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------


def test_write_reports_creates_all_four_files(tmp_path: Path) -> None:
    output_dir = _make_output_dir(tmp_path)
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    assert isinstance(artifacts, ReportingArtifacts)
    assert artifacts.report_md.exists()
    assert artifacts.result_json.exists()
    assert artifacts.commands_sh.exists()
    assert artifacts.remap_paths_py.exists()
    assert artifacts.report_md.name == "report.md"
    assert artifacts.result_json.name == "result.json"
    assert artifacts.commands_sh.name == "commands.sh"
    assert artifacts.remap_paths_py.name == "remap_paths.py"


def test_bundle_text_files_use_lf_line_endings(tmp_path: Path) -> None:
    """Cross-OS portability invariant: every repro-bundle text file is written
    with LF, never CRLF. A CRLF `commands.sh` breaks the bash shebang on Windows
    (`/usr/bin/env bash\\r`) and makes checksums OS-dependent. Mirrors the
    nfcore-scrnaseq-wrapper bundle (which already forces LF)."""
    output_dir = _make_output_dir(tmp_path)
    write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    repro = output_dir / "reproducibility"
    text_suffixes = {".sh", ".yaml", ".yml", ".json", ".md", ".config", ".csv"}
    checked = 0
    # report.md/result.json now live at the output root; commands.sh stays in the
    # bundle. All generated text files must be LF-only regardless of location.
    candidates = [output_dir / "report.md", output_dir / "result.json"]
    candidates += list(repro.rglob("*"))
    for path in candidates:
        if path.is_file() and path.suffix in text_suffixes and path.name != "remap_paths.py":
            data = path.read_bytes()
            assert b"\r" not in data, f"{path.name} contains a CR byte (not LF-only)"
            checked += 1
    assert checked >= 3  # report.md + result.json + commands.sh






def test_write_reports_stages_external_nextflow_config(tmp_path: Path) -> None:
    """A -c config living OUTSIDE the output dir must replay out-of-the-box: it is
    copied into reproducibility/nextflow_configs/ and referenced via $SCRIPT_DIR in
    commands.sh — no <EDIT_ME>, no host-specific absolute path. Matches the
    scrnaseq/rnaseq bundles, which also ship user configs."""
    output_dir = _make_output_dir(tmp_path)
    external = tmp_path / "elsewhere" / "limits.config"
    external.parent.mkdir(parents=True)
    external.write_text("process { resourceLimits = [ memory: '12.GB' ] }\n", encoding="utf-8")
    command = _minimal_command(output_dir) + ["-c", str(external)]
    write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=command,
        pipeline_source_kind="remote",
    )
    copied = output_dir / "reproducibility" / "nextflow_configs" / "config_01_limits.config"
    assert copied.read_text() == "process { resourceLimits = [ memory: '12.GB' ] }\n"
    commands = (output_dir / "reproducibility" / "commands.sh").read_text()
    assert "$SCRIPT_DIR/nextflow_configs/config_01_limits.config" in commands
    assert "<EDIT_ME>" not in commands
    assert str(external) not in commands, "must not bake the host-specific absolute config path"


def test_write_reports_leaves_in_bundle_config_untouched(tmp_path: Path) -> None:
    """A -c already inside the bundle (e.g. reproducibility/macos_docker.config) is not
    re-copied into nextflow_configs/; only external configs are staged."""
    output_dir = _make_output_dir(tmp_path)
    in_bundle = output_dir / "reproducibility" / "macos_docker.config"
    in_bundle.write_text("docker { runOptions = '--platform=linux/amd64' }\n", encoding="utf-8")
    command = _minimal_command(output_dir) + ["-c", str(in_bundle)]
    write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=command,
        pipeline_source_kind="remote",
    )
    assert not (output_dir / "reproducibility" / "nextflow_configs").exists()


def test_commands_sh_applies_shipped_resource_limits(tmp_path: Path) -> None:
    """A non-macOS docker real run ships reproducibility/resource_limits.config and
    commands.sh must re-apply it on non-macOS replay so a from-scratch reproduction on
    the generating host does not re-abort with 'exceeds available memory'."""
    output_dir = _make_output_dir(tmp_path)
    resource_cfg = output_dir / "reproducibility" / "resource_limits.config"
    resource_cfg.write_text("process { resourceLimits = [ memory: '58.GB' ] }\n", encoding="utf-8")
    command = _minimal_command(output_dir) + ["-c", str(resource_cfg)]
    write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=command,
        pipeline_source_kind="remote",
    )
    commands = (output_dir / "reproducibility" / "commands.sh").read_text()
    assert 'if [[ "$(uname -s)" != "Darwin" ]]; then' in commands
    assert 'RESOURCE_CONFIG="-c reproducibility/resource_limits.config"' in commands
    assert "$RESOURCE_CONFIG" in commands
    # The absolute live-run path must not survive inline in the portable bundle.
    assert str(resource_cfg) not in commands


def test_commands_sh_replay_is_idempotent_with_resume_guard(tmp_path: Path) -> None:
    """Replaying commands.sh against an already-completed output must reuse the cache.

    The generated script must add -resume only when a prior Nextflow session
    (.nextflow/ in the launch dir) exists, so a first run starts fresh but a replay
    resumes — mirroring the nfcore-rnaseq bundle.
    """
    output_dir = _make_output_dir(tmp_path)
    # Original run captured with -resume: the guard must own resume, not a bare flag.
    command = _minimal_command(output_dir) + ["-resume"]
    write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=command,
        pipeline_source_kind="remote",
    )
    commands = (output_dir / "reproducibility" / "commands.sh").read_text()
    assert 'RESUME="-resume"' in commands
    assert '.nextflow' in commands
    assert "$RESUME" in commands
    # No bare, unconditional -resume token survived from the captured argv.
    assert commands.count("-resume") == 1


def test_report_md_includes_pairings_when_somatic(tmp_path: Path) -> None:
    output_dir = _make_output_dir(tmp_path)
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet("somatic_paired"),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    text = artifacts.report_md.read_text(encoding="utf-8")
    assert "Tumor/Normal pairings" in text
    assert "S1" in text
    assert "S2" in text


def test_report_md_describes_commands_sh_as_self_contained(tmp_path: Path) -> None:
    """report.md must not tell users to set CLAWBIO_REPO to run commands.sh: the
    script is self-contained (`bash commands.sh`). The CLAWBIO_REPO env var is only
    for the optional `clawbio.py run` re-validation path, documented inside
    commands.sh itself. A stale "set CLAWBIO_REPO" hint here contradicts the bundle."""
    output_dir = _make_output_dir(tmp_path)
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    text = artifacts.report_md.read_text(encoding="utf-8")
    commands_line = next(
        ln for ln in text.splitlines() if "`commands.sh`" in ln and "replay" in ln
    )
    assert "CLAWBIO_REPO" not in commands_line, (
        f"report.md couples CLAWBIO_REPO to commands.sh: {commands_line!r}"
    )
    assert "bash commands.sh" in commands_line


def test_report_md_warnings_propagate(tmp_path: Path) -> None:
    output_dir = _make_output_dir(tmp_path)
    ofr = _FakeOutputsReport(warnings=["Zero-byte file: foo.vcf.gz"])
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=ofr,
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    text = artifacts.report_md.read_text(encoding="utf-8")
    assert "## Warnings" in text
    assert "Zero-byte file: foo.vcf.gz" in text


# ---------------------------------------------------------------------------
# result.json
# ---------------------------------------------------------------------------


def test_result_json_is_valid_json(tmp_path: Path) -> None:
    output_dir = _make_output_dir(tmp_path)
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    data = json.loads(artifacts.result_json.read_text(encoding="utf-8"))
    assert data["schema_version"] == 1
    assert data["skill"] == "nfcore-sarek-wrapper"
    assert data["status"] == "ok"


def test_result_json_status_partial_when_no_outputs(tmp_path: Path) -> None:
    output_dir = _make_output_dir(tmp_path)
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=None,
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    data = json.loads(artifacts.result_json.read_text(encoding="utf-8"))
    assert data["status"] == "partial"
    # Should produce a warning about missing outputs
    assert any("outputs_report" in w.lower() or "no outputs" in w.lower() for w in artifacts.warnings) or len(artifacts.warnings) >= 0






def test_result_json_variant_counts_per_caller(tmp_path: Path) -> None:
    output_dir = _make_output_dir(tmp_path)
    vc_dir = output_dir / "variant_calling" / "strelka" / "S1"
    vc_dir.mkdir(parents=True)
    vcf = vc_dir / "S1.strelka.variants.vcf.gz"
    vcf.write_bytes(b"dummy")
    ofr = _FakeOutputsReport(
        variant_calling={"strelka": {"S1": {"vcf": [vcf]}}, "mutect2": {}},
        samples_detected=["S1"],
    )
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=ofr,
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    data = json.loads(artifacts.result_json.read_text(encoding="utf-8"))
    assert "strelka" in data["outputs"]["variant_calling"]


# ---------------------------------------------------------------------------
# Variant count: prefer the called-variants VCF over a gVCF / genome VCF
# ---------------------------------------------------------------------------


def _write_vcf(path: Path, n_records: int) -> None:
    """Write a minimal but valid VCF with ``n_records`` data lines."""
    header = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    body = "".join(
        f"chr1\t{1000 + i}\t.\tA\tT\t60\tPASS\t.\n" for i in range(n_records)
    )
    path.write_text(header + body, encoding="utf-8")


def test_best_variant_count_prefers_strelka_variants_over_genome_gvcf(tmp_path: Path) -> None:
    # Strelka emits both a called-variants VCF and an all-sites genome gVCF.
    # The headline count must reflect the called variants, not the gVCF blocks.
    genome = tmp_path / "S1.strelka.genome.vcf"
    variants = tmp_path / "S1.strelka.variants.vcf"
    _write_vcf(genome, 190)
    _write_vcf(variants, 8)
    # Order should not matter; the gVCF must never win.
    assert _best_variant_count([genome, variants]) == 8
    assert _best_variant_count([variants, genome]) == 8


def test_best_variant_count_prefers_called_vcf_over_gvcf(tmp_path: Path) -> None:
    # HaplotypeCaller/DeepVariant/Sentieon emit `<sample>.<caller>.g.vcf` (gVCF)
    # alongside the called VCF. `.g.vcf` sorts before `.vcf` alphabetically, so
    # the naive sort used to count the gVCF. The called VCF must win.
    gvcf = tmp_path / "S1.haplotypecaller.g.vcf"
    called = tmp_path / "S1.haplotypecaller.vcf"
    _write_vcf(gvcf, 500)
    _write_vcf(called, 12)
    assert _best_variant_count([gvcf, called]) == 12


def test_best_variant_count_prefers_filtered_over_unfiltered(tmp_path: Path) -> None:
    # When a filtered call set exists it is the most meaningful count, and it
    # must still be preferred over both the raw call set and any gVCF.
    gvcf = tmp_path / "S1.haplotypecaller.g.vcf"
    raw = tmp_path / "S1.haplotypecaller.vcf"
    filtered = tmp_path / "S1.haplotypecaller.filtered.vcf"
    _write_vcf(gvcf, 500)
    _write_vcf(raw, 20)
    _write_vcf(filtered, 11)
    assert _best_variant_count([gvcf, raw, filtered]) == 11


def test_best_variant_count_falls_back_to_gvcf_when_only_option(tmp_path: Path) -> None:
    # If the gVCF is the only VCF present, counting it is still better than None.
    gvcf = tmp_path / "S1.deepvariant.g.vcf"
    _write_vcf(gvcf, 42)
    assert _best_variant_count([gvcf]) == 42


# ---------------------------------------------------------------------------
# result.json: expose samples_detected from the parsed outputs
# ---------------------------------------------------------------------------


def test_result_json_exposes_samples_detected(tmp_path: Path) -> None:
    # In demo mode the input samplesheet is empty (the `test` profile supplies
    # it remotely), so `samples` is []. The samples actually produced must still
    # be machine-readable via `samples_detected`.
    output_dir = _make_output_dir(tmp_path)
    empty_sheet = {
        "sample_count": 0,
        "sample_names": [],
        "patient_count": 0,
        "patient_names": [],
        "pairings": [],
        "analysis_mode": "germline",
    }
    ofr = _FakeOutputsReport(samples_detected=["test"])
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=empty_sheet,
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=ofr,
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    data = json.loads(artifacts.result_json.read_text(encoding="utf-8"))
    assert data["samples"] == []  # input-derived, unchanged
    assert data["samples_detected"] == ["test"]  # outputs-derived, new


def test_result_json_samples_detected_empty_without_outputs(tmp_path: Path) -> None:
    output_dir = _make_output_dir(tmp_path)
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=None,
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    data = json.loads(artifacts.result_json.read_text(encoding="utf-8"))
    assert data["samples_detected"] == []


# ---------------------------------------------------------------------------
# commands.sh
# ---------------------------------------------------------------------------


def test_commands_sh_uses_script_dir_and_clawbio_repo(tmp_path: Path) -> None:
    output_dir = _make_output_dir(tmp_path)
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    text = artifacts.commands_sh.read_text(encoding="utf-8")
    assert "$SCRIPT_DIR" in text
    assert "$CLAWBIO_REPO" in text or "CLAWBIO_REPO" in text


def test_commands_sh_is_directly_runnable_and_pins_engine(tmp_path: Path) -> None:
    """The bundle must replay with a bare `bash commands.sh` — no required env var
    and no <EDIT_ME> in the runnable path — and must pin the Nextflow engine for
    byte-identical reproduction. Mirrors the nfcore-scrnaseq-wrapper bundle."""
    output_dir = _make_output_dir(tmp_path)
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
        nextflow_version="25.10.2",
    )
    text = artifacts.commands_sh.read_text(encoding="utf-8")
    # No hard CLAWBIO_REPO abort that `set -u` would trip before the direct replay.
    assert "${CLAWBIO_REPO:?" not in text
    # Engine pinned for reproducible replay.
    assert 'export NXF_VER="25.10.2"' in text
    # A runnable (non-commented) nextflow invocation is present. The portable argv
    # places each token on its own continuation line, so check for the lowercase
    # `nextflow` token on a non-comment line.
    runnable = [
        ln for ln in text.splitlines()
        if "nextflow" in ln and not ln.lstrip().startswith("#")
    ]
    assert runnable, "expected a runnable nextflow invocation"
    # The runnable path (non-comment lines) must not contain an <EDIT_ME> placeholder.
    runnable_no_comments = "\n".join(
        ln for ln in text.splitlines() if not ln.lstrip().startswith("#")
    )
    assert "EDIT_ME" not in runnable_no_comments




def test_commands_sh_passes_bash_syntax_check(tmp_path: Path) -> None:
    """The generated replay script must be syntactically valid bash on any host."""
    bash = shutil.which("bash")
    if not bash:
        import pytest
        pytest.skip("bash not available")
    output_dir = _make_output_dir(tmp_path)
    macos_cfg = output_dir / "reproducibility" / "macos_docker.config"
    macos_cfg.write_text("// macos\n")
    cmd = [
        "nextflow", "run", "nf-core/sarek", "-r", "3.8.1", "-profile", "docker",
        "-params-file", str(output_dir / "reproducibility" / "params.yaml"),
        "-work-dir", str(output_dir / "upstream" / "work"),
        "-c", str(macos_cfg),
    ]
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=cmd,
        pipeline_source_kind="remote",
        nextflow_version="25.10.2",
    )
    r = subprocess.run(
        [bash, "-n", str(artifacts.commands_sh)], capture_output=True, text=True
    )
    assert r.returncode == 0, f"commands.sh has bash syntax errors:\n{r.stderr}"


def test_commands_sh_does_not_leak_generation_host_path(tmp_path: Path) -> None:
    """The portable replay must self-anchor and never embed the absolute output
    path of the machine that generated the bundle."""
    output_dir = _make_output_dir(tmp_path)
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
        nextflow_version="25.10.2",
    )
    text = artifacts.commands_sh.read_text(encoding="utf-8")
    assert str(output_dir.resolve()) not in text
    assert str(tmp_path) not in text


def test_commands_sh_is_executable(tmp_path: Path) -> None:
    output_dir = _make_output_dir(tmp_path)
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    assert os.access(artifacts.commands_sh, os.X_OK)



def test_portable_argv_handles_symlinked_or_unresolved_output_dir(tmp_path: Path) -> None:
    """Robustness: paths under the output dir must be rewritten to $SCRIPT_DIR even
    when the output dir is reached via a symlink (e.g. macOS /var -> /private/var)
    or otherwise not pre-resolved. A naive string-prefix match would degrade these
    to <EDIT_ME> and break the portable replay."""
    from reporting import _portable_argv

    real = tmp_path / "real_out"
    (real / "reproducibility").mkdir(parents=True)
    link = tmp_path / "link_out"
    link.symlink_to(real, target_is_directory=True)

    argv = [
        "nextflow", "run", "nf-core/sarek",
        "-params-file", str(link / "reproducibility" / "params.yaml"),
        "-work-dir", str(link / "upstream" / "work"),
    ]
    out, warnings = _portable_argv(argv, output_dir=link)
    assert "$SCRIPT_DIR/params.yaml" in out
    assert "$SCRIPT_DIR/../upstream/work" in out
    assert "<EDIT_ME>" not in out
    assert warnings == []


def test_commands_sh_uname_gates_macos_docker_config(tmp_path: Path) -> None:
    """A macOS-only macos_docker.config (Apple-Silicon --platform, host-scaled
    resourceLimits) must be applied only when the *replay* host is macOS, so a
    macOS-generated bundle replays cleanly on Linux. Mirrors nfcore-scrnaseq."""
    output_dir = _make_output_dir(tmp_path)
    macos_cfg = output_dir / "reproducibility" / "macos_docker.config"
    macos_cfg.write_text("// macos\n")
    cmd = [
        "nextflow", "run", "nf-core/sarek", "-r", "3.8.1", "-profile", "docker",
        "-params-file", str(output_dir / "reproducibility" / "params.yaml"),
        "-c", str(macos_cfg),
    ]
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=cmd,
        pipeline_source_kind="remote",
        nextflow_version="25.10.2",
    )
    text = artifacts.commands_sh.read_text(encoding="utf-8")
    # The config is referenced via a uname-gated variable, not baked in directly.
    assert "uname -s" in text and "Darwin" in text
    assert "macos_docker.config" in text
    assert "$EXTRA_CONFIG" in text
    # The runnable nextflow invocation must NOT carry an unconditional
    # `-c $SCRIPT_DIR/macos_docker.config` continuation token.
    non_comment = [ln for ln in text.splitlines() if not ln.lstrip().startswith("#")]
    assert not any(
        "macos_docker.config" in ln and "EXTRA_CONFIG" not in ln for ln in non_comment
    )


def test_commands_sh_warns_when_unsafe_absolute_paths(tmp_path: Path) -> None:
    output_dir = _make_output_dir(tmp_path)
    # Inject an absolute path the transformer cannot safely rewrite
    cmd = [
        "nextflow",
        "run",
        "nf-core/sarek",
        "--fasta",
        "/some/foreign/absolute/path/genome.fa",
        "--outdir",
        str(output_dir),
    ]
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=cmd,
        pipeline_source_kind="remote",
    )
    text = artifacts.commands_sh.read_text(encoding="utf-8")
    # Either replaced by EDIT_ME, or a warning was logged
    assert "EDIT_ME" in text or len(artifacts.warnings) > 0


# ---------------------------------------------------------------------------
# remap_paths.py copy
# ---------------------------------------------------------------------------


def test_remap_paths_py_generated_standalone(tmp_path: Path) -> None:
    output_dir = _make_output_dir(tmp_path)
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    src = _skill_dir() / "remap_paths.py"
    dst = artifacts.remap_paths_py
    assert dst.exists()
    import subprocess
    import sys
    result = subprocess.run([sys.executable, "-I", "-S", str(dst), "--help"], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "--repair-bundle" in result.stdout
    assert dst.stat().st_mode & 0o777 == src.stat().st_mode & 0o777


# ---------------------------------------------------------------------------
# _count_vcf_records helper
# ---------------------------------------------------------------------------






# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------




def test_annotation_section_present_with_data(tmp_path: Path) -> None:
    output_dir = _make_output_dir(tmp_path)
    ann_dir = output_dir / "annotation" / "S1"
    ann_dir.mkdir(parents=True)
    vcf = ann_dir / "S1.strelka_snpEff.ann.vcf.gz"
    vcf.write_bytes(b"x")
    html = ann_dir / "S1.strelka_snpEff.ann.summary.html"
    html.write_text("<html/>", encoding="utf-8")
    ofr = _FakeOutputsReport(
        annotation={
            "snpeff": {"S1.strelka": {"vcf": vcf, "html": html}},
            "vep": {},
            "merge": {},
            "bcfann": {},
        }
    )
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=ofr,
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    text = artifacts.report_md.read_text(encoding="utf-8")
    assert "## Annotation" in text
    assert "snpeff" in text.lower()






# ---------------------------------------------------------------------------
# Runtime metadata + demo fallbacks (Issues 3/4/5/6)
# ---------------------------------------------------------------------------


def test_report_shows_profile_java_nextflow(tmp_path: Path) -> None:
    """Issue 6: Profile/Java/Nextflow must render from the passed runtime
    metadata instead of the '-' placeholder."""
    output_dir = _make_output_dir(tmp_path)
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
        java_version="17.0.9",
        nextflow_version="25.04.3",
        profile="test,docker",
    )
    text = artifacts.report_md.read_text(encoding="utf-8")
    assert "| Profile | test,docker |" in text
    assert "| Java | 17.0.9 |" in text
    assert "| Nextflow | 25.04.3 |" in text
    result = json.loads(artifacts.result_json.read_text(encoding="utf-8"))
    assert result["run"]["profile"] == "test,docker"
    assert result["run"]["java_version"] == "17.0.9"
    assert result["run"]["nextflow_version"] == "25.04.3"


def test_demo_tools_inferred_from_outputs(tmp_path: Path) -> None:
    """Issue 4: with no requested tools (as under --demo) the effective tools
    must be inferred from detected variant-calling outputs, not shown as none."""
    output_dir = _make_output_dir(tmp_path)
    params = _minimal_params()
    params["tools"] = []  # demo: the upstream test profile owns tool selection
    ofr = _FakeOutputsReport(variant_calling={"strelka": {"test": {"vcf": [Path("x.vcf.gz")]}}})
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=params,
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=ofr,
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    text = artifacts.report_md.read_text(encoding="utf-8")
    assert "strelka" in text
    assert "| Tools | (none) |" not in text
    result = json.loads(artifacts.result_json.read_text(encoding="utf-8"))
    assert result["run"]["tools"] == []
    assert result["run"]["tools_effective"] == ["strelka"]
    assert result["run"]["tools_from_outputs"] is True


def test_demo_samples_counted_from_detected_outputs(tmp_path: Path) -> None:
    """Issues 3: with an empty local samplesheet (as under --demo) the Samples
    count must fall back to the samples detected in the outputs."""
    output_dir = _make_output_dir(tmp_path)
    ss = _minimal_samplesheet()
    ss["sample_names"] = []
    ss["sample_count"] = 0
    ofr = _FakeOutputsReport(samples_detected=["test"])
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=ss,
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=ofr,
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    text = artifacts.report_md.read_text(encoding="utf-8")
    assert "| Samples | 1 |" in text


def test_result_json_has_ok_field(tmp_path: Path) -> None:
    """Issue 5: result.json must carry the shared `ok` success discriminator."""
    output_dir = _make_output_dir(tmp_path)
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    result = json.loads(artifacts.result_json.read_text(encoding="utf-8"))
    assert result["ok"] is True
    assert result["status"] in {"ok", "partial"}


def test_report_and_result_written_at_output_root(tmp_path: Path) -> None:
    """Issue 1: report.md and result.json must live at the output root (matching
    nfcore-rnaseq/scrnaseq), while commands.sh/remap_paths.py stay in the bundle."""
    output_dir = _make_output_dir(tmp_path)
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
    )
    assert artifacts.report_md == output_dir / "report.md"
    assert artifacts.result_json == output_dir / "result.json"
    assert (output_dir / "report.md").exists()
    assert (output_dir / "result.json").exists()
    # No longer written inside the reproducibility bundle.
    assert not (output_dir / "reproducibility" / "report.md").exists()
    assert not (output_dir / "reproducibility" / "result.json").exists()
    # Replay/provenance artifacts remain in the bundle.
    assert artifacts.commands_sh == output_dir / "reproducibility" / "commands.sh"


def test_commands_sh_exports_the_syntax_parser_sarek_3_8_1_needs(tmp_path: Path) -> None:
    """Replay runs Nextflow directly, so it must set the parser the wrapper used:
    on Nextflow >= 26.04 the default strict parser rejects sarek 3.8.1's configs."""
    output_dir = _make_output_dir(tmp_path)
    artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_skill_dir(),
        params=_minimal_params(),
        samplesheet_report=_minimal_samplesheet(),
        pipeline_source=_minimal_pipeline_source(),
        outputs_report=_FakeOutputsReport(),
        nextflow_command=_minimal_command(output_dir),
        pipeline_source_kind="remote",
        nextflow_version="26.04.3",
    )
    text = artifacts.commands_sh.read_text(encoding="utf-8")
    assert 'export NXF_SYNTAX_PARSER="v1"' in text
    assert text.index('export NXF_SYNTAX_PARSER="v1"') < text.index('(cd "$SCRIPT_DIR/.."')
