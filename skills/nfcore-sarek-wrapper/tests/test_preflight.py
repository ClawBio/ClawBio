"""Tests for nfcore-sarek-wrapper / preflight.py."""
from __future__ import annotations

import subprocess

import pytest

from errors import ErrorCode, SkillError
from preflight import (
    PreflightResult,
    REFERENCE_PATH_PARAMS,
    _check_annotation_cache,
    _check_email,
    _check_enums,
    _check_flag_compatibility,
    _check_genome_known,
    _check_output_dir,
    _check_pipeline_source,
    _check_profile_string,
    _check_reference_paths,
    _check_resume_drift,
    _check_step,
    _check_tools_against_pairing,
    run_preflight,
)

from preflight_fixtures import make_params, make_pipeline_source, make_samplesheet


# ---------------------------------------------------------------------------
# Stubs for external binaries (Java / Nextflow / Docker)
# ---------------------------------------------------------------------------


@pytest.fixture
def stub_binaries(monkeypatch):
    """Stub shutil.which + subprocess.run so external binaries appear present."""
    import preflight as p

    def fake_which(name: str):
        return f"/usr/bin/{name}"

    def fake_run(args, *a, **kw):
        binary = args[0] if isinstance(args, (list, tuple)) else args
        out = ""
        rc = 0
        if "java" in binary:
            out = 'openjdk version "17.0.10" 2024-01-16'
        elif "nextflow" in binary:
            out = "25.10.2 build 5912"
        elif "docker" in binary:
            out = "Server Version: 24.0"
        elif "podman" in binary:
            out = "Podman version 4.0"
        else:
            out = ""
        return subprocess.CompletedProcess(
            args=args, returncode=rc, stdout=out, stderr=""
        )

    monkeypatch.setattr(p.shutil, "which", fake_which)
    monkeypatch.setattr(p.subprocess, "run", fake_run)
    yield


@pytest.fixture
def stub_no_java(monkeypatch):
    import preflight as p

    monkeypatch.setattr(p.shutil, "which", lambda name: None if name == "java" else f"/usr/bin/{name}")


@pytest.fixture
def stub_old_java(monkeypatch):
    import preflight as p

    def fake_which(name):
        return f"/usr/bin/{name}"

    def fake_run(args, *a, **kw):
        binary = args[0]
        if "java" in binary:
            return subprocess.CompletedProcess(args, 0, "openjdk version \"11.0.20\"", "")
        if "nextflow" in binary:
            return subprocess.CompletedProcess(args, 0, "25.10.2", "")
        return subprocess.CompletedProcess(args, 0, "Server Version: 24.0", "")

    monkeypatch.setattr(p.shutil, "which", fake_which)
    monkeypatch.setattr(p.subprocess, "run", fake_run)


@pytest.fixture
def stub_old_nextflow(monkeypatch):
    import preflight as p

    def fake_which(name):
        return f"/usr/bin/{name}"

    def fake_run(args, *a, **kw):
        binary = args[0]
        if "java" in binary:
            return subprocess.CompletedProcess(args, 0, "openjdk version \"17.0.10\"", "")
        if "nextflow" in binary:
            return subprocess.CompletedProcess(args, 0, "24.10.0", "")
        return subprocess.CompletedProcess(args, 0, "Server Version: 24.0", "")

    monkeypatch.setattr(p.shutil, "which", fake_which)
    monkeypatch.setattr(p.subprocess, "run", fake_run)


@pytest.fixture
def stub_no_docker(monkeypatch):
    import preflight as p

    def fake_which(name):
        return None if name == "docker" else f"/usr/bin/{name}"

    def fake_run(args, *a, **kw):
        binary = args[0]
        if "java" in binary:
            return subprocess.CompletedProcess(args, 0, "openjdk version \"17.0.10\"", "")
        if "nextflow" in binary:
            return subprocess.CompletedProcess(args, 0, "25.10.2", "")
        return subprocess.CompletedProcess(args, 0, "", "")

    monkeypatch.setattr(p.shutil, "which", fake_which)
    monkeypatch.setattr(p.subprocess, "run", fake_run)


@pytest.fixture
def stub_docker_not_running(monkeypatch):
    import preflight as p

    def fake_which(name):
        return f"/usr/bin/{name}"

    def fake_run(args, *a, **kw):
        binary = args[0]
        if "java" in binary:
            return subprocess.CompletedProcess(args, 0, "openjdk version \"17.0.10\"", "")
        if "nextflow" in binary:
            return subprocess.CompletedProcess(args, 0, "25.10.2", "")
        if "docker" in binary:
            return subprocess.CompletedProcess(args, 1, "", "Cannot connect to the Docker daemon")
        return subprocess.CompletedProcess(args, 0, "", "")

    monkeypatch.setattr(p.shutil, "which", fake_which)
    monkeypatch.setattr(p.subprocess, "run", fake_run)


# ---------------------------------------------------------------------------
# Utility tests
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# §5.1 Java / Nextflow / profile / backend / output dir / email / pipeline
# ---------------------------------------------------------------------------


def test_java_missing(stub_no_java):
    from preflight import _check_java

    with pytest.raises(SkillError) as exc:
        _check_java()
    # A *missing* binary must be MISSING_JAVA, not JAVA_VERSION_TOO_OLD (the code must not
    # contradict the message). Aligned with nfcore-scrnaseq-wrapper.
    assert exc.value.error_code == ErrorCode.MISSING_JAVA


def test_java_version_unparseable_is_version_failure(monkeypatch):
    import preflight as p

    monkeypatch.setattr(p.shutil, "which", lambda n: "/usr/bin/java")
    monkeypatch.setattr(p, "_command_output", lambda args: "no version here")
    with pytest.raises(SkillError) as exc:
        p._check_java()
    # java IS present but its version can't be confirmed → version failure, not MISSING.
    assert exc.value.error_code == ErrorCode.JAVA_VERSION_TOO_OLD


def test_nextflow_missing(monkeypatch):
    import preflight as p

    monkeypatch.setattr(p.shutil, "which", lambda n: None if n == "nextflow" else "/usr/bin/x")
    with pytest.raises(SkillError) as exc:
        p._check_nextflow()
    assert exc.value.error_code == ErrorCode.MISSING_NEXTFLOW


def test_nextflow_version_unparseable_is_version_failure(monkeypatch):
    import preflight as p

    monkeypatch.setattr(p.shutil, "which", lambda n: "/usr/bin/nextflow")
    monkeypatch.setattr(p._runtime, "run_probe", lambda command, timeout_seconds: __import__("subprocess").CompletedProcess(command, 0, "no version here", ""))
    with pytest.raises(SkillError) as exc:
        p._check_nextflow()
    # nextflow IS present but its version can't be confirmed → version failure,
    # not NEXTFLOW_NOT_FOUND (which would wrongly imply the binary is absent).
    assert exc.value.error_code == ErrorCode.NEXTFLOW_VERSION_TOO_OLD


def test_nextflow_version_preserves_zero_padded_month(monkeypatch):
    """A detected version like '26.04.3' must be reported verbatim, never
    reconstructed from the comparison tuple as '26.4.3' (which is not a real
    Nextflow release and would break NXF_VER / conda pins on replay)."""
    import preflight as p

    monkeypatch.setattr(p.shutil, "which", lambda n: f"/usr/bin/{n}")
    monkeypatch.setattr(p._runtime, "run_probe", lambda command, timeout_seconds: __import__("subprocess").CompletedProcess(command, 0, "  N E X T F L O W\n  version 26.04.3 build 5928", ""))
    assert p._check_nextflow() == "26.04.3"


def test_java_version_preserves_exact_string(monkeypatch):
    import preflight as p

    monkeypatch.setattr(p.shutil, "which", lambda n: f"/usr/bin/{n}")
    monkeypatch.setattr(
        p, "_command_output", lambda args: 'openjdk version "17.0.10" 2024-01-16'
    )
    assert p._check_java() == "17.0.10"


def test_profile_empty():
    with pytest.raises(SkillError) as exc:
        _check_profile_string("")
    assert exc.value.error_code == ErrorCode.INVALID_PROFILE


def test_profile_custom_token_is_allowed_for_institutional_configs():
    _check_profile_string("docker,my_institutional_hpc")


def test_backend_docker_missing(stub_no_docker, monkeypatch):
    import preflight as p

    with pytest.raises(SkillError) as exc:
        p._check_backends_for_profile("docker")
    # Specific backend codes (aligned with nfcore-scrnaseq-wrapper).
    assert exc.value.error_code == ErrorCode.MISSING_DOCKER


def test_backend_specific_codes_per_runtime(monkeypatch):
    import preflight as p

    # Binary absent → MISSING_<runtime>; daemon present but down → <runtime>_NOT_RUNNING.
    monkeypatch.setattr(p.shutil, "which", lambda name: None)
    cases = {
        "podman": ErrorCode.MISSING_PODMAN,
        "singularity": ErrorCode.MISSING_SINGULARITY,
        "apptainer": ErrorCode.MISSING_SINGULARITY,
        "conda": ErrorCode.MISSING_CONDA,
        "mamba": ErrorCode.MISSING_CONDA,
        "shifter": ErrorCode.MISSING_HPC_RUNTIME,
        "charliecloud": ErrorCode.MISSING_HPC_RUNTIME,
    }
    for token, code in cases.items():
        with pytest.raises(SkillError) as exc:
            p._check_backends_for_profile(token)
        assert exc.value.error_code == code, token


def test_sentieon_tool_requires_license_secret(monkeypatch):
    import preflight as p

    monkeypatch.delenv("SENTIEON_LICENSE_BASE64", raising=False)
    monkeypatch.delenv("SENTIEON_LICENSE", raising=False)
    monkeypatch.setattr(p, "_nextflow_secret_exists", lambda name: False)

    with pytest.raises(SkillError) as exc:
        _check_tools_against_pairing(
            params={"step": "mapping", "tools": ["sentieon_haplotyper"]},
            tools=["sentieon_haplotyper"],
            present_modes={"germline"},
        )

    assert exc.value.error_code == ErrorCode.MISSING_SENTIEON_LICENSE


def test_sentieon_aligner_requires_license_secret(monkeypatch):
    import preflight as p

    monkeypatch.delenv("SENTIEON_LICENSE_BASE64", raising=False)
    monkeypatch.delenv("SENTIEON_LICENSE", raising=False)
    monkeypatch.setattr(p, "_nextflow_secret_exists", lambda name: False)

    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params={"step": "mapping", "aligner": "sentieon-bwamem"}, tools=[])

    assert exc.value.error_code == ErrorCode.MISSING_SENTIEON_LICENSE


def test_backend_docker_not_running(stub_docker_not_running):
    import preflight as p

    with pytest.raises(SkillError) as exc:
        p._check_backends_for_profile("docker")
    assert exc.value.error_code == ErrorCode.DOCKER_NOT_RUNNING


def test_output_dir_inside_repo(tmp_path):
    repo_root = tmp_path / "repo"
    repo_root.mkdir()
    out = repo_root / "results"
    with pytest.raises(SkillError) as exc:
        _check_output_dir(out, repo_root=repo_root, resume=False)
    assert exc.value.error_code == ErrorCode.OUTPUT_DIR_INSIDE_REPO


def test_output_dir_nonempty_no_resume(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    out = tmp_path / "out"
    out.mkdir()
    (out / "leftover.txt").write_text("x")
    with pytest.raises(SkillError) as exc:
        _check_output_dir(out, repo_root=repo, resume=False)
    assert exc.value.error_code == ErrorCode.OUTPUT_DIR_NOT_EMPTY


@pytest.mark.parametrize(
    "junk_name",
    [".DS_Store", ".gitkeep", ".gitignore", "Thumbs.db", "check_result.json"],
)
def test_output_dir_tolerates_ignored_root_files(tmp_path, junk_name):
    """Stray OS/VCS scratch files and a prior --check summary at the output root
    must not trip OUTPUT_DIR_NOT_EMPTY (parity with nfcore-rnaseq/scrnaseq). On
    macOS, Finder writes .DS_Store into any visited folder, so rejecting it would
    make sarek fail on directories the sibling wrappers accept."""
    repo = tmp_path / "repo"
    repo.mkdir()
    out = tmp_path / "out"
    out.mkdir()
    (out / junk_name).write_text("x")
    # Must not raise.
    _check_output_dir(out, repo_root=repo, resume=False)


def test_output_dir_check_then_run_same_dir(tmp_path):
    """A --check run writes check_result.json at the output root; a subsequent
    real run in the same directory must still pass the empty-dir gate."""
    repo = tmp_path / "repo"
    repo.mkdir()
    out = tmp_path / "out"
    out.mkdir()
    (out / "check_result.json").write_text("{}")
    (out / "reproducibility").mkdir()
    # reproducibility/ + check_result.json only → tolerated.
    _check_output_dir(out, repo_root=repo, resume=False)


def test_email_invalid():
    with pytest.raises(SkillError) as exc:
        _check_email("not_an_email")
    assert exc.value.error_code == ErrorCode.INVALID_EMAIL


def test_pipeline_source_remote_no_version():
    with pytest.raises(SkillError) as exc:
        _check_pipeline_source({"source_kind": "remote_repo", "resolved_version": ""})
    assert exc.value.error_code == ErrorCode.PIPELINE_SOURCE_INVALID


# ---------------------------------------------------------------------------
# §5.2 step enum
# ---------------------------------------------------------------------------


def test_step_unknown():
    with pytest.raises(SkillError) as exc:
        _check_step({"step": "no_such_step"})
    assert exc.value.error_code == ErrorCode.INVALID_STEP


def test_variant_calling_requires_tools_in_preflight(tmp_path, stub_binaries):
    with pytest.raises(SkillError) as exc:
        run_preflight(
            params=make_params(step="variant_calling", tools=[]),
            samplesheet=make_samplesheet(),
            pipeline_source=make_pipeline_source(),
            output_dir=tmp_path / "out",
            repo_root=tmp_path / "repo",
        )
    assert exc.value.error_code == ErrorCode.INVALID_TOOLS


def test_annotate_requires_tools_in_preflight(tmp_path, stub_binaries):
    with pytest.raises(SkillError) as exc:
        run_preflight(
            params=make_params(step="annotate", tools=[]),
            samplesheet=make_samplesheet(),
            pipeline_source=make_pipeline_source(),
            output_dir=tmp_path / "out",
            repo_root=tmp_path / "repo",
        )
    assert exc.value.error_code == ErrorCode.INVALID_TOOLS


def test_input_error_surfaces_before_backend_probe(tmp_path, stub_no_docker):
    # Cheap, environment-independent input validation must run before the backend
    # probe so a typo'd --step reports INVALID_STEP instead of being masked by a
    # missing/stopped container backend.
    with pytest.raises(SkillError) as exc:
        run_preflight(
            params=make_params(step="no_such_step"),
            samplesheet=make_samplesheet(),
            pipeline_source=make_pipeline_source(),
            output_dir=tmp_path / "out",
            repo_root=tmp_path / "repo",
        )
    assert exc.value.error_code == ErrorCode.INVALID_STEP


def test_unknown_tool_surfaces_before_backend_probe(tmp_path, stub_no_docker):
    with pytest.raises(SkillError) as exc:
        run_preflight(
            params=make_params(step="variant_calling", tools=["foobar"]),
            samplesheet=make_samplesheet(),
            pipeline_source=make_pipeline_source(),
            output_dir=tmp_path / "out",
            repo_root=tmp_path / "repo",
        )
    assert exc.value.error_code == ErrorCode.INVALID_TOOLS


def test_backend_still_checked_when_input_valid(tmp_path, stub_no_docker):
    # The reorder must NOT skip the backend probe: valid input + no backend still
    # raises the (now runtime-specific) docker-missing code.
    with pytest.raises(SkillError) as exc:
        run_preflight(
            params=make_params(step="variant_calling", tools=["strelka"]),
            samplesheet=make_samplesheet(),
            pipeline_source=make_pipeline_source(),
            output_dir=tmp_path / "out",
            repo_root=tmp_path / "repo",
        )
    assert exc.value.error_code == ErrorCode.MISSING_DOCKER


# ---------------------------------------------------------------------------
# §5.3 tools × pairing matrix
# ---------------------------------------------------------------------------


def _tools_check(
    *,
    tools: list[str],
    analysis_mode: str | None = "germline",
    present_modes: set[str] | None = None,
    params: dict | None = None,
):
    # Tests may pass a single analysis_mode (legacy) or an explicit set of
    # per-patient present_modes (for mixed-samplesheet scenarios).
    if present_modes is None:
        present_modes = {analysis_mode} if analysis_mode else {"germline"}
    return _check_tools_against_pairing(
        params=params or {},
        tools=tools,
        present_modes=present_modes,
    )


def test_tools_unknown():
    with pytest.raises(SkillError) as exc:
        _tools_check(tools=["nonexistenttool"])
    assert exc.value.error_code == ErrorCode.INVALID_TOOLS


# --- #22 merge applies snpeff + vep rules; #24 vep cache exceptions ---------

def test_merge_requires_snpeff_db_when_cache_set():
    # merge applies the snpEff rule set: an explicit --snpeff_cache without
    # --snpeff_db is the (only) hard error, matching sarek.
    with pytest.raises(SkillError) as exc:
        _tools_check(
            tools=["merge"],
            params={
                "step": "annotate", "snpeff_cache": "/cache",
                "vep_cache_version": "110", "vep_genome": "GRCh38", "vep_species": "homo_sapiens",
                "genome": "null", "igenomes_ignore": True,
            },
        )
    assert exc.value.error_code == ErrorCode.MISSING_SNPEFF_DB


def test_vep_missing_cache_errors_without_igenomes():
    with pytest.raises(SkillError) as exc:
        _tools_check(tools=["vep"], params={"step": "annotate", "igenomes_ignore": True})
    assert exc.value.error_code == ErrorCode.MISSING_VEP_CACHE


def test_strelka_tumor_only_rejected():
    # Official tool matrix: Strelka supports germline + somatic, NOT tumor-only.
    with pytest.raises(SkillError) as exc:
        _tools_check(tools=["strelka"], analysis_mode="tumor_only")
    assert exc.value.error_code == ErrorCode.INVALID_PAIRING


def test_mpileup_somatic_paired_runs_on_normal_member():
    assert _tools_check(tools=["mpileup"], analysis_mode="somatic_paired") == []


def test_indexcov_tumor_only_rejected():
    with pytest.raises(SkillError) as exc:
        _tools_check(tools=["indexcov"], analysis_mode="tumor_only")
    assert exc.value.error_code == ErrorCode.INVALID_PAIRING


def test_mutect2_tumor_only_no_pon_warns():
    warnings = _tools_check(
        tools=["mutect2"],
        analysis_mode="tumor_only",
        params={"pon": None},
    )
    assert any("Panel-of-Normals" in w for w in warnings)


def test_mutect2_default_pon_warns():
    warnings = _tools_check(
        tools=["mutect2"],
        analysis_mode="tumor_only",
        params={"pon": "/refs/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz",
                "germline_resource": "/g.vcf.gz"},
    )
    assert any("default GATK Panel-of-Normals" in w for w in warnings)


def test_ascat_requires_paired():
    with pytest.raises(SkillError) as exc:
        _tools_check(tools=["ascat"], analysis_mode="tumor_only")
    assert exc.value.error_code == ErrorCode.INVALID_PAIRING


def test_ascat_wes_missing_resources():
    with pytest.raises(SkillError) as exc:
        _tools_check(
            tools=["ascat"],
            analysis_mode="somatic_paired",
            params={"wes": True, "genome": "null", "igenomes_ignore": True},
        )
    assert exc.value.error_code == ErrorCode.INVALID_ASCAT_RESOURCES


def test_msisensor2_paired_rejected():
    with pytest.raises(SkillError) as exc:
        _tools_check(tools=["msisensor2"], analysis_mode="somatic_paired")
    assert exc.value.error_code == ErrorCode.INVALID_PAIRING


def test_msisensor2_tumor_only_missing_models():
    with pytest.raises(SkillError) as exc:
        _tools_check(
            tools=["msisensor2"],
            analysis_mode="tumor_only",
            params={"genome": "null", "igenomes_ignore": True},
        )
    assert exc.value.error_code == ErrorCode.MISSING_REFERENCE


def test_msisensorpro_tumor_only_rejected():
    with pytest.raises(SkillError) as exc:
        _tools_check(tools=["msisensorpro"], analysis_mode="tumor_only")
    assert exc.value.error_code == ErrorCode.INVALID_PAIRING


def test_msisensorpro_paired_missing_scan():
    # Sarek auto-generates the scan from the reference FASTA in PREPARE_GENOME,
    # so missing --msisensorpro_scan is a warning, not a hard error.
    warnings = _tools_check(tools=["msisensorpro"], analysis_mode="somatic_paired")
    assert any("msisensorpro_scan" in w for w in warnings)


def test_snpeff_without_db_or_cache_raises_when_no_igenomes():
    # Without --snpeff_db AND without an iGenomes genome that supplies SnpEff
    # config, Sarek will fail → the wrapper raises a hard error.
    with pytest.raises(SkillError) as exc:
        _tools_check(
            tools=["snpeff"],
            analysis_mode="germline",
            params={"step": "annotate", "genome": "null", "igenomes_ignore": True},
        )
    assert exc.value.error_code == ErrorCode.MISSING_SNPEFF_DB


def test_snpeff_with_cache_still_requires_db():
    with pytest.raises(SkillError) as exc:
        _tools_check(
            tools=["snpeff"],
            analysis_mode="germline",
            params={
                "step": "annotate", "snpeff_cache": "/cache",
                "genome": "null", "igenomes_ignore": True,
            },
        )
    assert exc.value.error_code == ErrorCode.MISSING_SNPEFF_DB


# nf-core/sarek 3.8.1 refuses these itself in subworkflows/local/samplesheet_to_channel
# when the run starts: --tools ascat/controlfreec/mutect2/msisensorpro with no tumour
# sample, ascat/deepvariant/haplotypecaller/msisensorpro with no normal sample, ascat/
# controlfreec/varlociraptor with a sample of sex NA, and varlociraptor with a tumour
# lacking contamination. Preflight keeps no copy of these rules.
_SAREK_OWN_TOOL_RULES = {
    "mutect2 without a tumour": dict(tools=["mutect2"], present_modes={"germline"}),
    "controlfreec without a tumour": dict(
        tools=["controlfreec"], present_modes={"germline"}, params={"mappability": "/m", "chr_dir": "/c"}),
    "haplotypecaller without a normal": dict(tools=["haplotypecaller"], present_modes={"tumor_only"}),
    "deepvariant without a normal": dict(tools=["deepvariant"], present_modes={"tumor_only"}),
}


@pytest.mark.parametrize("case", sorted(_SAREK_OWN_TOOL_RULES))
def test_tool_rules_sarek_enforces_reach_it(case):
    _tools_check(**_SAREK_OWN_TOOL_RULES[case])


@pytest.mark.parametrize("tools,mode,rows", [
    (["controlfreec"], "somatic_paired", [{"sample": "N1", "status": 0, "sex": "NA", "contamination": None, "line": 2},
                                          {"sample": "T1", "status": 1, "sex": "NA", "contamination": None, "line": 3}]),
    (["ascat"], "somatic_paired", [{"sample": "N1", "status": 0, "sex": "NA", "contamination": None, "line": 2},
                                   {"sample": "T1", "status": 1, "sex": "XX", "contamination": 0.1, "line": 3}]),
    (["mutect2", "varlociraptor"], "tumor_only", [{"sample": "T1", "status": 1, "sex": "NA", "contamination": None, "line": 2}]),
])
def test_sex_and_contamination_rules_are_left_to_sarek(tmp_path, stub_binaries, tools, mode, rows):
    """samplesheet_to_channel refuses sex NA for ascat/controlfreec/varlociraptor and a
    varlociraptor tumour without contamination; the full preflight passes them on."""
    (tmp_path / "repo").mkdir()
    pairings = [{"patient": "P1", "normal": "N1", "tumor": "T1", "mode": mode}]
    run_preflight(
        params=make_params(tools=tools),
        samplesheet=make_samplesheet(analysis_mode=mode, rows_by_patient={"P1": rows}, pairings=pairings),
        pipeline_source=make_pipeline_source(),
        output_dir=tmp_path / "results",
        repo_root=tmp_path / "repo",
    )


@pytest.mark.parametrize("tool,present_modes", [
    ("sentieon_haplotyper", {"tumor_only"}),
    ("sentieon_dnascope", {"tumor_only"}),
    ("muse", {"tumor_only"}),
    ("lofreq", {"somatic_paired"}),
    ("strelka", {"tumor_only"}),
])
def test_tool_rules_sarek_does_not_check_are_still_refused(tool, present_modes, monkeypatch):
    """Sarek's own tool checks cover only the tools named above; for these it would
    run and produce nothing for the samples given, so preflight still refuses."""
    import preflight as p
    monkeypatch.setattr(p, "_require_sentieon_license", lambda context: None)
    with pytest.raises(SkillError) as exc:
        _tools_check(tools=[tool], present_modes=present_modes)
    assert exc.value.error_code == ErrorCode.INVALID_PAIRING


def test_bcfann_requires_annotations_tbi_and_header_lines():
    with pytest.raises(SkillError) as exc:
        _tools_check(tools=["bcfann"], analysis_mode="germline", params={"step": "annotate"})
    assert exc.value.error_code == ErrorCode.MISSING_REFERENCE


def test_vep_missing_cache_params():
    with pytest.raises(SkillError) as exc:
        _tools_check(
            tools=["vep"],
            analysis_mode="germline",
            params={"genome": "null", "igenomes_ignore": True},
        )
    assert exc.value.error_code == ErrorCode.MISSING_VEP_CACHE


def test_sentieon_tnscope_germline_rejected(monkeypatch):
    # sentieon_tnscope is somatic-only; in a germline run it must be rejected as an
    # invalid pairing BEFORE a Sentieon license is required. Force the no-license
    # condition (as on CI) so this is deterministic regardless of the host
    # environment: the pairing error must win over MISSING_SENTIEON_LICENSE.
    import preflight as p
    monkeypatch.delenv("SENTIEON_LICENSE_BASE64", raising=False)
    monkeypatch.delenv("SENTIEON_LICENSE", raising=False)
    monkeypatch.setattr(p, "_nextflow_secret_exists", lambda name: False)
    with pytest.raises(SkillError) as exc:
        _tools_check(tools=["sentieon_tnscope"], analysis_mode="germline")
    assert exc.value.error_code == ErrorCode.INVALID_PAIRING


# ---------------------------------------------------------------------------
# §5.4 flag compatibility
# ---------------------------------------------------------------------------


def test_spark_with_save_output_as_bam():
    params = make_params(use_gatk_spark="markduplicates", save_output_as_bam=True, save_mapped=True)
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=[])
    assert exc.value.error_code == ErrorCode.SPARK_OUTPUT_INCOMPATIBLE


def test_markduplicates_spark_with_header_umi_rejected():
    # Sarek's exact rule: MarkDuplicatesSpark can't do UMI-based dedup with
    # umi_in_read_header or umi_location.
    for umi in ({"umi_in_read_header": True}, {"umi_location": "per_read", "umi_length": 8}):
        params = make_params(use_gatk_spark="markduplicates", **umi)
        with pytest.raises(SkillError) as exc:
            _check_flag_compatibility(params=params, tools=[])
        assert exc.value.error_code == ErrorCode.UMI_SPARK_INCOMPATIBLE


def test_joint_germline_without_caller():
    params = make_params(joint_germline=True, tools=["manta"])
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=["manta"])
    assert exc.value.error_code == ErrorCode.JOINT_GERMLINE_REQUIREMENT


def test_joint_mutect2_without_mutect2():
    params = make_params(joint_mutect2=True, tools=["strelka"])
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=["strelka"])
    assert exc.value.error_code == ErrorCode.JOINT_MUTECT2_REQUIREMENT


def test_genome_missing_without_igenomes_ignore():
    params = make_params(genome="", igenomes_ignore=False)
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=[])
    assert exc.value.error_code == ErrorCode.GENOME_REQUIRED


def test_genome_and_fasta_without_igenomes_ignore_allows_documented_override():
    # Sarek documents that an individual iGenomes reference may be overwritten
    # while the remaining resources remain sourced through --genome.
    params = make_params(genome="GATK.GRCh38", fasta="/refs/genome.fa", igenomes_ignore=False)
    _check_flag_compatibility(params=params, tools=[])


def test_build_only_index_with_input_real():
    params = make_params(build_only_index=True, input="samplesheet.csv")
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=[])
    assert exc.value.error_code == ErrorCode.BUILD_ONLY_INDEX_CONFLICT


def test_parabricks_with_conda_profile():
    params = make_params(aligner="parabricks", profile="conda")
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=[])
    assert exc.value.error_code == ErrorCode.PARABRICKS_CONDA_CONFLICT


def test_vep_dbnsfp_without_dbnsfp():
    params = make_params(vep_dbnsfp=True)
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=[])
    assert exc.value.error_code == ErrorCode.MISSING_DBNSFP


@pytest.mark.parametrize("flag,paths", [
    ("vep_condel", ("condel_config",)),
    ("vep_mastermind", ("mastermind_file",)),
    ("vep_spliceai", ("spliceai_snv", "spliceai_snv_tbi", "spliceai_indel", "spliceai_indel_tbi")),
])
def test_vep_plugins_require_documented_paths(flag, paths):
    params = make_params(**{flag: True})
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=["vep"])
    assert exc.value.error_code == ErrorCode.MISSING_REFERENCE
    for path_param in paths:
        assert path_param in exc.value.details["missing"]


def test_filter_vcfs_rejects_blank_expression():
    params = make_params(filter_vcfs=True, bcftools_filter_criteria="   ")
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=[])
    assert exc.value.error_code == ErrorCode.INVALID_BCFTOOLS_FILTER


def test_snv_consensus_too_few_callers():
    params = make_params(snv_consensus_calling=True, tools=["mutect2"])
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=["mutect2"])
    assert exc.value.error_code == ErrorCode.INVALID_CONSENSUS_CALLING


def test_snv_consensus_min_count_too_high():
    params = make_params(
        snv_consensus_calling=True,
        tools=["mutect2", "strelka"],
        consensus_min_count=5,
    )
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=["mutect2", "strelka"])
    assert exc.value.error_code == ErrorCode.INVALID_CONSENSUS_CALLING


def test_snpsift_missing_databases():
    params = make_params(tools=["snpsift"])
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=["snpsift"])
    assert exc.value.error_code == ErrorCode.MISSING_SNPSIFT_DATABASES


def test_no_intervals_with_intervals():
    params = make_params(no_intervals=True, intervals="/path/intervals.bed", genome="null", igenomes_ignore=True, fasta="/refs/genome.fa")
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=[])
    assert exc.value.error_code == ErrorCode.INVALID_INTERVALS


def test_wes_intervals_must_be_bed():
    params = make_params(wes=True, intervals="/refs/targets.interval_list", genome="null", igenomes_ignore=True, fasta="/refs/genome.fa")
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=[])
    assert exc.value.error_code == ErrorCode.INVALID_INTERVALS


def test_umi_location_requires_umi_length():
    params = make_params(umi_location="read1", umi_length=None)
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=[])
    assert exc.value.error_code == ErrorCode.INVALID_FLAG_COMBINATION


def test_umi_read_structure_conflicts_with_umi_location():
    params = make_params(umi_read_structure="+T +T", umi_location="read1", umi_length=8)
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=[])
    assert exc.value.error_code == ErrorCode.INVALID_FLAG_COMBINATION


def test_umi_in_read_header_requires_exact_read_structure_when_combined():
    params = make_params(umi_in_read_header=True, umi_read_structure="8M+T 8M+T")
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=[])
    assert exc.value.error_code == ErrorCode.INVALID_FLAG_COMBINATION


def test_parabricks_umi_flags_rejected():
    params = make_params(aligner="parabricks", umi_in_read_header=True)
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=[])
    assert exc.value.error_code == ErrorCode.INVALID_FLAG_COMBINATION


def test_sentieon_joint_germline_requires_gvcf_emit_mode():
    params = make_params(
        joint_germline=True,
        sentieon_haplotyper_emit_mode="variant",
    )
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=["sentieon_haplotyper"])
    assert exc.value.error_code == ErrorCode.JOINT_GERMLINE_REQUIREMENT


def test_bqsr_requires_dbsnp_or_known_indels_when_not_skipped():
    params = make_params(step="prepare_recalibration", skip_tools=[], genome="", igenomes_ignore=True, fasta="/ref.fa")
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=[])
    assert exc.value.error_code == ErrorCode.MISSING_REFERENCE


def test_bqsr_error_explains_baserecalibrator_runs_by_default():
    """A user who only requested --tools haplotypecaller may not know
    baserecalibrator runs by default in the mapping workflow. The error must
    say so, so the fix (--skip-tools baserecalibrator or supplying known sites)
    is understandable rather than surprising."""
    params = make_params(step="mapping", skip_tools=[], genome="", igenomes_ignore=True, fasta="/ref.fa")
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=params, tools=["haplotypecaller"])
    assert exc.value.error_code == ErrorCode.MISSING_REFERENCE
    assert "baserecalibrator runs by default" in exc.value.message.lower()


def test_bqsr_allows_upstream_test_profile_refs():
    params = make_params(step="mapping", skip_tools=[], genome="", fasta="", profile="test,docker")
    _check_flag_compatibility(params=params, tools=[])


# ---------------------------------------------------------------------------
# §5.5 reference path existence
# ---------------------------------------------------------------------------


def test_reference_path_uri_skipped(tmp_path):
    _check_reference_paths({"fasta": "s3://bucket/ref.fa"})


def test_reference_param_list_contains_expected():
    for key in ("fasta", "intervals", "vep_cache", "snpeff_cache", "pon", "ascat_alleles"):
        assert key in REFERENCE_PATH_PARAMS


# ---------------------------------------------------------------------------
# §5.6 annotation cache
# ---------------------------------------------------------------------------


def test_download_cache_unwritable(tmp_path, monkeypatch):
    import os as real_os
    out = tmp_path / "out"
    out.mkdir()
    # Force os.access to return False
    monkeypatch.setattr(real_os, "access", lambda p, m: False)
    import preflight as p

    monkeypatch.setattr(p.os, "access", lambda path, mode: False)
    with pytest.raises(SkillError) as exc:
        _check_annotation_cache(
            params={"download_cache": True, "outdir_cache": str(out)},
            output_dir=out,
        )
    assert exc.value.error_code == ErrorCode.MISSING_ANNOTATION_CACHE


def test_snpeff_cache_layout_warns(tmp_path):
    cache_root = tmp_path / "snpeff"
    cache_root.mkdir()
    warnings = _check_annotation_cache(
        params={"snpeff_cache": str(cache_root), "snpeff_db": "GRCh38.105"},
        output_dir=tmp_path / "out",
    )
    assert any("snpeff_db" in w.lower() or "snpeff" in w.lower() for w in warnings)


def test_vep_cache_layout_warns(tmp_path):
    cache_root = tmp_path / "vep"
    cache_root.mkdir()
    warnings = _check_annotation_cache(
        params={
            "vep_cache": str(cache_root),
            "vep_species": "homo_sapiens",
            "vep_cache_version": "111",
            "vep_genome": "GRCh38",
        },
        output_dir=tmp_path / "out",
    )
    assert any("vep" in w.lower() for w in warnings)


def test_umi_read_structure_invalid_rejected():
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params={"umi_read_structure": "8X nonsense"}, tools=[])
    assert exc.value.error_code == ErrorCode.INVALID_FLAG_COMBINATION


def test_dragmap_without_skip_baserecalibrator_warns():
    warnings = _check_flag_compatibility(
        params=make_params(aligner="dragmap", step="mapping"), tools=[]
    )
    assert any("DragMap" in w for w in warnings)


def test_invalid_enum_values_rejected():
    for params in (
        {"umi_location": "nope"},
        {"group_by_umi_strategy": "bad"},
        {"vep_out_format": "xml"},
        {"use_gatk_spark": "banana"},
    ):
        with pytest.raises(SkillError) as exc:
            _check_enums(params)
        assert exc.value.error_code == ErrorCode.INVALID_FLAG_COMBINATION


def test_invalid_sentieon_emit_mode_rejected():
    # Guards the dead-constant consolidation: SUPPORTED_SENTIEON_EMIT_MODE is the
    # single source of truth used by _check_enums for both emit-mode params.
    for name in ("sentieon_haplotyper_emit_mode", "sentieon_dnascope_emit_mode"):
        with pytest.raises(SkillError) as exc:
            _check_enums({name: "banana"})
        assert exc.value.error_code == ErrorCode.INVALID_FLAG_COMBINATION
        # All documented combinations remain accepted.
        for good in ("variant", "gvcf", "gvcf,variant", "variant,gvcf", "all"):
            _check_enums({name: good})


def test_invalid_sentieon_dnascope_pcr_indel_model_rejected():
    with pytest.raises(SkillError) as exc:
        _check_enums({"sentieon_dnascope_pcr_indel_model": "banana"})
    assert exc.value.error_code == ErrorCode.INVALID_FLAG_COMBINATION
    for good in ("NONE", "HOSTILE", "AGGRESSIVE", "CONSERVATIVE"):
        _check_enums({"sentieon_dnascope_pcr_indel_model": good})


# ---------------------------------------------------------------------------
# --genome iGenomes key validation
# Audit fix: SUPPORTED_IGENOMES_NAMES (the pinned conf/igenomes.config catalogue)
# and the INVALID_GENOME error code were both unused; the wrapper never caught a
# typo'd/invalid --genome key (e.g. bare "GRCh38" instead of "GATK.GRCh38") before
# launching Nextflow. These tests lock the new preflight gate.
# ---------------------------------------------------------------------------


def test_unknown_genome_rejected_with_default_igenomes_base():
    # Bare "GRCh38" is NOT a Sarek iGenomes key (the catalogue uses GATK.GRCh38,
    # Ensembl.GRCh37, NCBI.GRCh38, …) — it must fail fast with the default mirror.
    with pytest.raises(SkillError) as exc:
        _check_genome_known({"genome": "GRCh38"})
    assert exc.value.error_code == ErrorCode.INVALID_GENOME
    assert "GRCh38" in exc.value.message


def test_known_genome_keys_accepted():
    for key in ("GATK.GRCh38", "GATK.GRCh37", "Ensembl.GRCh37", "NCBI.GRCh38", "hg38", "mm10"):
        assert _check_genome_known({"genome": key}) == []


def test_absent_genome_uses_upstream_default_and_is_accepted():
    # No genome key → upstream default GATK.GRCh38; nothing to validate.
    assert _check_genome_known({}) == []
    assert _check_genome_known({"genome": ""}) == []


def test_genome_disable_sentinel_not_validated():
    # Custom-FASTA runs disable iGenomes — the genome label is then informational.
    assert _check_genome_known({"genome": "null"}) == []
    assert _check_genome_known({"genome": "GRCh38", "igenomes_ignore": True}) == []


def test_unknown_genome_with_custom_igenomes_base_warns_not_errors():
    # A custom iGenomes mirror may define other keys, so don't hard-fail — warn.
    warnings = _check_genome_known(
        {"genome": "MyCustomKey", "igenomes_base": "/data/custom_igenomes/"}
    )
    assert warnings and any("MyCustomKey" in w for w in warnings)


def test_testdata_genome_rejected_without_matching_igenomes_base():
    # nf-core/sarek's `testdata.nf-core.sarek` only resolves under the test-datasets
    # igenomes_base (conf/test.config), not the default s3://ngi-igenomes mirror. Using
    # it in a normal run without that override fails late at nf-schema validation, so
    # preflight must reject it early (finding #2).
    with pytest.raises(SkillError) as exc:
        _check_genome_known({"genome": "testdata.nf-core.sarek"})
    assert exc.value.error_code == ErrorCode.INVALID_GENOME
    assert "testdata.nf-core.sarek" in exc.value.message
    # The fix must point at both remedies.
    assert "igenomes-base" in exc.value.fix or "igenomes_base" in exc.value.fix
    assert "--demo" in exc.value.fix


def test_testdata_genome_accepted_with_test_profile():
    assert _check_genome_known({"genome": "testdata.nf-core.sarek", "profile": "test,docker"}) == []


def test_testdata_genome_accepted_with_test_datasets_igenomes_base():
    base = "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/"
    assert _check_genome_known({"genome": "testdata.nf-core.sarek", "igenomes_base": base}) == []


def test_macos_docker_tmp_warns(monkeypatch):
    """On macOS + Docker, an output dir under /tmp warns (Colima does not share
    /tmp into its VM). Mirrors nfcore-scrnaseq-wrapper's guard."""
    import preflight as p
    from pathlib import Path as _P
    monkeypatch.setattr(p.sys, "platform", "darwin")
    warns = p._check_macos_docker_tmp("docker", _P("/tmp/run_out"))
    assert warns and any("/tmp" in w for w in warns)
    # Non-docker backend, or output not under /tmp, does not warn.
    assert p._check_macos_docker_tmp("singularity", _P("/tmp/run_out")) == []


def test_macos_docker_tmp_skipped_off_darwin(monkeypatch):
    import preflight as p
    from pathlib import Path as _P
    monkeypatch.setattr(p.sys, "platform", "linux")
    assert p._check_macos_docker_tmp("docker", _P("/tmp/run_out")) == []


def test_unknown_genome_exempt_under_upstream_test_profile():
    # An upstream test/test_* profile supplies its own reference, mirroring the
    # GENOME_REQUIRED exemption, so it must not raise INVALID_GENOME.
    assert _check_genome_known({"genome": "GRCh38", "profile": "docker,test_full"}) == []
    assert _check_genome_known({"genome": "GRCh38", "profile": "test"}) == []


def test_unknown_genome_surfaces_in_run_preflight(tmp_path, stub_binaries):
    out = tmp_path / "out"
    out.mkdir()
    params = make_params(genome="GRCh38")
    params["profile"] = "docker"
    with pytest.raises(SkillError) as exc:
        run_preflight(
            params=params,
            samplesheet=make_samplesheet(),
            pipeline_source=make_pipeline_source(),
            output_dir=out,
            repo_root=tmp_path / "repo",
        )
    assert exc.value.error_code == ErrorCode.INVALID_GENOME


# ---------------------------------------------------------------------------
# Resume drift
# ---------------------------------------------------------------------------


def _resume_check(params, samplesheet, manifest, pipeline_source=None):
    return _check_resume_drift(
        params=params,
        samplesheet=samplesheet,
        pipeline_source=pipeline_source or {"source_kind": "remote_repo", "resolved_version": "3.8.1"},
        manifest=manifest,
    )


def _base_resume_state():
    params = {
        "step": "mapping",
        "aligner": "bwa-mem",
        "tools": ["mutect2"],
        "skip_tools": [],
        "joint_germline": False,
        "joint_mutect2": True,
        "wes": False,
        "profile": "docker",
        "params_checksum": "sha256:abc",
    }
    samplesheet = {
        "analysis_mode": "somatic_paired",
        "samplesheet_checksum": "sha256:xyz",
    }
    manifest = {
        "step": "mapping",
        "aligner": "bwa-mem",
        "tools": ["mutect2"],
        "skip_tools": [],
        "analysis_mode": "somatic_paired",
        "joint_germline": False,
        "joint_mutect2": True,
        "wes": False,
        "profile": "docker",
        "params_checksum": "sha256:abc",
        "samplesheet_checksum": "sha256:xyz",
        "pipeline_source": {"source_kind": "remote_repo", "resolved_version": "3.8.1"},
    }
    return params, samplesheet, manifest


@pytest.mark.parametrize(
    "field,new_value",
    [
        ("step", "variant_calling"),
        ("aligner", "bwa-mem2"),
        ("joint_germline", True),
        ("joint_mutect2", False),
        ("wes", True),
        # "arm" is omitted: it's always appended as the arm64 profile token,
        # so drift on "profile" covers it.  See _RESUME_TRACKED_FIELDS.
        ("profile", "singularity"),
        ("params_checksum", "sha256:zzz"),
    ],
)
def test_resume_drift_per_field(field, new_value):
    params, samplesheet, manifest = _base_resume_state()
    params[field] = new_value
    with pytest.raises(SkillError) as exc:
        _resume_check(params, samplesheet, manifest)
    assert exc.value.error_code == ErrorCode.RESUME_DRIFT
    assert field in exc.value.details["diffs"]


# ---------------------------------------------------------------------------
# Happy path: full run_preflight
# ---------------------------------------------------------------------------


def test_run_preflight_happy_path(tmp_path, stub_binaries):
    repo_root = tmp_path / "repo"
    repo_root.mkdir()
    output_dir = tmp_path / "results"

    params = make_params(
        step="mapping",
        aligner="bwa-mem",
        profile="docker",
        tools=["haplotypecaller"],
        genome="GATK.GRCh38",
    )
    samplesheet = make_samplesheet(analysis_mode="germline")
    pipeline_source = make_pipeline_source()

    result = run_preflight(
        params=params,
        samplesheet=samplesheet,
        pipeline_source=pipeline_source,
        output_dir=output_dir,
        repo_root=repo_root,
    )
    assert isinstance(result, PreflightResult)


# --- --demo offline gate (nf-core -profile test is remote by design) ---------

def test_check_demo_network_offline_test_profile_raises(monkeypatch):
    """--demo composes the upstream `-profile test`, whose inputs/references are
    remote GitHub URLs; under NXF_OFFLINE it must fail fast with a clear
    DEMO_REQUIRES_NETWORK error instead of a cryptic Nextflow abort."""
    import preflight as _pf

    monkeypatch.setenv("NXF_OFFLINE", "true")
    with pytest.raises(SkillError) as exc:
        _pf._check_demo_network({"profile": "docker,test"})
    assert exc.value.error_code == ErrorCode.DEMO_REQUIRES_NETWORK


def test_check_demo_network_offline_real_run_not_blocked(monkeypatch):
    """A real (non-test-profile) run is fully local; NXF_OFFLINE must not block it."""
    import preflight as _pf

    monkeypatch.setenv("NXF_OFFLINE", "true")
    _pf._check_demo_network({"profile": "docker"})


def test_check_demo_network_test_profile_online_not_blocked(monkeypatch):
    """A test-profile run with network available must not trip the guard."""
    import preflight as _pf

    monkeypatch.delenv("NXF_OFFLINE", raising=False)
    _pf._check_demo_network({"profile": "docker,test"})


def test_run_preflight_demo_offline_raises(tmp_path, stub_binaries, monkeypatch):
    repo_root = tmp_path / "repo"
    repo_root.mkdir()
    params = make_params(step="mapping", aligner="bwa-mem", profile="docker,test", tools=["haplotypecaller"])
    monkeypatch.setenv("NXF_OFFLINE", "true")
    with pytest.raises(SkillError) as exc:
        run_preflight(
            params=params,
            samplesheet=make_samplesheet(analysis_mode="germline"),
            pipeline_source=make_pipeline_source(),
            output_dir=tmp_path / "results",
            repo_root=repo_root,
        )
    assert exc.value.error_code == ErrorCode.DEMO_REQUIRES_NETWORK


# --- --allow-remote-inputs gate (local-first by default) ---------------------

def test_remote_samplesheet_input_rejected_by_default():
    import preflight as _pf
    with pytest.raises(SkillError) as exc:
        _pf._check_remote_inputs(
            {}, {"fastq_paths": ["s3://bucket/n_R1.fastq.gz"], "bam_paths": [], "cram_paths": []},
            allow_remote_inputs=False,
        )
    assert exc.value.error_code == ErrorCode.REMOTE_INPUT_NOT_ALLOWED


def test_remote_reference_param_rejected_by_default():
    import preflight as _pf
    with pytest.raises(SkillError) as exc:
        _pf._check_remote_inputs(
            {"fasta": "https://example.org/genome.fa"}, {}, allow_remote_inputs=False
        )
    assert exc.value.error_code == ErrorCode.REMOTE_INPUT_NOT_ALLOWED


def test_igenomes_base_remote_default_is_not_gated():
    """The public iGenomes mirror base (s3://ngi-igenomes/) is remote by design and
    must NOT trip the local-first gate."""
    import preflight as _pf
    _pf._check_remote_inputs(
        {"igenomes_base": "s3://ngi-igenomes/igenomes/"}, {}, allow_remote_inputs=False
    )


def test_remote_inputs_allowed_with_flag_warns(capsys):
    import preflight as _pf
    _pf._check_remote_inputs(
        {}, {"fastq_paths": ["gs://bucket/n_R1.fastq.gz"]}, allow_remote_inputs=True
    )
    assert "--allow-remote-inputs" in capsys.readouterr().err


def test_local_inputs_pass_the_remote_gate():
    import preflight as _pf
    _pf._check_remote_inputs({"fasta": "/refs/genome.fa"}, {"fastq_paths": ["/d/n_R1.fastq.gz"]}, allow_remote_inputs=False)


# nf-core/sarek 3.8.1's nextflow_schema.json declares these reference params
# `exists: true` (known_indels as a file-path-pattern that must match a file), most with
# a ^\S+ name pattern, and Sarek reads snpsift_databases through its own
# assets/schema_snpsift_databases.json and checkIfExists: nf-schema and the pipeline
# refuse a missing path, a wrong name or an intervals file that is neither .bed nor
# .interval_list when Nextflow starts, before any task runs.
_SCHEMA_CHECKED_PATHS = (
    "fasta", "fasta_fai", "dict", "bwa", "bwamem2", "dragmap", "dbsnp", "dbsnp_tbi",
    "known_indels", "known_indels_tbi", "known_snps", "known_snps_tbi", "pon", "pon_tbi",
    "intervals", "ascat_alleles", "ascat_loci", "ascat_loci_gc", "ascat_loci_rt", "chr_dir",
    "msisensor2_models", "msisensorpro_scan", "dbnsfp", "dbnsfp_tbi", "mastermind_file",
    "phenotypes_file_tbi", "spliceai_snv", "spliceai_snv_tbi", "spliceai_indel",
    "spliceai_indel_tbi", "bcftools_annotations", "bcftools_annotations_tbi", "condel_config",
    "cf_chrom_len", "bbsplit_fasta_list", "bbsplit_index", "snpsift_databases",
)


@pytest.mark.parametrize("name", _SCHEMA_CHECKED_PATHS)
def test_paths_the_schema_checks_are_left_to_nf_schema(tmp_path, name):
    _check_reference_paths({name: str(tmp_path / "missing" / f"{name}.x")})


def test_known_indels_pattern_matching_nothing_is_left_to_nf_schema(tmp_path):
    _check_reference_paths({"known_indels": str(tmp_path / "{a,b}.vcf.gz")})


@pytest.mark.parametrize("name", [
    "germline_resource", "germline_resource_tbi", "mappability", "ngscheckmate_bed",
    "phenotypes_file", "sentieon_dnascope_model", "snpeff_cache", "vep_cache", "cnvkit_reference",
    "bcftools_header_lines", "bcftools_columns", "varlociraptor_scenario_somatic", "multiqc_config",
])
def test_paths_the_schema_does_not_check_must_exist(tmp_path, name):
    """No `exists` in the 3.8.1 schema: a missing local path would fail inside a task."""
    assert name not in _SCHEMA_CHECKED_PATHS
    with pytest.raises(SkillError) as exc:
        _check_reference_paths({name: str(tmp_path / "missing" / name)})
    assert exc.value.error_code == ErrorCode.MISSING_REFERENCE


@pytest.mark.parametrize("name,refused", [
    ("germline_resource", True), ("germline_resource_tbi", True), ("intervals", True),
    ("fasta", False), ("dbsnp", False), ("pon", False),
])
def test_whitespace_is_refused_only_where_the_schema_pattern_allows_it(name, refused):
    """germline_resource(_tbi) and intervals have unanchored patterns (\\S+\\.vcf\\.gz$),
    which a path under a directory with spaces matches; the other patterns start ^\\S+."""
    import preflight as p
    value = {"germline_resource": "/my refs/g.vcf.gz", "germline_resource_tbi": "/my refs/g.vcf.gz.tbi",
             "intervals": "/my refs/t.bed", "fasta": "/my refs/g.fa", "dbsnp": "/my refs/d.vcf.gz",
             "pon": "/my refs/p.vcf.gz"}[name]
    if refused:
        with pytest.raises(SkillError) as exc:
            p._check_reference_path_whitespace({name: value})
        assert exc.value.error_code == ErrorCode.MISSING_REFERENCE
    else:
        p._check_reference_path_whitespace({name: value})


def test_wrong_intervals_extension_is_left_to_sarek():
    params = make_params(intervals="/refs/targets.txt", genome="null", igenomes_ignore=True,
                         fasta="/refs/genome.fa", skip_tools=["baserecalibrator"])
    _check_flag_compatibility(params=params, tools=[])


def test_snpsift_databases_csv_rules_are_sareks(tmp_path):
    """assets/schema_snpsift_databases.json, the fields-or-vardb error and checkIfExists
    in main.nf are Sarek's; the wrapper refuses only a vardb path with whitespace."""
    dbs = tmp_path / "snpsift_databases.csv"
    dbs.write_text("vcf,tbi,fields,prefix,vardb\n/db/dbsnp.txt,,,,\n", encoding="utf-8")
    _check_flag_compatibility(params=make_params(tools=["snpsift"], snpsift_databases=str(dbs)), tools=["snpsift"])
    dbs.write_text("vcf,tbi,fields,prefix,vardb\n/db/dbsnp.vcf.gz,,,,/my dbs/x.snpsift.vardb\n", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        _check_flag_compatibility(params=make_params(tools=["snpsift"], snpsift_databases=str(dbs)), tools=["snpsift"])
    assert exc.value.error_code == ErrorCode.MISSING_SNPSIFT_DATABASES


# The minimum versions are Nextflow's to enforce: each pinned release declares
# manifest.nextflowVersion with `!`, so Nextflow refuses an older engine when the run
# starts, and the Nextflow launcher refuses a Java outside 17-26 before that. The wrapper
# checks only that `java` and `nextflow` run and records their versions for the bundle.
@pytest.mark.parametrize("java_text,nextflow_text", [
    ('openjdk version "11.0.1"', "Nextflow version 24.10.0 build 5928"),
    ('openjdk version "16.0.2"', "Nextflow version 25.04.2 build 5947"),
])
def test_old_java_and_nextflow_are_left_to_nextflow(monkeypatch, java_text, nextflow_text):
    import preflight as p
    monkeypatch.setattr(p.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(p._runtime, "run_probe", lambda command, timeout_seconds: __import__("subprocess").CompletedProcess(
        command, 0, nextflow_text if command[0] == "nextflow" else "", java_text if command[0] == "java" else ""))
    assert p._check_java() == java_text.split('"')[1]
    assert p._check_nextflow() == nextflow_text.split()[2]


def test_a_launcher_refusal_is_named_when_nextflow_cannot_report_its_version(monkeypatch):
    """With an unsupported Java the launcher fails `nextflow -version` itself; the
    error must carry its reason instead of only saying the version is unknown."""
    import subprocess as sp
    import preflight as p
    from errors import SkillError as _SkillError
    refusal = ("ERROR: Cannot find Java or it's a wrong version -- please make sure that "
               "Java 17 or later (up to 26) is installed")
    monkeypatch.setattr(p.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(p._runtime, "run_probe", lambda command, timeout_seconds: sp.CompletedProcess(command, 1, "", refusal + "\n"))
    with pytest.raises(_SkillError) as exc:
        p._check_nextflow()
    assert refusal.removeprefix("ERROR: ") in exc.value.message
