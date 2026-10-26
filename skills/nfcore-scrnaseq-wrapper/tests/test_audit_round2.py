"""Regression tests for the second production-audit round.

Each test pins a finding so the corrected behaviour cannot silently regress:

* R2-1  ``--demo`` is hermetic: reference/index/protocol flags are never written
        into params.yaml (they would override the upstream ``test`` profile, since
        ``-params-file`` takes precedence over profile config) and the user is warned.
* R2-2  The Nextflow execution timeout is configurable via ``--timeout-hours``
        (0 disables it for long HPC/cloud runs); default preserves 12 h.
* R2-3  ``--expected-cells`` help documents the single-sample-only restriction.
* R2-4  The protocol matrix is cross-checked against a sibling ``assets/protocols.json``
        when a checkout is present (live drift guard, like the STAR args test).
* R2-5  The "under /tmp" predicate is a single shared helper reused by preflight
        and executor (no divergent copies).
"""

from __future__ import annotations

from argparse import Namespace
from pathlib import Path
import importlib.util
import json
import re
import sys

import pytest

_SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_SKILL_DIR))

import schemas
from params_builder import build_effective_params


def _load_skill_module():
    spec = importlib.util.spec_from_file_location(
        "nfcore_scrnaseq_wrapper", _SKILL_DIR / "nfcore_scrnaseq_wrapper.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _params_args(**overrides) -> Namespace:
    defaults = dict(
        demo=False,
        preset="standard",
        protocol="10XV3",
        email=None,
        email_on_fail=None,
        multiqc_title=None,
        multiqc_config=None,
        multiqc_logo=None,
        multiqc_methods_description=None,
        publish_dir_mode=None,
        trace_report_suffix=None,
        monochrome_logs=False,
        skip_cellbender=False,
        skip_fastqc=False,
        skip_emptydrops=False,
        skip_multiqc=False,
        skip_cellranger_renaming=False,
        skip_cellrangermulti_vdjref=False,
        fasta=None,
        gtf=None,
        transcript_fasta=None,
        txp2gene=None,
        simpleaf_index=None,
        kallisto_index=None,
        star_index=None,
        cellranger_index=None,
        barcode_whitelist=None,
        star_feature=None,
        star_ignore_sjdbgtf=False,
        seq_center=None,
        simpleaf_umi_resolution=None,
        kb_workflow=None,
        kb_t1c=None,
        kb_t2c=None,
        motifs=None,
        cellrangerarc_config=None,
        cellrangerarc_reference=None,
        cellranger_vdj_index=None,
        gex_frna_probe_set=None,
        gex_target_panel=None,
        gex_cmo_set=None,
        fb_reference=None,
        vdj_inner_enrichment_primers=None,
        gex_barcode_sample_assignment=None,
        cellranger_multi_barcodes=None,
        genome=None,
        igenomes_base=None,
        save_reference=False,
        save_align_intermeds=None,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


def _samplesheet(tmp_path: Path) -> Path:
    ss = tmp_path / "reproducibility" / "samplesheet.valid.csv"
    ss.parent.mkdir(parents=True, exist_ok=True)
    ss.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    return ss


# ── R2-1: --demo is hermetic (no reference/protocol leakage into params.yaml) ──


def test_demo_does_not_write_reference_or_protocol_params(tmp_path):
    ss = _samplesheet(tmp_path)
    fa = tmp_path / "g.fa"
    fa.write_text(">c\nA\n", encoding="utf-8")
    gtf = tmp_path / "g.gtf"
    gtf.write_text("x\n", encoding="utf-8")
    args = _params_args(
        demo=True,
        preset="star",
        protocol="10XV3",
        genome="GRCh38",
        fasta=str(fa),
        gtf=str(gtf),
        igenomes_base="/mnt/igenomes",
        barcode_whitelist=str(fa),
    )
    params = build_effective_params(
        args, normalized_samplesheet=ss, output_dir=tmp_path
    )
    # The test profile owns input + references + protocol; none may be overridden.
    for leaked in (
        "genome",
        "fasta",
        "gtf",
        "igenomes_base",
        "protocol",
        "barcode_whitelist",
    ):
        assert leaked not in params, (
            f"--demo must not write {leaked!r} (would override the test profile)"
        )
    # The hermetic demo guarantees stay intact.
    assert params["aligner"] == "star"
    assert params.get("igenomes_ignore") is True
    assert "input" not in params


def test_non_demo_still_writes_reference_and_protocol_params(tmp_path):
    ss = _samplesheet(tmp_path)
    fa = tmp_path / "g.fa"
    fa.write_text(">c\nA\n", encoding="utf-8")
    gtf = tmp_path / "g.gtf"
    gtf.write_text("x\n", encoding="utf-8")
    args = _params_args(
        demo=False, preset="star", protocol="10XV3", fasta=str(fa), gtf=str(gtf)
    )
    params = build_effective_params(
        args, normalized_samplesheet=ss, output_dir=tmp_path
    )
    assert params["protocol"] == "10XV3"
    assert params["fasta"] == fa.resolve().as_posix()
    assert params["gtf"] == gtf.resolve().as_posix()


def test_demo_ignored_flags_helper_lists_supplied_flags():
    module = _load_skill_module()
    args = _params_args(
        demo=True, preset="star", genome="GRCh38", fasta="/x/g.fa", protocol="10XV3"
    )
    ignored = module._demo_ignored_flags(args)
    assert "--genome" in ignored
    assert "--fasta" in ignored
    assert "--protocol" in ignored
    # nothing supplied → empty
    assert (
        module._demo_ignored_flags(
            _params_args(demo=True, preset="star", protocol=None)
        )
        == []
    )


# ── R2-2: configurable execution timeout ───────────────────────────────────────


def test_timeout_hours_default_preserves_twelve_hours():
    module = _load_skill_module()
    parser = module.build_parser()
    args = parser.parse_args(["--output", "/tmp/out", "--demo"])
    assert (
        module._resolve_timeout_seconds(args)
        == schemas.DEFAULT_TIMEOUT_SECONDS
        == 60 * 60 * 12
    )


def test_timeout_hours_zero_disables_the_cap():
    module = _load_skill_module()
    parser = module.build_parser()
    args = parser.parse_args(["--output", "/tmp/out", "--demo", "--timeout-hours", "0"])
    assert module._resolve_timeout_seconds(args) is None


def test_timeout_hours_custom_value_converts_to_seconds():
    module = _load_skill_module()
    parser = module.build_parser()
    args = parser.parse_args(
        ["--output", "/tmp/out", "--demo", "--timeout-hours", "1.5"]
    )
    assert module._resolve_timeout_seconds(args) == 5400


def test_timeout_hours_rejects_negative():
    module = _load_skill_module()
    parser = module.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["--output", "/tmp/out", "--demo", "--timeout-hours", "-3"])


def test_executor_without_a_timeout_never_times_out(tmp_path, monkeypatch):
    """timeout_seconds=None sets no wall-clock cap: the run completes normally
    instead of raising EXECUTION_FAILED, however long Nextflow keeps running."""
    import subprocess
    import executor

    class _Proc:
        pid = 4321

        def __init__(self):
            self.checks = 0

        def wait(self, timeout=None):
            # Still running for the first checks; without a cap none may end the run.
            self.checks += 1
            if self.checks < 5:
                raise subprocess.TimeoutExpired(["nextflow"], timeout)
            return 0

        def kill(self):  # pragma: no cover - not expected
            pass

    proc = _Proc()
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    result = executor.execute_nextflow(
        ["nextflow", "run", "test"],
        cwd=tmp_path,
        output_dir=tmp_path,
        timeout_seconds=None,
    )
    assert result["exit_code"] == 0
    assert proc.checks == 5


# ── R2-3: --expected-cells documents single-sample restriction ─────────────────


def test_expected_cells_help_documents_single_sample_restriction():
    module = _load_skill_module()
    parser = module.build_parser()
    action = next(a for a in parser._actions if "--expected-cells" in a.option_strings)
    assert re.search(
        r"single[- ]sample|one sample|per[- ]row", action.help, re.IGNORECASE
    ), (
        "--expected-cells help must document that the global override is single-sample only"
    )


# ── R2-4: protocol matrix cross-checked against a sibling protocols.json ───────


def _normalise_token(text: str) -> str:
    return re.sub(r"[\s_-]+", "", text.lower())


def test_protocol_matrix_matches_sibling_protocols_json_if_present():
    from nfcore_4_1_0_contract import PROTOCOLS_JSON_4_1_0

    candidates = [
        schemas.DEFAULT_LOCAL_PIPELINE_DIR / "assets" / "protocols.json",
        schemas.REPO_PARENT / "scrnaseq" / "assets" / "protocols.json",
    ]
    protocols_file = next((p for p in candidates if p.is_file()), None)
    if protocols_file is None:
        pytest.skip("no sibling nf-core/scrnaseq checkout available to cross-check")
    upstream = json.loads(protocols_file.read_text(encoding="utf-8"))
    # aligner key 'simpleaf' maps to the wrapper preset 'standard'.
    aligner_to_preset = {
        "simpleaf": "standard",
        "star": "star",
        "kallisto": "kallisto",
        "cellranger": "cellranger",
        "cellrangerarc": "cellrangerarc",
    }
    for aligner, preset in aligner_to_preset.items():
        assert aligner in upstream, f"{aligner} missing from upstream protocols.json"
        expected = {_normalise_token(key) for key in upstream[aligner].keys()}
        actual = set(PROTOCOLS_JSON_4_1_0[preset])
        assert actual == expected, (
            f"protocol matrix for {preset!r} drifted from upstream {aligner!r}: "
            f"wrapper={sorted(actual)} upstream={sorted(expected)}"
        )
    # cellrangermulti is samplesheet-driven and must stay out of the matrix.
    assert "cellrangermulti" not in PROTOCOLS_JSON_4_1_0


# ── R2-6: the ClawBio runner forwards every documented BooleanOptional negative ─


def _load_clawbio():
    """Return the ClawBio runner engine.

    The runner registry (``SKILLS``) and the ``run_skill`` / ``main`` callables
    live in :mod:`clawbio.cli`; the top-level ``clawbio.py`` is a thin shim that
    only re-exports and delegates to them. Tests inspect and monkeypatch the
    engine, so they must target ``clawbio.cli`` directly — patching a name on the
    shim would not change the binding ``main`` resolves inside ``clawbio.cli``.
    """
    import clawbio.cli as clawbio_cli

    return clawbio_cli


def test_runner_allowlist_includes_no_save_align_intermeds():
    """The runner silently drops flags absent from the allowlist; the documented
    --no-save-align-intermeds negative must be forwardable (it controls disk use
    on large runs) and treated as a value-less flag."""
    clawbio = _load_clawbio()
    entry = clawbio.SKILLS["scrnaseq-pipeline"]
    assert "--no-save-align-intermeds" in entry["allowed_extra_flags"]
    assert "--no-save-align-intermeds" in entry["allowed_extra_flags_without_values"]


def test_runner_help_exposes_no_save_align_intermeds_switch():
    """The public runner should expose the same two-way switch as the wrapper,
    instead of relying on parse_known_args to pass the negative form through."""
    import subprocess

    clawbio_path = _SKILL_DIR.parent.parent / "clawbio.py"
    result = subprocess.run(
        [sys.executable, str(clawbio_path), "run", "--help"],
        check=True,
        capture_output=True,
        text=True,
    )

    assert "--no-save-align-intermeds" in result.stdout


def test_runner_forwards_explicit_no_save_align_intermeds_and_new_safety_flags(
    monkeypatch,
):
    """Flags parsed by the top-level runner must still reach the wrapper; the
    security allowlist only sees unknown extras, so native flags need explicit
    forwarding."""
    clawbio = _load_clawbio()
    captured = {}

    def fake_run_skill(**kwargs):
        captured.update(kwargs)
        return {
            "success": True,
            "exit_code": 0,
            "duration_seconds": 0,
            "output_dir": None,
            "files": [],
            "stdout": "",
            "stderr": "",
        }

    monkeypatch.setattr(clawbio, "run_skill", fake_run_skill)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "clawbio.py",
            "run",
            "scrnaseq-pipeline",
            "--output",
            "out",
            "--no-save-align-intermeds",
            "--work-dir",
            "s3://bucket/scrnaseq/work",
            "--require-local-pipeline",
            "--allow-conda-cellranger",
        ],
    )

    with pytest.raises(SystemExit) as exc:
        clawbio.main()

    assert exc.value.code == 0
    assert captured["skill_name"] == "scrnaseq-pipeline"
    assert captured["extra_args"] == [
        "--work-dir",
        "s3://bucket/scrnaseq/work",
        "--require-local-pipeline",
        "--allow-conda-cellranger",
        "--no-save-align-intermeds",
    ]


def test_runner_allowlist_includes_nfcore_aligner_alias():
    clawbio = _load_clawbio()
    entry = clawbio.SKILLS["scrnaseq-pipeline"]
    assert "--aligner" in entry["allowed_extra_flags"]


def test_runner_value_less_flags_are_marked_without_values():
    """Every value-less wrapper flag in the allowlist must also be in
    allowed_extra_flags_without_values, or the runner consumes the next token as
    its value (corrupting the forwarded command)."""
    clawbio = _load_clawbio()
    module = _load_skill_module()
    no_value = {
        o
        for a in module.build_parser()._actions
        if a.nargs == 0
        for o in a.option_strings
    }
    entry = clawbio.SKILLS["scrnaseq-pipeline"]
    allow = set(entry.get("allowed_extra_flags", set()))
    without = set(entry.get("allowed_extra_flags_without_values", set()))
    missing = (allow & no_value) - without
    assert missing == set(), (
        f"value-less flags not marked without-values: {sorted(missing)}"
    )


def test_runner_allowlist_includes_config_short_flag():
    """README documents `-c/--config`; the runner must forward the short form too,
    or an HPC user's `-c site.config` is silently dropped before reaching Nextflow."""
    clawbio = _load_clawbio()
    entry = clawbio.SKILLS["scrnaseq-pipeline"]
    assert "-c" in entry["allowed_extra_flags"]


def test_every_wrapper_cli_flag_is_forwardable_or_native():
    """No wrapper CLI flag may be silently dropped by the runner: each must be in
    the allowlist or be one of the natively-handled flags (--input/--output/--demo).
    Guards the whole class of 'documented flag silently dropped' bugs."""
    clawbio = _load_clawbio()
    module = _load_skill_module()
    cli = {
        o
        for a in module.build_parser()._actions
        for o in a.option_strings
        if o.startswith("-")
    }
    cli -= {"--help", "-h"}
    allow = set(clawbio.SKILLS["scrnaseq-pipeline"].get("allowed_extra_flags", set()))
    native = {"--input", "--output", "--demo"}
    dropped = cli - allow - native
    assert dropped == set(), (
        f"runner would silently drop wrapper flags: {sorted(dropped)}"
    )


def test_runner_allowlist_forwards_every_wrapper_boolean_optional_negative():
    """Any --no-* negative registered by the wrapper CLI must be in the runner
    allowlist, or it would be silently dropped before reaching the wrapper."""
    clawbio = _load_clawbio()
    module = _load_skill_module()
    negatives = {
        o
        for a in module.build_parser()._actions
        for o in a.option_strings
        if o.startswith("--no-")
    }
    allow = set(clawbio.SKILLS["scrnaseq-pipeline"].get("allowed_extra_flags", set()))
    missing = negatives - allow
    assert missing == set(), (
        f"runner allowlist drops documented negatives: {sorted(missing)}"
    )


# ── R2-7: SKILL.md output tree matches where files are actually written ────────

_PROVENANCE_JSONS = (
    "inputs.json",
    "invocation.json",
    "preflight.json",
    "upstream.json",
    "outputs.json",
    "runtime.json",
    "skill.json",
)


def test_skill_output_tree_places_provenance_json_under_reproducibility():
    """provenance.write_provenance_bundle writes the 7 JSONs into reproducibility/
    (not a separate provenance/ tree). SKILL.md must document them there; the
    provenance/ directory holds only handoff.json (only with --run-downstream)."""
    text = (_SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
    repro_idx = text.index("reproducibility/")
    prov_idx = text.index(
        "provenance/                       # Written only with --run-downstream"
    )
    assert repro_idx < prov_idx, "reproducibility/ must precede provenance/ in the tree"
    for name in _PROVENANCE_JSONS:
        idx = text.index(name)
        assert repro_idx < idx < prov_idx, (
            f"{name} must be documented under reproducibility/, not provenance/"
        )
    # provenance/ documents only handoff.json
    assert text.index("handoff.json") > prov_idx


def test_provenance_bundle_writes_jsons_into_reproducibility_not_provenance(tmp_path):
    """Lock the runtime location the doc now describes."""
    import provenance

    args = Namespace(
        preset="star",
        profile="docker",
        pipeline_version="4.1.0",
        resume=False,
        demo=False,
        check=False,
        allow_dirty_pipeline=False,
        extra_config=[],
    )
    pipeline_source = {
        "source_kind": "remote_repo",
        "source_ref": "nf-core/scrnaseq",
        "resolved_version": "4.1.0",
        "branch": "",
        "dirty": False,
    }
    ss = tmp_path / "reproducibility" / "samplesheet.valid.csv"
    ss.parent.mkdir(parents=True, exist_ok=True)
    ss.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    params_path = tmp_path / "reproducibility" / "params.yaml"
    params_path.write_text("aligner: star\n", encoding="utf-8")
    logs = tmp_path / "logs"
    logs.mkdir()
    (logs / "stdout.txt").write_text("", encoding="utf-8")
    (logs / "stderr.txt").write_text("", encoding="utf-8")
    provenance.write_provenance_bundle(
        tmp_path,
        args=args,
        pipeline_source=pipeline_source,
        preflight_result={
            "java": {"version": "17"},
            "nextflow": {"version": "25.04.0"},
            "references": {},
        },
        params_path=params_path,
        params_payload={"aligner": "star"},
        normalized_samplesheet=ss,
        samplesheet_summary={
            "sample_count": 0,
            "fastq_paths": [],
            "unknown_columns": [],
        },
        parsed_outputs={"h5ad_candidates": [], "preferred_h5ad": ""},
        execution_result={
            "stdout_path": str(logs / "stdout.txt"),
            "stderr_path": str(logs / "stderr.txt"),
        },
        command_str="nextflow run ...",
    )
    for name in _PROVENANCE_JSONS:
        assert (tmp_path / "reproducibility" / name).exists(), (
            f"{name} must be in reproducibility/"
        )
    assert not (tmp_path / "provenance").exists(), (
        "no separate provenance/ dir without downstream handoff"
    )


# ── R2-8: preflight records only provided references (no empty-string noise) ───


def _ref_args(**over):
    base = dict(
        preset="star",
        genome=None,
        fasta=None,
        gtf=None,
        transcript_fasta=None,
        txp2gene=None,
        simpleaf_index=None,
        kallisto_index=None,
        star_index=None,
        cellranger_index=None,
        barcode_whitelist=None,
        cellrangerarc_reference=None,
        cellrangerarc_config=None,
        motifs=None,
        kb_t1c=None,
        kb_t2c=None,
        cellranger_vdj_index=None,
        gex_frna_probe_set=None,
        gex_target_panel=None,
        gex_cmo_set=None,
        fb_reference=None,
        vdj_inner_enrichment_primers=None,
        gex_barcode_sample_assignment=None,
        cellranger_multi_barcodes=None,
    )
    base.update(over)
    return Namespace(**base)


def test_check_references_records_only_provided_references(tmp_path):
    """preflight_result['references'] flows into inputs.json/preflight.json; it must
    list only references the user actually supplied, not 20+ empty-string fields."""
    import preflight

    fa = tmp_path / "g.fa"
    fa.write_text(">c\nA\n", encoding="utf-8")
    gtf = tmp_path / "g.gtf"
    gtf.write_text("x\n", encoding="utf-8")
    resolved = preflight._check_references(
        _ref_args(preset="star", fasta=str(fa), gtf=str(gtf))
    )
    assert set(resolved) == {"fasta", "gtf"}
    assert all(v for v in resolved.values()), (
        "no empty-string reference entries may be recorded"
    )


def test_check_references_preserves_genome_shortcut(tmp_path):
    import preflight

    resolved = preflight._check_references(_ref_args(preset="star", genome="GRCh38"))
    assert resolved == {"genome": "GRCh38"}


# ── R2-9: every error_code raised in the codebase is a defined ErrorCode ───────


def test_all_error_codes_used_in_source_are_defined_constants():
    """Guard against typo'd or undefined error codes: every literal error_code="X"
    and every ErrorCode.X used across the skill source must resolve to a constant."""
    import inspect
    import glob
    import re
    from errors import ErrorCode

    consts = {
        v
        for k, v in inspect.getmembers(ErrorCode)
        if not k.startswith("_") and isinstance(v, str)
    }
    const_names = {k for k, _ in inspect.getmembers(ErrorCode) if not k.startswith("_")}

    undefined = []
    for path in glob.glob(str(_SKILL_DIR / "*.py")):
        text = Path(path).read_text(encoding="utf-8")
        for m in re.finditer(r'error_code\s*=\s*"([A-Z_]+)"', text):
            if m.group(1) not in consts:
                undefined.append((Path(path).name, m.group(1)))
        for m in re.finditer(r"error_code\s*=\s*ErrorCode\.([A-Za-z_]+)", text):
            if m.group(1) not in const_names:
                undefined.append((Path(path).name, "ErrorCode." + m.group(1)))
    assert undefined == [], f"undefined error codes used in source: {undefined}"


def test_dynamic_executable_error_codes_are_defined():
    """_check_executable builds MISSING_<TOOL> dynamically; for every executable the
    wrapper actually checks, the produced code must be a defined ErrorCode constant."""
    from errors import ErrorCode

    for tool in ("java", "nextflow"):
        code = f"MISSING_{tool.upper().replace('-', '_')}"
        assert getattr(ErrorCode, code, None) == code, (
            f"{code} must be a defined ErrorCode"
        )


# ── R2-10: documentation ↔ code parity for parameter governance + metadata ─────


def _read_yaml_frontmatter():
    import yaml

    text = (_SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
    return yaml.safe_load(text.split("---")[1])


def test_unsupported_param_lists_in_docs_match_code():
    """SKILL.md and README 'Intentionally unsupported' lists must equal the code's
    INTENTIONALLY_UNSUPPORTED_PARAMS, and supported∪unsupported must partition the
    official params with no overlap."""
    from nfcore_4_1_0_contract import (
        INTENTIONALLY_UNSUPPORTED_PARAMS,
        WRAPPER_SUPPORTED_UPSTREAM_PARAMS,
        OFFICIAL_PARAMS,
    )

    code = set(INTENTIONALLY_UNSUPPORTED_PARAMS)
    for doc in ("SKILL.md", "README.md"):
        text = (_SKILL_DIR / doc).read_text(encoding="utf-8")
        m = re.search(
            r"Intentionally unsupported upstream parameters:\*\*\s*(.+)", text
        )
        assert m, f"{doc} missing the unsupported-params line"
        listed = set(re.findall(r"`([a-z_]+)`", m.group(1)))
        assert listed == code, (
            f"{doc} unsupported list drifted from code: {listed ^ code}"
        )
    assert set(OFFICIAL_PARAMS) == set(WRAPPER_SUPPORTED_UPSTREAM_PARAMS) | code
    assert not (set(WRAPPER_SUPPORTED_UPSTREAM_PARAMS) & code)


def test_catalog_entry_mirrors_skill_frontmatter():
    """catalog.json's scrnaseq entry must mirror the authored SKILL.md frontmatter
    (description, tags, version, trigger_keywords) so the machine index never drifts."""
    import json

    fm = _read_yaml_frontmatter()
    meta = fm.get("metadata", {})
    catalog = json.loads(
        (_SKILL_DIR.parent / "catalog.json").read_text(encoding="utf-8")
    )
    items = catalog if isinstance(catalog, list) else catalog.get("skills", [])
    entry = next(e for e in items if e.get("name") == "nfcore-scrnaseq-wrapper")
    # Frontmatter follows the canonical SKILL-TEMPLATE schema: description is
    # top-level; version/tags under metadata:; trigger_keywords under
    # metadata.openclaw: (mirrors nfcore-sarek/rnaseq test_catalog_consistency).
    assert entry["description"] == fm["description"]
    assert entry["version"] == meta["version"]
    assert entry["tags"] == meta["tags"]
    assert entry["trigger_keywords"] == meta["openclaw"]["trigger_keywords"]


def test_declared_packages_cover_external_imports():
    """If the skill imports an external package (e.g. yaml/PyYAML), it must be
    declared in the SKILL.md frontmatter packages list and catalog dependencies."""
    import glob
    import json

    sources = "\n".join(
        Path(p).read_text(encoding="utf-8") for p in glob.glob(str(_SKILL_DIR / "*.py"))
    )
    fm = _read_yaml_frontmatter()
    declared = {p.lower() for p in fm["metadata"]["dependencies"].get("packages", [])}
    if re.search(r"^\s*import yaml\b", sources, re.MULTILINE):
        assert "pyyaml" in declared, (
            "code imports yaml but pyyaml is not in frontmatter packages"
        )
        catalog = json.loads(
            (_SKILL_DIR.parent / "catalog.json").read_text(encoding="utf-8")
        )
        items = catalog if isinstance(catalog, list) else catalog.get("skills", [])
        entry = next(e for e in items if e.get("name") == "nfcore-scrnaseq-wrapper")
        assert "pyyaml" in {d.lower() for d in entry["dependencies"]}, (
            "catalog deps missing pyyaml"
        )


# ── R2-5: single shared "under /tmp" predicate ─────────────────────────────────


def test_is_under_tmp_helper_behaviour():
    from clawbio.common.nfcore.paths import is_under_tmp

    assert is_under_tmp(Path("/tmp/clawbio_run")) is True
    assert is_under_tmp(Path("/private/tmp/clawbio_run")) is True
    home_path = Path.home() / "clawbio_run"
    assert is_under_tmp(home_path) is False


def test_preflight_and_executor_reuse_the_shared_helper():
    import preflight
    from clawbio.common.nfcore import execution, paths

    # The pre-run WARNING and the post-failure hint must use one predicate.
    assert preflight.is_under_tmp is paths.is_under_tmp
    assert execution.is_under_tmp is paths.is_under_tmp
