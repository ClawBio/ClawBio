"""Regression tests for the production-audit fixes.

Each test pins a specific audit finding so the corrected behaviour cannot
silently regress:

* H-1  macOS+Docker resource ceilings (cpus/memory) are emitted ONLY for demo
       runs, never for real datasets.
* H-2  Auxiliary local files (barcode whitelist, CMO/probe/feature sets, ...)
       must NOT trigger ``igenomes_ignore`` and must NOT conflict with --genome.
* H-3  The reference-field categories are centralised in ``schemas`` and reused
       by both ``params_builder`` and ``preflight`` (no divergent copies).
* H-5  The STAR_ALIGN base ext.args is a single pinned constant; if a sibling
       pipeline checkout is present, it must match ``conf/modules.config``.
* H-7  ``--skip-downstream`` documents its force-skip precedence.
"""

from __future__ import annotations

from argparse import Namespace
from pathlib import Path
import importlib.util
import re
import sys

import pytest

_SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_SKILL_DIR))

import schemas
from params_builder import build_effective_params
from repro_commands import build_macos_docker_config
from nfcore_4_1_0_contract import STAR_ALIGN_BASE_EXT_ARGS


def _load_skill_module():
    spec = importlib.util.spec_from_file_location(
        "nfcore_scrnaseq_wrapper", _SKILL_DIR / "nfcore_scrnaseq_wrapper.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _params_args(tmp_path: Path, **overrides) -> Namespace:
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


# ── H-3: centralised field categories ─────────────────────────────────────────


def test_reference_field_categories_are_centralised_and_disjoint():
    genome = set(schemas.GENOME_REFERENCE_FIELDS)
    aux = set(schemas.AUXILIARY_PATH_FIELDS)
    assert genome.isdisjoint(aux), "genome and auxiliary field sets must not overlap"
    # barcode whitelist + CMO/probe/feature sets are auxiliary, never genome refs
    for name in (
        "barcode_whitelist",
        "gex_cmo_set",
        "fb_reference",
        "cellranger_multi_barcodes",
    ):
        assert name in aux and name not in genome
    # the index/fasta/gtf fields are genome refs
    for name in (
        "fasta",
        "gtf",
        "star_index",
        "simpleaf_index",
        "kallisto_index",
        "cellranger_index",
    ):
        assert name in genome


# ── H-2: auxiliary files must not trigger igenomes_ignore ──────────────────────


def test_genome_plus_barcode_whitelist_does_not_set_igenomes_ignore(tmp_path):
    ss = _samplesheet(tmp_path)
    wl = tmp_path / "whitelist.txt"
    wl.write_text("ACGT\n", encoding="utf-8")
    args = _params_args(
        tmp_path, preset="standard", genome="GRCh38", barcode_whitelist=str(wl)
    )
    params = build_effective_params(
        args, normalized_samplesheet=ss, output_dir=tmp_path
    )
    assert params.get("genome") == "GRCh38"
    assert params.get("barcode_whitelist") == wl.resolve().as_posix()
    assert "igenomes_ignore" not in params, (
        "a barcode whitelist is not a genome reference and must not suppress iGenomes"
    )


def test_genome_plus_cmo_set_does_not_set_igenomes_ignore(tmp_path):
    ss = _samplesheet(tmp_path)
    cmo = tmp_path / "cmo.csv"
    cmo.write_text("id,name\n", encoding="utf-8")
    args = _params_args(
        tmp_path,
        preset="cellrangermulti",
        protocol=None,
        genome="GRCh38",
        gex_cmo_set=str(cmo),
    )
    params = build_effective_params(
        args, normalized_samplesheet=ss, output_dir=tmp_path
    )
    assert params.get("genome") == "GRCh38"
    assert "igenomes_ignore" not in params


# ── H-2 (preflight side): genome + auxiliary file is not a conflict ────────────


def _preflight_ref_args(tmp_path: Path, **overrides) -> Namespace:
    base = dict(
        preset="standard",
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
    base.update(overrides)
    return Namespace(**base)


def test_preflight_genome_plus_barcode_whitelist_is_not_a_conflict(tmp_path):
    import preflight

    wl = tmp_path / "whitelist.txt"
    wl.write_text("ACGT\n", encoding="utf-8")
    args = _preflight_ref_args(
        tmp_path, preset="standard", genome="GRCh38", barcode_whitelist=str(wl)
    )
    resolved = preflight._check_references(args)  # must not raise
    assert resolved["genome"] == "GRCh38"
    assert resolved["barcode_whitelist"] == str(wl)


def test_preflight_genome_plus_star_index_still_conflicts(tmp_path):
    import preflight
    from errors import SkillError

    star_idx = tmp_path / "star_idx"
    star_idx.mkdir()
    args = _preflight_ref_args(
        tmp_path, preset="star", genome="GRCh38", star_index=str(star_idx)
    )
    with pytest.raises(SkillError) as exc:
        preflight._check_references(args)
    assert exc.value.error_code == "CONFLICTING_REFERENCES"


# ── H-1 / H-5: macOS docker config split + pinned STAR args ────────────────────


def test_macos_config_demo_has_resource_limits():
    content = build_macos_docker_config(demo=True)
    assert "resourceLimits" in content
    assert "stageInMode" in content and '"copy"' in content
    assert "--outTmpDir" in content
    assert "--platform linux/amd64" in content


def test_macos_config_real_run_has_no_resource_limits():
    content = build_macos_docker_config(demo=False)
    assert "resourceLimits" not in content, (
        "real (non-demo) runs must not be capped to test-profile cpus/memory"
    )
    # but the genuine workarounds must still be present
    assert "stageInMode" in content and '"copy"' in content
    assert "--outTmpDir" in content
    assert "--platform linux/amd64" in content
    assert "STAR_ALIGN" in content


def test_build_extra_configs_real_run_has_separate_runtime_cap_on_macos(tmp_path, monkeypatch):
    module = _load_skill_module()
    monkeypatch.setattr(module.platform, "system", lambda: "Darwin")
    args = Namespace(profile="docker", demo=False)
    configs = module._build_extra_nextflow_configs(args, tmp_path)
    assert len(configs) == 2
    assert "resourceLimits" not in configs[0].read_text(encoding="utf-8")
    assert "resourceLimits" in configs[1].read_text(encoding="utf-8")


def test_build_extra_configs_demo_capped_on_macos(tmp_path, monkeypatch):
    module = _load_skill_module()
    monkeypatch.setattr(module.platform, "system", lambda: "Darwin")
    args = Namespace(profile="docker", demo=True)
    configs = module._build_extra_nextflow_configs(args, tmp_path)
    assert len(configs) == 2
    assert "resourceLimits" in configs[0].read_text(encoding="utf-8")


def test_pinned_star_args_match_sibling_checkout_if_present():
    """If a real pipeline checkout is reachable, the pinned base must match it."""
    candidates = [
        schemas.DEFAULT_LOCAL_PIPELINE_DIR / "conf" / "modules.config",
        schemas.REPO_PARENT / "scrnaseq" / "conf" / "modules.config",
    ]
    modules_config = next((p for p in candidates if p.is_file()), None)
    if modules_config is None:
        pytest.skip("no sibling nf-core/scrnaseq checkout available to cross-check")
    text = modules_config.read_text(encoding="utf-8")
    m = re.search(
        r"withName:\s*STAR_ALIGN\s*\{\s*ext\.args\s*=\s*\{\s*\"([^\"]+)\"", text
    )
    assert m, "could not locate STAR_ALIGN ext.args in conf/modules.config"
    assert STAR_ALIGN_BASE_EXT_ARGS == m.group(1).strip(), (
        "pinned STAR_ALIGN_BASE_EXT_ARGS drifted from the pipeline checkout"
    )


# ── H-10: FASTQ remap must use boundary-aware prefix matching ──────────────────


def _write_ss(path: Path, rows: list[dict]) -> None:
    import csv

    fields = ["sample", "fastq_1", "fastq_2"]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def test_remap_csv_does_not_match_sibling_directory_prefix(tmp_path):
    import remap_paths

    ss = tmp_path / "samplesheet.valid.csv"
    _write_ss(
        ss,
        [
            {
                "sample": "S",
                "fastq_1": "/data/sample10/a_R1.fastq.gz",
                "fastq_2": "/data/sample10/a_R2.fastq.gz",
            }
        ],
    )
    # /data/sample1 is a sibling of /data/sample10 — must NOT match.
    changes = remap_paths.remap_csv(ss, "/data/sample1", "/NEW", dry_run=True)
    assert changes == [], f"sibling-dir prefix must not match, got {changes}"


def test_remap_csv_matches_exact_directory_prefix(tmp_path):
    import remap_paths

    ss = tmp_path / "samplesheet.valid.csv"
    _write_ss(
        ss,
        [
            {
                "sample": "S",
                "fastq_1": "/data/sample1/a_R1.fastq.gz",
                "fastq_2": "/data/sample1/a_R2.fastq.gz",
            }
        ],
    )
    changes = remap_paths.remap_csv(ss, "/data/sample1", "/new/s1", dry_run=True)
    assert {c[2] for c in changes} == {
        "/new/s1/a_R1.fastq.gz",
        "/new/s1/a_R2.fastq.gz",
    }


# ── H-16: schema-typed params survive YAML round-trip with the right type ──────


def test_star_ignore_sjdbgtf_serialized_as_string_not_bool(tmp_path):
    """The schema types star_ignore_sjdbgtf as 'string'. It must serialize as the
    quoted string 'true' so nf-schema does not see a YAML boolean and reject the
    params file. Guards against a serializer change silently breaking validation."""
    import yaml
    from params_builder import build_effective_params, serialize_params_yaml

    ss = _samplesheet(tmp_path)
    fa = tmp_path / "g.fa"
    fa.write_text(">c\nA\n", encoding="utf-8")
    gtf = tmp_path / "g.gtf"
    gtf.write_text("x\n", encoding="utf-8")
    args = _params_args(
        tmp_path,
        preset="star",
        fasta=str(fa),
        gtf=str(gtf),
        star_feature="Gene Velocyto",
        star_ignore_sjdbgtf=True,
    )
    params = build_effective_params(
        args, normalized_samplesheet=ss, output_dir=tmp_path
    )
    reloaded = yaml.safe_load(serialize_params_yaml(params))
    v = reloaded["star_ignore_sjdbgtf"]
    assert isinstance(v, str) and v == "true", (
        f"expected string 'true', got {v!r} ({type(v).__name__})"
    )


# ── H-12: pre-execution scaffolding is allowlisted consistently ────────────────


def test_macos_docker_config_does_not_block_rerun(tmp_path):
    """All pre-execution reproducibility artifacts (samplesheet, params.yaml,
    macos_docker.config) must be treated alike by the output-dir reuse policy."""
    import preflight

    repro = tmp_path / "reproducibility"
    repro.mkdir()
    (repro / "samplesheet.valid.csv").write_text(
        "sample,fastq_1,fastq_2\n", encoding="utf-8"
    )
    (repro / "params.yaml").write_text("aligner: star\n", encoding="utf-8")
    (repro / "macos_docker.config").write_text(
        'process { stageInMode = "copy" }\n', encoding="utf-8"
    )
    # must not raise — a failed macOS+docker run leaves exactly this scaffolding
    preflight.check_output_dir_available(tmp_path, resume=False)


def test_real_output_artifacts_still_block_rerun(tmp_path):
    """A genuine prior result (report.md / upstream/) must still block a non-resume re-run."""
    import preflight
    from errors import SkillError

    (tmp_path / "report.md").write_text("# prior run\n", encoding="utf-8")
    with pytest.raises(SkillError) as exc:
        preflight.check_output_dir_available(tmp_path, resume=False)
    assert exc.value.error_code == "OUTPUT_DIR_NOT_EMPTY"


def test_output_dir_inside_repo_is_rejected(tmp_path, monkeypatch):
    """Pipeline outputs must not be written inside the ClawBio source tree — parity
    with nfcore-rnaseq / nfcore-sarek, which both refuse this. The dedicated
    OUTPUT_DIR_INSIDE_REPO code is self-describing (vs the generic
    OUTPUT_DIR_NOT_WRITABLE). PROJECT_ROOT is monkeypatched to a temp dir so the
    check is hermetic and never materialises a directory inside the real repo."""
    import preflight
    from errors import SkillError

    monkeypatch.setattr(preflight, "PROJECT_ROOT", tmp_path, raising=False)
    inside = tmp_path / "skills" / "nfcore-scrnaseq-wrapper" / "tmp-test-output"
    with pytest.raises(SkillError) as exc:
        preflight.check_output_dir_available(inside, resume=False)
    assert exc.value.error_code == "OUTPUT_DIR_INSIDE_REPO"
    assert not inside.exists(), "guard must reject before creating the directory"


# ── H-17: detected versions must preserve their exact string (NXF_VER/conda) ───


def _patch_versions(monkeypatch, *, nextflow_text=None, java_text=None):
    import preflight

    monkeypatch.setattr(preflight.shutil, "which", lambda name: f"/usr/bin/{name}")

    def fake_output(args):
        if args and args[0] == "nextflow":
            return nextflow_text or ""
        if args and args[0] == "java":
            return java_text or ""
        return ""

    monkeypatch.setattr(preflight, "_command_output", fake_output)


def test_nextflow_version_preserves_zero_padded_month(monkeypatch):
    """Detected '25.04.0' must be stored verbatim, not reconstructed to '25.4.0'
    (which is not a real Nextflow release and would break NXF_VER / conda on replay)."""
    import preflight

    _patch_versions(
        monkeypatch, nextflow_text="  N E X T F L O W\n  version 25.04.0 build 5947"
    )
    assert preflight._check_nextflow()["version"] == "25.04.0"


def test_java_version_preserves_exact_string(monkeypatch):
    import preflight

    _patch_versions(monkeypatch, java_text='openjdk version "17.0.10" 2024-01-16')
    assert preflight._check_java()["version"] == "17.0.10"


# ── H-14: single source for the minimum Nextflow version ───────────────────────


def test_nextflow_min_display_matches_comparison_tuple():
    import preflight

    parsed = preflight._parse_version_tuple(schemas.NEXTFLOW_MIN_VERSION_DISPLAY)
    assert parsed == schemas.NEXTFLOW_MIN_VERSION, (
        "the human-readable Nextflow min version must parse to the comparison tuple"
    )


def test_preflight_messages_use_canonical_nextflow_version(monkeypatch):
    """User-facing fix messages must show the canonical '25.04.0', never the
    tuple rendering '25.4.0'."""
    import preflight
    from errors import SkillError

    assert (
        ".".join(map(str, schemas.NEXTFLOW_MIN_VERSION)) == "25.4.0"
    )  # the anti-pattern

    monkeypatch.setattr(preflight.shutil, "which", lambda name: f"/usr/bin/{name}")
    # unparseable version -> MISSING_NEXTFLOW fix message
    monkeypatch.setattr(preflight, "_command_output", lambda args: "no version here")
    try:
        preflight._check_nextflow()
        assert False, "expected SkillError"
    except SkillError as exc:
        assert schemas.NEXTFLOW_MIN_VERSION_DISPLAY in exc.fix
        assert "25.4.0" not in exc.fix

    # too-old version -> NEXTFLOW_VERSION_TOO_OLD fix message
    monkeypatch.setattr(
        preflight, "_command_output", lambda args: "version 24.10.0 build 1"
    )
    try:
        preflight._check_nextflow()
        assert False, "expected SkillError"
    except SkillError as exc:
        assert schemas.NEXTFLOW_MIN_VERSION_DISPLAY in exc.fix
        assert "25.4.0" not in exc.fix


def test_unparseable_version_is_version_failure_not_missing(monkeypatch):
    """A present binary whose version can't be parsed is a version-gate failure,
    not a MISSING_* binary (the code must not contradict the 'is installed'
    message). Aligned with nfcore-sarek-wrapper."""
    import preflight
    from errors import ErrorCode, SkillError

    monkeypatch.setattr(preflight.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(preflight, "_command_output", lambda args: "no version here")

    with pytest.raises(SkillError) as exc:
        preflight._check_java()
    assert exc.value.error_code == ErrorCode.JAVA_VERSION_TOO_OLD

    with pytest.raises(SkillError) as exc:
        preflight._check_nextflow()
    assert exc.value.error_code == ErrorCode.NEXTFLOW_VERSION_TOO_OLD


# ── H-7: --skip-downstream documents its precedence ────────────────────────────


def test_skip_downstream_help_documents_precedence():
    module = _load_skill_module()
    parser = module.build_parser()
    action = next(a for a in parser._actions if "--skip-downstream" in a.option_strings)
    assert re.search(r"even if|overrides|force", action.help, re.IGNORECASE), (
        "--skip-downstream help must document that it overrides --run-downstream"
    )


# ── F1: Smart-seq2 rejected; smartseq (Smart-seq3) gated to star/kallisto ──────


def _protocol_args(preset, protocol):
    return Namespace(demo=False, preset=preset, protocol=protocol)


@pytest.mark.parametrize(
    "preset", ["standard", "star", "kallisto", "cellranger", "cellrangermulti"]
)
def test_smartseq2_rejected_for_every_preset(preset):
    import preflight
    from errors import SkillError

    with pytest.raises(SkillError) as exc:
        preflight._check_protocol_compatibility(_protocol_args(preset, "smartseq2"))
    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert "Smart-seq2" in exc.value.message


@pytest.mark.parametrize("protocol", ["smartseq2", "Smart-seq2", "SMARTSEQ2"])
def test_smartseq2_rejected_case_and_separator_insensitive(protocol):
    import preflight
    from errors import SkillError

    with pytest.raises(SkillError):
        preflight._check_protocol_compatibility(_protocol_args("star", protocol))


def test_smartseq3_still_accepted_for_star_and_kallisto():
    import preflight

    for preset in ("star", "kallisto"):
        preflight._check_protocol_compatibility(
            _protocol_args(preset, "smartseq")
        )  # no raise


def test_smartseq3_rejected_for_standard():
    import preflight
    from errors import SkillError

    with pytest.raises(SkillError) as exc:
        preflight._check_protocol_compatibility(_protocol_args("standard", "smartseq"))
    assert "Smart-seq3" in exc.value.message


# ── F2: feature_type=cmo requires --cellranger-multi-barcodes ───────────────────


def _multi_args(**over):
    base = dict(
        preset="cellrangermulti",
        gex_frna_probe_set=None,
        cellranger_multi_barcodes=None,
        fb_reference=None,
        cellranger_vdj_index=None,
        skip_cellrangermulti_vdjref=False,
    )
    base.update(over)
    return Namespace(**base)


def test_cmo_without_barcodes_fails_preflight():
    import preflight
    from errors import SkillError

    with pytest.raises(SkillError) as exc:
        preflight._check_cellrangermulti_multiplexing_policy(
            _multi_args(), {"gex", "cmo"}
        )
    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["missing_field"] == "cellranger_multi_barcodes"
    assert exc.value.details["feature_type"] == "cmo"


def test_cmo_with_barcodes_passes(tmp_path):
    import preflight

    bc = tmp_path / "multi_barcodes.csv"
    bc.write_text("sample,cmo_ids\n", encoding="utf-8")
    preflight._check_cellrangermulti_multiplexing_policy(
        _multi_args(cellranger_multi_barcodes=str(bc)), {"gex", "cmo"}
    )  # no raise


def test_non_cmo_multi_does_not_require_barcodes():
    import preflight

    preflight._check_cellrangermulti_multiplexing_policy(
        _multi_args(), {"gex", "vdj"}
    )  # no raise


# ── F3: missing .h5ad is a required-output failure ──────────────────────────────


def test_missing_h5ad_is_required_output_failure():
    from outputs_parser import validate_expected_outputs

    parsed = {
        "top_level_entries": ["simpleaf", "multiqc", "pipeline_info"],
        "pipeline_info_dir": "/x/pipeline_info",
        "multiqc_report": "/x/multiqc/multiqc_report.html",
        "h5ad_candidates": [],
    }
    v = validate_expected_outputs(
        parsed, aligner="simpleaf", skip_multiqc=False, skip_fastqc=True
    )
    assert "simpleaf/mtx_conversions/*.h5ad" in v["missing_required"]
    assert v["missing_optional"] == []


def test_present_h5ad_yields_no_missing():
    from outputs_parser import validate_expected_outputs

    parsed = {
        "top_level_entries": ["star", "multiqc", "pipeline_info"],
        "pipeline_info_dir": "/x/pipeline_info",
        "multiqc_report": "/x/multiqc/multiqc_report.html",
        "h5ad_candidates": ["/x/star/mtx_conversions/combined_filtered_matrix.h5ad"],
    }
    v = validate_expected_outputs(
        parsed, aligner="star", skip_multiqc=False, skip_fastqc=True
    )
    assert v["missing_required"] == [] and v["missing_optional"] == []


# ── F9: downstream handoff records orchestrator path + checksum in provenance ──


def test_handoff_provenance_records_path_and_checksum(tmp_path):
    module = _load_skill_module()
    orch = tmp_path / "orch.py"
    orch.write_text("print('hi')\n", encoding="utf-8")
    import hashlib

    expected = hashlib.sha256(orch.read_bytes()).hexdigest()
    assert module._sha256_file_safe(orch) == expected
    assert module._sha256_file_safe(tmp_path / "missing.py") == ""
    record = {
        "orchestrator": orch.as_posix(),
        "orchestrator_sha256": expected,
        "status": "ok",
    }
    module._write_handoff_provenance(tmp_path, record)
    import json as _json

    written = _json.loads((tmp_path / "provenance" / "handoff.json").read_text())
    assert written["orchestrator_sha256"] == expected and written["status"] == "ok"


# ── F7: clawbio runner allowlist stays in sync with the wrapper CLI ─────────────


def test_clawbio_allowlist_is_subset_of_wrapper_flags():
    # The runner registry lives in clawbio.cli; clawbio.py is a thin shim that
    # delegates to it.
    import clawbio.cli as clawbio

    entry = clawbio.SKILLS["scrnaseq-pipeline"]
    allow = set(entry.get("allowed_extra_flags", set())) | set(
        entry.get("allowed_extra_flags_without_values", set())
    )

    module = _load_skill_module()
    # Consider every option string (short -c and long --config alike) so the
    # allowlist↔CLI contract covers short flags too, not only the long forms.
    wrapper_flags = {
        o
        for a in module.build_parser()._actions
        for o in a.option_strings
        if o.startswith("-")
    }
    drift = allow - wrapper_flags
    assert drift == set(), (
        f"clawbio allowlist flags not in wrapper CLI (drift): {sorted(drift)}"
    )
