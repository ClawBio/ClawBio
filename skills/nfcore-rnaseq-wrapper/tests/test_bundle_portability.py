"""Cross-OS portability guard for the reproducibility bundle.

Generates a full bundle via the real writers and asserts every text artifact is
LF-only (no CR), so a bundle is byte-identical and bash/checksum-safe regardless
of the OS that produced it. This is the regression guard for routing every
rnaseq bundle writer through the shared ``write_text_lf`` choke-point, matching
the nfcore-scrnaseq and nfcore-sarek wrappers.
"""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

_SKILL_DIR = Path(__file__).resolve().parent.parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

import params_builder
import provenance
import reporting
import remap_paths


# Text artifacts that must be LF-only on every OS.
_TEXT_ARTIFACTS = (
    "commands.sh",
    "params.yaml",
    "samplesheet.valid.csv",
    "checksums.sha256",
    "environment.yml",
    "manifest.json",
    "remap_paths.py",
    "runtime.json",
    "inputs.json",
    "outputs.json",
    "invocation.json",
    "skill.json",
    "preflight.json",
    "upstream.json",
)


def _args(out: Path, normalized: Path, fa: Path, gtf: Path) -> SimpleNamespace:
    return SimpleNamespace(
        aligner="star_salmon",
        profile="docker",
        demo=False,
        resume=False,
        check=False,
        input=str(normalized),
        output=str(out),
        pipeline_version="3.26.0",
        pipeline_local=None,
        pseudo_aligner=None,
        prokaryotic=False,
        arm=False,
        trimmer="trimgalore",
        genome=None,
        fasta=str(fa),
        gtf=str(gtf),
        nextflow_config=[],
        timeout_hours=12,
        run_downstream=False,
        skip_downstream=False,
        deseq2_vst=None,
    )


def _build_full_bundle(tmp_path: Path) -> Path:
    out = tmp_path / "run"
    (out / "reproducibility").mkdir(parents=True)
    (out / "upstream" / "results").mkdir(parents=True)
    (out / "logs").mkdir(parents=True)

    refs = tmp_path / "refs"
    refs.mkdir()
    fa = refs / "genome.fa"
    fa.write_text(">chr1\nACGT\n", encoding="utf-8")
    gtf = refs / "genes.gtf"
    gtf.write_text("# gtf\n", encoding="utf-8")

    normalized = out / "reproducibility" / "samplesheet.valid.csv"
    normalized.write_text(
        "sample,fastq_1,fastq_2,strandedness\nS1,/data/S1_R1.fastq.gz,/data/S1_R2.fastq.gz,auto\n",
        encoding="utf-8",
    )

    args = _args(out, normalized, fa, gtf)
    params_path, params = params_builder.build_params_file(
        args, normalized_samplesheet=normalized, output_dir=out
    )

    pipeline_source = {
        "source_kind": "remote_tag",
        "source_ref": "3.26.0",
        "resolved_version": "3.26.0",
    }
    preflight = {
        "java": {"version": "21.0.1", "path": "/usr/bin/java"},
        "nextflow": {"version": "25.04.3", "path": "/usr/bin/nextflow"},
        "references": {"fasta": str(fa), "gtf": str(gtf)},
        "samplesheet": {"strandedness_counts": {"auto": 1}},
        "warnings": [],
    }
    parsed_outputs = {
        "aligner_effective": "star_salmon",
        "preferred_counts_tsv": "",
        "rds_file": "",
        "tpm_tsv": "",
        "multiqc_report": "",
        "pipeline_info_dir": "",
        "samples_detected": 1,
        "handoff_available": False,
        "skip_quantification_merge": False,
        "hisat2_no_quant": False,
    }
    command_str = f"cd {out.as_posix()} && nextflow run nf-core/rnaseq -r 3.26.0"

    reporting.write_repro_commands(out, args=args)
    reporting.write_report(
        out,
        args=args,
        pipeline_source=pipeline_source,
        preflight_result=preflight,
        parsed_outputs=parsed_outputs,
        command_str=command_str,
    )
    reporting.write_result(
        out,
        args=args,
        pipeline_source=pipeline_source,
        parsed_outputs=parsed_outputs,
        command_str=command_str,
    )

    stdout = out / "logs" / "stdout.txt"
    stdout.write_text("ok\n", encoding="utf-8")
    stderr = out / "logs" / "stderr.txt"
    stderr.write_text("", encoding="utf-8")
    provenance.write_provenance_bundle(
        out,
        args=args,
        pipeline_source=pipeline_source,
        preflight_result=preflight,
        params_path=params_path,
        params_payload=params,
        normalized_samplesheet=normalized,
        samplesheet_summary={"sample_count": 1, "fastq_paths": []},
        parsed_outputs=parsed_outputs,
        execution_result={"stdout_path": str(stdout), "stderr_path": str(stderr)},
        command_str=command_str,
    )
    return out


def test_remap_write_text_lf_helper_normalizes_cr(tmp_path):
    """The standalone bundle copy cannot import clawbio, so it embeds the runtime's
    LF writer. It must normalise CRLF/CR to LF and write bytes (no OS newline
    translation)."""
    target = tmp_path / "out.txt"
    remap_paths._nfcore_write_text_lf(target, "a\r\nb\rc\n")
    assert target.read_bytes() == b"a\nb\nc\n"


def test_all_bundle_text_artifacts_are_lf_only(tmp_path):
    out = _build_full_bundle(tmp_path)
    repro = out / "reproducibility"
    offenders = []
    for name in _TEXT_ARTIFACTS:
        p = repro / name
        if p.exists() and b"\r" in p.read_bytes():
            offenders.append(name)
    for name in ("report.md", "result.json"):
        p = out / name
        if p.exists() and b"\r" in p.read_bytes():
            offenders.append(name)
    assert offenders == [], f"CRLF/CR found in bundle artifacts: {offenders}"


def test_commands_sh_is_pure_lf(tmp_path):
    out = _build_full_bundle(tmp_path)
    # A CR anywhere in a bash script can corrupt variables / line continuations.
    assert b"\r" not in (out / "reproducibility" / "commands.sh").read_bytes()


def test_params_yaml_is_pure_lf(tmp_path):
    out = _build_full_bundle(tmp_path)
    assert b"\r" not in (out / "reproducibility" / "params.yaml").read_bytes()
