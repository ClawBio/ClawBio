from __future__ import annotations

from argparse import Namespace
import json
from pathlib import Path
import sys

_SKILL_DIR = Path(__file__).resolve().parent.parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))


def _purge_foreign_bare_modules(*names: str) -> None:
    for name in names:
        module = sys.modules.get(name)
        module_file = Path(getattr(module, "__file__", "") or "")
        if module is not None and _SKILL_DIR not in module_file.parents and module_file != _SKILL_DIR / f"{name}.py":
            sys.modules.pop(name, None)


def _purge_local_modules(*names: str) -> None:
    for name in names:
        module = sys.modules.get(name)
        module_file = Path(getattr(module, "__file__", "") or "")
        if module is not None and (module_file == _SKILL_DIR / f"{name}.py" or _SKILL_DIR in module_file.parents):
            sys.modules.pop(name, None)


_purge_foreign_bare_modules("preflight", "provenance", "schemas")

from preflight import check_resume_params_checksum
import provenance as provenance_module
from provenance import (
    _reference_checksums,
    build_inputs_payload,
    build_runtime_payload,
    build_upstream_payload,
    write_provenance_bundle,
)

_purge_local_modules("preflight", "provenance", "schemas")
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))


def _args(**kwargs) -> Namespace:
    defaults = dict(
        aligner="star_salmon",
        profile="docker",
        resume=False,
        demo=False,
        check=False,
        pipeline_version="3.26.0",
        pseudo_aligner=None,
        prokaryotic=False,
    )
    defaults.update(kwargs)
    return Namespace(**defaults)


def _preflight(reference: Path | None = None) -> dict:
    refs = {"fasta": str(reference)} if reference else {}
    return {
        "java": {"version": "21.0.1", "path": "/usr/bin/java"},
        "nextflow": {"version": "25.04.3", "path": "/usr/bin/nextflow"},
        "references": refs,
        "warnings": ["low disk"],
    }


def test_build_upstream_payload_uses_rnaseq_pipeline_constant():
    payload = build_upstream_payload({"source_kind": "remote_tag", "source_ref": "3.26.0", "resolved_version": "3.26.0"})
    assert payload["pipeline"] == "nf-core/rnaseq"
    assert payload["version"] == "3.26.0"


def test_build_inputs_payload_includes_params_and_samplesheet_checksums(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    params = tmp_path / "params.yaml"
    samplesheet.write_text("sample,fastq_1,strandedness\nS1,a,auto\n", encoding="utf-8")
    params.write_text("aligner: star_salmon\n", encoding="utf-8")
    payload = build_inputs_payload(
        normalized_samplesheet=samplesheet,
        samplesheet_summary={"sample_count": 1, "fastq_paths": ["/fq/S1.fastq.gz"]},
        preflight_result=_preflight(),
        params_path=params,
    )
    assert payload["samplesheet_path"] == str(samplesheet)
    assert payload["samplesheet_checksum"]
    assert payload["params_checksum"]
    assert payload["samples_count"] == 1


def test_build_inputs_payload_preserves_remote_input_uris(tmp_path):
    samplesheet = tmp_path / "samplesheet.valid.csv"
    params = tmp_path / "params.yaml"
    samplesheet.write_text("sample,fastq_1,strandedness\nS1,s3://bucket/S1_R1.fastq.gz,auto\n", encoding="utf-8")
    params.write_text("aligner: star_salmon\n", encoding="utf-8")
    fastq_uri = "s3://bucket/S1_R1.fastq.gz"
    bam_uri = "s3://bucket/S1.genome.bam"
    payload = build_inputs_payload(
        normalized_samplesheet=samplesheet,
        samplesheet_summary={"sample_count": 1, "fastq_paths": [fastq_uri], "bam_paths": [bam_uri]},
        preflight_result=_preflight(),
        params_path=params,
    )
    assert payload["fastq_paths"] == [fastq_uri]
    assert payload["fastq_checksums"][fastq_uri] == "<remote-uri>"
    assert payload["bam_paths"] == [bam_uri]
    assert payload["bam_checksums"][bam_uri] == "<remote-uri>"


def test_build_inputs_payload_includes_reference_checksums(tmp_path):
    ref = tmp_path / "genome.fa"
    ref.write_text(">chr1\nAC\n", encoding="utf-8")
    samplesheet = tmp_path / "samplesheet.csv"
    params = tmp_path / "params.yaml"
    samplesheet.write_text("sample,fastq_1,strandedness\n", encoding="utf-8")
    params.write_text("{}", encoding="utf-8")
    payload = build_inputs_payload(
        normalized_samplesheet=samplesheet,
        samplesheet_summary={"sample_count": 0, "fastq_paths": []},
        preflight_result=_preflight(ref),
        params_path=params,
    )
    assert payload["reference_paths"]["fasta"] == str(ref)
    assert payload["reference_checksums"]["fasta"]


def test_build_runtime_payload_has_required_contract_keys(tmp_path):
    payload = build_runtime_payload(
        tmp_path,
        args=_args(resume=True),
        preflight_result=_preflight(),
        command_str="nextflow run nf-core/rnaseq",
        timestamp="2026-05-09T00:00:00+00:00",
    )
    for key in ("python_version", "java_version", "nextflow_version", "os", "arch", "host", "started_at", "finished_at", "duration_seconds", "cpu_count", "ram_total_gb"):
        assert key in payload
    assert payload["resume_used"] is True


def test_write_provenance_bundle_writes_seven_json_files(tmp_path):
    params = tmp_path / "reproducibility" / "params.yaml"
    samplesheet = tmp_path / "reproducibility" / "samplesheet.valid.csv"
    stdout = tmp_path / "logs" / "stdout.txt"
    stderr = tmp_path / "logs" / "stderr.txt"
    for path in (params, samplesheet, stdout, stderr):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x", encoding="utf-8")
    provenance_dir, manifest = write_provenance_bundle(
        tmp_path,
        args=_args(),
        pipeline_source={"source_kind": "local_checkout", "source_ref": "/repo/rnaseq", "resolved_version": "abc123", "git_commit": "abc123", "git_branch": "main", "git_dirty": False},
        preflight_result=_preflight(),
        params_path=params,
        params_payload={"aligner": "star_salmon"},
        normalized_samplesheet=samplesheet,
        samplesheet_summary={"sample_count": 2, "fastq_paths": []},
        parsed_outputs={"preferred_counts_tsv": "/counts.tsv", "handoff_available": True},
        execution_result={"stdout_path": str(stdout), "stderr_path": str(stderr)},
        command_str="nextflow run nf-core/rnaseq",
    )
    assert manifest.exists()
    assert sorted(p.name for p in provenance_dir.glob("*.json")) == [
        "inputs.json",
        "invocation.json",
        "outputs.json",
        "preflight.json",
        "runtime.json",
        "skill.json",
        "upstream.json",
    ]


def test_manifest_contains_checksums_and_pipeline_source(tmp_path):
    params = tmp_path / "reproducibility" / "params.yaml"
    samplesheet = tmp_path / "reproducibility" / "samplesheet.valid.csv"
    stdout = tmp_path / "logs" / "stdout.txt"
    stderr = tmp_path / "logs" / "stderr.txt"
    for path in (params, samplesheet, stdout, stderr):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x", encoding="utf-8")
    params_payload = {"aligner": "star_salmon", "pseudo_aligner": "salmon"}
    _, manifest_path = write_provenance_bundle(
        tmp_path,
        args=_args(pseudo_aligner="salmon", prokaryotic=True),
        pipeline_source={"source_kind": "remote_tag", "source_ref": "3.26.0", "resolved_version": "3.26.0"},
        preflight_result=_preflight(),
        params_path=params,
        params_payload=params_payload,
        normalized_samplesheet=samplesheet,
        samplesheet_summary={"sample_count": 1, "fastq_paths": []},
        parsed_outputs={},
        execution_result={"stdout_path": str(stdout), "stderr_path": str(stderr)},
        command_str="cmd",
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["skill_name"] == "nfcore-rnaseq-wrapper"
    assert manifest["params_checksum"] == provenance_module.params_payload_checksum(params_payload)
    assert manifest["params_file_sha256"]
    assert manifest["samplesheet_checksum"]
    assert manifest["pipeline_source"]["pipeline"] == "nf-core/rnaseq"
    assert manifest["pseudo_aligner"] == "salmon"
    assert manifest["prokaryotic"] is True
    assert manifest["checksums_file"].endswith("checksums.sha256")
    # `demo` is part of the resume compatibility contract: --demo swaps in the upstream
    # `test` profile's own samplesheet and bundled references, so a resume must not
    # cross the demo/real boundary. Recording it lets preflight reject that precisely.
    assert manifest["demo"] is False


def test_manifest_written_by_provenance_validates_for_resume(tmp_path):
    params = tmp_path / "reproducibility" / "params.yaml"
    samplesheet = tmp_path / "reproducibility" / "samplesheet.valid.csv"
    stdout = tmp_path / "logs" / "stdout.txt"
    stderr = tmp_path / "logs" / "stderr.txt"
    for path in (params, samplesheet, stdout, stderr):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x", encoding="utf-8")
    params_payload = {"outdir": str(tmp_path / "upstream" / "results"), "aligner": "star_salmon", "pseudo_aligner": "salmon"}
    write_provenance_bundle(
        tmp_path,
        args=_args(pseudo_aligner="salmon"),
        pipeline_source={"source_kind": "remote_tag", "source_ref": "3.26.0", "resolved_version": "3.26.0"},
        preflight_result=_preflight(),
        params_path=params,
        params_payload=params_payload,
        normalized_samplesheet=samplesheet,
        samplesheet_summary={"sample_count": 1, "fastq_paths": []},
        parsed_outputs={},
        execution_result={"stdout_path": str(stdout), "stderr_path": str(stderr)},
        command_str="cmd",
    )
    assert check_resume_params_checksum(params_payload, tmp_path) == provenance_module.params_payload_checksum(params_payload)


def test_runtime_payload_finished_at_reflects_duration():
    """finished_at must be started_at + duration_seconds, not a copy of the same timestamp."""
    from datetime import datetime

    duration = 3600.0
    timestamp = "2026-05-16T10:00:00+00:00"
    payload = build_runtime_payload(
        Path("/tmp"),
        args=_args(),
        preflight_result=_preflight(),
        command_str="nextflow run",
        timestamp=timestamp,
        duration_seconds=duration,
    )
    started = datetime.fromisoformat(payload["started_at"])
    finished = datetime.fromisoformat(payload["finished_at"])
    assert finished != started, "finished_at must differ from started_at when duration > 0"
    assert abs((finished - started).total_seconds() - duration) < 1e-3


def test_manifest_environment_yml_mode_exact_when_versions_known(tmp_path):
    """environment_yml_mode must be 'exact_versions' when java and nextflow versions are known."""
    params = tmp_path / "reproducibility" / "params.yaml"
    samplesheet = tmp_path / "reproducibility" / "samplesheet.valid.csv"
    stdout = tmp_path / "logs" / "stdout.txt"
    stderr = tmp_path / "logs" / "stderr.txt"
    for path in (params, samplesheet, stdout, stderr):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x", encoding="utf-8")
    _, manifest_path = write_provenance_bundle(
        tmp_path,
        args=_args(),
        pipeline_source={"source_kind": "remote_tag", "source_ref": "3.26.0", "resolved_version": "3.26.0"},
        preflight_result=_preflight(),  # _preflight() includes java + nextflow versions
        params_path=params,
        params_payload={},
        normalized_samplesheet=samplesheet,
        samplesheet_summary={"sample_count": 0, "fastq_paths": []},
        parsed_outputs={},
        execution_result={"stdout_path": str(stdout), "stderr_path": str(stderr)},
        command_str="cmd",
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["environment_yml_mode"] == "exact_versions", (
        "Should be 'exact_versions' when java and nextflow versions are both present in preflight"
    )


def test_reference_checksums_skips_empty_string_values(tmp_path):
    # Regression test: _reference_checksums must NOT treat "" as Path(".") and
    # attempt to checksum the CWD. Empty-string reference fields are "not provided"
    # by the user and must be omitted from the output dict entirely.
    real_ref = tmp_path / "genome.fa"
    real_ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    refs = {
        "genome": "",
        "fasta": str(real_ref),
        "gtf": "",
        "gff": "",
        "transcript_fasta": "",
        "star_index": "",
        "salmon_index": "",
    }
    result = _reference_checksums(refs)
    # Only the real file should appear in the output
    assert "fasta" in result
    assert result["fasta"] and result["fasta"] != "<missing>"
    # Empty-string keys must NOT be in the output (they are not user-provided paths)
    for key in ("genome", "gtf", "gff", "transcript_fasta", "star_index", "salmon_index"):
        assert key not in result, f"Empty-string key '{key}' should be omitted, not checksummed as CWD"


def test_write_provenance_bundle_with_real_run_references(tmp_path):
    # Regression: a custom (non-demo) run passes preflight_result["references"] that
    # includes all reference fields, most as empty strings. write_provenance_bundle
    # must complete and write all 7 JSON files without hanging or crashing.
    params = tmp_path / "reproducibility" / "params.yaml"
    samplesheet = tmp_path / "reproducibility" / "samplesheet.valid.csv"
    stdout = tmp_path / "logs" / "stdout.txt"
    stderr = tmp_path / "logs" / "stderr.txt"
    fasta = tmp_path / "genome.fa"
    gtf = tmp_path / "genome.gtf"
    for path in (params, samplesheet, stdout, stderr, fasta, gtf):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x", encoding="utf-8")
    # Simulate what _collect_reference_values returns for a --fasta + --gtf run
    real_run_preflight = {
        "java": {"version": "21.0.1", "path": "/usr/bin/java"},
        "nextflow": {"version": "25.04.3", "path": "/usr/bin/nextflow"},
        "references": {
            "genome": "",
            "fasta": str(fasta),
            "gtf": str(gtf),
            "gff": "",
            "transcript_fasta": "",
            "additional_fasta": "",
            "gene_bed": "",
            "splicesites": "",
            "star_index": "",
            "rsem_index": "",
            "hisat2_index": "",
            "bowtie2_index": "",
            "salmon_index": "",
            "kallisto_index": "",
        },
        "warnings": [],
    }
    provenance_dir, manifest = write_provenance_bundle(
        tmp_path,
        args=_args(),
        pipeline_source={"source_kind": "remote_repo", "source_ref": "nf-core/rnaseq", "resolved_version": "3.26.0"},
        preflight_result=real_run_preflight,
        params_path=params,
        params_payload={"aligner": "star_salmon"},
        normalized_samplesheet=samplesheet,
        samplesheet_summary={"sample_count": 2, "fastq_paths": ["/reads/R1.fastq.gz"]},
        parsed_outputs={"preferred_counts_tsv": "/counts.tsv"},
        execution_result={"stdout_path": str(stdout), "stderr_path": str(stderr)},
        command_str="nextflow run nf-core/rnaseq",
    )
    assert manifest.exists()
    assert sorted(p.name for p in provenance_dir.glob("*.json")) == [
        "inputs.json", "invocation.json", "outputs.json",
        "preflight.json", "runtime.json", "skill.json", "upstream.json",
    ]
    inputs = json.loads((provenance_dir / "inputs.json").read_text(encoding="utf-8"))
    # Empty-string refs must not appear in reference_checksums
    for key in ("genome", "gff", "transcript_fasta", "additional_fasta",
                "gene_bed", "splicesites", "star_index", "rsem_index",
                "hisat2_index", "bowtie2_index", "salmon_index", "kallisto_index"):
        assert key not in inputs["reference_checksums"], (
            f"Empty-string reference '{key}' must not appear in reference_checksums"
        )


def test_checksum_input_file_returns_remote_uri_sentinel_for_s3():
    """_checksum_input_file must return '<remote-uri>' for s3:// paths rather than
    '<missing>', which would be a false provenance failure for cloud-based FASTQ inputs."""
    result = provenance_module._checksum_input_file("s3://my-bucket/sample.fastq.gz")
    assert result == "<remote-uri>", (
        f"s3:// path must yield '<remote-uri>' sentinel, got {result!r}"
    )


def test_reference_checksums_returns_bash_variable_sentinel():
    """_reference_checksums must return '<bash-variable>' for $-prefixed reference paths
    instead of '<missing>'.  Commands.sh may embed paths like $REFS/hg38.fa that are
    environment-variable templates, not literal filesystem paths."""
    result = provenance_module._reference_checksums({
        "fasta": "$REFS/hg38.fa",
        "gtf": "${REFS}/hg38.gtf",
        "star_index": "/real/existing/path_that_may_not_exist",
    })
    assert result.get("fasta") == "<bash-variable>", (
        f"$-prefixed path must yield '<bash-variable>' sentinel, got fasta={result.get('fasta')!r}"
    )
    assert result.get("gtf") == "<bash-variable>", (
        f"${{}}-prefixed path must yield '<bash-variable>' sentinel, got gtf={result.get('gtf')!r}"
    )


def test_dir_checksum_checksums_small_directories(tmp_path):
    """_dir_checksum must return a SHA-256 hex string for directories under 5 GB."""
    small_dir = tmp_path / "ref"
    small_dir.mkdir()
    (small_dir / "genome.fa").write_text(">chr1\nACGT\n", encoding="utf-8")

    result = provenance_module._dir_checksum(small_dir)

    assert len(result) == 64, f"Expected 64-char hex SHA256, got {len(result)}: {result!r}"
    assert all(c in "0123456789abcdef" for c in result)


def test_dir_checksum_skips_large_directories(tmp_path, monkeypatch):
    """_dir_checksum must return '<skipped:X.Xgb>' for directories above 5 GB."""
    large_dir = tmp_path / "star_index"
    large_dir.mkdir()
    real_file = large_dir / "Genome"
    real_file.write_text("x", encoding="utf-8")

    # Make the file appear to be 6 GB by patching stat().st_size
    _BYTES_6GB = 6 * 1024 ** 3
    original_stat = real_file.stat()

    class FakeStat:
        def __init__(self):
            self.st_size = _BYTES_6GB
            # forward other attrs
            for attr in ("st_mode", "st_ino", "st_dev", "st_nlink",
                         "st_uid", "st_gid", "st_atime", "st_mtime", "st_ctime"):
                setattr(self, attr, getattr(original_stat, attr))

    # Patch stat() ONLY for the target "Genome" file (so it reports 6 GB). Every
    # other stat — including the directory stats pathlib.Path.rglob() performs while
    # walking the tree inside _dir_checksum — must use the real implementation. A
    # globally-mocked Path.stat returns a regular-file st_mode for directories too,
    # which breaks rglob's is-directory checks and can yield an empty traversal (and
    # thus a wrong hash) on some Python versions.
    _real_stat = Path.stat

    def _fake_stat(self, **kw):
        if self.name == "Genome":
            return FakeStat()
        return _real_stat(self, **kw)

    monkeypatch.setattr("pathlib.Path.stat", _fake_stat)

    result = provenance_module._dir_checksum(large_dir)

    assert result.startswith("<skipped:"), f"Expected skip sentinel, got: {result!r}"
    assert result.endswith("gb>"), f"Expected sentinel to end with 'gb>', got: {result!r}"
