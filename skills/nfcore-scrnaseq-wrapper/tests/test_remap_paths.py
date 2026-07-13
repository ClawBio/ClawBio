from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest

from remap_paths import cmd_output_dir_hint, cmd_remap, main, remap_csv, verify_paths


def test_cmd_output_dir_hint_is_accepted_noop(capsys):
    # CLI parity with the rnaseq/sarek wrappers: --output-dir is accepted, but the
    # scrnaseq bundle self-relocates (commands.sh anchors to its own location), so it
    # explains that no rewrite is needed and succeeds (finding #4).
    rc = cmd_output_dir_hint("/some/new/output")
    assert rc == 0
    out = capsys.readouterr().out.lower()
    assert "self" in out and "no" in out

_FASTQ_HEADER = "sample,fastq_1,fastq_2,expected_cells,seq_center\n"


def _write_samplesheet(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["sample", "fastq_1", "fastq_2", "expected_cells", "seq_center"],
        )
        writer.writeheader()
        writer.writerows(rows)


# ── find_samplesheet ──────────────────────────────────────────────────────────


# ── remap_csv ─────────────────────────────────────────────────────────────────


def test_remap_csv_replaces_matching_prefix(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(
        ss,
        [
            {
                "sample": "S1",
                "fastq_1": "/old/data/S1_R1.fastq.gz",
                "fastq_2": "/old/data/S1_R2.fastq.gz",
                "expected_cells": "",
                "seq_center": "",
            },
        ],
    )
    changes = remap_csv(ss, "/old/data", "/new/data", dry_run=False)
    assert len(changes) == 2
    assert changes[0] == (
        "fastq_1",
        "/old/data/S1_R1.fastq.gz",
        "/new/data/S1_R1.fastq.gz",
    )
    assert changes[1] == (
        "fastq_2",
        "/old/data/S1_R2.fastq.gz",
        "/new/data/S1_R2.fastq.gz",
    )
    rows = list(csv.DictReader(ss.read_text(encoding="utf-8").splitlines()))
    assert rows[0]["fastq_1"] == "/new/data/S1_R1.fastq.gz"
    assert rows[0]["fastq_2"] == "/new/data/S1_R2.fastq.gz"


def test_remap_csv_only_replaces_prefix_not_middle(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(
        ss,
        [
            {
                "sample": "S1",
                "fastq_1": "/data/old/S1_R1.fastq.gz",
                "fastq_2": "/data/old/S1_R2.fastq.gz",
                "expected_cells": "",
                "seq_center": "",
            },
        ],
    )
    changes = remap_csv(ss, "/old", "/new", dry_run=False)
    assert changes == [], "Should not replace /old in the middle of the path"


# ── verify_paths ──────────────────────────────────────────────────────────────


def test_verify_paths_returns_missing_paths(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(
        ss,
        [
            {
                "sample": "S1",
                "fastq_1": "/nonexistent/R1.fastq.gz",
                "fastq_2": "/nonexistent/R2.fastq.gz",
                "expected_cells": "",
                "seq_center": "",
            },
        ],
    )
    missing = verify_paths(ss)
    assert len(missing) == 2
    assert "/nonexistent/R1.fastq.gz" in missing
    assert "/nonexistent/R2.fastq.gz" in missing


# ── cmd_remap (integration) ───────────────────────────────────────────────────


def test_cmd_remap_succeeds_when_paths_are_remapped_to_existing_files(tmp_path):
    r1 = tmp_path / "fastqs" / "S1_R1.fastq.gz"
    r2 = tmp_path / "fastqs" / "S1_R2.fastq.gz"
    r1.parent.mkdir()
    r1.write_bytes(b"")
    r2.write_bytes(b"")
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(
        ss,
        [
            {
                "sample": "S1",
                "fastq_1": "/old/fastqs/S1_R1.fastq.gz",
                "fastq_2": "/old/fastqs/S1_R2.fastq.gz",
                "expected_cells": "",
                "seq_center": "",
            },
        ],
    )
    rc = cmd_remap(
        "/old/fastqs", str(tmp_path / "fastqs"), dry_run=False, bundle_dir=tmp_path
    )
    assert rc == 0
    rows = list(csv.DictReader(ss.read_text(encoding="utf-8").splitlines()))
    assert rows[0]["fastq_1"] == str(tmp_path / "fastqs" / "S1_R1.fastq.gz")


def test_cmd_remap_returns_nonzero_when_new_paths_missing(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(
        ss,
        [
            {
                "sample": "S1",
                "fastq_1": "/old/S1_R1.fastq.gz",
                "fastq_2": "/old/S1_R2.fastq.gz",
                "expected_cells": "",
                "seq_center": "",
            },
        ],
    )
    rc = cmd_remap("/old", "/nonexistent", dry_run=False, bundle_dir=tmp_path)
    assert rc != 0


# ── cmd_verify (integration) ──────────────────────────────────────────────────


# ── find_commands_sh ──────────────────────────────────────────────────────────


# ── params.yaml reference remapping ──────────────────────────────────────────


def _write_bundle(tmp_path):
    bundle = tmp_path / "reproducibility"
    bundle.mkdir(parents=True, exist_ok=True)
    return bundle


def test_remap_params_references_rewrites_local_paths(tmp_path):
    import remap_paths

    bundle = _write_bundle(tmp_path)
    params = bundle / "params.yaml"
    params.write_text(
        "outdir: upstream/results\n"
        "fasta: /old/refs/genome.fa\n"
        "gtf: /old/refs/genes.gtf\n"
        "star_index: s3://bucket/star\n",
        encoding="utf-8",
    )

    changes = remap_paths.remap_params_references(
        params, "/old/refs", "/new/refs", dry_run=False
    )

    text = params.read_text(encoding="utf-8")
    assert "fasta: /new/refs/genome.fa" in text
    assert "gtf: /new/refs/genes.gtf" in text
    assert "star_index: s3://bucket/star" in text
    assert {c[0] for c in changes} == {"fasta", "gtf"}


def test_verify_params_references_skips_uris_and_missing(tmp_path):
    import remap_paths

    bundle = _write_bundle(tmp_path)
    real = bundle / "genome.fa"
    real.write_text("x", encoding="utf-8")
    params = bundle / "params.yaml"
    params.write_text(
        f"fasta: {real}\n"
        "gtf: /does/not/exist.gtf\n"
        "transcript_fasta: https://example.org/tx.fa\n",
        encoding="utf-8",
    )

    missing = remap_paths.verify_params_references(params)
    assert missing == ["/does/not/exist.gtf"]


def test_cmd_remap_references_returns_nonzero_when_refs_missing(tmp_path):
    """A non-dry-run reference remap that leaves references missing must signal
    failure via a non-zero exit code, matching cmd_remap and cmd_verify (audit H-6)."""
    import remap_paths

    bundle = _write_bundle(tmp_path)
    params = bundle / "params.yaml"
    params.write_text("fasta: /old/refs/genome.fa\n", encoding="utf-8")

    rc = remap_paths.cmd_remap_references(
        "/old/refs", "/nonexistent/refs", dry_run=False, bundle_dir=bundle
    )

    assert rc != 0


def test_cmd_remap_references_succeeds_when_refs_exist(tmp_path):
    """When the remapped references resolve on disk, the exit code is 0 (audit H-6)."""
    import remap_paths

    bundle = _write_bundle(tmp_path)
    new_refs = tmp_path / "refs"
    new_refs.mkdir()
    (new_refs / "genome.fa").write_text("x", encoding="utf-8")
    params = bundle / "params.yaml"
    params.write_text("fasta: /old/refs/genome.fa\n", encoding="utf-8")

    rc = remap_paths.cmd_remap_references(
        "/old/refs", str(new_refs), dry_run=False, bundle_dir=bundle
    )

    assert rc == 0


def test_params_reference_keys_match_params_builder():
    import remap_paths
    import params_builder

    # Ordered equality (stronger than set): the bundle remapper must walk params.yaml
    # reference keys in the exact canonical order, not merely the same membership.
    assert tuple(remap_paths._PARAMS_REFERENCE_KEYS) == tuple(
        params_builder._REFERENCE_PATH_FIELDS
    )


def test_bundle_copy_of_real_script_is_canonical_and_functional(tmp_path):
    """The bundle copy of the real remap_paths.py exposes canonical reference keys
    and a working params.yaml regex."""
    import reporting
    from schemas import ALL_REFERENCE_PATH_FIELDS

    reporting._write_remap_script(tmp_path)
    copied = (tmp_path / "remap_paths.py").read_text(encoding="utf-8")

    namespace: dict = {"__file__": str(tmp_path / "remap_paths.py")}
    exec(compile(copied, "remap_paths_bundle", "exec"), namespace)
    assert namespace["_PARAMS_REFERENCE_KEYS"] == tuple(ALL_REFERENCE_PATH_FIELDS)
    assert namespace["_PARAMS_REF_RE"].search("fasta: /x/genome.fa\n") is not None


def test_bundle_copy_regenerates_keys_even_if_in_repo_literal_drifts(
    tmp_path, monkeypatch
):
    """By construction: even if the in-repo literal between the AUTO-GENERATED
    sentinels drifts, the bundle copy is regenerated from the canonical source so
    the standalone bundle can never go out of sync."""
    import reporting
    from schemas import ALL_REFERENCE_PATH_FIELDS

    drifted_source = (
        "from pathlib import Path\n"
        "import re\n"
        "_BUNDLE_DIR = Path(__file__).resolve().parent\n"
        "# >>> AUTO-GENERATED reference keys (do not edit by hand) >>>\n"
        '_PARAMS_REFERENCE_KEYS = ("WRONG_DRIFTED_KEY",)\n'
        "# <<< AUTO-GENERATED reference keys <<<\n"
        "_PARAMS_REF_RE = re.compile('|'.join(re.escape(k) for k in _PARAMS_REFERENCE_KEYS))\n"
    )
    fake_src = tmp_path / "src_remap_paths.py"
    fake_src.write_text(drifted_source, encoding="utf-8")
    monkeypatch.setattr(reporting, "_REMAP_SCRIPT_SRC", fake_src)

    out_dir = tmp_path / "bundle"
    reporting._write_remap_script(out_dir)
    copied = (out_dir / "remap_paths.py").read_text(encoding="utf-8")

    namespace: dict = {"__file__": str(out_dir / "remap_paths.py")}
    exec(compile(copied, "remap_paths_bundle", "exec"), namespace)
    assert "WRONG_DRIFTED_KEY" not in copied
    assert namespace["_PARAMS_REFERENCE_KEYS"] == tuple(ALL_REFERENCE_PATH_FIELDS)


# ── cmd_repair_bundle (crash-recovery parity with rnaseq/sarek) ────────────────


def _make_repairable_bundle(tmp_path: Path) -> Path:
    """Minimal scrnaseq bundle: reproducibility/ holds commands.sh + params.yaml +
    samplesheet; the output dir also has upstream/results and logs. Returns the
    reproducibility (bundle) directory. Callers add/remove the regenerable files."""
    output_dir = tmp_path / "run"
    bundle = output_dir / "reproducibility"
    bundle.mkdir(parents=True)
    (bundle / "commands.sh").write_text("#!/usr/bin/env bash\nnextflow run nf-core/scrnaseq\n", encoding="utf-8")
    (bundle / "params.yaml").write_text("outdir: upstream/results\naligner: simpleaf\n", encoding="utf-8")
    (bundle / "samplesheet.valid.csv").write_text(
        "sample,fastq_1,fastq_2\nS1,/data/a.fastq.gz,/data/b.fastq.gz\n", encoding="utf-8"
    )
    (output_dir / "upstream" / "results").mkdir(parents=True)
    (output_dir / "upstream" / "results" / "matrix.h5ad").write_bytes(b"h5ad-bytes")
    (output_dir / "logs").mkdir()
    (output_dir / "logs" / "nextflow_stdout.txt").write_text("nextflow log\n", encoding="utf-8")
    return bundle


def test_cmd_repair_bundle_noop_when_all_files_present(tmp_path):
    from remap_paths import cmd_repair_bundle

    bundle = _make_repairable_bundle(tmp_path)
    for name in ("manifest.json", "checksums.sha256", "environment.yml"):
        (bundle / name).write_text("placeholder\n", encoding="utf-8")
    assert cmd_repair_bundle(bundle_dir=bundle) == 0


def test_cmd_repair_bundle_regenerates_all_missing_files(tmp_path):
    from remap_paths import cmd_repair_bundle

    bundle = _make_repairable_bundle(tmp_path)
    result = cmd_repair_bundle(bundle_dir=bundle)
    assert result == 0
    for name in ("manifest.json", "checksums.sha256", "environment.yml"):
        assert (bundle / name).exists(), f"{name} must be regenerated"


def test_cmd_repair_bundle_returns_nonzero_when_commands_sh_missing(tmp_path):
    from remap_paths import cmd_repair_bundle

    bundle = _make_repairable_bundle(tmp_path)
    (bundle / "commands.sh").unlink()
    # Cannot repair without the original command record.
    assert cmd_repair_bundle(bundle_dir=bundle) == 1


def test_cmd_repair_bundle_stubs_marked_post_hoc(tmp_path):
    import json

    from remap_paths import cmd_repair_bundle

    bundle = _make_repairable_bundle(tmp_path)
    cmd_repair_bundle(bundle_dir=bundle)
    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    assert manifest.get("regenerated_post_hoc") is True
    assert "regenerated_post_hoc: true" in (bundle / "environment.yml").read_text(encoding="utf-8")
    # Environment stub must name the scrnaseq wrapper, not a sibling.
    assert "clawbio-nfcore-scrnaseq-wrapper" in (bundle / "environment.yml").read_text(encoding="utf-8")


def test_cmd_repair_bundle_is_idempotent(tmp_path):
    from remap_paths import cmd_repair_bundle

    bundle = _make_repairable_bundle(tmp_path)
    assert cmd_repair_bundle(bundle_dir=bundle) == 0
    assert cmd_repair_bundle(bundle_dir=bundle) == 0  # second run is a no-op


def test_regenerated_checksums_verify_against_files(tmp_path):
    """Safety guarantee: every entry in the regenerated checksums.sha256 must match
    the actual file digest, with labels resolving relative to output_dir — i.e.
    `sha256sum -c checksums.sha256` run from the output dir passes."""
    from remap_paths import _sha256_file, cmd_repair_bundle

    bundle = _make_repairable_bundle(tmp_path)
    cmd_repair_bundle(bundle_dir=bundle)
    output_dir = bundle.parent
    text = (bundle / "checksums.sha256").read_text(encoding="utf-8")
    assert text.strip(), "checksums.sha256 must not be empty"
    seen_labels = []
    for line in text.splitlines():
        digest, label = line.split("  ", 1)
        seen_labels.append(label)
        target = output_dir / label
        assert target.is_file(), f"checksum label does not resolve to a file: {label}"
        assert _sha256_file(target) == digest, f"stale/wrong digest for {label}"
    # The checksum file must never hash itself.
    assert "reproducibility/checksums.sha256" not in seen_labels
    # Must cover the regenerated params.yaml (the file remap mutates) and the output.
    assert "reproducibility/params.yaml" in seen_labels
    assert "upstream/results/matrix.h5ad" in seen_labels


def test_regenerated_checksums_are_lf_only(tmp_path):
    from remap_paths import cmd_repair_bundle

    bundle = _make_repairable_bundle(tmp_path)
    cmd_repair_bundle(bundle_dir=bundle)
    for name in ("checksums.sha256", "manifest.json", "environment.yml"):
        assert b"\r" not in (bundle / name).read_bytes(), f"{name} must be LF-only"


def _make_realistic_bundle(tmp_path: Path) -> Path:
    """A bundle laid out like a real post-run scrnaseq output dir, including files the
    original checksum manifest deliberately excludes (the reproducibility/*.json
    provenance tree, environment.yml/manifest.json, extra upstream outputs, extra
    logs). Returns the reproducibility (bundle) directory."""
    output_dir = tmp_path / "run"
    repro = output_dir / "reproducibility"
    repro.mkdir(parents=True)
    (repro / "commands.sh").write_text("#!/usr/bin/env bash\nnextflow run nf-core/scrnaseq\n", encoding="utf-8")
    (repro / "params.yaml").write_text("outdir: upstream/results\naligner: simpleaf\n", encoding="utf-8")
    (repro / "samplesheet.valid.csv").write_text(
        "sample,fastq_1,fastq_2\nS1,/data/a.fastq.gz,/data/b.fastq.gz\n", encoding="utf-8"
    )
    # Decoys under reproducibility/ that the original manifest never hashes.
    for name in ("runtime.json", "inputs.json", "invocation.json", "upstream.json",
                 "preflight.json", "outputs.json", "samplesheet.json"):
        (repro / name).write_text('{"decoy": true}\n', encoding="utf-8")

    results = output_dir / "upstream" / "results"
    mtx = results / "simpleaf" / "mtx_conversions"
    mtx.mkdir(parents=True)
    (mtx / "combined_filtered_matrix.h5ad").write_bytes(b"h5ad-combined")
    (mtx / "S1").mkdir()
    (mtx / "S1" / "S1_filtered_matrix.h5ad").write_bytes(b"h5ad-per-sample")
    multiqc = results / "multiqc"
    multiqc.mkdir(parents=True)
    (multiqc / "multiqc_report.html").write_text("<html>mqc</html>\n", encoding="utf-8")
    # Decoy upstream output that is not an h5ad candidate or the MultiQC report.
    pipeline_info = results / "pipeline_info"
    pipeline_info.mkdir()
    (pipeline_info / "execution_report.html").write_text("<html>exec</html>\n", encoding="utf-8")

    logs = output_dir / "logs"
    logs.mkdir()
    (logs / "stdout.txt").write_text("nextflow stdout\n", encoding="utf-8")
    (logs / "stderr.txt").write_text("nextflow stderr\n", encoding="utf-8")
    # Decoy log the original manifest never hashes.
    (logs / "nextflow.log").write_text("timestamped log line\n", encoding="utf-8")
    return repro


def test_repaired_checksums_match_freshly_generated_manifest(tmp_path):
    """Blocking-audit guard: the checksums.sha256 that --repair-bundle regenerates must
    be byte-identical to the one the wrapper itself would write, i.e. cover exactly the
    original allowlist and not the reproducibility/*.json tree, environment.yml/
    manifest.json, other upstream outputs, or extra logs."""
    from outputs_parser import find_h5ad_candidates, find_multiqc_report
    from provenance import write_reproducibility_checksums
    from remap_paths import _regenerate_checksums

    bundle = _make_realistic_bundle(tmp_path)
    output_dir = bundle.parent
    results = output_dir / "upstream" / "results"

    fresh_path = write_reproducibility_checksums(
        output_dir,
        normalized_samplesheet=bundle / "samplesheet.valid.csv",
        params_path=bundle / "params.yaml",
        preflight_result={},
        parsed_outputs={
            "h5ad_candidates": find_h5ad_candidates(results),
            "multiqc_report": find_multiqc_report(results),
        },
        execution_result={
            "stdout_path": str(output_dir / "logs" / "stdout.txt"),
            "stderr_path": str(output_dir / "logs" / "stderr.txt"),
        },
    )
    fresh_bytes = fresh_path.read_bytes()
    fresh_labels = {line.split("  ", 1)[1] for line in fresh_path.read_text().splitlines()}

    fresh_path.unlink()
    _regenerate_checksums(bundle)
    repaired_bytes = (bundle / "checksums.sha256").read_bytes()
    repaired_labels = {line.split("  ", 1)[1] for line in (bundle / "checksums.sha256").read_text().splitlines()}

    assert repaired_labels == fresh_labels
    assert repaired_bytes == fresh_bytes, "repaired manifest must be byte-identical to a fresh one"

    # The allowlist we expect, and the decoys that must never appear.
    assert "reproducibility/params.yaml" in repaired_labels
    assert "reproducibility/samplesheet.valid.csv" in repaired_labels
    assert "logs/stdout.txt" in repaired_labels
    assert "logs/stderr.txt" in repaired_labels
    assert any(lbl.endswith(".h5ad") for lbl in repaired_labels)
    assert "upstream/results/multiqc/multiqc_report.html" in repaired_labels
    for decoy in (
        "reproducibility/runtime.json",
        "reproducibility/inputs.json",
        "reproducibility/environment.yml",
        "reproducibility/manifest.json",
        "logs/nextflow.log",
        "upstream/results/pipeline_info/execution_report.html",
    ):
        assert decoy not in repaired_labels, f"{decoy} must be excluded from the manifest"


def test_repair_bundle_does_not_overwrite_existing_valid_checksums(tmp_path):
    """Secondary-audit guard: when only a stub (manifest.json/environment.yml) is
    missing, an existing valid checksums.sha256 is left untouched — no silent rewrite."""
    from remap_paths import cmd_repair_bundle

    bundle = _make_realistic_bundle(tmp_path)
    (bundle / "checksums.sha256").write_text("sentinel-original-manifest\n", encoding="utf-8")
    before = (bundle / "checksums.sha256").read_bytes()

    assert cmd_repair_bundle(bundle_dir=bundle) == 0  # regenerates the two missing stubs
    assert (bundle / "manifest.json").exists()
    assert (bundle / "environment.yml").exists()
    assert (bundle / "checksums.sha256").read_bytes() == before, (
        "an existing checksums.sha256 must not be overwritten when only stubs are missing"
    )


# ── main(): requested actions must compose, never be silently dropped ─────────


def _bundle_for_main(bundle: Path) -> Path:
    bundle.mkdir(parents=True, exist_ok=True)
    fastq = bundle / "S1.fastq.gz"
    fastq.write_bytes(b"")
    (bundle / "samplesheet.valid.csv").write_text(
        f"sample,fastq_1,fastq_2\nS1,{fastq},\n", encoding="utf-8"
    )
    return bundle


def test_main_output_dir_and_verify_compose(tmp_path, capsys):
    """`--output-dir X --verify` must print the relocation hint AND still verify.

    The dispatch used to be a chain of early returns with --output-dir first, so the
    hint short-circuited and --verify was silently discarded: the user asked for a
    readiness check and got none, with a zero exit code implying all was well. (The
    rnaseq/sarek bundles had the mirror image of this bug, dropping --output-dir.)
    """
    bundle = _bundle_for_main(tmp_path / "reproducibility")

    rc = main(["--output-dir", "/somewhere/new", "--verify"], bundle_dir=bundle)

    out = capsys.readouterr().out
    assert "self-relocating" in out, "the --output-dir hint should still be shown"
    assert "FASTQ paths" in out, "--verify was silently dropped by --output-dir"
    assert rc == 0


def test_main_rejects_half_a_prefix_pair(tmp_path):
    """--old without --new (or vice versa) must fail loudly, not be silently ignored."""
    bundle = _bundle_for_main(tmp_path / "reproducibility")
    for argv in (["--old", "/a"], ["--new", "/b"], ["--refs-old", "/a"], ["--refs-new", "/b"]):
        with pytest.raises(SystemExit) as exc:
            main(argv, bundle_dir=bundle)
        assert exc.value.code != 0
