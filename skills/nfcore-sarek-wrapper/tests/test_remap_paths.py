"""Tests for the standalone bundle-side remap_paths.py helper."""
from __future__ import annotations

from pathlib import Path

import remap_paths


def test_regenerate_checksums_mirrors_provenance_exclusions(tmp_path: Path):
    """The stdlib regenerator must hash the same file set as provenance.py.

    provenance.py excludes work/, .nextflow/, reproducibility/ and logs/ (and
    .log files). The standalone regenerator copied into each bundle must agree,
    otherwise ``sha256sum -c`` after a remap would not match the original.
    """
    output_dir = tmp_path / "out"
    bundle = output_dir / "reproducibility"
    (output_dir / "upstream" / "results").mkdir(parents=True)
    (output_dir / "upstream" / "results" / "var.vcf.gz").write_bytes(b"vcf")
    bundle.mkdir(parents=True)
    (bundle / "params.yaml").write_text("x", encoding="utf-8")
    (bundle / "logs").mkdir()
    (bundle / "logs" / "stdout.txt").write_text("log", encoding="utf-8")
    (output_dir / "work").mkdir()
    (output_dir / "work" / "junk").write_text("junk", encoding="utf-8")
    (output_dir / ".nextflow").mkdir()
    (output_dir / ".nextflow" / "history").write_text("h", encoding="utf-8")

    remap_paths._regenerate_checksums(bundle)
    text = (bundle / "checksums.sha256").read_text()

    assert "upstream/results/var.vcf.gz" in text
    assert "work/" not in text
    assert ".nextflow" not in text
    # The reproducibility/ bundle (incl. its logs/) is excluded, mirroring provenance.py.
    assert "reproducibility/" not in text
    assert "logs/" not in text
    # The checksum file never references itself.
    assert "checksums.sha256" not in text


def test_regenerate_checksums_excludes_nested_upstream_work(tmp_path: Path):
    """The default Nextflow work dir is <output>/upstream/work — nested one level
    below the output root. The standalone regenerator must exclude it there too,
    otherwise a remapped bundle's manifest would list ephemeral scratch files."""
    output_dir = tmp_path / "out"
    bundle = output_dir / "reproducibility"
    (output_dir / "upstream" / "results").mkdir(parents=True)
    (output_dir / "upstream" / "results" / "var.vcf.gz").write_bytes(b"vcf")
    bundle.mkdir(parents=True)
    (bundle / "params.yaml").write_text("x", encoding="utf-8")
    (output_dir / "upstream" / "work" / "ab" / "cd12").mkdir(parents=True)
    (output_dir / "upstream" / "work" / "ab" / "cd12" / "task.bam").write_text(
        "ephemeral", encoding="utf-8"
    )

    remap_paths._regenerate_checksums(bundle)
    text = (bundle / "checksums.sha256").read_text()

    assert "upstream/results/var.vcf.gz" in text
    assert "upstream/work/" not in text
    assert "task.bam" not in text


# ---------------------------------------------------------------------------
# params.yaml reference remapping (full cross-OS portability)
# ---------------------------------------------------------------------------


def _write_params(tmp_path: Path, body: str) -> Path:
    repro = tmp_path / "reproducibility"
    repro.mkdir(parents=True, exist_ok=True)
    p = repro / "params.yaml"
    p.write_text(body, encoding="utf-8")
    return p


def test_remap_params_references_rewrites_reference_paths(tmp_path: Path):
    p = _write_params(
        tmp_path,
        "# header comment with /orig/refs that must NOT change\n"
        "input: reproducibility/samplesheet.valid.csv\n"
        "fasta: /orig/refs/genome.fa\n"
        "dbsnp: /orig/refs/dbsnp.vcf.gz\n"
        "pon: false\n"
        "vep_cache: s3://annotation-cache/vep_cache/\n"
        "bcftools_columns: CHROM,POS,REF,ALT\n",
    )
    changes = remap_paths.remap_params_references(p, "/orig/refs", "/new/refs", dry_run=False)
    keys = {k for k, _old, _new in changes}
    assert keys == {"fasta", "dbsnp"}
    text = p.read_text()
    assert "fasta: /new/refs/genome.fa" in text
    assert "dbsnp: /new/refs/dbsnp.vcf.gz" in text
    assert "pon: false" in text                      # disable sentinel untouched
    assert "vep_cache: s3://annotation-cache/vep_cache/" in text  # URI untouched
    assert "bcftools_columns: CHROM,POS,REF,ALT" in text          # non-path untouched
    assert "# header comment with /orig/refs" in text            # comments untouched
    assert p.with_suffix(".yaml.bak").exists()


def test_remap_params_references_dry_run_does_not_write(tmp_path: Path):
    p = _write_params(tmp_path, "fasta: /orig/g.fa\n")
    before = p.read_text()
    changes = remap_paths.remap_params_references(p, "/orig", "/new", dry_run=True)
    assert len(changes) == 1
    assert p.read_text() == before
    assert not p.with_suffix(".yaml.bak").exists()


def test_remap_params_references_preserves_quotes(tmp_path: Path):
    p = _write_params(tmp_path, "fasta: '/orig/my genome.fa'\n")
    remap_paths.remap_params_references(p, "/orig", "/new", dry_run=False)
    assert "fasta: '/new/my genome.fa'" in p.read_text()


def test_verify_params_references_reports_missing(tmp_path: Path):
    p = _write_params(
        tmp_path,
        "fasta: /no/such/genome.fa\n"
        "vep_cache: s3://x/\n"
        "pon: false\n",
    )
    missing = remap_paths.verify_params_references(p)
    assert missing == ["/no/such/genome.fa"]


def test_params_reference_keys_match_provenance():
    # The standalone key list must stay in sync with the authoritative set so a
    # newly-added reference param is not silently skipped during remap.
    from provenance import REFERENCE_PATH_PARAMS

    assert set(remap_paths._PARAMS_REFERENCE_KEYS) == set(REFERENCE_PATH_PARAMS)


def test_cmd_remap_references_rewrites_params_yaml(tmp_path: Path):
    bundle = tmp_path / "reproducibility"
    bundle.mkdir(parents=True)
    newrefs = tmp_path / "refs"
    newrefs.mkdir()
    (newrefs / "genome.fa").write_text("x", encoding="utf-8")  # so post-remap verify passes
    (bundle / "params.yaml").write_text("fasta: /orig/genome.fa\n", encoding="utf-8")
    (bundle / "commands.sh").write_text("#!/usr/bin/env bash\nnextflow run nf-core/sarek\n", encoding="utf-8")
    rc = remap_paths.cmd_remap_references("/orig", str(newrefs), dry_run=False, bundle_dir=bundle)
    assert rc == 0
    assert f"fasta: {newrefs.as_posix()}/genome.fa" in (bundle / "params.yaml").read_text()


def test_replay_hint_is_self_contained_no_required_clawbio_repo():
    """commands.sh replays with a bare ``bash commands.sh`` — no required env var
    (enforced by test_reporting.test_commands_sh_is_directly_runnable_and_pins_engine).

    remap_paths.py's own replay hints must match, or the two halves of the same
    bundle contradict each other: telling the user ``CLAWBIO_REPO=... bash
    commands.sh`` is required is stale and wrong. ``CLAWBIO_REPO`` belongs only to
    the optional ``clawbio.py run`` re-validation path, which commands.sh documents
    itself — remap_paths.py must not couple it to ``bash commands.sh``.
    """
    src = (Path(__file__).resolve().parent.parent / "remap_paths.py").read_text(
        encoding="utf-8"
    )
    offenders = [
        f"{i}: {ln.strip()}"
        for i, ln in enumerate(src.splitlines(), 1)
        if "CLAWBIO_REPO" in ln and "commands.sh" in ln and "bash" in ln
    ]
    assert not offenders, (
        "remap_paths.py couples CLAWBIO_REPO to `bash commands.sh`, but commands.sh is "
        "self-contained. Drop the CLAWBIO_REPO prefix from the replay hint:\n"
        + "\n".join(offenders)
    )
    assert "CLAWBIO_REPO is always required" not in src
    assert "CLAWBIO_REPO must always be set" not in src
    # A self-contained replay hint must still be present.
    assert "bash commands.sh" in src


def test_bundle_remap_roundtrip_simulates_another_machine(tmp_path: Path):
    """End-to-end cross-machine portability: move the FASTQs to a NEW prefix that
    contains a space and a non-ASCII char (as a real replay machine well might),
    remap the bundle samplesheet, and confirm every path resolves and the file
    stays LF-only. Mirrors nfcore-scrnaseq-wrapper's bundle roundtrip so both skills
    prove the same guarantee, and exercises spaces/unicode in relocated paths."""
    data = tmp_path / "data"
    data.mkdir()
    r1 = data / "s1_R1.fastq.gz"
    r1.write_bytes(b"\x1f\x8b")
    r2 = data / "s1_R2.fastq.gz"
    r2.write_bytes(b"\x1f\x8b")

    repro = tmp_path / "reproducibility"
    repro.mkdir()
    ss = repro / "samplesheet.valid.csv"
    ss.write_text(
        "patient,sample,fastq_1,fastq_2\n"
        f"P1,S1,{r1.as_posix()},{r2.as_posix()}\n",
        encoding="utf-8",
    )

    old_prefix = data.resolve().as_posix()
    # Destination prefix with a space AND a non-ASCII character.
    moved = tmp_path / "relocated café data"
    data.rename(moved)
    new_prefix = moved.resolve().as_posix()

    # Before remap: the moved paths must be reported missing.
    assert remap_paths.verify_paths(ss), "paths should be missing after the move"

    changes = remap_paths.remap_csv(ss, old_prefix, new_prefix, dry_run=False)
    assert changes, "remap should have rewritten the FASTQ paths"

    # After remap: every path resolves on the 'new machine', incl. space+unicode.
    assert remap_paths.verify_paths(ss) == [], "all paths must resolve after remap"
    # The rewrite preserved LF endings (cross-OS byte stability — no CRLF).
    assert b"\r" not in ss.read_bytes()


def test_module_docstring_does_not_claim_all_inputs_are_absolute():
    """The remap_paths.py header must not state unconditionally that every
    samplesheet input path is stored as an absolute path. nf-core/sarek accepts
    remote data URLs (passed through unchanged), which this helper's verify step
    already skips — so the blanket claim is false and self-contradictory with the
    header's own note that remote URIs are skipped. Mirrors the rnaseq fix.
    """
    doc = remap_paths.__doc__ or ""
    # The absolute-path statement must be scoped to local inputs, not unconditional.
    assert "Local samplesheet input paths" in doc
    # The remote-input case must be acknowledged in the header intro itself.
    assert "passed through unchanged" in doc


# ── main(): requested actions must compose, never be silently dropped ─────────


def _bundle_for_main(bundle: Path, *, old_output: str) -> None:
    """A bundle whose commands.sh points --output at ``old_output``."""
    bundle.mkdir(parents=True, exist_ok=True)
    fastq = bundle / "S1.fastq.gz"
    fastq.write_bytes(b"")
    (bundle / "samplesheet.valid.csv").write_text(
        "patient,sample,lane,fastq_1,fastq_2\n"
        f"P1,S1,1,{fastq},\n",
        encoding="utf-8",
    )
    (bundle / "commands.sh").write_text(
        "nextflow run nf-core/sarek \\\n"
        f"    --output {old_output} \\\n",
        encoding="utf-8",
    )


def test_main_output_dir_and_verify_compose(tmp_path: Path):
    """`--output-dir X --verify` must rewrite --output AND then verify.

    The dispatch used to be a chain of early returns with --verify ahead of
    --output-dir, so combining them silently discarded the rewrite and then reported
    the bundle as ready to replay while commands.sh still pointed at the OLD output
    directory — a confident green light on a stale bundle.
    """
    bundle = tmp_path / "reproducibility"
    new_output = tmp_path / "relocated_run"
    _bundle_for_main(bundle, old_output="/original/output")

    remap_paths.main(["--output-dir", str(new_output), "--verify"], bundle_dir=bundle)

    content = (bundle / "commands.sh").read_text(encoding="utf-8")
    assert str(new_output) in content, "--output-dir rewrite was silently dropped by --verify"
    assert "/original/output" not in content


def test_main_verify_reports_on_the_rewritten_bundle(tmp_path: Path, capsys):
    """--verify must run AFTER the rewrite, so it never reports the stale output dir."""
    bundle = tmp_path / "reproducibility"
    new_output = tmp_path / "relocated_run"
    _bundle_for_main(bundle, old_output="/original/output")

    remap_paths.main(["--output-dir", str(new_output), "--verify"], bundle_dir=bundle)

    out = capsys.readouterr().out
    assert "/original/output" not in out, "--verify reported on the stale, pre-rewrite bundle"


def test_main_rejects_half_a_prefix_pair(tmp_path: Path):
    """--old without --new (or vice versa) must fail loudly, not be silently ignored."""
    import pytest

    bundle = tmp_path / "reproducibility"
    _bundle_for_main(bundle, old_output="/original/output")
    for argv in (["--old", "/a"], ["--new", "/b"], ["--refs-old", "/a"], ["--refs-new", "/b"]):
        with pytest.raises(SystemExit) as exc:
            remap_paths.main(argv, bundle_dir=bundle)
        assert exc.value.code != 0
