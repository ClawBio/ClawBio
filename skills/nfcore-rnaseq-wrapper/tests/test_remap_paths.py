from __future__ import annotations

import csv
from pathlib import Path
import sys

import pytest

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


_purge_foreign_bare_modules("remap_paths")

from remap_paths import (
    _format_output_value,
    cmd_remap,
    cmd_remap_references,
    cmd_repair_bundle,
    cmd_update_output,
    cmd_verify,
    find_samplesheet,
    main,
    remap_commands_references,
    remap_csv,
    update_commands_output,
    verify_paths,
    verify_reference_paths,
)

# Capture the module docstring before it is purged from sys.modules below, so the
# header-accuracy test can assert on it without re-importing.
_REMAP_PATHS_DOC = sys.modules["remap_paths"].__doc__ or ""

_purge_local_modules("remap_paths")
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))


def _write_minimal_bundle_files(bundle_dir: Path) -> None:
    """Write the three required bundle files with minimal stub content.

    Use this in tests that need --verify to pass but don't care about
    the accuracy of manifest/checksums/environment content.
    """
    (bundle_dir / "manifest.json").write_text('{"stub": true}', encoding="utf-8")
    (bundle_dir / "checksums.sha256").write_text("", encoding="utf-8")
    (bundle_dir / "environment.yml").write_text("name: test\n", encoding="utf-8")


def _write_samplesheet(path: Path, fastq_1: str, fastq_2: str = "") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample", "fastq_1", "fastq_2", "strandedness"])
        writer.writeheader()
        writer.writerow({"sample": "S1", "fastq_1": fastq_1, "fastq_2": fastq_2, "strandedness": "auto"})
    return path


def _write_multi_sample_sheet(path: Path, rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample", "fastq_1", "fastq_2", "strandedness"])
        writer.writeheader()
        writer.writerows(rows)
    return path


# ── find_samplesheet ────────────────────────────────────────────────────────


def test_find_samplesheet_returns_valid_csv_when_present(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    ss.write_text("sample,fastq_1,strandedness\n", encoding="utf-8")
    assert find_samplesheet(bundle_dir=tmp_path) == ss


# ── remap_csv ───────────────────────────────────────────────────────────────


def test_remap_csv_replaces_matching_prefix(tmp_path):
    ss = _write_samplesheet(tmp_path / "samplesheet.valid.csv", "/old/data/S1_R1.fastq.gz", "/old/data/S1_R2.fastq.gz")
    changes = remap_csv(ss, "/old/data", "/new/data", dry_run=False)
    assert len(changes) == 2
    assert changes[0] == ("fastq_1", "/old/data/S1_R1.fastq.gz", "/new/data/S1_R1.fastq.gz")
    rows = list(csv.DictReader(ss.read_text(encoding="utf-8").splitlines()))
    assert rows[0]["fastq_1"] == "/new/data/S1_R1.fastq.gz"
    assert rows[0]["fastq_2"] == "/new/data/S1_R2.fastq.gz"


# ── verify_paths ────────────────────────────────────────────────────────────


def test_verify_paths_returns_missing_paths(tmp_path):
    ss = _write_samplesheet(tmp_path / "samplesheet.valid.csv", "/nonexistent/R1.fastq.gz")
    missing = verify_paths(ss)
    assert "/nonexistent/R1.fastq.gz" in missing


# ── cmd_remap ───────────────────────────────────────────────────────────────


def test_cmd_remap_returns_nonzero_when_no_samplesheet(tmp_path):
    assert cmd_remap("/old", "/new", dry_run=False, bundle_dir=tmp_path) != 0


def test_cmd_remap_exit_code_consistent_across_reruns_when_target_missing(tmp_path):
    """Re-running the same remap must give the SAME exit code for the SAME readiness
    state. A first run that remaps onto missing targets is not replay-ready (non-zero);
    a second, idempotent run (nothing left to change) must NOT silently return 0 while
    the bundle is still not replay-ready — that inconsistency breaks scripted relocation
    and re-execution."""
    _write_samplesheet(tmp_path / "samplesheet.valid.csv", "/old/data/R1.fastq.gz")
    first = cmd_remap("/old/data", "/missing/data", dry_run=False, bundle_dir=tmp_path)
    second = cmd_remap("/old/data", "/missing/data", dry_run=False, bundle_dir=tmp_path)
    assert first != 0
    assert second != 0, "idempotent re-run must report the same not-ready state, not exit 0"


def test_cmd_remap_returns_zero_when_remapped_paths_exist(tmp_path):
    """When the remapped FASTQ paths exist on disk, the bundle is replay-ready → exit 0."""
    real = tmp_path / "newdata"
    real.mkdir()
    (real / "R1.fastq.gz").write_bytes(b"")
    _write_samplesheet(tmp_path / "samplesheet.valid.csv", "/old/data/R1.fastq.gz")
    rc = cmd_remap("/old/data", real.as_posix(), dry_run=False, bundle_dir=tmp_path)
    assert rc == 0


def test_cmd_remap_dry_run_does_not_gate_on_missing_paths(tmp_path):
    """--dry-run is a preview: it must never fail on missing targets (files unchanged)."""
    _write_samplesheet(tmp_path / "samplesheet.valid.csv", "/old/data/R1.fastq.gz")
    assert cmd_remap("/old/data", "/missing/data", dry_run=True, bundle_dir=tmp_path) == 0


# ── cmd_verify ──────────────────────────────────────────────────────────────


def test_cmd_verify_returns_nonzero_when_paths_missing(tmp_path):
    _write_samplesheet(tmp_path / "samplesheet.valid.csv", "/nonexistent/R1.fastq.gz")
    assert cmd_verify(bundle_dir=tmp_path) != 0


def test_cmd_verify_passes_when_output_dir_does_not_exist(tmp_path):
    """A fresh replay produces a non-existent output dir (the wrapper creates it).
    --verify must not fail solely because the output dir is missing."""
    fastq = tmp_path / "S1.fastq.gz"
    fastq.write_bytes(b"")
    _write_samplesheet(tmp_path / "samplesheet.valid.csv", str(fastq))
    commands_sh = tmp_path / "commands.sh"
    commands_sh.write_text(
        "python skill.py \\\n"
        f"    --output {tmp_path / 'fresh_output'} \\\n",
        encoding="utf-8",
    )
    _write_minimal_bundle_files(tmp_path)
    result = cmd_verify(bundle_dir=tmp_path)
    assert result == 0, "--verify must pass (exit 0) when output dir doesn't exist yet"


def test_cmd_verify_fails_when_bundle_files_missing(tmp_path, capsys):
    # --verify must treat an incomplete bundle (missing manifest.json,
    # checksums.sha256, environment.yml) as a verification failure (exit != 0),
    # not a mere warning. A missing bundle file means the audit contract is broken.
    fastq = tmp_path / "S1.fastq.gz"
    fastq.write_bytes(b"")
    _write_samplesheet(tmp_path / "samplesheet.valid.csv", str(fastq))
    commands_sh = tmp_path / "commands.sh"
    commands_sh.write_text(
        "python skill.py \\\n"
        "    --output /some/new/output \\\n",
        encoding="utf-8",
    )
    # Do NOT write manifest.json, checksums.sha256, or environment.yml
    result = cmd_verify(bundle_dir=tmp_path)
    captured = capsys.readouterr()
    assert result != 0, "--verify must fail (non-zero) when bundle files are missing"
    assert "manifest" in captured.out.lower() or "incomplete" in captured.out.lower(), (
        "Expected a message about missing manifest.json in an incomplete bundle"
    )


# ── cmd_repair_bundle ───────────────────────────────────────────────────────


def test_cmd_repair_bundle_noop_when_all_files_present(tmp_path):
    """When all required bundle files already exist, --repair-bundle is a no-op
    and returns 0 immediately without modifying any files."""
    commands_sh = tmp_path / "commands.sh"
    commands_sh.write_text(
        "python skill.py \\\n"
        f"    --output {tmp_path} \\\n",
        encoding="utf-8",
    )
    (tmp_path / "manifest.json").write_text('{"complete": true}', encoding="utf-8")
    (tmp_path / "checksums.sha256").write_text("abc123  some_file\n", encoding="utf-8")
    (tmp_path / "environment.yml").write_text("name: test\n", encoding="utf-8")

    result = cmd_repair_bundle(bundle_dir=tmp_path)
    assert result == 0, "--repair-bundle must return 0 when all required files already exist"
    # Files must not be modified when bundle is already complete
    assert (tmp_path / "manifest.json").read_text(encoding="utf-8") == '{"complete": true}'


def test_cmd_repair_bundle_returns_zero_and_creates_all_files(tmp_path):
    """When bundle files are missing, --repair-bundle must create all three and return 0."""
    (tmp_path / "commands.sh").write_text(
        "python skill.py \\\n"
        f"    --output {tmp_path} \\\n",
        encoding="utf-8",
    )
    # No bundle files present
    assert not (tmp_path / "manifest.json").exists()
    assert not (tmp_path / "checksums.sha256").exists()
    assert not (tmp_path / "environment.yml").exists()

    result = cmd_repair_bundle(bundle_dir=tmp_path)
    assert result == 0, "--repair-bundle must return 0 after successful repair"
    assert (tmp_path / "manifest.json").exists(), "manifest.json must be created"
    assert (tmp_path / "checksums.sha256").exists(), "checksums.sha256 must be created"
    assert (tmp_path / "environment.yml").exists(), "environment.yml must be created"


def test_cmd_repair_bundle_manifest_stub_has_post_hoc_marker(tmp_path):
    """Regenerated manifest.json must be marked regenerated_post_hoc: true so
    audit consumers can distinguish it from the original."""
    import json

    (tmp_path / "commands.sh").write_text(
        "python skill.py \\\n    --output /some/output \\\n",
        encoding="utf-8",
    )
    cmd_repair_bundle(bundle_dir=tmp_path)

    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    assert manifest.get("regenerated_post_hoc") is True, (
        "manifest.json must contain regenerated_post_hoc: true"
    )
    assert "commands_sh" in manifest, "manifest.json must reference commands_sh for traceability"


def test_cmd_repair_bundle_checksums_updated_when_only_manifest_missing(tmp_path):
    """When only manifest.json is missing, repair must recompute checksums.sha256
    so the checksum file covers the newly-created manifest stub.

    Regression: previously checksums were only recomputed if checksums.sha256
    itself was in the missing set, leaving the checksum file stale after a
    partial repair."""
    (tmp_path / "commands.sh").write_text(
        "python skill.py \\\n    --output /some/output \\\n",
        encoding="utf-8",
    )
    (tmp_path / "params.yaml").write_text("genome: GRCh38\n", encoding="utf-8")

    # Pre-create environment.yml and checksums.sha256 — only manifest.json is absent
    (tmp_path / "environment.yml").write_text(
        "name: test\nregenerated_post_hoc: true\n", encoding="utf-8"
    )
    # Checksums pre-repair: covers only params.yaml and environment.yml, NOT manifest
    import hashlib

    def _sha256(p):
        h = hashlib.sha256()
        h.update(p.read_bytes())
        return h.hexdigest()

    old_checksums = (
        f"{_sha256(tmp_path / 'params.yaml')}  params.yaml\n"
        f"{_sha256(tmp_path / 'environment.yml')}  environment.yml\n"
    )
    (tmp_path / "checksums.sha256").write_text(old_checksums, encoding="utf-8")

    result = cmd_repair_bundle(bundle_dir=tmp_path)
    assert result == 0

    new_checksums = (tmp_path / "checksums.sha256").read_text(encoding="utf-8")
    assert "manifest.json" in new_checksums, (
        "checksums.sha256 must be recomputed after repair to include the new manifest.json"
    )


def test_cmd_repair_bundle_checksums_covers_bundle_contents(tmp_path):
    """Regenerated checksums.sha256 must include entries for files present
    in the bundle directory."""
    commands_sh = tmp_path / "commands.sh"
    commands_sh.write_text(
        "python skill.py \\\n    --output /some/output \\\n",
        encoding="utf-8",
    )
    # Put a content file in the bundle dir (simulates params.yaml, samplesheet, etc.)
    (tmp_path / "params.yaml").write_text("genome: GRCh38\n", encoding="utf-8")

    cmd_repair_bundle(bundle_dir=tmp_path)

    checksums_text = (tmp_path / "checksums.sha256").read_text(encoding="utf-8")
    # checksums.sha256 must not list itself
    assert "checksums.sha256" not in checksums_text
    # params.yaml lives in the bundle dir and must be covered
    assert "params.yaml" in checksums_text, (
        "checksums.sha256 must include an entry for params.yaml which is present in the bundle"
    )


def test_cmd_repair_bundle_is_idempotent(tmp_path):
    """Running --repair-bundle twice must be a no-op on the second call
    (returns 0 and does not raise). This allows automation to run it safely."""
    (tmp_path / "commands.sh").write_text(
        "python skill.py \\\n    --output /some/output \\\n",
        encoding="utf-8",
    )
    first = cmd_repair_bundle(bundle_dir=tmp_path)
    assert first == 0, "First repair must succeed"

    second = cmd_repair_bundle(bundle_dir=tmp_path)
    assert second == 0, "Second repair (idempotent no-op) must also return 0"


def test_cmd_repair_bundle_returns_nonzero_when_commands_sh_missing(tmp_path):
    """--repair-bundle must fail (non-zero) when commands.sh is absent —
    it cannot know what to repair without the original command."""
    result = cmd_repair_bundle(bundle_dir=tmp_path)
    assert result != 0, "--repair-bundle must fail when commands.sh is missing"


def test_cmd_verify_mentions_repair_bundle_when_bundle_files_missing(tmp_path, capsys):
    """When --verify finds missing bundle files it must suggest --repair-bundle
    so users know how to fix the problem, and must return non-zero."""
    fastq = tmp_path / "S1.fastq.gz"
    fastq.write_bytes(b"")
    _write_samplesheet(tmp_path / "samplesheet.valid.csv", str(fastq))
    commands_sh = tmp_path / "commands.sh"
    commands_sh.write_text(
        "python skill.py \\\n"
        "    --output /some/output \\\n",
        encoding="utf-8",
    )
    # Do NOT write manifest.json, checksums.sha256, or environment.yml
    result = cmd_verify(bundle_dir=tmp_path)
    captured = capsys.readouterr()
    assert result != 0, "--verify must fail (non-zero) when bundle files are missing"
    assert "--repair-bundle" in captured.out, (
        "--verify must mention --repair-bundle when bundle files are missing"
    )


# ── integration smoke ───────────────────────────────────────────────────────


def test_remap_csv_backup_preserves_full_filename(tmp_path):
    """After remap_csv the backup must keep the full original filename plus .bak
    (e.g. samplesheet.valid.csv.bak), not lose the .csv extension
    (which with_suffix('.bak') would produce samplesheet.valid.bak)."""
    fastq_new = tmp_path / "new" / "S1.fastq.gz"
    fastq_new.parent.mkdir()
    fastq_new.write_bytes(b"")
    samplesheet = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(samplesheet, "/old/S1.fastq.gz")

    remap_csv(samplesheet, old_prefix="/old", new_prefix=str(tmp_path / "new"), dry_run=False)

    expected_backup = tmp_path / "samplesheet.valid.csv.bak"
    wrong_backup = tmp_path / "samplesheet.valid.bak"
    assert expected_backup.exists(), (
        f"Backup must be named '{expected_backup.name}' to preserve the .csv extension; "
        f"found wrong_backup_exists={wrong_backup.exists()}"
    )
    assert not wrong_backup.exists(), (
        "with_suffix('.bak') creates a misleading backup that drops .csv — must use "
        "with_name(name + '.bak') instead."
    )


def test_verify_reference_paths_skips_unbraced_bash_variable(tmp_path):
    """verify_reference_paths must skip paths starting with '$' (not just '${...}').
    An unbraced variable like $REFS/genome.fa is valid shell syntax and cannot be
    checked with a filesystem stat — it should not be reported as missing."""
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    commands_sh = tmp_path / "commands.sh"
    commands_sh.write_text(
        f"python skill.py \\\n"
        f"    --fasta $REFS/genome.fa \\\n"         # unbraced — should be skipped
        f"    --gtf ${{REFS}}/genes.gtf \\\n"       # braced — already skipped
        f"    --star-index {fasta} \\\n"             # real path — must exist
        f"    --output /tmp/out \\\n",
        encoding="utf-8",
    )
    missing = verify_reference_paths(commands_sh)
    assert "$REFS/genome.fa" not in missing, (
        "Unbraced bash variable $REFS/genome.fa must be skipped, not reported as missing; "
        f"got missing={missing!r}"
    )
    assert "${REFS}/genes.gtf" not in missing, (
        "Braced bash variable ${REFS}/genes.gtf must be skipped (already passing); "
        f"got missing={missing!r}"
    )
    assert str(fasta) not in missing, (
        f"Existing path {fasta} must not be reported as missing; got missing={missing!r}"
    )


# ── find_samplesheet priority order ─────────────────────────────────────────

def test_find_samplesheet_prefers_valid_over_demo(tmp_path):
    """valid > demo: when both exist, samplesheet.valid.csv must be returned."""
    valid = tmp_path / "samplesheet.valid.csv"
    demo = tmp_path / "samplesheet.demo.csv"
    valid.write_text("sample,fastq_1,fastq_2,strandedness\n", encoding="utf-8")
    demo.write_text("sample,fastq_1,fastq_2,strandedness\n", encoding="utf-8")
    result = find_samplesheet(bundle_dir=tmp_path)
    assert result == valid, f"Expected valid > demo; got {result}"


# ── remap_commands_references ────────────────────────────────────────────────

def test_remap_commands_references_replaces_prefix(tmp_path):
    """remap_commands_references must rewrite reference paths whose prefix matches."""
    commands_sh = tmp_path / "commands.sh"
    commands_sh.write_text(
        "python skill.py \\\n"
        "    --fasta /old/refs/genome.fa \\\n"
        "    --gtf /old/refs/genes.gtf \\\n"
        "    --output /tmp/out \\\n",
        encoding="utf-8",
    )
    changes = remap_commands_references(commands_sh, "/old/refs", "/new/refs", dry_run=False)
    assert len(changes) == 2, f"Expected 2 changes; got {changes!r}"
    text = commands_sh.read_text(encoding="utf-8")
    assert "/new/refs/genome.fa" in text
    assert "/new/refs/genes.gtf" in text
    assert "/old/refs" not in text
    assert (tmp_path / "commands.sh.bak").exists(), "Backup commands.sh.bak must be created"


# ── verify_reference_paths ───────────────────────────────────────────────────

def test_verify_reference_paths_skips_remote_uri(tmp_path):
    """s3://, https:// and other remote URIs must be skipped (unreachable via stat)."""
    commands_sh = tmp_path / "commands.sh"
    commands_sh.write_text(
        "python skill.py \\\n"
        "    --fasta s3://my-bucket/genome.fa \\\n"
        "    --gtf https://example.com/genes.gtf \\\n"
        "    --output /tmp/out \\\n",
        encoding="utf-8",
    )
    missing = verify_reference_paths(commands_sh)
    assert not any("s3://" in m for m in missing), f"s3:// URIs must be skipped; got {missing!r}"
    assert not any("https://" in m for m in missing), f"https:// URIs must be skipped; got {missing!r}"


def test_verify_reference_paths_reports_nonexistent_local_path(tmp_path):
    """A local path that does not exist on disk must be reported as missing."""
    commands_sh = tmp_path / "commands.sh"
    commands_sh.write_text(
        "python skill.py \\\n"
        "    --fasta /does/not/exist/genome.fa \\\n"
        "    --output /tmp/out \\\n",
        encoding="utf-8",
    )
    missing = verify_reference_paths(commands_sh)
    assert "/does/not/exist/genome.fa" in missing, (
        f"Non-existent local path must be reported; got {missing!r}"
    )


# ── cmd_remap_references ─────────────────────────────────────────────────────

def test_cmd_remap_references_returns_nonzero_when_no_commands_sh(tmp_path):
    """cmd_remap_references must return 1 if commands.sh is absent."""
    result = cmd_remap_references("/old", "/new", dry_run=False, bundle_dir=tmp_path)
    assert result == 1


# ── cmd_update_output ────────────────────────────────────────────────────────

def test_cmd_update_output_returns_nonzero_when_no_commands_sh(tmp_path):
    """cmd_update_output must return 1 when commands.sh is absent."""
    result = cmd_update_output("/new/out", dry_run=False, bundle_dir=tmp_path)
    assert result == 1


def test_cmd_update_output_replaces_output_path(tmp_path):
    """cmd_update_output must rewrite the --output value in commands.sh."""
    commands_sh = tmp_path / "commands.sh"
    commands_sh.write_text(
        "python skill.py \\\n    --output /old/out \\\n",
        encoding="utf-8",
    )
    result = cmd_update_output("/new/out", dry_run=False, bundle_dir=tmp_path)
    assert result == 0
    assert "/new/out" in commands_sh.read_text(encoding="utf-8")
    assert (tmp_path / "commands.sh.bak").exists()


def test_verify_paths_skips_remote_uri_fastq(tmp_path):
    """verify_paths must not report remote FASTQ/BAM URIs (s3://, https://, etc.) as
    missing — they cannot be stat-checked locally and are valid nf-core samplesheet
    entries for cloud execution."""
    samplesheet = tmp_path / "samplesheet.valid.csv"
    samplesheet.write_text(
        "sample,fastq_1,fastq_2,strandedness\n"
        "S1,s3://bucket/S1_R1.fastq.gz,s3://bucket/S1_R2.fastq.gz,auto\n"
        "S2,https://example.com/S2.fastq.gz,,auto\n",
        encoding="utf-8",
    )
    missing = verify_paths(samplesheet)
    assert not any("s3://" in m for m in missing), (
        f"s3:// URIs must be skipped in verify_paths; got missing={missing!r}"
    )
    assert not any("https://" in m for m in missing), (
        f"https:// URIs must be skipped in verify_paths; got missing={missing!r}"
    )


def test_remap_paths_integration_smoke(tmp_path):
    """Full portability workflow: remap paths → repair bundle → verify ready to replay."""
    fastq = tmp_path / "fastqs" / "S1.fastq.gz"
    fastq.parent.mkdir()
    fastq.write_bytes(b"")
    samplesheet = tmp_path / "samplesheet.valid.csv"
    _write_samplesheet(samplesheet, "/old/S1.fastq.gz")
    # Step 1: remap FASTQ paths to the new machine location
    assert cmd_remap("/old", str(tmp_path / "fastqs"), dry_run=False, bundle_dir=tmp_path) == 0
    assert str(fastq) in samplesheet.read_text(encoding="utf-8")
    # Step 2: repair the bundle (creates checksums + stubs for missing provenance files)
    (tmp_path / "commands.sh").write_text(
        "python skill.py \\\n    --output /some/output \\\n",
        encoding="utf-8",
    )
    assert cmd_repair_bundle(bundle_dir=tmp_path) == 0
    # Step 3: verify the complete bundle is ready to replay
    assert cmd_verify(bundle_dir=tmp_path) == 0


# ── _format_output_value ────────────────────────────────────────────────────


def test_format_output_value_dollar_gets_double_quoted():
    """Paths containing $ must be double-quoted to prevent shell expansion."""
    result = _format_output_value("/data/$USER/out", old_value="/old/out")
    assert result.startswith('"') and result.endswith('"')
    assert "\\$USER" in result


def test_format_output_value_plain_path_unchanged():
    """Clean paths with no special characters are returned as-is."""
    result = _format_output_value("/data/clean/out", old_value="/old/out")
    assert result == "/data/clean/out"


# ── Audit follow-up: remap_csv must preserve the wrapper's '\n' line terminator ──
# (samplesheet_builder writes LF; the default csv.writer would silently rewrite the
# file with CRLF on the documented remap step, changing every line's bytes).


def test_remap_csv_preserves_unix_line_endings(tmp_path):
    ss = tmp_path / "samplesheet.valid.csv"
    ss.write_text(
        "sample,fastq_1,fastq_2,strandedness\n"
        "s1,/old/data/s1_R1.fastq.gz,,auto\n",
        encoding="utf-8",
    )
    changes = remap_csv(ss, "/old/data", "/new/data", dry_run=False)
    assert changes  # the row was remapped
    raw = ss.read_bytes()
    assert b"\r\n" not in raw
    assert raw.endswith(b"\n")
    assert b"/new/data/s1_R1.fastq.gz" in raw


def test_update_commands_output_rewrites_resume_guard_manifest(tmp_path):
    """--output-dir must relocate BOTH the --output argument and the idempotent-resume
    guard's manifest path, so a remapped bundle checks the NEW output directory for a
    completed run (not the original one). Otherwise a same-machine remap could
    wrongly resume against a stale manifest."""
    commands = tmp_path / "commands.sh"
    commands.write_text(
        'CLAWBIO_RESUME=""\n'
        'if [[ -f "/old/out/reproducibility/manifest.json" ]]; then\n'
        '  CLAWBIO_RESUME="--resume"\n'
        'fi\n'
        '"${PYTHON:-python3}" "$SKILL_SCRIPT" $CLAWBIO_RESUME \\\n'
        "    --output /old/out \\\n"
        "    --profile docker\n",
        encoding="utf-8",
    )
    changed = update_commands_output(commands, "/new/out")
    assert changed
    updated = commands.read_text(encoding="utf-8")
    assert 'if [[ -f "/new/out/reproducibility/manifest.json" ]]; then' in updated
    assert "--output /new/out" in updated
    assert "/old/out" not in updated


def test_module_docstring_does_not_claim_all_inputs_are_absolute():
    """The remap_paths.py header must not state unconditionally that every FASTQ/BAM
    and reference path is stored as an absolute path. A run started with
    --allow-remote-inputs stores remote URIs (which the remap/verify helpers here
    already skip), so the blanket claim is false and self-contradictory with the
    header's own note that remote URIs are skipped. Mirrors the fix to the
    commands.sh portability notice and the scrnaseq remap_paths header.
    """
    doc = _REMAP_PATHS_DOC
    assert "are stored as\nabsolute paths (required by Nextflow)" not in doc
    # The remote-input case must be acknowledged in the header itself.
    assert "--allow-remote-inputs" in doc
    assert "remote" in doc.lower()


# ── main(): requested actions must compose, never be silently dropped ─────────


def _bundle_for_main(tmp_path: Path, old_output: str) -> Path:
    """A complete bundle whose commands.sh points --output at ``old_output``."""
    fastq = tmp_path / "S1.fastq.gz"
    fastq.write_bytes(b"")
    _write_samplesheet(tmp_path / "samplesheet.valid.csv", str(fastq))
    (tmp_path / "commands.sh").write_text(
        "python skill.py \\\n"
        f"    --output {old_output} \\\n",
        encoding="utf-8",
    )
    _write_minimal_bundle_files(tmp_path)
    return tmp_path


def test_main_output_dir_and_verify_compose(tmp_path):
    """`--output-dir X --verify` must rewrite --output AND then verify.

    The dispatch used to be a chain of early returns with --verify ahead of
    --output-dir, so combining them silently discarded the rewrite and then printed
    "ready to replay" over a bundle that still pointed at the OLD output directory —
    the worst possible outcome: a confident green light on a stale bundle.
    """
    new_output = tmp_path / "relocated_run"
    _bundle_for_main(tmp_path, old_output="/original/output")

    rc = main(["--output-dir", str(new_output), "--verify"], bundle_dir=tmp_path)

    content = (tmp_path / "commands.sh").read_text(encoding="utf-8")
    assert str(new_output) in content, "--output-dir rewrite was silently dropped by --verify"
    assert "/original/output" not in content
    assert rc == 0


def test_main_verify_runs_after_the_rewrite_not_before(tmp_path, capsys):
    """--verify must report on the POST-rewrite bundle, so its output names the new dir."""
    new_output = tmp_path / "relocated_run"
    _bundle_for_main(tmp_path, old_output="/original/output")

    main(["--output-dir", str(new_output), "--verify"], bundle_dir=tmp_path)

    out = capsys.readouterr().out
    assert "/original/output" not in out, "--verify reported on the stale, pre-rewrite bundle"


def test_main_remap_and_verify_compose(tmp_path):
    """The same silent-drop bug applied to --old/--new + --verify."""
    real = tmp_path / "newdata"
    real.mkdir()
    (real / "S1.fastq.gz").write_bytes(b"")
    _write_samplesheet(tmp_path / "samplesheet.valid.csv", "/olddata/S1.fastq.gz")
    (tmp_path / "commands.sh").write_text(
        "python skill.py \\\n"
        f"    --output {tmp_path / 'out'} \\\n",
        encoding="utf-8",
    )
    _write_minimal_bundle_files(tmp_path)

    rc = main(["--old", "/olddata", "--new", str(real), "--verify"], bundle_dir=tmp_path)

    sheet = (tmp_path / "samplesheet.valid.csv").read_text(encoding="utf-8")
    assert str(real) in sheet, "--old/--new remap was silently dropped by --verify"
    assert rc == 0


def test_main_rejects_half_a_prefix_pair(tmp_path):
    """--old without --new (or vice versa) must fail loudly, not be silently ignored."""
    _bundle_for_main(tmp_path, old_output="/original/output")
    for argv in (["--old", "/a"], ["--new", "/b"], ["--refs-old", "/a"], ["--refs-new", "/b"]):
        with pytest.raises(SystemExit) as exc:
            main(argv, bundle_dir=tmp_path)
        assert exc.value.code != 0
