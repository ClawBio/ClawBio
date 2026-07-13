#!/usr/bin/env python3
"""
remap_paths.py — Make this reproducibility bundle portable across machines.

Local FASTQ/BAM paths (fastq_1, fastq_2, genome_bam, transcriptome_bam) and local
reference/index paths (--fasta, --gtf, --star-index, etc.) are stored as absolute
paths (required by Nextflow). A run started with --allow-remote-inputs instead
stores remote URIs (s3://, https://, ...), which resolve from anywhere and are
skipped by the remap and verify steps below. Before replaying on a different machine:

  1. Remap FASTQ/BAM paths in the samplesheet (if data moved):
       python3 remap_paths.py --old /original/data/dir --new /new/data/dir

  2. Remap reference/index paths in commands.sh (if references moved):
       python3 remap_paths.py --refs-old /original/refs --refs-new /new/refs

  3. Update the --output path in commands.sh (if output dir changed):
       python3 remap_paths.py --output-dir /new/output/dir

  4. Verify everything is ready:
       python3 remap_paths.py --verify

  5. Regenerate missing bundle files (manifest.json, checksums.sha256, environment.yml):
       python3 remap_paths.py --repair-bundle

  Preview any change without modifying files by adding --dry-run.

  Remote URIs (s3://, https://, etc.) are recognized and skipped automatically
  during path verification and checksumming — they are never resolved locally.

  CLAWBIO_REPO must always be set to replay:
       CLAWBIO_REPO=/path/to/ClawBio bash commands.sh
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shlex
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

_BUNDLE_DIR = Path(__file__).resolve().parent


def _write_text_lf(path: Path, text: str) -> None:
    """Write ``text`` with LF line endings on every OS (self-contained).

    This script is shipped inside the reproducibility bundle and runs standalone
    at replay (no ClawBio package available), so it cannot import the shared
    ``clawbio.common.textio`` helper. CRLF/CR are normalised to LF and the result
    is written as bytes (bypassing text-mode newline translation) so a rewrite on
    Windows never reintroduces CRLF into commands.sh, params.yaml, or checksums.
    Mirrors the nfcore-sarek-wrapper helper.
    """
    normalised = text.replace("\r\n", "\n").replace("\r", "\n")
    Path(path).write_bytes(normalised.encode("utf-8"))


_FASTQ_COLUMNS = ("fastq_1", "fastq_2")
_BAM_COLUMNS = ("genome_bam", "transcriptome_bam")
_REMAPPABLE_COLUMNS = _FASTQ_COLUMNS + _BAM_COLUMNS
_REQUIRED_BUNDLE_FILES = ("manifest.json", "checksums.sha256", "environment.yml")

# All flags in commands.sh that hold local absolute paths that may need
# remapping when the bundle is replayed on a different machine.
_REFERENCE_FLAGS = (
    "fasta",
    "gtf",
    "gff",
    "transcript-fasta",
    "additional-fasta",
    "gene-bed",
    "splicesites",
    "star-index",
    "rsem-index",
    "hisat2-index",
    "bowtie2-index",
    "salmon-index",
    "kallisto-index",
    "sortmerna-index",
    "kraken-db",
    "bbsplit-fasta-list",
    "bbsplit-index",
    "sylph-db",
    "sylph-taxonomy",
    "multiqc-config",
    "multiqc-logo",
    "multiqc-methods-description",
    "igenomes-base",
    "pipeline-local",
    "nextflow-config",
    "ribo-database-manifest",
)

# Matches `    --output /path \` lines — requires line to start with whitespace
# (not #) so comment lines containing --output are never modified. The value
# may be unquoted, single-quoted, or double-quoted.
_OUTPUT_FLAG_RE = re.compile(
    r"""^([ \t]+--output[ \t]+)(?P<value>"[^"\n]*"|'[^'\n]*'|\S+)([ \t]*(?:\\[ \t]*)?)$""",
    re.MULTILINE,
)

# Matches the idempotent-replay guard line written by reporting.py:
#   if [[ -f "<output_dir>/reproducibility/manifest.json" ]]; then
# The captured value is the output directory, kept in sync with the --output flag
# so a relocated bundle checks the NEW output dir (not the original) for a prior run.
_RESUME_GUARD_RE = re.compile(
    r"""^(if \[\[ -f ")(?P<value>[^"\n]*)(/reproducibility/manifest\.json" \]\]; then)$""",
    re.MULTILINE,
)

# Matches any of the reference flags in commands.sh, same quoting rules.
_REFERENCE_FLAG_RE = re.compile(
    r"""^([ \t]+--(?:"""
    + "|".join(re.escape(f) for f in _REFERENCE_FLAGS)
    + r""")[ \t]+)(?P<value>"[^"\n]*"|'[^'\n]*'|\S+)([ \t]*(?:\\[ \t]*)?)$""",
    re.MULTILINE,
)


def find_samplesheet(bundle_dir: Path | None = None) -> Path | None:
    search_dir = bundle_dir or _BUNDLE_DIR
    for name in ("samplesheet.valid.csv", "samplesheet.demo.csv", "samplesheet.noinput.csv"):
        p = search_dir / name
        if p.exists():
            return p
    return None


def find_commands_sh(bundle_dir: Path | None = None) -> Path | None:
    search_dir = bundle_dir or _BUNDLE_DIR
    p = search_dir / "commands.sh"
    return p if p.exists() else None


def _prefix_matches(path: str, prefix: str) -> bool:
    if path == prefix:
        return True
    sep = "" if prefix.endswith("/") else "/"
    return path.startswith(prefix + sep)


def remap_csv(
    samplesheet: Path,
    old_prefix: str,
    new_prefix: str,
    *,
    dry_run: bool,
) -> list[tuple[str, str, str]]:
    """Return list of (column, old_path, new_path) for every changed cell."""
    text = samplesheet.read_text(encoding="utf-8")
    fieldnames = list(csv.DictReader(text.splitlines()).fieldnames or [])
    rows = list(csv.DictReader(text.splitlines()))
    changes: list[tuple[str, str, str]] = []

    for row in rows:
        for col in _REMAPPABLE_COLUMNS:
            if col in row and row[col] and _prefix_matches(row[col], old_prefix):
                new_val = new_prefix + row[col][len(old_prefix):]
                changes.append((col, row[col], new_val))
                if not dry_run:
                    row[col] = new_val

    if not dry_run and changes:
        backup = samplesheet.with_name(samplesheet.name + ".bak")
        shutil.copy2(samplesheet, backup)
        with samplesheet.open("w", newline="", encoding="utf-8") as fh:
            # lineterminator="\n" matches how the wrapper first wrote the samplesheet
            # (samplesheet_builder._write_normalized_samplesheet); the csv default of
            # "\r\n" would otherwise flip every line ending on this remap step.
            writer = csv.DictWriter(fh, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    return changes


def remap_commands_references(
    commands_sh: Path,
    old_prefix: str,
    new_prefix: str,
    *,
    dry_run: bool,
) -> list[tuple[str, str, str]]:
    """Remap reference/index paths in commands.sh. Return (flag, old, new) triples."""
    content = commands_sh.read_text(encoding="utf-8")
    changes: list[tuple[str, str, str]] = []

    def _replace(m: re.Match) -> str:
        raw = m.group("value")
        value = _unquote_output_value(raw)
        if not _prefix_matches(value, old_prefix):
            return m.group(0)
        new_value = new_prefix + value[len(old_prefix):]
        flag = m.group(1).strip()
        changes.append((flag, value, new_value))
        return f"{m.group(1)}{_format_output_value(new_value, old_value=raw)}{m.group(3)}"

    updated = _REFERENCE_FLAG_RE.sub(_replace, content)

    if not dry_run and changes:
        backup = commands_sh.with_suffix(".sh.bak")
        shutil.copy2(commands_sh, backup)
        _write_text_lf(commands_sh, updated)

    return changes


def verify_paths(samplesheet: Path) -> list[str]:
    """Return FASTQ/BAM paths in the samplesheet that don't exist on disk.

    Remote URIs (s3://, https://, etc.) and bash variable references ($VAR/…
    or ${VAR}/…) are skipped — they cannot be stat-checked locally and are
    valid entries for cloud or replay-time execution.
    """
    missing: list[str] = []
    for row in csv.DictReader(samplesheet.read_text(encoding="utf-8").splitlines()):
        for col in _REMAPPABLE_COLUMNS:
            val = row.get(col, "")
            if not val or "://" in val or val.startswith("$"):
                continue
            if not Path(val).exists():
                missing.append(val)
    return missing


def verify_reference_paths(commands_sh: Path) -> list[str]:
    """Return local reference paths in commands.sh that don't exist on disk.

    Remote URIs (s3://, https://, etc.) and bash variable references are skipped
    because they cannot be verified with a simple filesystem check.
    """
    content = commands_sh.read_text(encoding="utf-8")
    missing: list[str] = []
    for m in _REFERENCE_FLAG_RE.finditer(content):
        value = _unquote_output_value(m.group("value"))
        if not value or value.startswith("$") or "://" in value:
            continue
        if not Path(value).exists():
            missing.append(value)
    return missing


def update_commands_output(commands_sh: Path, new_output_dir: str) -> bool:
    """Replace the --output value in commands.sh. Return True if changed."""
    original = commands_sh.read_text(encoding="utf-8")
    if not _OUTPUT_FLAG_RE.search(original):
        return False
    updated = _OUTPUT_FLAG_RE.sub(
        lambda m: f"{m.group(1)}{_format_output_value(new_output_dir, old_value=m.group('value'))}{m.group(3)}",
        original,
    )
    # Keep the idempotent-replay guard's manifest path in sync with the --output flag.
    updated = _RESUME_GUARD_RE.sub(
        lambda m: f"{m.group(1)}{new_output_dir}{m.group(3)}",
        updated,
    )
    if updated == original:
        return False
    backup = commands_sh.with_suffix(".sh.bak")
    shutil.copy2(commands_sh, backup)
    _write_text_lf(commands_sh, updated)
    return True


def cmd_remap(
    old_prefix: str,
    new_prefix: str,
    *,
    dry_run: bool,
    bundle_dir: Path | None = None,
) -> int:
    samplesheet = find_samplesheet(bundle_dir=bundle_dir)
    if samplesheet is None:
        print("ERROR: No samplesheet found in this bundle directory.", file=sys.stderr)
        return 1

    label = "[DRY RUN] " if dry_run else ""
    print(f"{label}Remapping FASTQ/BAM paths in: {samplesheet.name}")

    changes = remap_csv(samplesheet, old_prefix, new_prefix, dry_run=dry_run)

    if not changes:
        print(f"No FASTQ/BAM paths start with {old_prefix!r} — nothing to change.")
        # In a real run the exit code reflects replay-readiness, not whether THIS
        # invocation rewrote anything — otherwise an idempotent re-run (nothing left
        # to change) would exit 0 while the bundle is still not replay-ready, an
        # inconsistency that breaks scripted relocation. --dry-run is a preview and
        # never gates on readiness.
        return 0 if dry_run else _report_samplesheet_readiness(samplesheet)

    verb = "Would change" if dry_run else "Changed"
    print(f"\n{verb} {len(changes)} path(s):")
    for col, old_val, new_val in changes:
        print(f"  [{col}]")
        print(f"    - {old_val}")
        print(f"    + {new_val}")

    if dry_run:
        print("\nRe-run without --dry-run to apply these changes.")
        return 0

    print(f"\nBackup saved: {samplesheet.name + '.bak'}")
    return _report_samplesheet_readiness(samplesheet)


def _report_samplesheet_readiness(samplesheet: Path) -> int:
    """Verify every local FASTQ/BAM path resolves; return 0 only when replay-ready.

    Consistent across invocations (idempotent): the exit code depends solely on the
    samplesheet's current readiness, never on whether the calling command changed
    anything. Remote URIs and ``$VAR`` references are skipped by verify_paths.
    """
    missing = verify_paths(samplesheet)
    if missing:
        print(f"\nWARNING: {len(missing)} path(s) do not exist on this machine:")
        for m in missing:
            print(f"  {m}")
        print(
            "\nStage the data at these paths (or re-run with the correct --new prefix), "
            "then run `python3 remap_paths.py --verify` to confirm the bundle is replay-ready."
        )
        return 1

    print("\nAll FASTQ/BAM paths verified.")
    return 0


def cmd_remap_references(
    old_prefix: str,
    new_prefix: str,
    *,
    dry_run: bool,
    bundle_dir: Path | None = None,
) -> int:
    commands_sh = find_commands_sh(bundle_dir=bundle_dir)
    if commands_sh is None:
        print("ERROR: commands.sh not found in this bundle directory.", file=sys.stderr)
        return 1

    label = "[DRY RUN] " if dry_run else ""
    print(f"{label}Remapping reference paths in: commands.sh")

    changes = remap_commands_references(commands_sh, old_prefix, new_prefix, dry_run=dry_run)

    if not changes:
        print(f"No reference paths start with {old_prefix!r} — nothing to change.")
        return 0

    verb = "Would change" if dry_run else "Changed"
    print(f"\n{verb} {len(changes)} reference path(s):")
    for flag, old_val, new_val in changes:
        print(f"  [{flag}]")
        print(f"    - {old_val}")
        print(f"    + {new_val}")

    if dry_run:
        print("\nRe-run without --dry-run to apply these changes.")
        return 0

    print(f"\nBackup saved: {commands_sh.with_suffix('.sh.bak').name}")
    return 0


def cmd_update_output(new_output_dir: str, *, dry_run: bool, bundle_dir: Path | None = None) -> int:
    commands_sh = find_commands_sh(bundle_dir=bundle_dir)
    if commands_sh is None:
        print("ERROR: commands.sh not found in this bundle directory.", file=sys.stderr)
        return 1

    if dry_run:
        content = commands_sh.read_text(encoding="utf-8")
        m = _OUTPUT_FLAG_RE.search(content)
        if not m:
            print("No --output flag found in commands.sh — nothing to change.")
            return 0
        print("[DRY RUN] Would change --output in commands.sh:")
        print(f"    - {_unquote_output_value(m.group('value'))}")
        print(f"    + {new_output_dir}")
        return 0

    changed = update_commands_output(commands_sh, new_output_dir)
    if not changed:
        print("No --output flag found in commands.sh — nothing to change.")
        return 0

    print(f"Updated --output in commands.sh → {new_output_dir}")
    print(f"Backup saved: {commands_sh.with_suffix('.sh.bak').name}")
    return 0


def cmd_verify(bundle_dir: Path | None = None) -> int:
    ok = True
    warned = False

    samplesheet = find_samplesheet(bundle_dir=bundle_dir)
    if samplesheet is None:
        print("ERROR: No samplesheet found in this bundle directory.", file=sys.stderr)
        return 1

    missing_reads = verify_paths(samplesheet)
    if not missing_reads:
        print(f"FASTQ/BAM paths: all exist in {samplesheet.name}")
    else:
        ok = False
        print(f"FASTQ/BAM paths: {len(missing_reads)} missing in {samplesheet.name}:")
        for m in missing_reads:
            print(f"  {m}")
        print("  → fix: python3 remap_paths.py --old <old_prefix> --new <new_prefix>")

    commands_sh = find_commands_sh(bundle_dir=bundle_dir)
    if commands_sh is not None:
        content = commands_sh.read_text(encoding="utf-8")

        m = _OUTPUT_FLAG_RE.search(content)
        if m:
            output_path = _unquote_output_value(m.group("value"))
            out = Path(output_path)
            if out.exists():
                print(f"Output dir:  exists ({output_path})")
                # Non-blocking warning: if the dir has content outside reproducibility/
                # and --resume is absent, re-running would fail with "already exists" errors.
                has_content = any(
                    p for p in out.iterdir()
                    if p.name != "reproducibility"
                )
                has_resume = "--resume" in content
                if has_content and not has_resume:
                    warned = True
                    print(
                        "  WARNING: output dir is non-empty and commands.sh has no --resume.\n"
                        "  A fresh replay may fail because output already exists.\n"
                        "  Add --resume to replay into the existing run, or use a new --output path."
                    )
            else:
                # Not an error — the wrapper creates the output dir on a fresh replay.
                print(f"Output dir:  will be created on replay ({output_path})")
                print("  → to change it: python3 remap_paths.py --output-dir <new_output_dir>")

        missing_refs = verify_reference_paths(commands_sh)
        if not missing_refs:
            print("Reference paths: all exist in commands.sh")
        else:
            ok = False
            print(f"Reference paths: {len(missing_refs)} missing in commands.sh:")
            for r in missing_refs:
                print(f"  {r}")
            print("  → fix: python3 remap_paths.py --refs-old <old_prefix> --refs-new <new_prefix>")

    bd = bundle_dir or _BUNDLE_DIR
    missing_bundle = [f for f in _REQUIRED_BUNDLE_FILES if not (bd / f).exists()]
    if missing_bundle:
        ok = False  # incomplete bundle is a verification failure, not just a warning
        print(
            "\nBundle is incomplete — the following files are missing:\n"
            + "".join(f"  {bd / f}\n" for f in missing_bundle)
            + "  This usually means the wrapper crashed during post-processing.\n"
            + "  Run --repair-bundle to regenerate them:\n"
            + "    python3 remap_paths.py --repair-bundle"
        )

    if ok:
        replay_cmd = f"  CLAWBIO_REPO=/path/to/ClawBio bash {shlex.quote(str(bd / 'commands.sh'))}"
        if warned:
            print("\nAll checks passed (with warnings — review above before replaying):")
        else:
            print("\nAll checks passed — ready to replay:")
        print(replay_cmd)
    return 0 if ok else 1


def _sha256_file(path: Path) -> str:
    """Compute SHA-256 hex digest of a single file (stdlib-only, portable)."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _regenerate_checksums(bundle_dir: Path) -> None:
    """Write checksums.sha256 over all bundle content (stdlib-only, portable).

    Mirrors the logic of write_reproducibility_checksums in provenance.py but
    requires only stdlib so it can run on machines without ClawBio installed.

    bundle_dir is the reproducibility/ directory (where this script lives).
    Its parent is the output_dir that the wrapper used. Paths in the checksums
    file are relative to that parent so they are location-independent.
    """
    output_dir = bundle_dir.parent
    checksum_path = bundle_dir / "checksums.sha256"
    roots = [
        output_dir / "upstream" / "results",
        bundle_dir,  # the reproducibility/ dir itself
        output_dir / "provenance",
        output_dir / "logs",
    ]
    lines: list[str] = []
    for root in roots:
        if not root.exists():
            continue
        for f in sorted(root.rglob("*")):
            if not f.is_file():
                continue
            if f.name == "checksums.sha256":
                continue  # never hash the checksum file itself
            try:
                rel = f.relative_to(output_dir).as_posix()
            except ValueError:
                rel = f.name
            lines.append(f"{_sha256_file(f)}  {rel}")
    _write_text_lf(checksum_path, "\n".join(lines) + ("\n" if lines else ""))


def _regenerate_manifest_stub(bundle_dir: Path, commands_sh: Path) -> None:
    """Write a post-hoc manifest.json stub.

    The original manifest requires runtime metadata (args, params checksum,
    java/nextflow versions detected at run time) that cannot be reconstructed
    from the bundle alone.  The stub is clearly marked so audit consumers can
    distinguish it from an original.
    """
    manifest = {
        "regenerated_post_hoc": True,
        "note": (
            "Original manifest was lost because the wrapper crashed during "
            "post-processing. This stub was written by remap_paths.py --repair-bundle. "
            "For a complete audit bundle with original metadata, re-run the wrapper."
        ),
        "commands_sh": str(commands_sh),
        "regenerated_at": datetime.now(timezone.utc).isoformat(),
    }
    _write_text_lf(bundle_dir / "manifest.json", json.dumps(manifest, indent=2))


def _regenerate_environment_stub(bundle_dir: Path) -> None:
    """Write a post-hoc environment.yml stub.

    The original environment.yml captures exact Java/Nextflow versions detected
    at run time.  Without the wrapper context those cannot be reconstructed, so
    this stub records only a named placeholder that auditors can recognise.
    """
    _write_text_lf(
        bundle_dir / "environment.yml",
        "# regenerated post-hoc by remap_paths.py --repair-bundle\n"
        "# Original environment snapshot was lost because the wrapper crashed during\n"
        "# post-processing. For the original snapshot, re-run the wrapper.\n"
        "name: clawbio-nfcore-rnaseq-wrapper\n"
        "regenerated_post_hoc: true\n",
    )


def cmd_repair_bundle(bundle_dir: Path | None = None) -> int:
    """Regenerate missing bundle files from existing bundle contents.

    checksums.sha256 is recomputed from whatever files are present — fully
    accurate.  manifest.json and environment.yml are written as post-hoc stubs
    (marked regenerated_post_hoc: true) because the original runtime metadata
    cannot be reconstructed without re-running the wrapper.

    Returns 0 when all required files exist after repair, 1 on failure.
    Idempotent: re-running on a complete bundle is a no-op (returns 0).
    """
    bd = bundle_dir or _BUNDLE_DIR
    commands_sh = find_commands_sh(bundle_dir=bd)
    if commands_sh is None:
        print(
            "ERROR: commands.sh not found in bundle directory.\n"
            "  Cannot repair bundle without the original command record.",
            file=sys.stderr,
        )
        return 1

    missing = [f for f in _REQUIRED_BUNDLE_FILES if not (bd / f).exists()]
    if not missing:
        print("Bundle is complete - all required files are present. Nothing to repair.")
        return 0

    print(f"Repairing bundle - regenerating {len(missing)} missing file(s):")
    try:
        wrote_new_files = False
        if "manifest.json" in missing:
            _regenerate_manifest_stub(bd, commands_sh)
            print("  [OK] manifest.json  (post-hoc stub - original metadata unavailable)")
            wrote_new_files = True
        if "environment.yml" in missing:
            _regenerate_environment_stub(bd)
            print("  [OK] environment.yml  (post-hoc stub - original snapshot unavailable)")
            wrote_new_files = True
        # Always recompute checksums when any file was written (or checksums itself was
        # missing), so the checksum file covers the full post-repair bundle state.
        if wrote_new_files or "checksums.sha256" in missing:
            _regenerate_checksums(bd)
            print("  [OK] checksums.sha256  (recomputed from current bundle contents)")
    except Exception as exc:
        print(f"ERROR: Repair failed: {exc}", file=sys.stderr)
        return 1

    still_missing = [f for f in _REQUIRED_BUNDLE_FILES if not (bd / f).exists()]
    if still_missing:
        print(
            f"ERROR: Repair incomplete - still missing: {', '.join(still_missing)}",
            file=sys.stderr,
        )
        return 1

    stubs_regenerated = [f for f in ("manifest.json", "environment.yml") if f in missing]
    if stubs_regenerated:
        stub_list = " and ".join(stubs_regenerated)
        stub_word = "stubs" if len(stubs_regenerated) > 1 else "stub"
        are_or_is = "are" if len(stubs_regenerated) > 1 else "is a"
        print(
            f"\nRepair complete. Note: {stub_list} {are_or_is} post-hoc {stub_word}\n"
            "(marked regenerated_post_hoc: true). For a fully-original audit bundle,\n"
            "re-run the wrapper instead."
        )
    else:
        print("\nRepair complete.")
    return 0


def main(argv: list[str] | None = None, *, bundle_dir: Path | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Make this reproducibility bundle portable across machines.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  Remap FASTQ/BAM paths in the samplesheet:
    python3 remap_paths.py --old /Users/alice/fastqs --new /home/bob/fastqs

  Remap reference/index paths in commands.sh:
    python3 remap_paths.py --refs-old /Users/alice/refs --refs-new /home/bob/refs

  Update the --output directory in commands.sh:
    python3 remap_paths.py --output-dir /home/bob/my_run

  Preview any change without modifying files:
    python3 remap_paths.py --old /Users/alice/fastqs --new /home/bob/fastqs --dry-run

  Verify everything is ready to replay:
    python3 remap_paths.py --verify

  Regenerate missing bundle files (manifest.json, checksums.sha256, environment.yml):
    python3 remap_paths.py --repair-bundle

  Replay (CLAWBIO_REPO is always required):
    CLAWBIO_REPO=/path/to/ClawBio bash commands.sh
""",
    )
    parser.add_argument("--old", metavar="PREFIX", help="Original FASTQ/BAM path prefix to replace in samplesheet")
    parser.add_argument("--new", metavar="PREFIX", help="New FASTQ/BAM path prefix for this machine")
    parser.add_argument("--refs-old", metavar="PREFIX", help="Original reference/index path prefix to replace in commands.sh")
    parser.add_argument("--refs-new", metavar="PREFIX", help="New reference/index path prefix for this machine")
    parser.add_argument("--output-dir", metavar="PATH", help="New --output directory for commands.sh")
    parser.add_argument("--dry-run", action="store_true", help="Preview changes without modifying files")
    parser.add_argument("--verify", action="store_true", help="Check all paths exist on this machine")
    parser.add_argument("--repair-bundle", action="store_true", help="Regenerate missing bundle files (manifest.json, checksums.sha256, environment.yml)")
    args = parser.parse_args(argv)

    # A prefix remap needs BOTH halves. Accepting one silently (the old dispatch fell
    # through to print_help) hides a typo behind a "nothing to do" exit.
    if (args.old is None) != (args.new is None):
        parser.error("--old and --new must be given together")
    if (args.refs_old is None) != (args.refs_new is None):
        parser.error("--refs-old and --refs-new must be given together")

    requested = (
        args.repair_bundle
        or args.verify
        or args.output_dir is not None
        or args.old is not None
        or args.refs_old is not None
    )
    if not requested:
        parser.print_help()
        return 1

    # Every requested action runs, in the order the module docstring documents:
    # repair the bundle, apply the remaps, then verify LAST so --verify reports on the
    # bundle as it now stands. The previous dispatch was a chain of early returns with
    # --verify ahead of --output-dir, so `--output-dir X --verify` silently dropped the
    # rewrite and then printed "ready to replay" over a bundle still pointing at the old
    # output directory — a green light on a stale bundle. Actions compose; none is
    # discarded. The first non-zero step decides the exit code, but later steps still run
    # so the user sees every problem in one pass.
    exit_code = 0

    def _record(rc: int) -> None:
        nonlocal exit_code
        if exit_code == 0 and rc != 0:
            exit_code = rc

    if args.repair_bundle:
        _record(cmd_repair_bundle(bundle_dir=bundle_dir))
    if args.old is not None:
        _record(cmd_remap(args.old, args.new, dry_run=args.dry_run, bundle_dir=bundle_dir))
    if args.refs_old is not None:
        _record(cmd_remap_references(args.refs_old, args.refs_new, dry_run=args.dry_run, bundle_dir=bundle_dir))
    if args.output_dir is not None:
        _record(cmd_update_output(args.output_dir, dry_run=args.dry_run, bundle_dir=bundle_dir))
    if args.verify:
        _record(cmd_verify(bundle_dir=bundle_dir))
    return exit_code


def _format_output_value(new_output_dir: str, *, old_value: str) -> str:
    if old_value.startswith('"') and old_value.endswith('"'):
        return f'"{_escape_double_quoted(new_output_dir)}"'
    if old_value.startswith("'") and old_value.endswith("'"):
        if "'" not in new_output_dir:
            return f"'{new_output_dir}'"
        return shlex.quote(new_output_dir)
    if re.search(r"\s", new_output_dir):
        return f'"{_escape_double_quoted(new_output_dir)}"'
    if re.search(r"""[$'`\\]""", new_output_dir):
        return f'"{_escape_double_quoted(new_output_dir)}"'
    return new_output_dir


def _escape_double_quoted(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"').replace("$", "\\$").replace("`", "\\`")


def _unquote_output_value(value: str) -> str:
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        return value[1:-1]
    return value


if __name__ == "__main__":
    raise SystemExit(main())
