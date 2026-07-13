#!/usr/bin/env python3
"""
remap_paths.py — Make this reproducibility bundle portable across machines.

The replay command (commands.sh) self-anchors to this bundle's location, so the
output directory never needs patching. Only the machine-specific data paths —
FASTQ paths in the samplesheet and reference/index paths in params.yaml — need
remapping before replaying on a different machine:

  1. Update FASTQ paths in the samplesheet:
       python3 remap_paths.py --old /original/data/dir --new /new/data/dir

  2. Update reference/index paths in params.yaml (fasta, gtf, *_index, ...):
       python3 remap_paths.py --refs-old /old/refs --refs-new /new/refs

  3. Verify everything is ready:
       python3 remap_paths.py --verify

  4. Regenerate missing bundle files after a crash (manifest.json, checksums.sha256,
     environment.yml):
       python3 remap_paths.py --repair-bundle

  Preview any change without modifying files by adding --dry-run.
  Remote URIs (s3://, https://, ...) and $VAR references are always skipped.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

_BUNDLE_DIR = Path(__file__).resolve().parent
# Regenerable bundle files (crash recovery). commands.sh is the required input and
# is never regenerated. Matches nfcore-rnaseq/nfcore-sarek.
_REQUIRED_BUNDLE_FILES = ("manifest.json", "checksums.sha256", "environment.yml")


def _write_text_lf(path: Path, text: str) -> None:
    """Write ``text`` with LF line endings on every OS (self-contained).

    This script is shipped inside the reproducibility bundle and runs standalone at
    replay (no ClawBio package available), so it cannot import the shared
    ``clawbio.common.textio`` helper. CRLF/CR are normalised to LF and the result is
    written as bytes (bypassing text-mode newline translation) so a rewrite on
    Windows never reintroduces CRLF into checksums.sha256, manifest.json, or
    environment.yml. Mirrors the nfcore-rnaseq/nfcore-sarek helper.
    """
    normalised = text.replace("\r\n", "\n").replace("\r", "\n")
    Path(path).write_bytes(normalised.encode("utf-8"))


_FASTQ_COLUMNS = ("fastq_1", "fastq_2", "fastq_barcode")
# Reference path KEYS as they appear in params.yaml (underscore form).
# The block between the AUTO-GENERATED sentinels mirrors
# schemas.ALL_REFERENCE_PATH_FIELDS and is regenerated verbatim into every bundle
# copy by reporting._write_remap_script, so the standalone bundle can never drift
# from the canonical list even if this in-repo literal is edited. A parity test
# (test_params_reference_keys_match_params_builder) guards the in-repo copy too.
# >>> AUTO-GENERATED reference keys (do not edit by hand) >>>
_PARAMS_REFERENCE_KEYS = (
    "fasta",
    "gtf",
    "transcript_fasta",
    "txp2gene",
    "simpleaf_index",
    "kallisto_index",
    "star_index",
    "cellranger_index",
    "barcode_whitelist",
    "kb_t1c",
    "kb_t2c",
    "motifs",
    "cellrangerarc_config",
    "cellranger_vdj_index",
    "gex_frna_probe_set",
    "gex_target_panel",
    "gex_cmo_set",
    "fb_reference",
    "vdj_inner_enrichment_primers",
    "gex_barcode_sample_assignment",
    "cellranger_multi_barcodes",
)
# <<< AUTO-GENERATED reference keys <<<

# Matches a top-level `key: value` line in params.yaml for any reference key.
# Longest keys first so e.g. `cellranger_vdj_index` isn't shadowed by a prefix.
_PARAMS_REF_RE = re.compile(
    r"""^(?P<prefix>[ \t]*(?:"""
    + "|".join(
        re.escape(k) for k in sorted(_PARAMS_REFERENCE_KEYS, key=len, reverse=True)
    )
    + r""")[ \t]*:[ \t]+)(?P<q>["']?)(?P<value>.*?)(?P=q)(?P<trail>[ \t]*)$""",
    re.MULTILINE,
)


def find_samplesheet(bundle_dir: Path | None = None) -> Path | None:
    search_dir = bundle_dir or _BUNDLE_DIR
    for name in ("samplesheet.valid.csv", "samplesheet.demo.csv"):
        p = search_dir / name
        if p.exists():
            return p
    return None


def find_commands_sh(bundle_dir: Path | None = None) -> Path | None:
    search_dir = bundle_dir or _BUNDLE_DIR
    p = search_dir / "commands.sh"
    return p if p.exists() else None


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
        for col in _FASTQ_COLUMNS:
            # Boundary-aware match (same semantics as the params.yaml remapper):
            # a prefix of "/data/sample1" must NOT match "/data/sample10/...".
            if col in row and row[col] and _prefix_matches(row[col], old_prefix):
                new_val = new_prefix + row[col][len(old_prefix) :]
                changes.append((col, row[col], new_val))
                if not dry_run:
                    row[col] = new_val

    if not dry_run and changes:
        backup = samplesheet.with_suffix(".bak")
        shutil.copy2(samplesheet, backup)
        # lineterminator="\n": preserve the bundle's LF line endings (csv defaults
        # to CRLF) so remapping does not change the file's line-ending convention.
        with samplesheet.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    return changes


def verify_paths(samplesheet: Path) -> list[str]:
    """Return FASTQ paths in the samplesheet that don't exist on disk.

    Remote URIs (s3://, https://, ...) and $VAR references are skipped.
    """
    missing: list[str] = []
    for row in csv.DictReader(samplesheet.read_text(encoding="utf-8").splitlines()):
        for col in _FASTQ_COLUMNS:
            val = row.get(col, "")
            if not val or "://" in val or val.startswith("$"):
                continue
            if not Path(val).exists():
                missing.append(val)
    return missing


def _prefix_matches(path: str, prefix: str) -> bool:
    if path == prefix:
        return True
    sep = "" if prefix.endswith("/") else "/"
    return path.startswith(prefix + sep)


def find_params(bundle_dir: Path | None = None) -> Path | None:
    search_dir = bundle_dir or _BUNDLE_DIR
    p = search_dir / "params.yaml"
    return p if p.exists() else None


def remap_params_references(
    params_yaml: Path,
    old_prefix: str,
    new_prefix: str,
    *,
    dry_run: bool,
) -> list[tuple[str, str, str]]:
    """Rewrite reference/index paths stored in params.yaml. Return (key, old, new)."""
    content = params_yaml.read_text(encoding="utf-8")
    changes: list[tuple[str, str, str]] = []

    def _replace(m: "re.Match[str]") -> str:
        value = m.group("value")
        if "://" in value or not _prefix_matches(value, old_prefix):
            return m.group(0)
        new_value = new_prefix + value[len(old_prefix) :]
        key = m.group("prefix").strip().rstrip(":").strip()
        changes.append((key, value, new_value))
        return f"{m.group('prefix')}{m.group('q')}{new_value}{m.group('q')}{m.group('trail')}"

    updated = _PARAMS_REF_RE.sub(_replace, content)

    if not dry_run and changes:
        backup = params_yaml.with_suffix(".yaml.bak")
        shutil.copy2(params_yaml, backup)
        # write bytes with LF normalised: keep params.yaml byte-stable across OS
        # (Path.write_text would emit CRLF on Windows and change the checksum).
        params_yaml.write_bytes(
            updated.replace("\r\n", "\n").replace("\r", "\n").encode("utf-8")
        )

    return changes


def verify_params_references(params_yaml: Path) -> list[str]:
    """Return local reference paths in params.yaml that don't exist on disk.

    URIs, $VAR references, the ``false`` sentinel and glob patterns are skipped.
    """
    content = params_yaml.read_text(encoding="utf-8")
    missing: list[str] = []
    for m in _PARAMS_REF_RE.finditer(content):
        value = m.group("value").strip()
        if (
            not value
            or value.lower() in ("false", "true")
            or "://" in value
            or value.startswith("$")
            or any(ch in value for ch in "*{}")
        ):
            continue
        if not Path(value).exists():
            missing.append(value)
    return missing


def cmd_remap_references(
    old_prefix: str,
    new_prefix: str,
    *,
    dry_run: bool,
    bundle_dir: Path | None = None,
) -> int:
    params_yaml = find_params(bundle_dir=bundle_dir)
    if params_yaml is None:
        print("ERROR: params.yaml not found in this bundle directory.", file=sys.stderr)
        return 1

    label = "[DRY RUN] " if dry_run else ""
    print(f"{label}Remapping reference paths in: params.yaml")
    changes = remap_params_references(
        params_yaml, old_prefix, new_prefix, dry_run=dry_run
    )

    if not changes:
        print(f"No reference paths start with {old_prefix!r} — nothing to change.")
        return 0

    verb = "Would change" if dry_run else "Changed"
    print(f"\n{verb} {len(changes)} reference path(s):")
    for key, old_val, new_val in changes:
        print(f"  [{key}]")
        print(f"    - {old_val}")
        print(f"    + {new_val}")

    if dry_run:
        print("\nRe-run without --dry-run to apply these changes.")
        return 0

    print(f"\nBackup saved: {params_yaml.with_suffix('.yaml.bak').name}")
    missing = verify_params_references(params_yaml)
    if missing:
        print(
            f"\nWARNING: {len(missing)} reference path(s) do not exist on this machine:"
        )
        for m in missing:
            print(f"  {m}")
        print(
            "\nCheck the new prefix is correct, or run --verify after the files are in place."
        )
        # Signal failure so automated portability workflows can detect a broken
        # remap, matching cmd_remap (FASTQs) and cmd_verify (audit H-6).
        return 1
    print("\nAll reference paths verified.")
    return 0


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
    print(f"{label}Remapping FASTQ paths in: {samplesheet.name}")

    changes = remap_csv(samplesheet, old_prefix, new_prefix, dry_run=dry_run)

    if not changes:
        print(f"No FASTQ paths start with {old_prefix!r} — nothing to change.")
        return 0

    verb = "Would change" if dry_run else "Changed"
    print(f"\n{verb} {len(changes)} path(s):")
    for col, old_val, new_val in changes:
        print(f"  [{col}]")
        print(f"    - {old_val}")
        print(f"    + {new_val}")

    if dry_run:
        print("\nRe-run without --dry-run to apply these changes.")
        return 0

    print(f"\nBackup saved: {samplesheet.with_suffix('.bak').name}")

    missing = verify_paths(samplesheet)
    if missing:
        print(f"\nWARNING: {len(missing)} path(s) do not exist on this machine:")
        for m in missing:
            print(f"  {m}")
        print(
            "\nCorrect the paths and run again, or verify the FASTQ files are accessible."
        )
        return 1

    print("\nAll FASTQ paths verified — ready to replay:")
    print(f"  bash {samplesheet.parent / 'commands.sh'}")
    return 0


def cmd_verify(bundle_dir: Path | None = None) -> int:
    samplesheet = find_samplesheet(bundle_dir=bundle_dir)
    if samplesheet is None:
        print("ERROR: No samplesheet found in this bundle directory.", file=sys.stderr)
        return 1

    ok = True
    missing = verify_paths(samplesheet)
    if not missing:
        print(f"FASTQ paths: all exist in {samplesheet.name}")
    else:
        ok = False
        print(f"FASTQ paths: {len(missing)} missing in {samplesheet.name}:")
        for m in missing:
            print(f"  {m}")
        print("  → fix: python3 remap_paths.py --old <old_prefix> --new <new_prefix>")

    params_yaml = find_params(bundle_dir=bundle_dir)
    if params_yaml is not None:
        missing_refs = verify_params_references(params_yaml)
        if not missing_refs:
            print("Reference paths: all exist in params.yaml")
        else:
            ok = False
            print(f"Reference paths: {len(missing_refs)} missing in params.yaml:")
            for m in missing_refs:
                print(f"  {m}")
            print(
                "  → fix: python3 remap_paths.py --refs-old <old_prefix> --refs-new <new_prefix>"
            )

    if ok:
        print("\nAll checks passed — ready to replay:")
        print(f"  bash {(bundle_dir or _BUNDLE_DIR) / 'commands.sh'}")
    return 0 if ok else 1


def _sha256_file(path: Path) -> str:
    """Compute the SHA-256 hex digest of a single file (stdlib-only, portable)."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _find_h5ad_candidates(upstream_dir: Path) -> list[Path]:
    """Stdlib-only mirror of outputs_parser.find_h5ad_candidates.

    The canonical mtx_conversions matrices at both nesting depths, else a full-tree
    fallback. Ordered by string exactly as the wrapper orders them, so the resulting
    manifest matches a freshly generated one line-for-line.
    """
    if not upstream_dir.exists():
        return []
    canonical = sorted(
        {
            str(path)
            for pattern in ("*/mtx_conversions/*.h5ad", "*/mtx_conversions/*/*.h5ad")
            for path in upstream_dir.glob(pattern)
        }
    )
    if canonical:
        return [Path(p) for p in canonical]
    return [Path(p) for p in sorted(str(path) for path in upstream_dir.rglob("*.h5ad"))]


def _find_multiqc_report(upstream_dir: Path) -> Path | None:
    """Stdlib-only mirror of outputs_parser.find_multiqc_report (first sorted match)."""
    if not upstream_dir.exists():
        return None
    for pattern in ("multiqc/**/multiqc_report.html", "**/multiqc_report.html"):
        matches = sorted(upstream_dir.glob(pattern))
        if matches:
            return matches[0]
    return None


def _bundle_checksum_targets(output_dir: Path) -> list[Path]:
    """Reconstruct, from the bundle on disk, the exact allowlist that the wrapper's
    provenance.write_reproducibility_checksums hashes (stdlib-only).

    Membership and order mirror the wrapper: the normalized samplesheet, params.yaml,
    the stdout/stderr logs, the h5ad candidates, then the MultiQC report — each entry
    only when it exists inside output_dir, deduped preserving first-seen order. The
    reproducibility/*.json provenance tree, commands.sh, environment.yml, manifest.json,
    the copied policy files, and any other upstream/results or logs/ files are
    deliberately excluded, matching the original manifest. Because params.yaml is the
    only in-manifest file a later --refs-new remap mutates — exactly as in a fresh
    bundle — the regenerated manifest stays reproducible and self-consistent.
    """
    repro = output_dir / "reproducibility"
    upstream = output_dir / "upstream" / "results"
    logs = output_dir / "logs"

    # Only one samplesheet name exists per run (valid for a real run, demo for --demo);
    # listing both and dropping the absent one keeps the surviving order identical to
    # the wrapper's, which places the samplesheet first.
    ordered: list[Path] = [
        repro / "samplesheet.valid.csv",
        repro / "samplesheet.demo.csv",
        repro / "params.yaml",
        logs / "stdout.txt",
        logs / "stderr.txt",
    ]
    ordered.extend(_find_h5ad_candidates(upstream))
    multiqc = _find_multiqc_report(upstream)
    if multiqc is not None:
        ordered.append(multiqc)

    resolved_output_dir = output_dir.resolve()
    seen: set[Path] = set()
    targets: list[Path] = []
    for path in ordered:
        if path in seen:
            continue
        seen.add(path)
        if not path.is_file():
            continue
        try:
            path.resolve().relative_to(resolved_output_dir)
        except ValueError:
            continue  # defensive: never emit a bare-basename label that breaks -c
        targets.append(path)
    return targets


def _regenerate_checksums(bundle_dir: Path) -> None:
    """Write checksums.sha256 over the bundle's original checksum allowlist (stdlib-only).

    bundle_dir is the reproducibility/ directory (where this script lives); its parent
    is the output_dir the wrapper used. Labels are relative to output_dir so
    ``sha256sum -c`` succeeds when run from there, on any OS. The file-set matches what
    provenance.write_reproducibility_checksums originally wrote (see
    _bundle_checksum_targets), so a repaired manifest is byte-identical to a fresh one.
    checksums.sha256 never appears in the allowlist, so it never hashes itself.
    """
    output_dir = bundle_dir.parent
    checksum_path = bundle_dir / "checksums.sha256"
    lines = [
        f"{_sha256_file(f)}  {f.relative_to(output_dir).as_posix()}"
        for f in _bundle_checksum_targets(output_dir)
    ]
    _write_text_lf(checksum_path, "\n".join(lines) + ("\n" if lines else ""))


def _regenerate_manifest_stub(bundle_dir: Path, commands_sh: Path) -> None:
    """Write a post-hoc manifest.json stub, clearly marked so audit consumers can
    tell it from an original. The original manifest needs runtime metadata (args,
    params checksum, Java/Nextflow versions detected at run time) that cannot be
    reconstructed from the bundle alone."""
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
    """Write a post-hoc environment.yml stub. The original captures exact
    Java/Nextflow versions detected at run time, which cannot be reconstructed
    without the wrapper context, so this records only a named placeholder."""
    _write_text_lf(
        bundle_dir / "environment.yml",
        "# regenerated post-hoc by remap_paths.py --repair-bundle\n"
        "# Original environment snapshot was lost because the wrapper crashed during\n"
        "# post-processing. For the original snapshot, re-run the wrapper.\n"
        "name: clawbio-nfcore-scrnaseq-wrapper\n"
        "regenerated_post_hoc: true\n",
    )


def cmd_repair_bundle(bundle_dir: Path | None = None) -> int:
    """Regenerate missing bundle files from existing bundle contents.

    checksums.sha256 is recomputed from whatever files are present — fully accurate.
    manifest.json and environment.yml are written as post-hoc stubs (marked
    regenerated_post_hoc: true) because the original runtime metadata cannot be
    reconstructed without re-running the wrapper. Returns 0 when all required files
    exist after repair, 1 on failure. Idempotent: a complete bundle is a no-op.
    Parity with nfcore-rnaseq/nfcore-sarek.
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
        if "manifest.json" in missing:
            _regenerate_manifest_stub(bd, commands_sh)
            print("  [OK] manifest.json  (post-hoc stub - original metadata unavailable)")
        if "environment.yml" in missing:
            _regenerate_environment_stub(bd)
            print("  [OK] environment.yml  (post-hoc stub - original snapshot unavailable)")
        # checksums.sha256 lists only the wrapper's original allowlist (samplesheet,
        # params.yaml, stdout, stderr, h5ad candidates, MultiQC report); it never
        # covers the manifest.json / environment.yml stubs, so regenerating those
        # cannot change it. Recompute only when checksums.sha256 itself is missing —
        # this never silently overwrites an existing, valid manifest.
        if "checksums.sha256" in missing:
            _regenerate_checksums(bd)
            print("  [OK] checksums.sha256  (recomputed from the original checksum allowlist)")
    except Exception as exc:  # pragma: no cover - defensive
        print(f"ERROR: Repair failed: {exc}", file=sys.stderr)
        return 1

    still_missing = [f for f in _REQUIRED_BUNDLE_FILES if not (bd / f).exists()]
    if still_missing:
        print(
            f"ERROR: Repair incomplete - still missing: {', '.join(still_missing)}",
            file=sys.stderr,
        )
        return 1

    stubs = [f for f in ("manifest.json", "environment.yml") if f in missing]
    if stubs:
        stub_list = " and ".join(stubs)
        word = "stubs" if len(stubs) > 1 else "stub"
        are_is = "are" if len(stubs) > 1 else "is a"
        print(
            f"\nRepair complete. Note: {stub_list} {are_is} post-hoc {word}\n"
            "(marked regenerated_post_hoc: true). For a fully-original audit bundle,\n"
            "re-run the wrapper instead."
        )
    else:
        print("\nRepair complete.")
    return 0


def cmd_output_dir_hint(new_output_dir: str) -> int:
    """Handle ``--output-dir`` for CLI parity with the rnaseq/sarek wrappers.

    Those wrappers bake an absolute ``--output`` into ``commands.sh`` and expose
    ``--output-dir`` to rewrite it. The scrnaseq bundle instead self-relocates:
    ``commands.sh`` resolves its output directory from its own location and
    ``params.yaml`` stores output-relative paths, so moving the output directory needs
    no rewrite. This command accepts the flag and explains that, so a user reaching for
    ``--output-dir`` out of habit gets clear guidance instead of an argument error.
    """
    print(
        "The nfcore-scrnaseq-wrapper reproducibility bundle is self-relocating:\n"
        "  • commands.sh resolves its output directory from its own location, and\n"
        "  • params.yaml stores paths relative to that directory.\n"
        f"So there is nothing to rewrite for a new output directory ({new_output_dir}).\n"
        "Just move the whole output directory and run `bash reproducibility/commands.sh`.\n"
        "(If your FASTQ/reference paths also changed, use --old/--new and --refs-old/--refs-new.)"
    )
    return 0


def main(argv: list[str] | None = None, *, bundle_dir: Path | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Make this reproducibility bundle portable across machines.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  Remap FASTQ paths (required when FASTQs live at a different prefix):
    python3 remap_paths.py --old /Users/alice/fastqs --new /home/bob/fastqs

  Remap reference/index paths in params.yaml:
    python3 remap_paths.py --refs-old /Users/alice/refs --refs-new /home/bob/refs

  Preview any change without modifying files:
    python3 remap_paths.py --old /Users/alice/fastqs --new /home/bob/fastqs --dry-run

  Verify everything is ready to replay:
    python3 remap_paths.py --verify

  Regenerate missing bundle files (manifest.json, checksums.sha256, environment.yml):
    python3 remap_paths.py --repair-bundle
""",
    )
    parser.add_argument(
        "--old", metavar="PREFIX", help="Original FASTQ path prefix to replace"
    )
    parser.add_argument(
        "--new", metavar="PREFIX", help="New FASTQ path prefix for this machine"
    )
    parser.add_argument(
        "--refs-old",
        metavar="PREFIX",
        help="Original reference/index path prefix to replace in params.yaml",
    )
    parser.add_argument(
        "--refs-new",
        metavar="PREFIX",
        help="New reference/index path prefix for this machine",
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Preview changes without modifying files"
    )
    parser.add_argument(
        "--verify", action="store_true", help="Check all paths exist on this machine"
    )
    parser.add_argument(
        "--repair-bundle",
        action="store_true",
        help="Regenerate missing bundle files (manifest.json, checksums.sha256, environment.yml)",
    )
    parser.add_argument(
        "--output-dir",
        metavar="PATH",
        help=(
            "Accepted for parity with the rnaseq/sarek wrappers. This bundle "
            "self-relocates, so no rewrite is needed — just move the output directory."
        ),
    )
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

    # Every requested action runs, in the documented order: repair the bundle, apply the
    # remaps, show the --output-dir relocation hint, then verify LAST so --verify reports
    # on the bundle as it now stands. The previous dispatch was a chain of early returns
    # with --output-dir first, so `--output-dir X --verify` short-circuited on the hint and
    # silently discarded --verify: the user asked for a readiness check, got none, and saw
    # a zero exit code implying all was well. (rnaseq/sarek had the mirror image, dropping
    # --output-dir instead.) Actions compose; none is discarded. The first non-zero step
    # decides the exit code, but later steps still run so every problem shows in one pass.
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
        _record(cmd_output_dir_hint(args.output_dir))
    if args.verify:
        _record(cmd_verify(bundle_dir=bundle_dir))
    return exit_code
    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
