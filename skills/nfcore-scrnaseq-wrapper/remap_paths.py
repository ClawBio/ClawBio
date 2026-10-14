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
import re
import sys
from pathlib import Path

_BUNDLE_DIR = Path(__file__).resolve().parent

# >>> SHARED NFCORE RUNTIME >>>
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))
from clawbio.common.nfcore.bundle_runtime import (
    missing_csv_paths as _nfcore_missing_csv_paths,
    missing_yaml_references as _nfcore_missing_yaml_references,
    remap_csv as _nfcore_remap_csv,
    remap_yaml_references as _nfcore_remap_yaml_references,
    repair_bundle as _nfcore_repair_bundle,
    run_actions as _nfcore_run_actions,
    sha256_file as _nfcore_sha256_file,
    write_environment_stub as _nfcore_write_environment_stub,
    write_manifest_stub as _nfcore_write_manifest_stub,
    write_text_lf as _nfcore_write_text_lf,
)
# <<< SHARED NFCORE RUNTIME <<<

# Regenerable bundle files (crash recovery). commands.sh is the required input and
# is never regenerated. Matches nfcore-rnaseq/nfcore-sarek.
_REQUIRED_BUNDLE_FILES = ("manifest.json", "checksums.sha256", "environment.yml")


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
    return _nfcore_remap_csv(
        samplesheet,
        old_prefix,
        new_prefix,
        columns=_FASTQ_COLUMNS,
        dry_run=dry_run,
        backup_path=samplesheet.with_suffix(".bak"),
    )


def verify_paths(samplesheet: Path) -> list[str]:
    """Return FASTQ paths in the samplesheet that don't exist on disk."""
    return _nfcore_missing_csv_paths(samplesheet, columns=_FASTQ_COLUMNS)


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
    """Rewrite reference/index paths stored in params.yaml. Return (key, old, new).

    Remote URIs are never rewritten: remote references are gated behind
    --allow-remote-inputs and staged by Nextflow.
    """
    return _nfcore_remap_yaml_references(
        params_yaml,
        old_prefix,
        new_prefix,
        dry_run=dry_run,
        pattern=_PARAMS_REF_RE,
        remap_uris=False,
    )


def verify_params_references(params_yaml: Path) -> list[str]:
    """Return local reference paths in params.yaml that don't exist on disk."""
    return _nfcore_missing_yaml_references(params_yaml, pattern=_PARAMS_REF_RE)


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
        f"{_nfcore_sha256_file(f)}  {f.relative_to(output_dir).as_posix()}"
        for f in _bundle_checksum_targets(output_dir)
    ]
    _nfcore_write_text_lf(checksum_path, "\n".join(lines) + ("\n" if lines else ""))


def _regenerate_environment_stub(bundle_dir: Path) -> None:
    _nfcore_write_environment_stub(bundle_dir, name="clawbio-nfcore-scrnaseq-wrapper")


def cmd_repair_bundle(bundle_dir: Path | None = None) -> int:
    # checksums.sha256 lists only the wrapper's original allowlist (samplesheet,
    # params.yaml, stdout, stderr, h5ad candidates, MultiQC report), never the
    # manifest/environment stubs, so writing a stub cannot make it stale.
    return _nfcore_repair_bundle(
        bundle_dir or _BUNDLE_DIR,
        required_files=_REQUIRED_BUNDLE_FILES,
        commands_finder=find_commands_sh,
        regenerate_manifest=_nfcore_write_manifest_stub,
        regenerate_environment=_regenerate_environment_stub,
        regenerate_checksums=_regenerate_checksums,
        checksums_cover_stubs=False,
    )


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
    actions = []

    if args.repair_bundle:
        actions.append(lambda: cmd_repair_bundle(bundle_dir=bundle_dir))
    if args.old is not None:
        actions.append(lambda: cmd_remap(args.old, args.new, dry_run=args.dry_run, bundle_dir=bundle_dir))
    if args.refs_old is not None:
        actions.append(lambda: cmd_remap_references(args.refs_old, args.refs_new, dry_run=args.dry_run, bundle_dir=bundle_dir))
    if args.output_dir is not None:
        actions.append(lambda: cmd_output_dir_hint(args.output_dir))
    if args.verify:
        actions.append(lambda: cmd_verify(bundle_dir=bundle_dir))
    return _nfcore_run_actions(actions)


if __name__ == "__main__":
    raise SystemExit(main())
