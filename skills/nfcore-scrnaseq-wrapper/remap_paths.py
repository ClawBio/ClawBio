#!/usr/bin/env python3
"""
remap_paths.py — Update absolute FASTQ paths in this reproducibility bundle.

FASTQ paths in the samplesheet are stored as absolute paths (required by
Nextflow). Before replaying this run on a different machine, update them:

    python remap_paths.py --old /original/data/dir --new /new/data/dir

Preview changes without modifying files:
    python remap_paths.py --old /original/data/dir --new /new/data/dir --dry-run

Verify all current paths exist on this machine:
    python remap_paths.py --verify
"""
from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

_BUNDLE_DIR = Path(__file__).resolve().parent
_FASTQ_COLUMNS = ("fastq_1", "fastq_2", "fastq_barcode")


def find_samplesheet(bundle_dir: Path | None = None) -> Path | None:
    search_dir = bundle_dir or _BUNDLE_DIR
    for name in ("samplesheet.valid.csv", "samplesheet.demo.csv"):
        p = search_dir / name
        if p.exists():
            return p
    return None


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
            if col in row and row[col] and row[col].startswith(old_prefix):
                new_val = new_prefix + row[col][len(old_prefix):]
                changes.append((col, row[col], new_val))
                if not dry_run:
                    row[col] = new_val

    if not dry_run and changes:
        backup = samplesheet.with_suffix(".bak")
        shutil.copy2(samplesheet, backup)
        with samplesheet.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    return changes


def verify_paths(samplesheet: Path) -> list[str]:
    """Return FASTQ paths in the samplesheet that don't exist on disk."""
    missing: list[str] = []
    for row in csv.DictReader(samplesheet.read_text(encoding="utf-8").splitlines()):
        for col in _FASTQ_COLUMNS:
            if col in row and row[col] and not Path(row[col]).exists():
                missing.append(row[col])
    return missing


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
        print("\nCorrect the paths and run again, or verify the FASTQ files are accessible.")
        return 1

    print("\nAll paths verified — ready to replay:")
    print(f"  bash {samplesheet.parent / 'commands.sh'}")
    return 0


def cmd_verify(bundle_dir: Path | None = None) -> int:
    samplesheet = find_samplesheet(bundle_dir=bundle_dir)
    if samplesheet is None:
        print("ERROR: No samplesheet found in this bundle directory.", file=sys.stderr)
        return 1

    missing = verify_paths(samplesheet)
    if not missing:
        print(f"All FASTQ paths in {samplesheet.name} exist on this machine.")
        print(f"Ready to replay: bash {samplesheet.parent / 'commands.sh'}")
        return 0

    print(f"Missing {len(missing)} FASTQ path(s) in {samplesheet.name}:")
    for m in missing:
        print(f"  {m}")
    print("\nTo remap paths: python remap_paths.py --old <old_prefix> --new <new_prefix>")
    return 1


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Update absolute FASTQ paths in the reproducibility bundle.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  Remap from the original machine to this one:
    python remap_paths.py --old /Users/alice/fastqs --new /home/bob/fastqs

  Preview changes (no files modified):
    python remap_paths.py --old /Users/alice/fastqs --new /home/bob/fastqs --dry-run

  Verify all paths exist on this machine:
    python remap_paths.py --verify
""",
    )
    parser.add_argument("--old", metavar="PREFIX", help="Original path prefix to replace")
    parser.add_argument("--new", metavar="PREFIX", help="New path prefix for this machine")
    parser.add_argument("--dry-run", action="store_true", help="Preview changes without modifying files")
    parser.add_argument("--verify", action="store_true", help="Check that all FASTQ paths exist")
    args = parser.parse_args()

    if args.verify:
        return cmd_verify()
    if args.old is not None and args.new is not None:
        return cmd_remap(args.old, args.new, dry_run=args.dry_run)
    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
