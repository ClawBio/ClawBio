#!/usr/bin/env python3
"""
remap_paths.py — Make this reproducibility bundle portable across machines.

FASTQ paths and the --output directory in commands.sh are stored as absolute
paths (required by Nextflow). Before replaying on a different machine:

  1. Update FASTQ paths in the samplesheet:
       python remap_paths.py --old /original/data/dir --new /new/data/dir

  2. Update the --output path in commands.sh (if output dir changed):
       python remap_paths.py --output-dir /new/output/dir

  3. Verify everything is ready:
       python remap_paths.py --verify

  Preview any change without modifying files by adding --dry-run.

  On a machine where ClawBio is installed at a non-standard location, set:
       CLAWBIO_REPO=/path/to/ClawBio bash commands.sh
"""
from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
from pathlib import Path

_BUNDLE_DIR = Path(__file__).resolve().parent
_FASTQ_COLUMNS = ("fastq_1", "fastq_2", "fastq_barcode")
# Matches `    --output /path \` lines — requires line to start with whitespace
# (not #) so comment lines containing --output are never modified.
_OUTPUT_FLAG_RE = re.compile(r"^([ \t]+--output[ \t]+)(\S+)([ \t]*(?:\\[ \t]*)?)$", re.MULTILINE)


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


def update_commands_output(commands_sh: Path, new_output_dir: str) -> bool:
    """Replace the --output value in commands.sh. Return True if changed."""
    original = commands_sh.read_text(encoding="utf-8")
    if not _OUTPUT_FLAG_RE.search(original):
        return False
    updated = _OUTPUT_FLAG_RE.sub(
        lambda m: f"{m.group(1)}{new_output_dir}{m.group(3)}",
        original,
    )
    if updated == original:
        return False
    backup = commands_sh.with_suffix(".sh.bak")
    shutil.copy2(commands_sh, backup)
    commands_sh.write_text(updated, encoding="utf-8")
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

    print("\nAll FASTQ paths verified — ready to replay:")
    print(f"  bash {samplesheet.parent / 'commands.sh'}")
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
        print(f"[DRY RUN] Would change --output in commands.sh:")
        print(f"    - {m.group(2)}")
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
        print("  → fix: python remap_paths.py --old <old_prefix> --new <new_prefix>")

    commands_sh = find_commands_sh(bundle_dir=bundle_dir)
    if commands_sh is not None:
        content = commands_sh.read_text(encoding="utf-8")
        m = _OUTPUT_FLAG_RE.search(content)
        if m:
            output_path = m.group(2)
            if Path(output_path).exists():
                print(f"Output dir:  exists ({output_path})")
            else:
                ok = False
                print(f"Output dir:  missing ({output_path})")
                print("  → fix: python remap_paths.py --output-dir <new_output_dir>")

    if ok:
        print(f"\nAll checks passed — ready to replay:")
        print(f"  bash {(bundle_dir or _BUNDLE_DIR) / 'commands.sh'}")
    return 0 if ok else 1


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Make this reproducibility bundle portable across machines.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  Remap FASTQ paths (required when FASTQs live at a different prefix):
    python remap_paths.py --old /Users/alice/fastqs --new /home/bob/fastqs

  Update the --output directory in commands.sh:
    python remap_paths.py --output-dir /home/bob/my_run

  Preview any change without modifying files:
    python remap_paths.py --old /Users/alice/fastqs --new /home/bob/fastqs --dry-run

  Verify everything is ready to replay:
    python remap_paths.py --verify

  If ClawBio is installed at a non-standard path on this machine:
    CLAWBIO_REPO=/path/to/ClawBio bash commands.sh
""",
    )
    parser.add_argument("--old", metavar="PREFIX", help="Original FASTQ path prefix to replace")
    parser.add_argument("--new", metavar="PREFIX", help="New FASTQ path prefix for this machine")
    parser.add_argument("--output-dir", metavar="PATH", help="New --output directory for commands.sh")
    parser.add_argument("--dry-run", action="store_true", help="Preview changes without modifying files")
    parser.add_argument("--verify", action="store_true", help="Check all paths exist on this machine")
    args = parser.parse_args()

    if args.verify:
        return cmd_verify()
    if args.output_dir is not None:
        return cmd_update_output(args.output_dir, dry_run=args.dry_run)
    if args.old is not None and args.new is not None:
        return cmd_remap(args.old, args.new, dry_run=args.dry_run)
    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
