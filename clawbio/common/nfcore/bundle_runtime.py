"""Standard-library relocation mechanics embedded into standalone bundle scripts."""
import csv
import shutil
import re
import shlex
import sys
from pathlib import Path
from typing import Callable, Sequence


def write_text_lf(path: Path, text: str) -> None:
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

def prefix_matches(path: str, prefix: str) -> bool:
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
    columns: tuple[str, ...] | list[str],
    backup_path: Path | None = None,
) -> list[tuple[str, str, str]]:
    """Return list of (column, old_path, new_path) for every changed cell."""
    text = samplesheet.read_text(encoding="utf-8")
    fieldnames = list(csv.DictReader(text.splitlines()).fieldnames or [])
    rows = list(csv.DictReader(text.splitlines()))
    changes: list[tuple[str, str, str]] = []

    for row in rows:
        for col in columns:
            if col in row and row[col] and prefix_matches(row[col], old_prefix):
                new_val = new_prefix + row[col][len(old_prefix):]
                changes.append((col, row[col], new_val))
                if not dry_run:
                    row[col] = new_val

    if not dry_run and changes:
        backup = backup_path if backup_path is not None else samplesheet.with_name(samplesheet.name + ".bak")
        shutil.copy2(samplesheet, backup)
        with samplesheet.open("w", newline="", encoding="utf-8") as fh:
            # lineterminator="\n" matches how the wrapper first wrote the samplesheet
            # (samplesheet_builder._write_normalized_samplesheet); the csv default of
            # "\r\n" would otherwise flip every line ending on this remap step.
            writer = csv.DictWriter(fh, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    return changes

def run_actions(actions: Sequence[Callable[[], int]]) -> int:
    exit_code = 0
    for action in actions:
        rc = action()
        if exit_code == 0 and rc != 0:
            exit_code = rc
    return exit_code


def escape_double_quoted(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"').replace("$", "\\$").replace("`", "\\`")

def format_output_value(new_output_dir: str, *, old_value: str) -> str:
    if old_value.startswith('"') and old_value.endswith('"'):
        return f'"{escape_double_quoted(new_output_dir)}"'
    if old_value.startswith("'") and old_value.endswith("'"):
        if "'" not in new_output_dir:
            return f"'{new_output_dir}'"
        return shlex.quote(new_output_dir)
    if re.search(r"\s", new_output_dir):
        return f'"{escape_double_quoted(new_output_dir)}"'
    if re.search(r"""[$'`\\]""", new_output_dir):
        return f'"{escape_double_quoted(new_output_dir)}"'
    return new_output_dir

def repair_bundle(bundle_dir: Path, *, required_files: tuple[str, ...], commands_finder: Callable, regenerate_manifest: Callable, regenerate_environment: Callable, regenerate_checksums: Callable) -> int:
    """Regenerate missing bundle files from existing bundle contents.

    checksums.sha256 is recomputed from whatever files are present — fully
    accurate.  manifest.json and environment.yml are written as post-hoc stubs
    (marked regenerated_post_hoc: true) because the original runtime metadata
    cannot be reconstructed without re-running the wrapper.

    Returns 0 when all required files exist after repair, 1 on failure.
    Idempotent: re-running on a complete bundle is a no-op (returns 0).
    """
    bd = bundle_dir
    commands_sh = commands_finder(bundle_dir=bd)
    if commands_sh is None:
        print(
            "ERROR: commands.sh not found in bundle directory.\n"
            "  Cannot repair bundle without the original command record.",
            file=sys.stderr,
        )
        return 1

    missing = [f for f in required_files if not (bd / f).exists()]
    if not missing:
        print("Bundle is complete - all required files are present. Nothing to repair.")
        return 0

    print(f"Repairing bundle - regenerating {len(missing)} missing file(s):")
    try:
        wrote_new_files = False
        if "manifest.json" in missing:
            regenerate_manifest(bd, commands_sh)
            print("  [OK] manifest.json  (post-hoc stub - original metadata unavailable)")
            wrote_new_files = True
        if "environment.yml" in missing:
            regenerate_environment(bd)
            print("  [OK] environment.yml  (post-hoc stub - original snapshot unavailable)")
            wrote_new_files = True
        # Always recompute checksums when any file was written (or checksums itself was
        # missing), so the checksum file covers the full post-repair bundle state.
        if wrote_new_files or "checksums.sha256" in missing:
            regenerate_checksums(bd)
            print("  [OK] checksums.sha256  (recomputed from current bundle contents)")
    except Exception as exc:
        print(f"ERROR: Repair failed: {exc}", file=sys.stderr)
        return 1

    still_missing = [f for f in required_files if not (bd / f).exists()]
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
