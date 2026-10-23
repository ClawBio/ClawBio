"""Standard-library relocation mechanics embedded into standalone bundle scripts."""

import csv
import hashlib
import json
import re
import shlex
import shutil
import sys
from collections.abc import Callable, Sequence
from datetime import datetime, timezone
from pathlib import Path


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


def sha256_file(path: Path) -> str:
    """SHA-256 hex digest of one file, streamed."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def prefix_matches(path: str, prefix: str) -> bool:
    if path == prefix:
        return True
    sep = "" if prefix.endswith("/") else "/"
    return path.startswith(prefix + sep)


def remap_prefix(path: str, old_prefix: str, new_prefix: str) -> str | None:
    """``path`` moved from under ``old_prefix`` to ``new_prefix``, or None if it is
    not under it. A trailing slash on either prefix makes no difference."""
    if not prefix_matches(path, old_prefix):
        return None
    old, new = old_prefix.rstrip("/"), new_prefix.rstrip("/")
    return (new + path[len(old):]) or "/"


def backup_once(path: Path, backup: Path) -> None:
    """Keep the file as the bundle produced it: a later relocation must not
    overwrite the first backup with already-relocated paths."""
    if not backup.exists():
        shutil.copy2(path, backup)


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
    with samplesheet.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    changes: list[tuple[str, str, str]] = []

    for row in rows:
        for col in columns:
            new_val = remap_prefix(row[col], old_prefix, new_prefix) if row.get(col) else None
            if new_val is not None:
                changes.append((col, row[col], new_val))
                if not dry_run:
                    row[col] = new_val

    if not dry_run and changes:
        backup = (
            backup_path
            if backup_path is not None
            else samplesheet.with_name(samplesheet.name + ".bak")
        )
        backup_once(samplesheet, backup)
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


# One bash word on a single line, as format_output_value writes it and as shlex.quote
# writes the original commands.sh: double-quoted runs (with backslash escapes),
# single-quoted runs, backslash escapes and plain characters, e.g. '/a b/it'"'"'s'.
SHELL_WORD = r"""(?:"(?:[^"\\\n]|\\.)*"|'[^'\n]*'|\\.|[^\s"'\\])+"""
# Characters shlex.quote leaves unquoted; anything else in a path needs quoting.
_SHELL_SAFE = re.compile(r"[\w@%+=:,./-]+", re.ASCII)


def escape_double_quoted(value: str) -> str:
    return (
        value.replace("\\", "\\\\")
        .replace('"', '\\"')
        .replace("$", "\\$")
        .replace("`", "\\`")
    )


def unquote_shell_word(word: str) -> str:
    """The literal string bash assigns to one SHELL_WORD, without expansions."""
    out: list[str] = []
    i = 0
    while i < len(word):
        char = word[i]
        if char == "'":
            end = word.index("'", i + 1)
            out.append(word[i + 1 : end])
            i = end + 1
        elif char == '"':
            i += 1
            while word[i] != '"':
                # Inside double quotes a backslash escapes only $ ` " and itself.
                if word[i] == "\\" and word[i + 1] in '$`"\\':
                    i += 1
                out.append(word[i])
                i += 1
            i += 1
        elif char == "\\":
            out.append(word[i + 1])
            i += 2
        else:
            out.append(char)
            i += 1
    return "".join(out)


def format_output_value(new_output_dir: str, *, old_value: str) -> str:
    """Quote a path for commands.sh, keeping the quoting style of the value it replaces.

    The result is always one SHELL_WORD that bash reads back as exactly
    ``new_output_dir``; an unquoted value stays bare only when no character needs
    quoting.
    """
    if "\n" in new_output_dir or "\r" in new_output_dir:
        raise ValueError("A path containing a line break cannot be written to commands.sh.")
    if old_value.startswith("'") and old_value.endswith("'") and len(old_value) > 1:
        if "'" not in new_output_dir:
            return f"'{new_output_dir}'"
        return shlex.quote(new_output_dir)
    if old_value.startswith('"') and old_value.endswith('"') and len(old_value) > 1:
        return f'"{escape_double_quoted(new_output_dir)}"'
    if _SHELL_SAFE.fullmatch(new_output_dir):
        return new_output_dir
    return f'"{escape_double_quoted(new_output_dir)}"'


# Absolute paths YAML reads back as the same plain string.
_YAML_PLAIN_SAFE = re.compile(r"/[\w@%+=,./~-]*", re.ASCII)
# YAML's double-quoted escapes (PyYAML reads all of these; yaml.dump writes \xXX).
_YAML_ESCAPES = {
    "0": "\0", "a": "\a", "b": "\b", "t": "\t", "\t": "\t", "n": "\n", "v": "\v",
    "f": "\f", "r": "\r", "e": "\x1b", " ": " ", '"': '"', "\\": "\\", "/": "/",
    "N": "\x85", "_": "\xa0", "L": "\u2028", "P": "\u2029",
}
_YAML_FOLDED_BREAK_RE = re.compile("[ \t]*([\x85\u2028\u2029])[ \t]*")
_YAML_ESCAPE_RE = re.compile(r"\\(x[0-9A-Fa-f]{2}|u[0-9A-Fa-f]{4}|U[0-9A-Fa-f]{8}|.)")


def _yaml_single_quotable(char: str) -> bool:
    # Printable and not a YAML line break (NEL, U+2028 and U+2029 fold like \n).
    code = ord(char)
    return (
        0x20 <= code <= 0x7E
        or (0xA0 <= code <= 0xD7FF and char not in "\u2028\u2029")
        or (0xE000 <= code <= 0xFFFD and char != "\ufeff")
        or code >= 0x10000
    )


def _yaml_escape(char: str) -> str:
    if char in '"\\':
        return "\\" + char
    if _yaml_single_quotable(char):
        return char
    code = ord(char)
    if code < 0x100:
        return f"\\x{code:02X}"
    return f"\\u{code:04X}" if code < 0x10000 else f"\\U{code:08X}"


def yaml_scalar_value(quote: str, raw: str) -> str:
    """The string a single-line YAML scalar holds, given its quote and inner text."""
    if quote != '"':
        # NEL, U+2028 and U+2029 are YAML line breaks that \n-based lines keep
        # inline; YAML drops the blanks around them (yaml.dump indents after one).
        raw = _YAML_FOLDED_BREAK_RE.sub(r"\1", raw)
    if quote == "'":
        return raw.replace("''", "'")
    if quote == '"':

        def unescape(m: re.Match[str]) -> str:
            code = m.group(1)
            if len(code) > 1:
                return chr(int(code[1:], 16))
            return _YAML_ESCAPES.get(code, m.group(0))

        return _YAML_ESCAPE_RE.sub(unescape, raw)
    return raw.strip(" \t")  # YAML trims only spaces and tabs, not U+00A0.


def format_yaml_scalar(value: str) -> str:
    """``value`` as a single-line YAML string that loads back unchanged.

    Plain when unambiguous, single-quoted when every character is printable, and
    double-quoted with escapes otherwise.
    """
    if _YAML_PLAIN_SAFE.fullmatch(value):
        return value
    if all(_yaml_single_quotable(char) for char in value):
        return "'" + value.replace("'", "''") + "'"
    return '"' + "".join(_yaml_escape(char) for char in value) + '"'


def remap_yaml_references(
    params_yaml: Path,
    old_prefix: str,
    new_prefix: str,
    *,
    dry_run: bool,
    pattern: re.Pattern[str],
    remap_uris: bool = True,
) -> list[tuple[str, str, str]]:
    """Rewrite reference paths matched by ``pattern`` in params.yaml.

    ``pattern`` must capture the key prefix, quote, value and trailing blanks of
    one ``key: value`` line. Return (key, old, new) for each changed value.
    """
    content = params_yaml.read_text(encoding="utf-8")
    changes: list[tuple[str, str, str]] = []

    def replace(m: re.Match[str]) -> str:
        value = yaml_scalar_value(m.group("q"), m.group("value"))
        new_value = None if not remap_uris and "://" in value else remap_prefix(value, old_prefix, new_prefix)
        if new_value is None:
            return m.group(0)
        key = m.group("prefix").strip().rstrip(":").strip()
        changes.append((key, value, new_value))
        return f"{m.group('prefix')}{format_yaml_scalar(new_value)}{m.group('trail')}"

    updated = pattern.sub(replace, content)
    if not dry_run and changes:
        backup_once(params_yaml, params_yaml.with_suffix(".yaml.bak"))
        write_text_lf(params_yaml, updated)
    return changes


def missing_yaml_references(params_yaml: Path, *, pattern: re.Pattern[str]) -> list[str]:
    """Local reference paths in params.yaml that do not exist on this machine.

    URIs, $VAR references, true/false sentinels and glob patterns are skipped:
    none of them can be checked with a filesystem lookup.
    """
    missing: list[str] = []
    for m in pattern.finditer(params_yaml.read_text(encoding="utf-8")):
        value = yaml_scalar_value(m.group("q"), m.group("value"))
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


def remap_flag_references(
    commands_sh: Path,
    old_prefix: str,
    new_prefix: str,
    *,
    dry_run: bool,
    pattern: re.Pattern[str],
) -> list[tuple[str, str, str]]:
    """Rewrite reference flags matched by ``pattern`` in commands.sh.

    ``pattern`` captures the flag (group 1), a SHELL_WORD named ``value`` and the
    line continuation (group 3). Return (flag, old, new) for each changed value.
    """
    content = commands_sh.read_text(encoding="utf-8")
    changes: list[tuple[str, str, str]] = []

    def replace(m: re.Match[str]) -> str:
        raw = m.group("value")
        value = unquote_shell_word(raw)
        new_value = remap_prefix(value, old_prefix, new_prefix)
        if new_value is None:
            return m.group(0)
        changes.append((m.group(1).strip(), value, new_value))
        return f"{m.group(1)}{format_output_value(new_value, old_value=raw)}{m.group(3)}"

    updated = pattern.sub(replace, content)
    if not dry_run and changes:
        backup_once(commands_sh, commands_sh.with_suffix(".sh.bak"))
        write_text_lf(commands_sh, updated)
    return changes


def missing_flag_references(commands_sh: Path, *, pattern: re.Pattern[str]) -> list[str]:
    """Local reference paths in commands.sh that do not exist on this machine.

    Remote URIs and bash variable references are skipped: neither can be checked
    with a filesystem lookup.
    """
    missing: list[str] = []
    for m in pattern.finditer(commands_sh.read_text(encoding="utf-8")):
        value = unquote_shell_word(m.group("value"))
        if not value or value.startswith("$") or "://" in value:
            continue
        if not Path(value).exists():
            missing.append(value)
    return missing


def missing_csv_paths(samplesheet: Path, *, columns: Sequence[str]) -> list[str]:
    """Local paths in ``columns`` of the samplesheet that do not exist on disk.

    Remote URIs and bash variable references ($VAR/... or ${VAR}/...) are skipped:
    they are valid for cloud or replay-time execution but cannot be stat-checked.
    """
    missing: list[str] = []
    with samplesheet.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for col in columns:
            value = row.get(col, "")
            if not value or "://" in value or value.startswith("$"):
                continue
            if not Path(value).exists():
                missing.append(value)
    return missing


def write_manifest_stub(bundle_dir: Path, commands_sh: Path) -> None:
    """Write a post-hoc manifest.json stub, marked so auditors can tell it apart.

    The original manifest records runtime metadata (arguments, params checksum,
    Java/Nextflow versions detected at run time) that the bundle alone cannot
    reconstruct.
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
    write_text_lf(bundle_dir / "manifest.json", json.dumps(manifest, indent=2))


def write_environment_stub(bundle_dir: Path, *, name: str) -> None:
    """Write a post-hoc environment.yml stub named like the wrapper's original.

    The original records the exact Java/Nextflow versions detected at run time,
    which the bundle alone cannot reconstruct.
    """
    write_text_lf(
        bundle_dir / "environment.yml",
        "# regenerated post-hoc by remap_paths.py --repair-bundle\n"
        "# Original environment snapshot was lost because the wrapper crashed during\n"
        "# post-processing. For the original snapshot, re-run the wrapper.\n"
        f"name: {name}\n"
        "regenerated_post_hoc: true\n",
    )


def repair_bundle(
    bundle_dir: Path,
    *,
    required_files: tuple[str, ...],
    commands_finder: Callable,
    regenerate_manifest: Callable,
    regenerate_environment: Callable,
    regenerate_checksums: Callable,
    checksums_cover_stubs: bool = True,
) -> int:
    """Regenerate missing bundle files from existing bundle contents.

    manifest.json and environment.yml are written as post-hoc stubs (marked
    regenerated_post_hoc: true) because the original runtime metadata cannot be
    reconstructed without re-running the wrapper. checksums.sha256 is recomputed
    when it is missing, and also after writing a stub when the checksum file
    covers the stubs; an allowlist that never lists them is left alone.

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
            print(
                "  [OK] manifest.json  (post-hoc stub - original metadata unavailable)"
            )
            wrote_new_files = True
        if "environment.yml" in missing:
            regenerate_environment(bd)
            print(
                "  [OK] environment.yml  (post-hoc stub - original snapshot unavailable)"
            )
            wrote_new_files = True
        if "checksums.sha256" in missing or (wrote_new_files and checksums_cover_stubs):
            regenerate_checksums(bd)
            print("  [OK] checksums.sha256  (recomputed with the bundle's checksum policy)")
    except Exception as exc:  # noqa: BLE001 - preserve best-effort repair diagnostics
        print(f"ERROR: Repair failed: {exc}", file=sys.stderr)
        return 1

    still_missing = [f for f in required_files if not (bd / f).exists()]
    if still_missing:
        print(
            f"ERROR: Repair incomplete - still missing: {', '.join(still_missing)}",
            file=sys.stderr,
        )
        return 1

    stubs_regenerated = [
        f for f in ("manifest.json", "environment.yml") if f in missing
    ]
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
