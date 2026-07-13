#!/usr/bin/env python3
"""
remap_paths.py — Make this reproducibility bundle portable across machines.

Local samplesheet input paths (fastq_1, fastq_2, spring_1, spring_2, bam, bai,
cram, crai, vcf, table) are stored as absolute paths (required by Nextflow);
remote data URLs (s3://, https://, ...) are passed through unchanged and need no
remapping. Reference/index paths are stored in params.yaml, not in commands.sh.
Before replaying on a different machine:

  1. Remap input data paths in the samplesheet (if data moved):
       python3 remap_paths.py --old /original/data/dir --new /new/data/dir

  2. Remap reference/index paths (--fasta, --dbsnp, …). For sarek these live in
     reproducibility/params.yaml; this remaps params.yaml (and commands.sh, if
     references were added there manually):
       python3 remap_paths.py --refs-old /original/refs --refs-new /new/refs

  3. Update the --output path in commands.sh (if output dir changed):
       python3 remap_paths.py --output-dir /new/output/dir

  4. Verify everything is ready:
       python3 remap_paths.py --verify

  Preview any change without modifying files by adding --dry-run.

  Remote URIs (s3://, https://, etc.) are recognized and skipped automatically
  during path verification and checksumming — they are never resolved locally.

  Replay is self-contained — no environment variable is required:
       bash commands.sh
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
    """
    normalised = text.replace("\r\n", "\n").replace("\r", "\n")
    Path(path).write_bytes(normalised.encode("utf-8"))


_FASTQ_COLUMNS = ("fastq_1", "fastq_2")
# All samplesheet columns that hold local file paths for nf-core/sarek.
_SAMPLESHEET_PATH_COLUMNS = (
    "fastq_1", "fastq_2",
    "spring_1", "spring_2",
    "bam", "bai",
    "cram", "crai",
    "vcf",
    "table",
)
_REMAPPABLE_COLUMNS = _SAMPLESHEET_PATH_COLUMNS
_REQUIRED_BUNDLE_FILES = ("manifest.json", "checksums.sha256", "environment.yml")

# Reference flags from the nf-core/sarek 3.8.1 schema (hyphenated CLI form).
# NOTE: for sarek the wrapper always uses -params-file, so these flags do NOT
# appear on the nextflow command line in commands.sh — they live in params.yaml.
# This list exists for documentation and in case users add flags manually.
_REFERENCE_FLAGS = (
    "fasta", "fasta-fai", "dict",
    "bwa", "bwamem2", "dragmap",
    "dbsnp", "dbsnp-tbi",
    "known-indels", "known-indels-tbi",
    "known-snps", "known-snps-tbi",
    "germline-resource", "germline-resource-tbi",
    "pon", "pon-tbi",
    "intervals",
    "ascat-alleles", "ascat-loci", "ascat-loci-gc", "ascat-loci-rt",
    "chr-dir", "mappability",
    "msisensor2-models", "msisensorpro-scan",
    "ngscheckmate-bed",
    "sentieon-dnascope-model",
    "snpeff-cache", "vep-cache",
    "dbnsfp", "dbnsfp-tbi",
    "mastermind-file",
    "phenotypes-file", "phenotypes-file-tbi",
    "spliceai-snv", "spliceai-snv-tbi",
    "spliceai-indel", "spliceai-indel-tbi",
    "bcftools-annotations", "bcftools-annotations-tbi",
    "bcftools-header-lines", "bcftools-columns",
    "condel-config", "cnvkit-reference", "cf-chrom-len",
    "bbsplit-fasta-list", "bbsplit-index",
    "snpsift-databases",
    "multiqc-config", "multiqc-logo", "multiqc-methods-description",
    "igenomes-base",
    "pipeline-local", "nextflow-config",
)

# Reference path KEYS as they appear in params.yaml (underscore form). For
# nf-core/sarek the wrapper always uses -params-file, so these — not the
# commands.sh flags above — are what actually need remapping across machines.
# This set MUST stay in sync with provenance.REFERENCE_PATH_PARAMS (enforced by
# test_params_reference_keys_match_provenance); it is hardcoded here because this
# script is stdlib-only and runs on machines without ClawBio installed.
_PARAMS_REFERENCE_KEYS = (
    "fasta", "fasta_fai", "dict",
    "bwa", "bwamem2", "dragmap",
    "dbsnp", "dbsnp_tbi",
    "known_indels", "known_indels_tbi",
    "known_snps", "known_snps_tbi",
    "germline_resource", "germline_resource_tbi",
    "pon", "pon_tbi",
    "intervals",
    "ascat_alleles", "ascat_loci", "ascat_loci_gc", "ascat_loci_rt",
    "chr_dir", "mappability",
    "msisensor2_models", "msisensorpro_scan",
    "ngscheckmate_bed",
    "sentieon_dnascope_model",
    "snpeff_cache", "vep_cache",
    "dbnsfp", "dbnsfp_tbi",
    "mastermind_file",
    "phenotypes_file", "phenotypes_file_tbi",
    "spliceai_snv", "spliceai_snv_tbi", "spliceai_indel", "spliceai_indel_tbi",
    "bcftools_annotations", "bcftools_annotations_tbi", "bcftools_header_lines",
    "condel_config", "cnvkit_reference", "cf_chrom_len",
    "bbsplit_fasta_list", "bbsplit_index",
    "bcftools_columns", "snpsift_databases",
    "varlociraptor_scenario_tumor_only",
    "varlociraptor_scenario_somatic",
    "varlociraptor_scenario_germline",
    "multiqc_config", "multiqc_logo", "multiqc_methods_description",
    "igenomes_base",
)

# Matches a top-level `key: value` line in params.yaml for any reference key.
# Longest keys first so e.g. `dbsnp_tbi` is not shadowed by `dbsnp`. The value
# may be unquoted or single/double quoted; surrounding quotes are preserved.
_PARAMS_REF_RE = re.compile(
    r"""^(?P<prefix>[ \t]*(?:"""
    + "|".join(re.escape(k) for k in sorted(_PARAMS_REFERENCE_KEYS, key=len, reverse=True))
    + r""")[ \t]*:[ \t]+)(?P<q>["']?)(?P<value>.*?)(?P=q)(?P<trail>[ \t]*)$""",
    re.MULTILINE,
)

# Matches `    --output /path \` lines — requires line to start with whitespace
# (not #) so comment lines containing --output are never modified. The value
# may be unquoted, single-quoted, or double-quoted.
_OUTPUT_FLAG_RE = re.compile(
    r"""^([ \t]+--output[ \t]+)(?P<value>"[^"\n]*"|'[^'\n]*'|\S+)([ \t]*(?:\\[ \t]*)?)$""",
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


def find_params(bundle_dir: Path | None = None) -> Path | None:
    search_dir = bundle_dir or _BUNDLE_DIR
    p = search_dir / "params.yaml"
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
        # lineterminator="\n": csv defaults to CRLF; force LF so a remap preserves
        # the LF endings written by samplesheet_builder and keeps the bundle
        # byte-stable across operating systems.
        with samplesheet.open("w", newline="", encoding="utf-8") as fh:
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


def remap_params_references(
    params_yaml: Path,
    old_prefix: str,
    new_prefix: str,
    *,
    dry_run: bool,
) -> list[tuple[str, str, str]]:
    """Remap reference/index paths stored in params.yaml. Return (key, old, new).

    For nf-core/sarek the wrapper passes references via -params-file, so this is
    the function that makes a bundle fully portable. URIs (s3://, https://, …),
    the ``false`` disable sentinel, and non-path values are left untouched
    because they do not start with a filesystem prefix.
    """
    content = params_yaml.read_text(encoding="utf-8")
    changes: list[tuple[str, str, str]] = []

    def _replace(m: re.Match) -> str:
        value = m.group("value")
        if not _prefix_matches(value, old_prefix):
            return m.group(0)
        new_value = new_prefix + value[len(old_prefix):]
        key = m.group("prefix").strip().rstrip(":").strip()
        changes.append((key, value, new_value))
        return f"{m.group('prefix')}{m.group('q')}{new_value}{m.group('q')}{m.group('trail')}"

    updated = _PARAMS_REF_RE.sub(_replace, content)

    if not dry_run and changes:
        backup = params_yaml.with_suffix(".yaml.bak")
        shutil.copy2(params_yaml, backup)
        _write_text_lf(params_yaml, updated)

    return changes


def verify_params_references(params_yaml: Path) -> list[str]:
    """Return local reference paths in params.yaml that don't exist on disk.

    URIs, bash-variable references, the ``false`` sentinel and glob patterns
    (``*``/``{}``) are skipped — they cannot be stat-checked with a simple
    filesystem lookup.
    """
    content = params_yaml.read_text(encoding="utf-8")
    missing: list[str] = []
    for m in _PARAMS_REF_RE.finditer(content):
        value = m.group("value").strip()
        if (
            not value
            or value.lower() == "false"
            or "://" in value
            or value.startswith("$")
            or any(ch in value for ch in "*{}")
        ):
            continue
        if not Path(value).exists():
            missing.append(value)
    return missing


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

    print(f"\nBackup saved: {samplesheet.name + '.bak'}")

    missing = verify_paths(samplesheet)
    if missing:
        print(f"\nWARNING: {len(missing)} path(s) do not exist on this machine:")
        for m in missing:
            print(f"  {m}")
        print("\nCorrect the paths and run again, or verify the FASTQ files are accessible.")
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
    params_yaml = find_params(bundle_dir=bundle_dir)
    commands_sh = find_commands_sh(bundle_dir=bundle_dir)
    if params_yaml is None and commands_sh is None:
        print("ERROR: neither params.yaml nor commands.sh found in this bundle directory.", file=sys.stderr)
        return 1

    label = "[DRY RUN] " if dry_run else ""
    changes: list[tuple[str, str, str]] = []

    # params.yaml is where nf-core/sarek reference paths actually live.
    if params_yaml is not None:
        print(f"{label}Remapping reference paths in: params.yaml")
        changes.extend(remap_params_references(params_yaml, old_prefix, new_prefix, dry_run=dry_run))
    # commands.sh only carries references if a user added flags manually.
    if commands_sh is not None:
        print(f"{label}Remapping reference paths in: commands.sh")
        changes.extend(remap_commands_references(commands_sh, old_prefix, new_prefix, dry_run=dry_run))

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

    backups = []
    if params_yaml is not None and params_yaml.with_suffix(".yaml.bak").exists():
        backups.append(params_yaml.with_suffix(".yaml.bak").name)
    if commands_sh is not None and commands_sh.with_suffix(".sh.bak").exists():
        backups.append(commands_sh.with_suffix(".sh.bak").name)
    if backups:
        print(f"\nBackup saved: {', '.join(backups)}")

    if params_yaml is not None:
        missing = verify_params_references(params_yaml)
        if missing:
            print(f"\nWARNING: {len(missing)} reference path(s) do not exist on this machine:")
            for ref in missing:
                print(f"  {ref}")
            print("\nCorrect the paths and run again, or verify the reference files are accessible.")
            return 1
        print("\nAll reference paths verified.")
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

        match_obj = _OUTPUT_FLAG_RE.search(content)
        if match_obj:
            output_path = _unquote_output_value(match_obj.group("value"))
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
        if missing_refs:
            ok = False
            print(f"Reference paths (commands.sh): {len(missing_refs)} missing:")
            for r in missing_refs:
                print(f"  {r}")
            print("  → fix: python3 remap_paths.py --refs-old <old_prefix> --refs-new <new_prefix>")

    # nf-core/sarek reference paths live in params.yaml (passed via -params-file).
    params_yaml = find_params(bundle_dir=bundle_dir)
    if params_yaml is not None:
        missing_params_refs = verify_params_references(params_yaml)
        if not missing_params_refs:
            print("Reference paths: all exist in params.yaml")
        else:
            ok = False
            print(f"Reference paths: {len(missing_params_refs)} missing in params.yaml:")
            for r in missing_params_refs:
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
        replay_cmd = f"  bash {shlex.quote(str(bd / 'commands.sh'))}"
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

    Exclusions mirror provenance.py._iter_checksum_paths exactly: the work/,
    .nextflow/, reproducibility/ (which also holds logs/) and logs/ trees — matched
    at any depth, so the default nested work dir <output>/upstream/work is caught —
    plus any .log file, are never hashed. Keeping the two generators in lockstep
    means ``sha256sum -c`` still verifies after a remap.
    """
    output_dir = bundle_dir.parent
    checksum_path = bundle_dir / "checksums.sha256"
    exclude_dirs = {"work", ".nextflow", "reproducibility", "logs"}
    lines: list[str] = []
    for f in sorted(output_dir.rglob("*")):
        if not f.is_file():
            continue
        try:
            rel = f.relative_to(output_dir)
        except ValueError:
            continue
        # Match an excluded tree at any depth, not only the top level: the
        # Nextflow work dir defaults to <output>/upstream/work, so a parts[0]-only
        # check would leak that nested scratch tree. Mirrors
        # provenance.py._iter_checksum_paths exactly.
        if any(part in exclude_dirs for part in rel.parts[:-1]):
            continue
        if f.suffix == ".log":
            continue
        lines.append(f"{_sha256_file(f)}  {rel.as_posix()}")
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
        "name: claw-sarek\n"
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

  Remap reference/index paths in params.yaml (and commands.sh if present):
    python3 remap_paths.py --refs-old /Users/alice/refs --refs-new /home/bob/refs

  Update the --output directory in commands.sh:
    python3 remap_paths.py --output-dir /home/bob/my_run

  Preview any change without modifying files:
    python3 remap_paths.py --old /Users/alice/fastqs --new /home/bob/fastqs --dry-run

  Verify everything is ready to replay:
    python3 remap_paths.py --verify

  Replay the run (self-contained — no environment variable required):
    bash commands.sh
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

    # Every requested action runs, in the order the module docstring documents: repair the
    # bundle, apply the remaps, then verify LAST so --verify reports on the bundle as it
    # now stands. The previous dispatch was a chain of early returns with --verify ahead of
    # --output-dir, so `--output-dir X --verify` silently dropped the rewrite and then
    # reported the bundle ready to replay while commands.sh still pointed at the old output
    # directory — a green light on a stale bundle. Actions compose; none is discarded. The
    # first non-zero step decides the exit code, but later steps still run so the user sees
    # every problem in one pass.
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
