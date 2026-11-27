#!/usr/bin/env python3
"""
Download and prepare Corpas 30x WGS reference genome data from Zenodo.

This script manages the full lifecycle of the reference genome dataset:
downloading VCF files, generating lightweight subsets for demos and CI,
computing QC baselines for regression tests, and verifying checksums.

Requirements:
  - Python 3.10+
  - bcftools >= 1.17 (for --subsets and --baselines)
  - Internet access (for --download only)

Usage:
  python scripts/prepare_corpas_30x.py --download     # fetch VCFs from Zenodo
  python scripts/prepare_corpas_30x.py --subsets       # generate skill subsets
  python scripts/prepare_corpas_30x.py --baselines     # compute QC baselines
  python scripts/prepare_corpas_30x.py --verify        # check SHA-256 checksums
  python scripts/prepare_corpas_30x.py --all           # do everything

Zenodo record: https://zenodo.org/records/19297389
DOI: 10.5281/zenodo.19297389
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import shutil
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
CLAWBIO_DIR = SCRIPT_DIR.parent
CORPAS_DIR = CLAWBIO_DIR / "corpas-30x"
FULL_DIR = CORPAS_DIR / "full"
SUBSETS_DIR = CORPAS_DIR / "subsets"
BASELINES_DIR = CORPAS_DIR / "baselines"
REGIONS_DIR = CORPAS_DIR / "regions"
MANIFEST_PATH = CORPAS_DIR / "manifest.json"

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("prepare_corpas_30x")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

BCFTOOLS_MIN_VERSION = (1, 17)


def _sha256(path: Path) -> str:
    """Compute SHA-256 hex digest of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _md5(path: Path) -> str:
    """Compute MD5 hex digest of a file."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_manifest() -> dict:
    """Load manifest.json."""
    if not MANIFEST_PATH.exists():
        log.error("manifest.json not found at %s", MANIFEST_PATH)
        sys.exit(1)
    return json.loads(MANIFEST_PATH.read_text())


def _save_manifest(manifest: dict) -> None:
    """Write manifest.json with consistent formatting."""
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")


def _check_bcftools() -> bool:
    """Check that bcftools is installed and meets the minimum version."""
    bcftools = shutil.which("bcftools")
    if not bcftools:
        log.error(
            "bcftools not found on PATH. Install bcftools >= %d.%d",
            *BCFTOOLS_MIN_VERSION,
        )
        return False

    try:
        result = subprocess.run(
            ["bcftools", "--version"],
            capture_output=True,
            text=True,
            check=True,
        )
        first_line = result.stdout.strip().split("\n")[0]
        # e.g. "bcftools 1.20"
        version_str = first_line.split()[-1]
        parts = version_str.split(".")
        major, minor = int(parts[0]), int(parts[1]) if len(parts) > 1 else 0
        if (major, minor) < BCFTOOLS_MIN_VERSION:
            log.error(
                "bcftools %d.%d found, but >= %d.%d required",
                major,
                minor,
                *BCFTOOLS_MIN_VERSION,
            )
            return False
        log.info("bcftools %s found at %s", version_str, bcftools)
        return True
    except (subprocess.CalledProcessError, ValueError, IndexError) as exc:
        log.error("Failed to determine bcftools version: %s", exc)
        return False


def _run_bcftools(args: list[str], output_path: Path | None = None) -> str | None:
    """Run a bcftools command, optionally writing stdout to a file."""
    cmd = ["bcftools"] + args
    log.info("Running: %s", " ".join(cmd))
    if output_path:
        with open(output_path, "wb") as out:
            result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE)
    else:
        result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        stderr = (
            result.stderr if isinstance(result.stderr, str) else result.stderr.decode()
        )
        log.error("bcftools failed (exit %d): %s", result.returncode, stderr.strip())
        return None
    if not output_path:
        return result.stdout
    return ""


def _download_file(url: str, dest: Path, expected_md5: str | None = None) -> bool:
    """Download a file with progress reporting and optional MD5 verification.

    Tries the direct HTTPS URL first. If that fails, falls back to the Zenodo
    API content endpoint.
    """
    log.info("Downloading %s ...", dest.name)

    # Build fallback URL from direct URL
    fallback_url = None
    if "/files/" in url and "?download=" in url:
        filename = url.split("/files/")[1].split("?")[0]
        fallback_url = (
            f"https://zenodo.org/api/records/19297389/files/{filename}/content"
        )

    for attempt_url in [url, fallback_url]:
        if attempt_url is None:
            continue
        try:
            req = urllib.request.Request(attempt_url)
            with urllib.request.urlopen(req, timeout=300) as resp:
                total = int(resp.headers.get("Content-Length", 0))
                downloaded = 0
                with open(dest, "wb") as f:
                    while True:
                        chunk = resp.read(1 << 20)  # 1 MB
                        if not chunk:
                            break
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total > 0:
                            pct = downloaded * 100 // total
                            mb = downloaded / (1 << 20)
                            print(
                                f"\r  {dest.name}: {mb:.1f} MB ({pct}%)",
                                end="",
                                flush=True,
                            )
                print()  # newline after progress

            # Verify MD5 if provided
            if expected_md5:
                actual_md5 = _md5(dest)
                if actual_md5 != expected_md5:
                    log.error(
                        "MD5 mismatch for %s: expected %s, got %s",
                        dest.name,
                        expected_md5,
                        actual_md5,
                    )
                    dest.unlink(missing_ok=True)
                    return False
                log.info("MD5 verified: %s", dest.name)

            return True

        except (urllib.error.URLError, urllib.error.HTTPError, OSError) as exc:
            log.warning("Download failed from %s: %s", attempt_url, exc)
            dest.unlink(missing_ok=True)
            continue

    log.error("All download attempts failed for %s", dest.name)
    return False


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------


def cmd_download() -> int:
    """Download full VCFs from Zenodo."""
    manifest = _load_manifest()
    full_files = manifest["files"]["full"]

    FULL_DIR.mkdir(parents=True, exist_ok=True)

    success_count = 0
    for local_name, info in full_files.items():
        dest = FULL_DIR / local_name
        if dest.exists():
            log.info("Already present: %s (skipping)", local_name)
            info["status"] = "present"
            success_count += 1
            continue

        url = info.get("zenodo_download_url")
        if not url:
            log.warning("No download URL for %s, skipping", local_name)
            continue

        if _download_file(url, dest, expected_md5=info.get("md5")):
            info["sha256"] = _sha256(dest)
            info["status"] = "present"
            success_count += 1
        else:
            info["status"] = "download_required"

    _save_manifest(manifest)
    log.info("Downloaded %d/%d files", success_count, len(full_files))

    # Index all .vcf.gz files for bcftools region queries
    if _check_bcftools():
        for local_name in full_files:
            if local_name.endswith(".vcf.gz"):
                vcf_path = FULL_DIR / local_name
                if vcf_path.exists():
                    log.info("Indexing %s ...", local_name)
                    _run_bcftools(["index", "--force", str(vcf_path)])

    return 0


def _symlink_to_tmpdir(src: Path, tmpdir: Path) -> Path:
    """Create a symlink in tmpdir pointing to src (avoids spaces in paths)."""
    link = tmpdir / src.name
    link.symlink_to(src.resolve())
    return link


def cmd_subsets() -> int:
    """Generate lightweight VCF subsets for demos and CI.

    Uses a temp directory for all bcftools operations to avoid issues with
    spaces in file paths (macOS iCloud paths contain spaces).
    """
    if not _check_bcftools():
        return 1

    manifest = _load_manifest()

    # Check that full VCFs are present
    snps_vcf = FULL_DIR / "snps.vcf.gz"
    indels_vcf = FULL_DIR / "indels.vcf.gz"
    sv_vcf = FULL_DIR / "sv.pass.vcf.gz"
    cnv_vcf = FULL_DIR / "cnv.vcf.gz"

    if not snps_vcf.exists() or not indels_vcf.exists():
        log.info(
            "Full VCFs not found in %s. Run with --download first.",
            FULL_DIR,
        )
        return 0  # not an error, just not ready yet

    SUBSETS_DIR.mkdir(parents=True, exist_ok=True)
    errors = 0

    with tempfile.TemporaryDirectory(prefix="clawbio_") as tmpdir:
        tmp = Path(tmpdir)

        # Symlink source VCFs + indices into temp dir (no spaces)
        for vcf in [snps_vcf, indels_vcf, sv_vcf, cnv_vcf]:
            if vcf.exists():
                _symlink_to_tmpdir(vcf, tmp)
                # Also symlink the index if it exists
                for ext in [".csi", ".tbi"]:
                    idx = vcf.parent / (vcf.name + ext)
                    if idx.exists():
                        _symlink_to_tmpdir(idx, tmp)

        t_snps = tmp / "snps.vcf.gz"
        t_indels = tmp / "indels.vcf.gz"

        # Ensure indices exist
        if not (tmp / "snps.vcf.gz.csi").exists() and not (tmp / "snps.vcf.gz.tbi").exists():
            log.info("Indexing snps.vcf.gz ...")
            _run_bcftools(["index", str(t_snps)])
        if not (tmp / "indels.vcf.gz.csi").exists() and not (tmp / "indels.vcf.gz.tbi").exists():
            log.info("Indexing indels.vcf.gz ...")
            _run_bcftools(["index", str(t_indels)])

        # --- chr20 SNPs + indels ---
        chr20_out = SUBSETS_DIR / "chr20_snps_indels.vcf.gz"
        log.info("Generating chr20 subset ...")
        chr20_snps_tmp = tmp / "_chr20_snps.vcf.gz"
        chr20_indels_tmp = tmp / "_chr20_indels.vcf.gz"
        chr20_out_tmp = tmp / "chr20_snps_indels.vcf.gz"

        # Detect chromosome naming convention (chr20 vs 20)
        detect_result = _run_bcftools(["index", "--stats", str(t_snps)])
        uses_chr_prefix = False
        if detect_result and "chr20" in detect_result:
            uses_chr_prefix = True
        chr20_name = "chr20" if uses_chr_prefix else "20"
        log.info("Chromosome naming: %s (prefix=%s)", chr20_name, uses_chr_prefix)

        _run_bcftools(
            ["view", "-r", chr20_name, str(t_snps), "-O", "z", "-o", str(chr20_snps_tmp)]
        )
        _run_bcftools(
            ["view", "-r", chr20_name, str(t_indels), "-O", "z", "-o", str(chr20_indels_tmp)]
        )

        # Index and concat
        if chr20_snps_tmp.exists() and chr20_indels_tmp.exists():
            _run_bcftools(["index", str(chr20_snps_tmp)])
            _run_bcftools(["index", str(chr20_indels_tmp)])
            _run_bcftools(
                [
                    "concat",
                    "-a",
                    str(chr20_snps_tmp),
                    str(chr20_indels_tmp),
                    "-O",
                    "z",
                    "-o",
                    str(chr20_out_tmp),
                ]
            )
        elif chr20_snps_tmp.exists():
            chr20_snps_tmp.rename(chr20_out_tmp)

        if chr20_out_tmp.exists() and chr20_out_tmp.stat().st_size > 0:
            shutil.copy2(chr20_out_tmp, chr20_out)
            manifest["files"]["subsets"]["chr20_snps_indels.vcf.gz"]["sha256"] = _sha256(chr20_out)
            manifest["files"]["subsets"]["chr20_snps_indels.vcf.gz"]["status"] = "present"
            log.info("Created %s", chr20_out.name)
        else:
            log.error("Failed to create chr20 subset")
            errors += 1

        # --- Position-based subsets (PGx, NutriGx) ---
        # The WGS VCF has no rsID annotations (ID column is "."), so we
        # extract by genomic position. Coordinates come from the ClawBio
        # demo patient files (GRCh37). The VCF uses chr-prefixed names.
        pgx_demo = CLAWBIO_DIR / "skills" / "pharmgx-reporter" / "demo_patient.txt"
        nutrigx_demo = CLAWBIO_DIR / "skills" / "nutrigx" / "tests" / "synthetic_patient.csv"

        for subset_name, coord_source, source_type in [
            ("pgx_loci.vcf.gz", pgx_demo, "demo_tsv"),
            ("nutrigx_loci.vcf.gz", nutrigx_demo, "demo_tsv"),
        ]:
            out_path = SUBSETS_DIR / subset_name
            if not coord_source.exists():
                log.warning("Coordinate source not found: %s", coord_source)
                errors += 1
                continue

            # Build regions list (tab-delimited chr\tpos for bcftools -R)
            regions: list[str] = []
            if source_type == "demo_tsv":
                for line in coord_source.read_text().split("\n"):
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        chrom = parts[1] if parts[1].startswith("chr") else f"chr{parts[1]}"
                        pos = parts[2]
                        regions.append(f"{chrom}\t{pos}")
            if not regions:
                log.warning("No regions extracted for %s", subset_name)
                errors += 1
                continue

            log.info("Extracting %d positions for %s ...", len(regions), subset_name)

            # Write regions file for bcftools
            tmp_regions = tmp / f"{subset_name}.regions"
            tmp_regions.write_text("\n".join(regions) + "\n")
            tmp_out = tmp / subset_name

            _run_bcftools(
                [
                    "view",
                    "-R",
                    str(tmp_regions),
                    str(t_snps),
                    "-O",
                    "z",
                    "-o",
                    str(tmp_out),
                ]
            )

            if tmp_out.exists() and tmp_out.stat().st_size > 0:
                # Count extracted variants
                count_result = _run_bcftools(["view", "-H", str(tmp_out)])
                variant_count = 0
                if count_result:
                    variant_count = len([l for l in count_result.strip().split("\n") if l])
                log.info(
                    "%s: extracted %d variants from %d target positions",
                    subset_name,
                    variant_count,
                    len(regions),
                )
                if variant_count < len(regions):
                    log.warning(
                        "%s: %d/%d positions had no variant call (ref-homozygous or filtered)",
                        subset_name,
                        len(regions) - variant_count,
                        len(regions),
                    )

                shutil.copy2(tmp_out, out_path)
                manifest["files"]["subsets"][subset_name]["sha256"] = _sha256(out_path)
                manifest["files"]["subsets"][subset_name]["status"] = "present"
                log.info("Created %s", subset_name)
            else:
                log.error("Failed to create %s (empty or missing output)", subset_name)
                out_path.unlink(missing_ok=True)
                errors += 1

    # --- SV and CNV copies (already small, no bcftools needed) ---
    for src_name, dest_name in [
        ("sv.pass.vcf.gz", "sv_calls.vcf.gz"),
        ("cnv.vcf.gz", "cnv_calls.vcf.gz"),
    ]:
        src = FULL_DIR / src_name
        dest = SUBSETS_DIR / dest_name
        if src.exists():
            shutil.copy2(src, dest)
            manifest["files"]["subsets"][dest_name]["sha256"] = _sha256(dest)
            manifest["files"]["subsets"][dest_name]["status"] = "present"
            log.info("Copied %s -> %s", src_name, dest_name)
        else:
            log.info("%s not found, skipping %s", src_name, dest_name)

    _save_manifest(manifest)
    return 1 if errors > 0 else 0


def cmd_baselines() -> int:
    """Compute QC baselines for regression tests.

    Uses a temp directory for bcftools operations to handle spaces in paths.
    """
    if not _check_bcftools():
        return 1

    snps_vcf = FULL_DIR / "snps.vcf.gz"
    indels_vcf = FULL_DIR / "indels.vcf.gz"
    sv_vcf = FULL_DIR / "sv.pass.vcf.gz"
    cnv_vcf = FULL_DIR / "cnv.vcf.gz"

    if not snps_vcf.exists():
        log.info("Full VCFs not found. Run with --download first.")
        return 0

    BASELINES_DIR.mkdir(parents=True, exist_ok=True)
    baselines: dict = {
        "description": "Frozen QC metrics for Corpas 30x WGS regression tests",
        "genome_build": "GRCh37",
    }

    with tempfile.TemporaryDirectory(prefix="clawbio_bl_") as tmpdir:
        tmp = Path(tmpdir)

        # Symlink VCFs + indices into temp dir
        for vcf in [snps_vcf, indels_vcf, sv_vcf, cnv_vcf]:
            if vcf.exists():
                _symlink_to_tmpdir(vcf, tmp)
                for ext in [".csi", ".tbi"]:
                    idx = vcf.parent / (vcf.name + ext)
                    if idx.exists():
                        _symlink_to_tmpdir(idx, tmp)

        t_snps = tmp / "snps.vcf.gz"
        t_indels = tmp / "indels.vcf.gz"
        t_sv = tmp / "sv.pass.vcf.gz"
        t_cnv = tmp / "cnv.vcf.gz"

        # --- bcftools stats on SNPs ---
        log.info("Computing SNP stats ...")
        stats_output = _run_bcftools(["stats", str(t_snps)])
        if stats_output:
            baselines["snp_stats"] = _parse_bcftools_stats(stats_output)

        # --- Het/Hom ratio via genotype counting ---
        log.info("Computing Het/Hom ratio ...")
        gt_output = _run_bcftools(["query", "-f", "[%GT]\\n", str(t_snps)])
        if gt_output:
            het = 0
            hom_alt = 0
            for gt in gt_output.strip().split("\n"):
                gt = gt.strip()
                if gt in ("0/1", "0|1", "1|0"):
                    het += 1
                elif gt in ("1/1", "1|1"):
                    hom_alt += 1
            if hom_alt > 0:
                het_hom = round(het / hom_alt, 4)
                baselines.setdefault("snp_stats", {})["heterozygous"] = het
                baselines["snp_stats"]["homozygous_alt"] = hom_alt
                baselines["snp_stats"]["het_hom_ratio"] = het_hom
                log.info("Het=%d, Hom=%d, Het/Hom=%.4f", het, hom_alt, het_hom)

        # --- bcftools stats on indels ---
        if t_indels.exists():
            log.info("Computing indel stats ...")
            indel_stats = _run_bcftools(["stats", str(t_indels)])
            if indel_stats:
                baselines["indel_stats"] = _parse_bcftools_stats(indel_stats)

        # --- SV counts by type ---
        if t_sv.exists():
            log.info("Computing SV counts ...")
            sv_output = _run_bcftools(["query", "-f", "%INFO/SVTYPE\n", str(t_sv)])
            if sv_output:
                sv_counts: dict[str, int] = {}
                for line in sv_output.strip().split("\n"):
                    svtype = line.strip()
                    if svtype:
                        sv_counts[svtype] = sv_counts.get(svtype, 0) + 1
                baselines["sv_counts"] = sv_counts
                baselines["sv_total"] = sum(sv_counts.values())

        # --- CNV counts ---
        if t_cnv.exists():
            log.info("Computing CNV counts ...")
            cnv_output = _run_bcftools(["view", "-H", str(t_cnv)])
            if cnv_output:
                baselines["cnv_total"] = len(
                    [l for l in cnv_output.strip().split("\n") if l]
                )

    # Write baselines
    out_path = BASELINES_DIR / "qc_summary.json"
    out_path.write_text(json.dumps(baselines, indent=2, ensure_ascii=False) + "\n")
    log.info("Baselines written to %s", out_path)

    # Update manifest
    manifest = _load_manifest()
    manifest["files"]["baselines"]["qc_summary.json"]["sha256"] = _sha256(out_path)
    manifest["files"]["baselines"]["qc_summary.json"]["status"] = "present"
    _save_manifest(manifest)

    # Print summary
    if "snp_stats" in baselines:
        s = baselines["snp_stats"]
        log.info(
            "SNPs: %s total, Ti/Tv=%.3f, Het/Hom=%.3f",
            s.get("total_snps", "?"),
            s.get("ts_tv_ratio", 0),
            s.get("het_hom_ratio", 0),
        )
    if "sv_counts" in baselines:
        log.info("SV counts: %s", baselines["sv_counts"])

    return 0


def _parse_bcftools_stats(output: str) -> dict:
    """Parse key metrics from bcftools stats output."""
    result: dict = {}
    for line in output.split("\n"):
        # SN = summary numbers
        if line.startswith("SN\t"):
            parts = line.split("\t")
            if len(parts) >= 4:
                key = parts[2].strip().rstrip(":")
                value = parts[3].strip()
                if key == "number of SNPs":
                    result["total_snps"] = int(value)
                elif key == "number of indels":
                    result["total_indels"] = int(value)
                elif key == "number of records":
                    result["total_records"] = int(value)
                elif key == "number of samples":
                    result["samples"] = int(value)
        # TSTV = transition/transversion
        elif line.startswith("TSTV\t"):
            parts = line.split("\t")
            if len(parts) >= 5:
                try:
                    result["transitions"] = int(parts[2])
                    result["transversions"] = int(parts[3])
                    result["ts_tv_ratio"] = float(parts[4])
                except (ValueError, IndexError):
                    pass
        # PSC = per-sample counts (for het/hom)
        elif line.startswith("PSC\t"):
            parts = line.split("\t")
            if len(parts) >= 6:
                try:
                    # PSC fields: id, sample, hom, het, ...
                    hom = int(parts[4]) if parts[4].strip() else 0
                    het = int(parts[3]) if parts[3].strip() else 0
                    if hom > 0:
                        result["heterozygous"] = het
                        result["homozygous_alt"] = hom
                        result["het_hom_ratio"] = round(het / hom, 4)
                except (ValueError, IndexError, ZeroDivisionError):
                    pass
    return result


def cmd_verify() -> int:
    """Verify SHA-256 checksums against manifest."""
    manifest = _load_manifest()
    errors = 0
    checked = 0

    for section in ["full", "subsets", "baselines"]:
        files = manifest["files"].get(section, {})
        for filename, info in files.items():
            status = info.get("status", "missing")
            expected_sha = info.get("sha256")
            if status not in ("present", "verified") or not expected_sha:
                continue

            path = CORPAS_DIR / section / filename
            if not path.exists():
                log.error("MISSING: %s/%s (manifest says %s)", section, filename, status)
                errors += 1
                continue

            actual_sha = _sha256(path)
            checked += 1
            if actual_sha == expected_sha:
                info["status"] = "verified"
                log.info("OK: %s/%s", section, filename)
            else:
                log.error(
                    "MISMATCH: %s/%s expected=%s actual=%s",
                    section,
                    filename,
                    expected_sha[:16] + "...",
                    actual_sha[:16] + "...",
                )
                errors += 1

    _save_manifest(manifest)
    log.info("Verified %d files, %d errors", checked, errors)
    return 1 if errors > 0 else 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare Corpas 30x WGS reference genome data for ClawBio",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--download", action="store_true", help="Download full VCFs from Zenodo"
    )
    parser.add_argument(
        "--subsets",
        action="store_true",
        help="Generate lightweight VCF subsets (requires bcftools >= 1.17)",
    )
    parser.add_argument(
        "--baselines",
        action="store_true",
        help="Compute QC baselines for regression tests (requires bcftools >= 1.17)",
    )
    parser.add_argument(
        "--verify", action="store_true", help="Verify SHA-256 checksums"
    )
    parser.add_argument("--all", action="store_true", help="Run all steps in order")

    args = parser.parse_args()

    if not any([args.download, args.subsets, args.baselines, args.verify, args.all]):
        parser.print_help()
        return 0

    exit_code = 0

    if args.all or args.download:
        rc = cmd_download()
        exit_code = max(exit_code, rc)

    if args.all or args.subsets:
        rc = cmd_subsets()
        exit_code = max(exit_code, rc)

    if args.all or args.baselines:
        rc = cmd_baselines()
        exit_code = max(exit_code, rc)

    if args.all or args.verify:
        rc = cmd_verify()
        exit_code = max(exit_code, rc)

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
