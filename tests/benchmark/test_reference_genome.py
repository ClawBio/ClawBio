"""
Regression benchmarks for the Corpas 30x WGS reference genome.

Structural tests (manifest, citation, rsID lists) always run.
Data-dependent tests skip gracefully when files are absent, so CI
stays green on fresh clones and becomes stronger once data is present.

Run with:
    python -m pytest tests/benchmark/test_reference_genome.py -v
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

CLAWBIO_DIR = Path(__file__).resolve().parent.parent.parent
CORPAS_DIR = CLAWBIO_DIR / "corpas-30x"
MANIFEST_PATH = CORPAS_DIR / "manifest.json"
BASELINES_DIR = CORPAS_DIR / "baselines"
SUBSETS_DIR = CORPAS_DIR / "subsets"
REGIONS_DIR = CORPAS_DIR / "regions"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load_manifest() -> dict:
    assert MANIFEST_PATH.exists(), f"manifest.json not found at {MANIFEST_PATH}"
    return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))


def _file_is_present(manifest: dict, section: str, filename: str) -> Path | None:
    """Return the file path if it exists and manifest says present/verified, else None."""
    info = manifest.get("files", {}).get(section, {}).get(filename, {})
    status = info.get("status", "missing")
    if status not in ("present", "verified"):
        return None
    path = CORPAS_DIR / section / filename
    if not path.exists():
        return None
    return path


# ---------------------------------------------------------------------------
# Structural tests (always run)
# ---------------------------------------------------------------------------


class TestManifest:
    """Tests for manifest.json structure and contents."""

    def test_manifest_exists(self):
        assert MANIFEST_PATH.exists(), "corpas-30x/manifest.json must exist"

    def test_manifest_valid_json(self):
        manifest = _load_manifest()
        assert isinstance(manifest, dict)

    def test_manifest_has_required_keys(self):
        manifest = _load_manifest()
        for key in [
            "version",
            "genome_build",
            "source_doi",
            "all_versions_doi",
            "license",
            "files",
        ]:
            assert key in manifest, f"manifest.json missing key: {key}"

    def test_manifest_doi_values(self):
        manifest = _load_manifest()
        assert manifest["source_doi"] == "10.5281/zenodo.19297389"
        assert manifest["all_versions_doi"] == "10.5281/zenodo.19285820"

    def test_manifest_file_sections(self):
        manifest = _load_manifest()
        files = manifest["files"]
        for section in ["full", "subsets", "baselines"]:
            assert section in files, f"manifest.json missing files section: {section}"

    def test_manifest_full_files_have_download_urls(self):
        manifest = _load_manifest()
        for filename, info in manifest["files"]["full"].items():
            assert "zenodo_download_url" in info, (
                f"full/{filename} missing zenodo_download_url"
            )
            assert "status" in info, f"full/{filename} missing status"

    def test_manifest_all_files_have_status(self):
        manifest = _load_manifest()
        valid_statuses = {
            "download_required",
            "pending_generation",
            "present",
            "verified",
        }
        for section in ["full", "subsets", "baselines"]:
            for filename, info in manifest["files"][section].items():
                status = info.get("status")
                assert status in valid_statuses, (
                    f"{section}/{filename} has invalid status: {status}"
                )


class TestCitation:
    """Tests for CITATION.cff structure."""

    def test_citation_exists(self):
        path = CORPAS_DIR / "CITATION.cff"
        assert path.exists(), "corpas-30x/CITATION.cff must exist"

    def test_citation_valid_yaml(self):
        try:
            import yaml
        except ImportError:
            pytest.skip("pyyaml not installed")
        path = CORPAS_DIR / "CITATION.cff"
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
        assert isinstance(data, dict)

    def test_citation_has_doi(self):
        try:
            import yaml
        except ImportError:
            pytest.skip("pyyaml not installed")
        path = CORPAS_DIR / "CITATION.cff"
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
        assert data.get("doi") == "10.5281/zenodo.19297389"

    def test_citation_has_all_versions_doi(self):
        try:
            import yaml
        except ImportError:
            pytest.skip("pyyaml not installed")
        path = CORPAS_DIR / "CITATION.cff"
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
        identifiers = data.get("identifiers", [])
        all_version_dois = [
            i["value"] for i in identifiers if i.get("type") == "doi"
        ]
        assert "10.5281/zenodo.19285820" in all_version_dois


class TestReadme:
    """Tests for corpas-30x/README.md."""

    def test_readme_exists(self):
        path = CORPAS_DIR / "README.md"
        assert path.exists(), "corpas-30x/README.md must exist"

    def test_readme_contains_doi(self):
        path = CORPAS_DIR / "README.md"
        text = path.read_text(encoding="utf-8")
        assert "10.5281/zenodo.19297389" in text
        assert "10.5281/zenodo.19285820" in text

    def test_readme_contains_research_caveat(self):
        path = CORPAS_DIR / "README.md"
        text = path.read_text(encoding="utf-8")
        assert "research and educational purposes only" in text


class TestRsidLists:
    """Tests for region rsID lists."""

    def test_pgx_rsids_exists(self):
        path = REGIONS_DIR / "pgx_rsids.json"
        assert path.exists(), "corpas-30x/regions/pgx_rsids.json must exist"

    def test_pgx_rsids_valid(self):
        path = REGIONS_DIR / "pgx_rsids.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        rsids = data.get("rsids", [])
        assert len(rsids) >= 25, f"Expected >= 25 PGx rsIDs, got {len(rsids)}"
        for rsid in rsids:
            assert rsid.startswith("rs"), f"Invalid rsID format: {rsid}"

    def test_nutrigx_rsids_exists(self):
        path = REGIONS_DIR / "nutrigx_rsids.json"
        assert path.exists(), "corpas-30x/regions/nutrigx_rsids.json must exist"

    def test_nutrigx_rsids_valid(self):
        path = REGIONS_DIR / "nutrigx_rsids.json"
        data = json.loads(path.read_text(encoding="utf-8"))
        rsids = data.get("rsids", [])
        assert len(rsids) >= 25, f"Expected >= 25 NutriGx rsIDs, got {len(rsids)}"
        for rsid in rsids:
            assert rsid.startswith("rs"), f"Invalid rsID format: {rsid}"

    def test_pgx_rsids_match_pharmgx_reporter(self):
        """Verify PGx rsIDs align with pharmgx-reporter skill."""
        pgx_path = REGIONS_DIR / "pgx_rsids.json"
        data = json.loads(pgx_path.read_text(encoding="utf-8"))
        rsids = set(data["rsids"])

        # Import PGX_SNPS from pharmgx_reporter if available
        reporter_path = (
            CLAWBIO_DIR / "skills" / "pharmgx-reporter" / "pharmgx_reporter.py"
        )
        if not reporter_path.exists():
            pytest.skip("pharmgx_reporter.py not found")

        import importlib.util

        spec = importlib.util.spec_from_file_location(
            "pharmgx_reporter", reporter_path
        )
        mod = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(mod)
        except Exception:
            pytest.skip("Could not import pharmgx_reporter")

        expected = set(mod.PGX_SNPS.keys())
        missing = expected - rsids
        assert not missing, f"PGx rsIDs missing from pgx_rsids.json: {missing}"


# ---------------------------------------------------------------------------
# Data-dependent tests (skip when files absent)
# ---------------------------------------------------------------------------


class TestSubsets:
    """Tests for generated VCF subsets. Skip if not yet generated."""

    def test_chr20_subset(self):
        manifest = _load_manifest()
        path = _file_is_present(manifest, "subsets", "chr20_snps_indels.vcf.gz")
        if path is None:
            pytest.skip("chr20 subset not yet generated")
        assert path.stat().st_size > 0, "chr20 subset is empty"

    def test_pgx_subset(self):
        manifest = _load_manifest()
        path = _file_is_present(manifest, "subsets", "pgx_loci.vcf.gz")
        if path is None:
            pytest.skip("PGx subset not yet generated")
        assert path.stat().st_size > 0, "PGx subset is empty"

    def test_nutrigx_subset(self):
        manifest = _load_manifest()
        path = _file_is_present(manifest, "subsets", "nutrigx_loci.vcf.gz")
        if path is None:
            pytest.skip("NutriGx subset not yet generated")
        assert path.stat().st_size > 0, "NutriGx subset is empty"

    def test_sv_subset(self):
        manifest = _load_manifest()
        path = _file_is_present(manifest, "subsets", "sv_calls.vcf.gz")
        if path is None:
            pytest.skip("SV subset not yet generated")
        assert path.stat().st_size > 0, "SV subset is empty"

    def test_cnv_subset(self):
        manifest = _load_manifest()
        path = _file_is_present(manifest, "subsets", "cnv_calls.vcf.gz")
        if path is None:
            pytest.skip("CNV subset not yet generated")
        assert path.stat().st_size > 0, "CNV subset is empty"


class TestBaselines:
    """Tests for computed QC baselines. Skip if not yet generated."""

    def test_qc_summary_exists(self):
        manifest = _load_manifest()
        path = _file_is_present(manifest, "baselines", "qc_summary.json")
        if path is None:
            pytest.skip("QC baselines not yet computed")
        data = json.loads(path.read_text(encoding="utf-8"))
        assert "genome_build" in data
        assert data["genome_build"] == "GRCh37"

    def test_snp_ti_tv_ratio(self):
        manifest = _load_manifest()
        path = _file_is_present(manifest, "baselines", "qc_summary.json")
        if path is None:
            pytest.skip("QC baselines not yet computed")
        data = json.loads(path.read_text(encoding="utf-8"))
        snp_stats = data.get("snp_stats", {})
        ts_tv = snp_stats.get("ts_tv_ratio")
        if ts_tv is None:
            pytest.skip("Ti/Tv ratio not in baselines")
        assert 1.9 <= ts_tv <= 2.2, f"Ti/Tv ratio {ts_tv} outside expected range [1.9, 2.2]"

    def test_snp_het_hom_ratio(self):
        manifest = _load_manifest()
        path = _file_is_present(manifest, "baselines", "qc_summary.json")
        if path is None:
            pytest.skip("QC baselines not yet computed")
        data = json.loads(path.read_text(encoding="utf-8"))
        snp_stats = data.get("snp_stats", {})
        het_hom = snp_stats.get("het_hom_ratio")
        if het_hom is None:
            pytest.skip("Het/Hom ratio not in baselines")
        assert 1.3 <= het_hom <= 1.9, f"Het/Hom ratio {het_hom} outside expected range [1.3, 1.9]"

    def test_sv_counts_present(self):
        manifest = _load_manifest()
        path = _file_is_present(manifest, "baselines", "qc_summary.json")
        if path is None:
            pytest.skip("QC baselines not yet computed")
        data = json.loads(path.read_text(encoding="utf-8"))
        sv_counts = data.get("sv_counts", {})
        if not sv_counts:
            pytest.skip("SV counts not in baselines")
        # Should have at least DEL and DUP
        assert len(sv_counts) >= 2, f"Expected >= 2 SV types, got {len(sv_counts)}"
