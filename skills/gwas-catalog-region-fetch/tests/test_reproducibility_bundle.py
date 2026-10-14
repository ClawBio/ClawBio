"""--demo writes a reproducibility bundle whose checksums verify from output_dir."""

from __future__ import annotations

import hashlib
from pathlib import Path

import gwas_catalog_region_fetch as grf


def _fake_result(cfg):
    return grf.RegionResult(
        accession=cfg["accession"], chromosome=str(cfg["chromosome"]),
        region_start_bp=int(cfg["start_bp"]), region_end_bp=int(cfg["end_bp"]),
        n_variants=1,
        variants=[grf.RegionVariant(
            variant_id="1_109274968_G_T", chromosome="1", position=109274968,
            ref="G", alt="T", beta=-0.27, se=0.01, p_value=1e-50,
            odds_ratio=None, effect_allele_frequency=0.2, raw={},
        )],
        release=grf.GWASCatalogRelease(
            accession=cfg["accession"], harmonised_url="https://example.invalid/x.h.tsv.gz",
            fetched_at_utc="2026-01-01T00:00:00Z",
        ),
        notes=[],
    )


def test_demo_writes_verifiable_reproducibility_bundle(tmp_path, monkeypatch):
    monkeypatch.setattr(grf, "_fetch_with_cache",
                        lambda *, client, cfg, cache_dir: _fake_result(cfg))
    out = tmp_path / "out"
    assert grf.main(["--demo", "--output", str(out)]) == 0

    repro = out / "reproducibility"
    commands = (repro / "commands.sh").read_text()
    assert "python skills/gwas-catalog-region-fetch/gwas_catalog_region_fetch.py --demo" in commands
    assert str(Path(grf.__file__).resolve()) not in commands
    env = (repro / "environment.yml").read_text()
    assert "pysam" in env and "requests" in env

    lines = (repro / "checksums.sha256").read_text().splitlines()
    labels = {line.split("  ", 1)[1] for line in lines}
    assert {"report.md", "variants.tsv", "reproducibility/environment.yml"} <= labels
    for line in lines:
        digest, label = line.split("  ", 1)
        assert hashlib.sha256((out / label).read_bytes()).hexdigest() == digest, label
