"""Unit tests for the on-demand 1000G LD client (plink 1.9 subprocess).

We mock the region fetch (requests / pysam) and the plink subprocess so
these tests are offline and don't require plink, pysam, or network. Live
smoke (real plink 1.9 against a tiny 1000G window) lives in
test_live_ld_1000g_region_compute.py.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from ondemand_client import (
    OnDemand1000GLDClient,
    OnDemandLDError,
    _alpha_plink_id,
    _parse_ld,
)
import ondemand_client as oc_module


# -----------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------


def _patch_plink_version(monkeypatch, version: str = "PLINK v1.90b6.27 64-bit (2023-05-09)") -> None:
    monkeypatch.setattr(oc_module, "_detect_plink_version", lambda _bin: version)


def _patch_samples(monkeypatch, samples: list[str] | None = None) -> None:
    samples = samples or ["NA12878", "NA12879", "NA12891"]
    monkeypatch.setattr(
        oc_module, "_resolve_super_pop_samples",
        lambda _super_pop, _cache_dir: list(samples),
    )


def _patch_region_fetch(monkeypatch, tmp_path: Path) -> Path:
    """Replace the tabix fetch with a no-op that returns an empty .vcf.gz path."""
    fake_vcf = tmp_path / "fake_region.vcf.gz"
    fake_vcf.write_bytes(b"")
    monkeypatch.setattr(
        oc_module, "_fetch_region_vcf",
        lambda chrom, start, end, cache_dir: fake_vcf,
    )
    return fake_vcf


def _write_ld(out_prefix: Path, lead: str, partner_pairs: list[tuple[str, float]]) -> None:
    """Emit a plink 1.9 .ld file at <out_prefix>.ld.

    plink 1.9 emits whitespace-separated columns:
    `CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2`.
    """

    def _to_panel(ot_id: str) -> str:
        return ot_id.replace("_", ":", 3)

    rows = [" CHR_A         BP_A        SNP_A  CHR_B         BP_B        SNP_B           R2"]
    for partner, r2 in partner_pairs:
        rows.append(
            f"   2     36910110  {_to_panel(lead)}    "
            f"   2     36932656  {_to_panel(partner)}    {r2:.6f}"
        )
    Path(f"{out_prefix}.ld").write_text("\n".join(rows) + "\n")


# -----------------------------------------------------------------
# Construction / validation
# -----------------------------------------------------------------


def test_client_validates_plink_present(monkeypatch, tmp_path):
    _patch_plink_version(monkeypatch)
    c = OnDemand1000GLDClient(super_pop="EUR", plink_bin="plink", cache_dir=tmp_path)
    assert c.plink_version.startswith("PLINK v1.90")
    assert c.super_pop == "EUR"


def test_client_raises_when_plink_missing(monkeypatch, tmp_path):
    def boom(_bin):
        raise OnDemandLDError("plink not found")
    monkeypatch.setattr(oc_module, "_detect_plink_version", boom)
    with pytest.raises(OnDemandLDError, match="plink not found"):
        OnDemand1000GLDClient(super_pop="EUR", plink_bin="plink", cache_dir=tmp_path)


# -----------------------------------------------------------------
# r2_with_lead happy path
# -----------------------------------------------------------------


def test_r2_with_lead_parses_ld_output(monkeypatch, tmp_path):
    _patch_plink_version(monkeypatch)
    _patch_samples(monkeypatch)
    _patch_region_fetch(monkeypatch, tmp_path)

    captured_cmd: dict = {}

    def fake_run(cmd, capture_output, text, check, timeout):
        captured_cmd["cmd"] = cmd
        out_prefix = Path(cmd[cmd.index("--out") + 1])
        _write_ld(out_prefix, lead="2_36910110_C_T", partner_pairs=[
            ("2_36932656_A_G", 0.94),
            ("2_36905984_C_T", 0.81),
            ("2_36897612_T_C", 0.40),
        ])
        return MagicMock(returncode=0, stdout="", stderr="")

    monkeypatch.setattr("subprocess.run", fake_run)
    client = OnDemand1000GLDClient(super_pop="EUR", plink_bin="plink", cache_dir=tmp_path)
    result = client.r2_with_lead(
        lead="2_36910110_C_T",
        partners=["2_36932656_A_G", "2_36905984_C_T", "2_36897612_T_C"],
        chromosome="2",
        window_bp=1_000_000,
    )

    assert result.lead_variant_id == "2_36910110_C_T"
    assert result.super_pop == "EUR"
    assert result.n_partners_requested == 3
    assert result.n_partners_returned == 3
    r2_by_id = {p.partner_variant_id: p.r2 for p in result.pairs}
    assert r2_by_id["2_36932656_A_G"] == pytest.approx(0.94)
    assert r2_by_id["2_36905984_C_T"] == pytest.approx(0.81)
    assert r2_by_id["2_36897612_T_C"] == pytest.approx(0.40)

    # Verify cmd shape uses plink 1.9 flags.
    cmd = captured_cmd["cmd"]
    assert "--r2" in cmd
    assert "--ld-snp" in cmd
    assert "--ld-window-r2" in cmd
    assert "--ld-window-kb" in cmd
    assert "--set-missing-var-ids" in cmd
    # plink 1.9 should NOT see plink2's matrix-only flag.
    assert "--r2-unphased" not in cmd


def test_r2_with_lead_handles_empty_partners(monkeypatch, tmp_path):
    _patch_plink_version(monkeypatch)
    client = OnDemand1000GLDClient(super_pop="EUR", plink_bin="plink", cache_dir=tmp_path)
    result = client.r2_with_lead(
        lead="2_1_C_T", partners=[], chromosome="2", window_bp=1_000,
    )
    assert result.n_partners_requested == 0
    assert result.pairs == []
    assert any("no partners requested" in n for n in result.notes)


def test_r2_with_lead_raises_on_plink_nonzero_exit(monkeypatch, tmp_path):
    _patch_plink_version(monkeypatch)
    _patch_samples(monkeypatch)
    _patch_region_fetch(monkeypatch, tmp_path)

    def boom(cmd, capture_output, text, check, timeout):
        return MagicMock(returncode=1, stdout="", stderr="missing variant id")
    monkeypatch.setattr("subprocess.run", boom)

    client = OnDemand1000GLDClient(super_pop="EUR", plink_bin="plink", cache_dir=tmp_path)
    with pytest.raises(OnDemandLDError, match="exited with code 1"):
        client.r2_with_lead(
            lead="2_1_C_T", partners=["2_2_A_G"], chromosome="2", window_bp=1_000,
        )


def test_r2_with_lead_raises_when_ld_absent(monkeypatch, tmp_path):
    _patch_plink_version(monkeypatch)
    _patch_samples(monkeypatch)
    _patch_region_fetch(monkeypatch, tmp_path)

    def succeeds_but_writes_nothing(cmd, capture_output, text, check, timeout):
        return MagicMock(returncode=0, stdout="ok", stderr="")
    monkeypatch.setattr("subprocess.run", succeeds_but_writes_nothing)

    client = OnDemand1000GLDClient(super_pop="EUR", plink_bin="plink", cache_dir=tmp_path)
    with pytest.raises(OnDemandLDError, match=r"no \.ld output"):
        client.r2_with_lead(
            lead="2_1_C_T", partners=["2_2_A_G"], chromosome="2", window_bp=1_000,
        )


# -----------------------------------------------------------------
# Parser-only tests
# -----------------------------------------------------------------


def test_parse_ld_skips_rows_without_lead_match(tmp_path):
    """An .ld row whose SNP_A and SNP_B are both unrelated to the lead is ignored."""
    path = tmp_path / "ld_out.ld"
    rows = [
        " CHR_A         BP_A        SNP_A  CHR_B         BP_B        SNP_B           R2",
        "   2            1   2:x:C:T      2            2   2:y:A:G        0.500000",
        "   2            1         lead    2            2   2:z:A:G        0.700000",
    ]
    path.write_text("\n".join(rows) + "\n")

    notes: list[str] = []
    pairs = _parse_ld(path, "lead", notes)
    assert len(pairs) == 1
    assert pairs[0].partner_variant_id == "2:z:A:G"
    assert pairs[0].r2 == pytest.approx(0.7)


def test_alpha_plink_id_sorts_alleles_alphabetically():
    """plink 1.9's --set-missing-var-ids '@:#:$1:$2' uses alphabetical-order
    alleles. Confirm `_alpha_plink_id` matches that convention so callers
    can pass OT-form ids in either allele order and still find the variant.
    """
    # ALT < REF alphabetically: must swap.
    assert _alpha_plink_id("1_109330668_G_A") == "1:109330668:A:G"
    assert _alpha_plink_id("1_109462434_T_A") == "1:109462434:A:T"
    # REF < ALT alphabetically: no swap (identity transformation modulo separator).
    assert _alpha_plink_id("1_109274968_G_T") == "1:109274968:G:T"
    assert _alpha_plink_id("2_36932656_A_G") == "2:36932656:A:G"
    # Same allele on both sides (degenerate; plink would never emit this
    # but the helper must not crash).
    assert _alpha_plink_id("1_100_A_A") == "1:100:A:A"
    # Unparseable input -> naive `:`-swap fallback, never raise.
    assert _alpha_plink_id("not_a_variant") == "not:a:variant"


def test_r2_with_lead_maps_alpha_form_back_to_ot_form(monkeypatch, tmp_path):
    """End-to-end mock: pass partners in OT-form where ALT < REF, confirm
    the client writes alphabetical-allele ids into `--extract`/`--ld-snp`
    and maps results back to the caller's OT-form.
    """
    _patch_plink_version(monkeypatch)
    _patch_samples(monkeypatch)
    _patch_region_fetch(monkeypatch, tmp_path)

    captured_cmd: dict = {}
    captured_extract: dict = {}

    def fake_run(cmd, capture_output, text, check, timeout):
        captured_cmd["cmd"] = cmd
        # Capture the extract list plink saw.
        extract_path = Path(cmd[cmd.index("--extract") + 1])
        captured_extract["lines"] = extract_path.read_text().splitlines()
        out_prefix = Path(cmd[cmd.index("--out") + 1])
        # plink emits SNP_A / SNP_B in its internal alphabetical-allele form.
        _write_ld(out_prefix, lead="1:109274968:G:T", partner_pairs=[
            # Partner the caller passed as "1_109330668_G_A"; plink sees "1:109330668:A:G".
            ("1:109330668:A:G", 0.85),
            # Partner the caller passed as "1_109462434_T_A"; plink sees "1:109462434:A:T".
            ("1:109462434:A:T", 0.40),
        ])
        return MagicMock(returncode=0, stdout="", stderr="")

    monkeypatch.setattr("subprocess.run", fake_run)
    client = OnDemand1000GLDClient(super_pop="EUR", plink_bin="plink", cache_dir=tmp_path)
    result = client.r2_with_lead(
        lead="1_109274968_G_T",
        partners=["1_109330668_G_A", "1_109462434_T_A"],
        chromosome="1",
        window_bp=1_000_000,
    )

    # --extract / --ld-snp must use plink's alphabetical-allele form.
    assert captured_extract["lines"] == [
        "1:109274968:G:T", "1:109330668:A:G", "1:109462434:A:T",
    ]
    assert "1:109274968:G:T" == captured_cmd["cmd"][captured_cmd["cmd"].index("--ld-snp") + 1]

    # Returned partner ids must be the caller's OT-form (not plink's alpha-form).
    returned_ids = sorted(p.partner_variant_id for p in result.pairs)
    assert returned_ids == ["1_109330668_G_A", "1_109462434_T_A"]


def test_r2_with_lead_returns_empty_pairs_when_lead_not_in_panel(monkeypatch, tmp_path):
    """plink errors with 'No valid variants specified by --ld-snp' when the
    lead variant is not in the 1000G panel (rare, indel, post-1000G). The
    client must degrade gracefully (empty pairs + explanatory note),
    NOT raise, so a regional plot renders without LD colouring instead of
    failing.
    """
    _patch_plink_version(monkeypatch)
    _patch_samples(monkeypatch)
    _patch_region_fetch(monkeypatch, tmp_path)

    def fake_run(cmd, capture_output, text, check, timeout):
        return MagicMock(
            returncode=5, stdout="",
            stderr="Error: No valid variants specified by --ld-snp/--ld-snps/--ld-snp-list.",
        )

    monkeypatch.setattr("subprocess.run", fake_run)
    client = OnDemand1000GLDClient(super_pop="EUR", plink_bin="plink", cache_dir=tmp_path)
    result = client.r2_with_lead(
        lead="1_999999_X_Y",
        partners=["1_999998_A_G", "1_999997_C_T"],
        chromosome="1",
        window_bp=1_000_000,
    )

    assert result.pairs == []
    assert result.n_partners_returned == 0
    assert result.n_partners_requested == 2
    assert any("not in 1000G" in n for n in result.notes)


def test_r2_with_lead_finds_a_lead_whose_alt_sorts_before_its_ref(monkeypatch, tmp_path):
    """A lead with REF=G, ALT=A is `1:67240275:A:G` in plink 1.9's id space, so the
    client must ask for that id. Asking for `1:67240275:G:A` makes plink fail with
    "No valid variants specified by --ld-snp" for every lead in that orientation
    (12,991 of the 27,603 rows in the 1 Mb window around rs11209026 at IL23R)."""
    _patch_plink_version(monkeypatch)
    _patch_samples(monkeypatch)
    _patch_region_fetch(monkeypatch, tmp_path)

    seen: dict = {}

    def fake_run(cmd, capture_output, text, check, timeout):
        seen["ld_snp"] = cmd[cmd.index("--ld-snp") + 1]
        seen["extract"] = Path(cmd[cmd.index("--extract") + 1]).read_text().splitlines()
        out_prefix = Path(cmd[cmd.index("--out") + 1])
        _write_ld(out_prefix, lead="1:67240275:A:G", partner_pairs=[("1:67242007:A:G", 0.9)])
        return MagicMock(returncode=0, stdout="", stderr="")

    monkeypatch.setattr("subprocess.run", fake_run)
    client = OnDemand1000GLDClient(super_pop="EUR", plink_bin="plink", cache_dir=tmp_path)
    result = client.r2_with_lead(
        lead="1_67240275_G_A",
        partners=["1_67242007_G_A"],
        chromosome="1",
        window_bp=1_000_000,
    )

    assert seen["ld_snp"] == "1:67240275:A:G"
    assert seen["extract"][0] == "1:67240275:A:G"
    assert [(p.partner_variant_id, p.r2) for p in result.pairs] == [("1_67242007_G_A", 0.9)]
