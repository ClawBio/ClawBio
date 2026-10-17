"""Render-time lead fallback: `lead.candidates` in the config.

The configured lead can be absent from one dataset or from the LD panel; a
credible-set member present in all three then carries the plot, and the swap is
stated. Every test goes through the real render path and observes the input the
figure was drawn from (`render_full_locuscompare` is wrapped, not replaced).

Verified 2026-10-13 that (a)-(d) FAIL when `_choose_lead` is made a no-op (it
returns the configured lead with no rows): (a) on its rows, (b)-(d) on the lead
drawn / the rows / the sentence. (e) passes without the rule by construction: it
pins that a config without candidates is rendered exactly as before.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pytest  # noqa: E402
import yaml  # noqa: E402

SKILL_DIR = Path(__file__).resolve().parent.parent
for _p in (
    SKILL_DIR,
    SKILL_DIR.parent / "eqtl-catalogue-region-fetch",
    SKILL_DIR.parent / "gwas-catalog-region-fetch",
    SKILL_DIR.parent / "ld-1000g-region-compute",
):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import cli  # noqa: E402
import locuscompare_region_render as lcr  # noqa: E402
import regional_plot  # noqa: E402
from _prefetched import PrefetchedEQTLClient, PrefetchedGWASClient  # noqa: E402
from eqtl_catalogue_region_fetch import RegionVariant as EQTLRegionVariant  # noqa: E402
from gwas_catalog_region_fetch import RegionVariant as GWASRegionVariant  # noqa: E402
from ld_1000g_region_compute import LDComputeError  # noqa: E402
from ondemand_client import OnDemandLDPair, OnDemandLDResult  # noqa: E402

CHROM = "1"
CENTRE = 1_000_000
WINDOW = 100_000  # bound = 10 kb
LEAD = "1_1000000_A_G"
C1 = "1_1001000_C_T"   # PIP 0.30
C2 = "1_1002000_G_A"   # PIP 0.20
C3 = "1_999000_T_C"    # PIP 0.10
FILL = ["1_1003000_A_C", "1_998000_C_A", "1_1004000_G_T"]
CANDIDATES = [(C3, 0.10), (LEAD, 0.40), (C1, 0.30), (C2, 0.20)]  # deliberately unsorted


def _vid_parts(vid):
    c, p, r, a = vid.split("_")
    return c, int(p), r, a


def _exposure(ids):
    out = []
    for i, vid in enumerate(ids):
        c, p, r, a = _vid_parts(vid)
        out.append(EQTLRegionVariant(
            variant_id=vid, chromosome=c, position=p, ref=r, alt=a,
            beta=0.2 + 0.01 * i, se=0.05, p_value=10 ** -(5 + i), maf=0.3,
            effect_allele_frequency=0.3,
        ))
    return out


def _outcome(ids):
    out = []
    for i, vid in enumerate(ids):
        c, p, r, a = _vid_parts(vid)
        out.append(GWASRegionVariant(
            variant_id=vid, chromosome=c, position=p, ref=r, alt=a,
            beta=0.05 + 0.01 * i, se=0.01, p_value=10 ** -(4 + i),
            odds_ratio=None, effect_allele_frequency=0.3,
        ))
    return out


ALL = [LEAD, C1, C2, C3, *FILL]


class FakeLD:
    """LD pairs for every partner when the lead is in `panel` (in either allele
    order, as the real client canonicalises), else none."""

    def __init__(self, panel, raise_on=()):
        self.panel = {lcr._variant_locus_key(v) for v in panel}
        self.raise_on = set(raise_on)
        self.calls: list[str] = []

    def r2_with_lead(self, lead, partners, chromosome, window_bp):
        self.calls.append(lead)
        if lead in self.raise_on:
            raise LDComputeError("plink failed")
        partners = list(partners)
        pairs = [OnDemandLDPair(partner_variant_id=p, r2=0.5) for p in partners] \
            if lcr._variant_locus_key(lead) in self.panel else []
        return OnDemandLDResult(
            panel_id="1000g_test", panel_version="v5b", super_pop="EUR",
            plink_version="1.9", chromosome=chromosome, lead_variant_id=lead,
            window_bp=window_bp, n_partners_requested=len(partners),
            n_partners_returned=len(pairs), pairs=pairs, fetched_at_utc="t",
        )


@pytest.fixture
def drawn(monkeypatch):
    """The input the figure was actually drawn from."""
    seen: dict = {}
    original = lcr.render_full_locuscompare

    def wrapper(inp, out_path):
        seen["inp"] = inp
        return original(inp, out_path)

    monkeypatch.setattr(lcr, "render_full_locuscompare", wrapper)
    return seen


def _render(tmp_path, *, exposure_ids=ALL, outcome_ids=ALL, ld=None, candidates=CANDIDATES,
            exposure_rows=None):
    spec = lcr.LocusCompareSpec(
        lead_variant_id=LEAD, chromosome=CHROM, lead_position_bp=CENTRE, window_bp=WINDOW,
        eqtl_dataset_id="QTD000271", molecular_trait_id=None, gwas_accession="GCST1",
        prefetched_gene_track=[], lead_candidates=list(candidates),
    )
    return lcr.render_locuscompare_for_lead(
        spec,
        eqtl_client=PrefetchedEQTLClient(variants=exposure_rows or _exposure(exposure_ids)),
        gwas_client=PrefetchedGWASClient(variants=_outcome(outcome_ids), accession="GCST1"),
        ld_client=ld,
        out_path=tmp_path / "lc.png",
    )


def _status(result):
    return {r["variant_id"]: r["status"] for r in result.manifest_block["lead_candidates"]}


def _assert_partition(result):
    rows = result.manifest_block["lead_candidates"]
    ids = [r["variant_id"] for r in rows]
    assert len(ids) == len(set(ids)) == 4  # configured lead listed once, with its PIP
    assert sum(r["status"] == "used" for r in rows) == 1
    assert {r["status"] for r in rows} <= set(lcr.LEAD_CANDIDATE_STATUSES)


# ------------------------------------------------------------------- (a)-(e)
def test_a_usable_lead_is_kept_with_one_ld_call(tmp_path, drawn):
    ld = FakeLD(panel=ALL)
    result = _render(tmp_path, ld=ld)
    assert ld.calls == [LEAD]
    assert drawn["inp"].lead_variant_id == LEAD
    assert result.manifest_block["lead_fallback"] is None
    assert _status(result) == {LEAD: "used", C1: "not_evaluated", C2: "not_evaluated",
                               C3: "not_evaluated"}
    assert result.manifest_block["lead_candidates"][0]["pip"] == 0.40
    _assert_partition(result)


def test_b_lead_absent_from_ld_panel_swaps_to_the_pip_best_member(tmp_path, drawn):
    ld = FakeLD(panel=[C1, C2, C3])
    result = _render(tmp_path, ld=ld)
    assert ld.calls == [LEAD, C1]
    inp = drawn["inp"]
    # The used lead drives the diamond, the labels and the WR slope.
    assert inp.lead_variant_id == C1
    assert C1 in inp.window_label
    assert inp.r2_by_variant[C1] == 1.0 and inp.r2_by_variant[LEAD] == 0.5
    assert result.manifest_block["lead_variant_id"] == C1
    assert result.manifest_block["lead_requested_variant_id"] == LEAD
    fb = result.manifest_block["lead_fallback"]
    assert fb["requested"] == LEAD and fb["used"] == C1
    assert "not in the LD reference panel" in fb["reason"]
    sentence = next(c for c in inp.caveats if c.startswith("The plot's lead is"))
    assert C1 in sentence and LEAD in sentence and "LD reference panel" in sentence
    assert sentence in result.manifest_block["ancestry_caveats"]
    assert _status(result) == {LEAD: "not_in_ld_panel", C1: "used", C2: "not_evaluated",
                               C3: "not_evaluated"}
    _assert_partition(result)


def test_c_lead_absent_from_exposure_swaps(tmp_path, drawn):
    ld = FakeLD(panel=ALL)
    result = _render(tmp_path, exposure_ids=[v for v in ALL if v != LEAD], ld=ld)
    assert ld.calls == [C1]
    assert drawn["inp"].lead_variant_id == C1
    assert "not in the exposure data" in result.manifest_block["lead_fallback"]["reason"]
    assert _status(result)[LEAD] == "not_in_exposure"
    _assert_partition(result)


def test_d_nothing_usable_keeps_the_configured_lead_and_says_so(tmp_path, drawn):
    ld = FakeLD(panel=ALL)
    exposure = [v for v in ALL if v not in (LEAD, C1)]
    outcome = [v for v in ALL if v not in (C2,)]
    # C3 made palindromic on both sides: present, but its alleles cannot be aligned.
    rows = _exposure(exposure)
    for r in rows:
        if r.variant_id == C3:
            r.ref, r.alt = "A", "T"
    result = _render(tmp_path, exposure_rows=rows, outcome_ids=outcome, ld=ld,
                     candidates=[(LEAD, 0.4), (C1, 0.3), (C2, 0.2), (C3, 0.1)])
    assert drawn["inp"].lead_variant_id == LEAD
    assert result.manifest_block["lead_fallback"] is None
    st = _status(result)
    assert st[LEAD] == "used"
    assert st[C1] == "not_in_exposure" and st[C2] == "not_in_outcome"
    assert st[C3] in ("palindromic", "alleles_irreconcilable")
    assert any("none of the 3 other credible-set candidates is usable" in c
               for c in drawn["inp"].caveats)
    _assert_partition(result)


def test_e_no_candidates_renders_as_before(tmp_path, drawn):
    ld = FakeLD(panel=[C1])
    result = _render(tmp_path, ld=ld, candidates=[])
    assert ld.calls == [LEAD]
    assert drawn["inp"].lead_variant_id == LEAD
    assert result.manifest_block["lead_fallback"] is None
    assert result.manifest_block["lead_candidates"] == []
    assert not any(c.startswith("The plot's lead is") for c in drawn["inp"].caveats)


# ------------------------------------------------------------------- edges
def test_no_ld_client_means_ld_is_not_a_criterion(tmp_path, drawn):
    result = _render(tmp_path, ld=None)
    assert drawn["inp"].lead_variant_id == LEAD
    assert result.manifest_block["lead_fallback"] is None


def test_in_both_datasets_beats_in_ld_when_no_member_is_in_all_three(tmp_path, drawn):
    ld = FakeLD(panel=[])
    result = _render(tmp_path, exposure_ids=[v for v in ALL if v != LEAD], ld=ld)
    # LEAD fails the datasets; C1, C2, C3 are in both but none in the panel: C1 wins.
    assert ld.calls == [C1, C2, C3]
    assert drawn["inp"].lead_variant_id == C1
    assert "both the exposure and the outcome data" in next(
        c for c in drawn["inp"].caveats if c.startswith("The plot's lead is"))
    assert _status(result) == {LEAD: "not_in_exposure", C1: "used",
                               C2: "not_in_ld_panel", C3: "not_in_ld_panel"}


def test_configured_lead_wins_ties(tmp_path, drawn):
    ld = FakeLD(panel=[])
    result = _render(tmp_path, ld=ld)
    assert drawn["inp"].lead_variant_id == LEAD
    assert result.manifest_block["lead_fallback"] is None
    assert _status(result)[LEAD] == "used"


def test_ld_failure_ends_the_loop_with_the_existing_caveat(tmp_path, drawn):
    ld = FakeLD(panel=[C2], raise_on=[C1])
    result = _render(tmp_path, ld=ld)
    assert ld.calls == [LEAD, C1]  # stops at the failure; C2 never tried
    assert drawn["inp"].lead_variant_id == LEAD
    # LEAD's own LD call succeeded (no pairs: not in the panel), so the caption says
    # that, not "compute failed", which was C1's failure.
    assert not any("compute failed" in c for c in drawn["inp"].caveats)
    assert any("not in the LD reference panel" in c for c in drawn["inp"].caveats)
    # The candidate whose LD run failed is recorded as such, and the configured lead
    # kept in its place keeps its own failure as data, not only in prose.
    assert _status(result)[C1] == "ld_failed"
    rows = {r["variant_id"]: r for r in result.manifest_block["lead_candidates"]}
    assert rows[LEAD]["status"] == "used" and rows[LEAD]["used_despite"] == "not_in_ld_panel"


def test_member_beyond_the_window_bound_is_not_used(tmp_path, drawn):
    far = "1_1020000_A_G"  # 20 kb from the centre; bound is 10 kb
    ld = FakeLD(panel=[far, C1])
    result = _render(tmp_path, exposure_ids=ALL + [far], outcome_ids=ALL + [far], ld=ld,
                     candidates=[(far, 0.9), (C1, 0.3)])
    assert drawn["inp"].lead_variant_id == C1
    assert _status(result)[far] == "beyond_window_bound"


def test_allele_flipped_candidate_still_matches(tmp_path, drawn):
    ld = FakeLD(panel=[C1])
    flipped_c1 = "1_1001000_T_C"
    result = _render(tmp_path, ld=ld, candidates=[(flipped_c1, 0.3)])
    # The configured lead is not in the panel; C1 is listed in the other allele order.
    # The row records the candidate as listed; the plot draws it under the id the
    # joined data carries, so the diamond, label and slope find it.
    assert _status(result)[flipped_c1] == "used"
    inp = drawn["inp"]
    assert inp.lead_variant_id == C1
    assert any(p.variant_id == C1 and p.r2_with_lead == 1.0 for p in inp.pairs)


def test_broken_partition_raises(monkeypatch, tmp_path):
    monkeypatch.setattr(lcr, "LEAD_CANDIDATE_STATUSES", ("used",))
    with pytest.raises(lcr.LeadChoiceError):
        _render(tmp_path, ld=FakeLD(panel=ALL))


# ------------------------------------------------------------------- config + CLI
def _write_inputs(tmp_path, exposure_ids, outcome_ids):
    head = "variant_id\tchromosome\tposition_bp\tallele_a\tallele_b\tbeta\tse\tp\n"
    for name, ids in (("exposure.tsv", exposure_ids), ("outcome.tsv", outcome_ids)):
        lines = []
        for i, vid in enumerate(ids):
            c, p, r, a = _vid_parts(vid)
            lines.append(f"{vid}\t{c}\t{p}\t{r}\t{a}\t{0.1 + 0.01 * i}\t0.05\t{10 ** -(4 + i)}\n")
        (tmp_path / name).write_text(head + "".join(lines))


def _config(candidates):
    return {
        "lead": {"variant_id": LEAD, "rs_id": "rs1", "chromosome": CHROM,
                 "position_bp": CENTRE, "window_bp": WINDOW, "candidates": candidates},
        "exposure": {"sumstats_path": "exposure.tsv", "trait_label": "x"},
        "outcome": {"sumstats_path": "outcome.tsv", "trait_label": "y"},
        "ld": {"source": "none"},
        "gene_track": {"source": "none"},
    }


@pytest.mark.parametrize("bad", [
    [{"variant_id": "rs123", "pip": 0.2}],
    [{"variant_id": "2_1001000_C_T", "pip": 0.2}],
    [{"variant_id": C1, "pip": "high"}],
    [{"variant_id": C1, "pip": 1.5}],
    [C1],
    {"variant_id": C1},
])
def test_malformed_or_other_chromosome_candidates_exit_2(tmp_path, bad, capsys):
    _write_inputs(tmp_path, ALL, ALL)
    (tmp_path / "config.json").write_text(json.dumps(_config(bad)))
    rc = cli.main(["--input", str(tmp_path / "config.json"), "--output", str(tmp_path / "out")])
    assert rc == 2
    assert "lead.candidates" in capsys.readouterr().err


def test_cli_manifest_report_and_file_name_on_a_swap(tmp_path, capsys):
    _write_inputs(tmp_path, [v for v in ALL if v != LEAD], ALL)
    cands = [{"variant_id": v, "pip": p} for v, p in CANDIDATES]
    (tmp_path / "config.json").write_text(json.dumps(_config(cands)))
    out = tmp_path / "out"
    assert cli.main(["--input", str(tmp_path / "config.json"), "--output", str(out)]) == 0
    m = yaml.safe_load((out / "manifest.yaml").read_text())
    assert m["lead_variant_id"] == C1
    assert m["lead_requested_variant_id"] == LEAD
    assert m["lead_rs_id"] is None  # rs1 names the requested lead, not the drawn one
    assert m["lead_fallback"]["used"] == C1
    assert len(m["lead_candidates"]) == 4
    # The file name stays keyed to the requested lead (the agent runner builds it).
    assert m["plot_path"] == f"{LEAD}_full_locuscompare.png"
    report = (out / "report.md").read_text()
    assert f"**Lead variant:** `{C1}`" in report and f"**Requested lead:** `{LEAD}` (rs1)" in report
    stdout = capsys.readouterr().out
    assert f"The plot's lead is {C1}" in stdout


def test_cli_without_candidates_keeps_the_manifest_lead(tmp_path):
    _write_inputs(tmp_path, ALL, ALL)
    cfg = _config(None)
    del cfg["lead"]["candidates"]
    (tmp_path / "config.json").write_text(json.dumps(cfg))
    out = tmp_path / "out"
    assert cli.main(["--input", str(tmp_path / "config.json"), "--output", str(out)]) == 0
    m = yaml.safe_load((out / "manifest.yaml").read_text())
    assert m["lead_variant_id"] == LEAD and m["lead_rs_id"] == "rs1"
    assert m["lead_fallback"] is None and m["lead_candidates"] == []


def test_example_09_config_carries_candidates_that_parse():
    cfg = yaml.safe_load(
        (Path(__file__).resolve().parent.parent / "examples"
         / "09_open_targets_credible_set_many_qtls" / "config.yaml").read_text())
    cands, err = cli._parse_lead_candidates(cfg["lead"]["candidates"], str(cfg["lead"]["chromosome"]))
    assert err is None and len(cands) == 5
    assert [p for _, p in cands] == sorted((p for _, p in cands), reverse=True)
    assert cfg["lead"]["variant_id"] == "1_55025188_T_C"


# --------------------------------------------------------------------------- review rework
P_AT = "1_1001500_A_T"  # strand-ambiguous on BOTH sides, PIP between LEAD and C1


def test_palindromic_member_is_passed_over_and_named(tmp_path, drawn):
    """A/T in both datasets: joined, but palindromic, so it must not carry the plot even
    though it outranks C1; the swap sentence must name it and say why."""
    ids = [*ALL, P_AT]
    result = _render(tmp_path, exposure_ids=ids, outcome_ids=ids,
                     ld=FakeLD(panel=[P_AT, C1]),
                     candidates=[*CANDIDATES, (P_AT, 0.35)])
    st = _status(result)
    assert st[P_AT] == "palindromic"
    assert st[C1] == "used" and drawn["inp"].lead_variant_id == C1
    sentence = next(c for c in drawn["inp"].caveats if c.startswith("The plot's lead is"))
    assert f"{P_AT} (PIP 0.350) is strand-ambiguous" in sentence
    assert "is not strand-ambiguous" in sentence  # the rule, stated, not "highest present"


def test_swap_keeps_the_x_range_and_title_on_the_fetched_window(tmp_path, drawn):
    _render(tmp_path, ld=FakeLD(panel=[C1]))
    inp = drawn["inp"]
    assert inp.lead_variant_id == C1
    assert inp.window_centre_variant_id == LEAD
    assert f"requested lead {LEAD} (window as fetched); lead shown {C1}" in inp.window_label
    assert regional_plot._compute_xlim_bp(inp.pairs, inp.window_centre_variant_id, inp.window_bp) == (
        CENTRE - WINDOW // 2, CENTRE + WINDOW // 2)


def test_ld_tries_are_capped(tmp_path):
    many = [f"1_{1_000_100 + 100 * i}_A_C" for i in range(15)]
    ids = [LEAD, *many, *FILL]
    ld = FakeLD(panel=[])
    result = _render(tmp_path, exposure_ids=ids, outcome_ids=ids, ld=ld,
                     candidates=[(v, 0.5 - 0.01 * i) for i, v in enumerate(many)])
    assert len(ld.calls) == lcr.MAX_LEAD_LD_ATTEMPTS
    assert sum(1 for s in _status(result).values() if s == "not_evaluated") == 16 - lcr.MAX_LEAD_LD_ATTEMPTS


def test_candidates_with_a_synthetic_ld_matrix_exit_2(tmp_path, capsys):
    _write_inputs(tmp_path, ALL, ALL)
    cfg = _config([{"variant_id": C1, "pip": 0.3}])
    cfg["ld"] = {"source": "synthetic", "ld_matrix_path": "ld.tsv"}
    (tmp_path / "config.json").write_text(json.dumps(cfg))
    assert cli.main(["--input", str(tmp_path / "config.json"), "--output", str(tmp_path / "out")]) == 2
    assert "synthetic" in capsys.readouterr().err


def test_a_long_caption_grows_the_figure_instead_of_overlapping(tmp_path, monkeypatch):
    from matplotlib.figure import Figure
    heights: list[float] = []
    real = Figure.savefig

    def spy(self, *a, **k):
        heights.append(self.get_figheight())
        return real(self, *a, **k)

    monkeypatch.setattr(Figure, "savefig", spy)
    spec_kw = dict(lead_variant_id=LEAD, chromosome=CHROM, lead_position_bp=CENTRE, window_bp=WINDOW,
                   eqtl_dataset_id="QTD000271", molecular_trait_id=None, gwas_accession="GCST1",
                   prefetched_gene_track=[])
    for notes in ([], ["a caveat " * 40] * 4):
        lcr.render_locuscompare_for_lead(
            lcr.LocusCompareSpec(**spec_kw, notes=notes),
            eqtl_client=PrefetchedEQTLClient(variants=_exposure(ALL)),
            gwas_client=PrefetchedGWASClient(variants=_outcome(ALL), accession="GCST1"),
            ld_client=None, out_path=tmp_path / f"lc{len(notes)}.png")
    assert heights[1] > heights[0] + 0.5


def _render_spy(monkeypatch):
    from matplotlib.figure import Figure
    seen: dict = {}
    real = Figure.savefig

    def spy(self, *a, **k):
        seen["xlim"] = self.axes[0].get_xlim()
        seen["title"] = self._suptitle.get_text() if self._suptitle else ""
        return real(self, *a, **k)

    monkeypatch.setattr(Figure, "savefig", spy)
    return seen


def test_rendered_axes_and_title_stay_on_the_fetched_window(tmp_path, monkeypatch):
    seen = _render_spy(monkeypatch)
    _render(tmp_path, ld=FakeLD(panel=[C1]))
    lo, hi = seen["xlim"]
    assert (lo + hi) / 2 == pytest.approx(CENTRE, abs=1)  # not C1's position
    assert f"of {LEAD}; lead shown {C1}" in seen["title"]


def test_the_swap_sentence_names_at_most_three_passed_over_members(tmp_path, drawn):
    absent = [f"1_{1_000_200 + 50 * i}_A_C" for i in range(6)]  # in no dataset
    cands = [(v, 0.9 - 0.01 * i) for i, v in enumerate(absent)] + [(C1, 0.3)]
    _render(tmp_path, ld=FakeLD(panel=[C1]), candidates=cands)
    sentence = next(c for c in drawn["inp"].caveats if c.startswith("The plot's lead is"))
    assert sum(v in sentence for v in absent) == 3
    assert "3 more higher-PIP members were passed over" in sentence


def test_the_on_demand_clients_own_error_degrades_like_an_ld_failure(tmp_path, drawn):
    """The on-demand client raises OnDemandLDError, which is not an LDComputeError."""

    class OnDemandFailingLD(FakeLD):
        def r2_with_lead(self, lead, partners, chromosome, window_bp):
            self.calls.append(lead)
            raise lcr.OnDemandLDError("plink not found")

    ld = OnDemandFailingLD(panel=[])
    result = _render(tmp_path, ld=ld)
    assert ld.calls == [LEAD]
    assert _status(result)[LEAD] == "used"
    assert any("compute failed" in c for c in drawn["inp"].caveats)
