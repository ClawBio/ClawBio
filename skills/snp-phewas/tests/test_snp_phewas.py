"""Unit tests (offline): catalog, tabix reader, orientation, lookup, partition, outages, CLI."""

from __future__ import annotations

import gzip
import os
import struct
import zlib

import pytest
import requests

import snp_phewas_schema as S
from snp_phewas_catalog import (
    UnknownDatasetError,
    default_dataset_ids,
    load_datasets,
    load_outcome_catalog,
    resolve_dataset_ids,
)
from snp_phewas_pheweb import FinnGenPhewebClient
from snp_phewas_tabix import (
    SourceFileMissing,
    SourceUnreachable,
    fetch_position,
    read_header,
)
from snp_phewas_lookup import lookup_variant
from snp_phewas_orient import orient
from snp_phewas_demo import (
    FG_DIR,
    FG_OUTCOMES,
    PAN_DIR,
    PAN_OUTCOMES,
    DemoFinnGenApi,
    demo_map,
    demo_rsid,
    demo_readers,
)

LDL = "biomarkers-30780-both_sexes-irnt"


def _row(result, ds, oid):
    (r,) = [x for x in result.rows if x.dataset_id == ds and x.outcome_id == oid]
    return r


# --- catalog ---------------------------------------------------------------------------

def test_catalog_lists_the_curated_datasets_with_license_and_attribution():
    ds = load_datasets()
    assert set(default_dataset_ids()) == {"pan_ukbb_eur", "finngen_r12"}
    assert {"pan_ukbb_afr", "pan_ukbb_amr", "pan_ukbb_csa", "pan_ukbb_eas", "pan_ukbb_mid"} <= set(ds)
    assert ds["pan_ukbb_eur"].license == "CC-BY-4.0" and ds["pan_ukbb_eur"].build == "GRCh37"
    assert "acknowledge the participants and investigators of the FinnGen study" in ds["finngen_r12"].attribution
    assert resolve_dataset_ids("all") == list(ds)
    # a repeated id would repeat every row
    assert resolve_dataset_ids("finngen_r12,finngen_r12") == ["finngen_r12"]
    with pytest.raises(UnknownDatasetError):
        resolve_dataset_ids("mvp_gwphewas")


def test_outcome_catalog_sizes_are_the_published_manifests():
    # Pan-UKBB manifest 2024-06-26: 7,223 phenotypes, 7,160 with a European arm.
    # FinnGen R12 manifest: 2,469 endpoints.
    assert len(load_outcome_catalog("pan_ukbb_eur")) == 7160
    assert len(load_outcome_catalog("pan_ukbb_afr")) == 2445
    fg = load_outcome_catalog("finngen_r12")
    assert len(fg) == 2469
    ldl = {o.outcome_id: o for o in load_outcome_catalog("pan_ukbb_eur")}[LDL]
    assert ldl.outcome_kind == "continuous" and ldl.n_cases == 400223 and ldl.n_controls is None
    chd = {o.outcome_id: o for o in fg}["I9_CHD"]
    assert (chd.n_cases, chd.n_controls, chd.outcome_kind) == (56650, 443698, "binary")


# --- tabix reader ----------------------------------------------------------------------

def test_reader_finds_the_variant_in_both_fixture_layouts():
    rows = fetch_position(PAN_DIR + LDL + ".tsv.bgz", PAN_DIR + LDL + ".tsv.bgz.tbi", "1", 109817590)
    assert len(rows) == 1 and rows[0][:4] == ["1", "109817590", "G", "T"]
    assert read_header(PAN_DIR + LDL + ".tsv.bgz")[:4] == ["chr", "pos", "ref", "alt"]
    fg = FG_DIR + "finngen_R12_I9_CHD.gz"
    rows = fetch_position(fg, fg + ".tbi", "1", 109274968)
    assert rows[0][4] == "rs12740374"
    assert read_header(fg)[0] == "chrom"


def test_reader_reports_an_absent_position_as_empty_after_reading_past_it():
    assert fetch_position(PAN_DIR + LDL + ".tsv.bgz", PAN_DIR + LDL + ".tsv.bgz.tbi", "1", 109817591) == []
    assert fetch_position(PAN_DIR + LDL + ".tsv.bgz", PAN_DIR + LDL + ".tsv.bgz.tbi", "2", 5) == []


def _bgzf_block(payload: bytes) -> bytes:
    comp = zlib.compressobj(6, zlib.DEFLATED, -15)
    data = comp.compress(payload) + comp.flush()
    bsize = 18 + len(data) + 8 - 1
    header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", bsize)
    return header + data + struct.pack("<II", zlib.crc32(payload) & 0xFFFFFFFF, len(payload))


def _tbi(names: list[str], linear: list[int]) -> bytes:
    nm = b"".join(n.encode() + b"\0" for n in names)
    body = b"TBI\x01" + struct.pack("<8i", len(names), 0, 1, 2, 2, ord("#"), 1, len(nm)) + nm
    body += struct.pack("<i", 0) + struct.pack("<i", len(linear)) + struct.pack(f"<{len(linear)}Q", *linear)
    return gzip.compress(body)


def test_reader_survives_an_index_offset_into_an_empty_block(tmp_path):
    """The Pan-UKBB index defect, reproduced: the linear offset names an EMPTY bgzip block
    with an in-block offset that is only meaningful in the next block. A reader that seeks
    into the empty block (htslib) returns nothing; this reader must find the record."""
    header = b"chr\tpos\tref\talt\n"
    block_a = header + b"1\t100\tA\tG\n"
    tile1 = b"1\t16500\tC\tT\n1\t16600\tG\tA\n"
    pre = b"1\t16400\tA\tC\n"  # first bytes of block B, before tile 1
    a = _bgzf_block(block_a)
    empty = _bgzf_block(b"")
    b = _bgzf_block(pre + tile1)
    eof = _bgzf_block(b"")
    data = tmp_path / "x.tsv.bgz"
    data.write_bytes(a + empty + b + eof)
    # Both tiles' offsets name the EMPTY block, so neither tile can rescue the other and only
    # inflating ACROSS the block boundary finds the record (the review showed tile 0's valid
    # offset used to rescue an htslib-style single-block read). tile 1 (16385..32768):
    # uoffset = len(pre), meaningful in block B only.
    voff_tile1 = (len(a) << 16) | len(pre)
    voff_tile0 = (len(a) << 16) | 1  # also past the empty block's end: no single-block rescue
    idx = tmp_path / "x.tsv.bgz.tbi"
    idx.write_bytes(_tbi(["1"], [voff_tile0, voff_tile1]))
    rows = fetch_position(str(data), str(idx), "1", 16600)
    assert rows == [["1", "16600", "G", "A"]]
    assert fetch_position(str(data), str(idx), "1", 16601) == []


def _single_block_file(tmp_path, lines: list[bytes], linear: list[int], name="y"):
    header = b"chr\tpos\tref\talt\n"
    payload = header + b"".join(lines)
    data = tmp_path / f"{name}.tsv.bgz"
    data.write_bytes(_bgzf_block(payload) + _bgzf_block(b""))
    idx = tmp_path / f"{name}.tsv.bgz.tbi"
    idx.write_bytes(_tbi(["1"], linear))
    return str(data), str(idx), header


def test_reader_starts_one_tile_early_against_an_off_by_one_index(tmp_path):
    """The record at 16390 sits in tile 1 but BEFORE the offset the index gives tile 1.
    Starting at tile 1's own offset reads past it; the one-tile slack finds it."""
    before = b"1\t100\tA\tG\n"
    target = b"1\t16390\tA\tG\n"
    after = b"1\t16500\tC\tT\n"
    data, idx, header = _single_block_file(tmp_path, [before, target, after], [0, 0])
    lin = [len(header), len(header) + len(before) + len(target)]  # tile 1 skips 16390
    with open(idx, "wb") as fh:
        fh.write(_tbi(["1"], lin))
    assert fetch_position(data, idx, "1", 16390) == [["1", "16390", "A", "G"]]


def _incompressible_lines(n: int, start: int = 1000) -> list[bytes]:
    import random

    rnd = random.Random(7)
    return [f"1\t{start + i}\tA\tG\t{rnd.getrandbits(256):064x}\n".encode() for i in range(n)]


def _multi_block_file(tmp_path, lines: list[bytes], linear: list[int] | None = None):
    """bgzip file of `lines` in ~60 KB blocks; returns (data, idx, [(coffset, first_pos)])."""
    header = b"chr\tpos\tref\talt\tpad\n"
    payloads, cur = [], header
    for ln in lines:
        if len(cur) + len(ln) > 60000:
            payloads.append(cur)
            cur = b""
        cur += ln
    payloads.append(cur)
    blocks, table, off = [], [], 0
    for pl in payloads:
        first = next(l for l in pl.split(b"\n") if l and l[:1].isdigit())
        table.append((off, int(first.split(b"\t")[1])))
        b = _bgzf_block(pl)
        blocks.append(b)
        off += len(b)
    data = tmp_path / "big.tsv.bgz"
    data.write_bytes(b"".join(blocks) + _bgzf_block(b""))
    idx = tmp_path / "big.tsv.bgz.tbi"
    idx.write_bytes(_tbi(["1"], linear or [len(header)]))
    return str(data), str(idx), table


def test_reader_falls_back_to_the_nearest_earlier_tile_when_both_are_empty(tmp_path, monkeypatch):
    """Tiles 2 and 3 have no linear entry. The scan must start at tile 1's offset, the
    nearest earlier one, not at byte 0: from byte 0 the read cap is hit and a real answer
    would be reported as an outage. The cap is set between the two distances."""
    import snp_phewas_tabix as ct

    # Dense in tile 0 (many bytes before tile 1), sparse after, so the two start points are
    # far apart in compressed bytes.
    positions = list(range(1000, 16000, 2)) + list(range(16400, 52000, 40))
    lines = [f"1\t{q}\tA\tG\t{h}\n".encode() for q, h in zip(positions, _hex_pad(len(positions)))]
    _, _, table = _multi_block_file(tmp_path, lines)
    def block_of(q):  # compressed offset of the block holding position q
        return max(c for c, first in table if first <= q)

    target = next(q for q in positions if q > 3 * 16384)
    tile1, target_block = block_of(16400), block_of(target)
    data, idx, _ = _multi_block_file(tmp_path, lines, [len(b"chr\tpos\tref\talt\tpad\n"),
                                                       tile1 << 16, 0, 0])
    assert target_block - tile1 < target_block  # the fallback is strictly closer
    monkeypatch.setattr(ct, "_INITIAL_READ", 4 * 1024)
    cap = 4 * 1024
    while cap < (target_block - tile1) + 64 * 1024:
        cap *= 2
    assert cap < target_block, "fixture too small to separate the two starting points"
    monkeypatch.setattr(ct, "_MAX_READ", cap)
    assert [r[1] for r in fetch_position(data, idx, "1", target)] == [str(target)]


def _hex_pad(n: int) -> list[str]:
    import random

    rnd = random.Random(11)
    return [f"{rnd.getrandbits(256):064x}" for _ in range(n)]


def test_reader_doubles_its_read_until_it_passes_the_position(tmp_path, monkeypatch):
    """A target further from the index offset than the first read: the reader must keep
    reading (and not answer from the first window). Read sizes are shrunk so the fixture
    stays small; the code path is the one a real 128 KB first read takes."""
    import snp_phewas_tabix as ct

    lines = _incompressible_lines(3000)  # ~230 KB of near-incompressible text
    data, idx, _ = _multi_block_file(tmp_path, lines)
    monkeypatch.setattr(ct, "_INITIAL_READ", 16 * 1024)
    monkeypatch.setattr(ct, "_MAX_READ", 1024 * 1024)
    pos = 1000 + 2900
    rows = fetch_position(data, idx, "1", pos)
    assert [r[1] for r in rows] == [str(pos)]
    assert fetch_position(data, idx, "1", 1000 + 3000) == []  # read to EOF, honestly absent


def test_reader_refuses_to_call_a_variant_absent_when_the_read_cap_is_hit(tmp_path, monkeypatch):
    import snp_phewas_tabix as ct

    lines = _incompressible_lines(3000)
    data, idx, _ = _multi_block_file(tmp_path, lines)
    monkeypatch.setattr(ct, "_INITIAL_READ", 16 * 1024)
    monkeypatch.setattr(ct, "_MAX_READ", 32 * 1024)
    with pytest.raises(SourceUnreachable, match="refusing to report the variant as absent"):
        fetch_position(data, idx, "1", 1000 + 2900)


_FG_HEADER = b"chrom\tpos\tref\talt\tbeta\tsebeta\tmlogp\taf_alt\n"


def _fg_file(tmp_path, record: bytes):
    data = tmp_path / "fg.gz"
    data.write_bytes(_bgzf_block(_FG_HEADER + record) + _bgzf_block(b""))
    idx = tmp_path / "fg.gz.tbi"
    idx.write_bytes(_tbi(["1"], [len(_FG_HEADER)]))
    return str(data), str(idx)


def _read_fg(data, idx):
    from snp_phewas_sources import FINNGEN_COLUMNS, read_outcome_file
    from snp_phewas_variant import Placement

    tgt = Placement("GRCh38", "1", 109274968, "G", "T")
    return read_outcome_file(data, idx, FINNGEN_COLUMNS, tgt, index_cache=None)


def test_a_record_with_the_wrong_field_count_is_unparsable_not_absent(tmp_path):
    good = _read_fg(*_fg_file(tmp_path, b"1\t109274968\tG\tT\t-0.08\t0.0086\t19.7\t0.2\n"))
    assert good.status == S.ST_OK
    rec = _read_fg(*_fg_file(tmp_path, b"1\t109274968\tG\tT\t-0.08\t0.0086\n"))
    assert rec.status == S.ST_RECORD_UNPARSABLE, rec


def test_a_corrupt_file_is_one_unreadable_outcome_not_a_failed_lookup(tmp_path):
    data, idx = _fg_file(tmp_path, b"1\t109274968\tG\tT\t-0.08\t0.0086\t19.7\t0.2\n")
    raw = bytearray(open(data, "rb").read())
    raw[20:50] = b"\xff" * 30  # damage the deflate stream of the first block
    open(data, "wb").write(bytes(raw))
    from snp_phewas_sources import _header
    _header.cache_clear()
    rec = _read_fg(data, idx)
    assert rec.status == S.ST_RECORD_UNPARSABLE, rec


def test_index_cache_is_keyed_by_url_and_scoped_by_release(tmp_path, monkeypatch):
    import snp_phewas_tabix as ct
    from snp_phewas_catalog import load_datasets
    from snp_phewas_sources import release_index_cache

    fetched = []
    monkeypatch.setattr(ct, "_http_get", lambda url, **k: fetched.append(url) or b"x")
    a = ct.ensure_index("https://h/one/x.tbi", str(tmp_path))
    b = ct.ensure_index("https://h/two/x.tbi", str(tmp_path))
    assert a != b and len(fetched) == 2  # same basename, two sources, two cache entries
    spec = load_datasets()["pan_ukbb_eur"]
    other = type(spec)(**{**spec.__dict__, "release": spec.release + "-next"})
    assert release_index_cache(spec, str(tmp_path)) != release_index_cache(other, str(tmp_path))


def test_reader_classifies_http_failures(monkeypatch):
    class R:
        def __init__(self, code):
            self.status_code, self.content = code, b""

    def get_404(*a, **k):
        return R(404)

    def get_503(*a, **k):
        return R(503)

    def boom(*a, **k):
        raise requests.ConnectionError("down")

    import snp_phewas_tabix as ct
    for fn, exc in ((get_404, SourceFileMissing), (get_503, SourceUnreachable), (boom, SourceUnreachable)):
        monkeypatch.setattr(ct._session(), "get", fn)
        with pytest.raises(exc):
            read_header("https://example.invalid/x.tsv.bgz")


# --- orientation -----------------------------------------------------------------------

def test_orientation_rules():
    o = orient("G", "T", "G", "T")
    assert o.matches and not o.flip and not o.palindromic
    o = orient("T", "G", "G", "T")
    assert o.matches and o.flip
    o = orient("A", "T", "A", "T")
    assert o.matches and o.palindromic and not o.flip and o.withheld_reason is None
    o = orient("T", "A", "A", "T")
    assert o.matches and o.palindromic and o.withheld_reason == S.WITHHELD_PALINDROMIC_ORDER
    assert not orient("G", "C", "G", "T").matches


# --- lookup ----------------------------------------------------------------------------

def test_lookup_reads_both_datasets_and_accounts_for_every_catalog_outcome(offline_kwargs):
    r = lookup_variant("rs12740374", **offline_kwargs)
    assert r.variant.variant_id_grch38 == "1_109274968_G_T"
    pan = r.dataset_rows("pan_ukbb_eur")
    assert len(pan) == 7160 and len(r.dataset_rows("finngen_r12")) == 2469
    assert r.status_counts("pan_ukbb_eur") == {"ok": 3, "not_queried": 7157}
    ldl = _row(r, "pan_ukbb_eur", LDL)
    assert ldl.variant_id_native == "1_109817590_G_T" and ldl.effect_allele == "T"
    assert ldl.eaf == pytest.approx(0.2215)
    e78 = _row(r, "pan_ukbb_eur", "icd10-E78-both_sexes")
    assert e78.eaf is not None  # case/control file: controls' frequency
    assert [u.outcome_id for u in r.uncatalogued] == ["K11_GINGIVITIS_PERIODONTAL"]
    assert r.uncatalogued[0].status == S.ST_NOT_IN_CATALOG
    assert any("files in" in n and "read from" in n for n in r.notes)
    assert any("only 3 of 7,160 outcomes were requested" in n for n in r.notes)


def test_lookup_flips_an_allele_swapped_query(offline_kwargs):
    ref = lookup_variant("rs12740374", **offline_kwargs)
    offline_kwargs["readers"] = demo_readers(finngen_route="tabix")
    sw = lookup_variant("1_109274968_T_G", **offline_kwargs)
    for ds, oid in (("pan_ukbb_eur", LDL), ("finngen_r12", "I9_CHD")):
        a, b = _row(ref, ds, oid), _row(sw, ds, oid)
        assert b.allele_flip_applied and not a.allele_flip_applied
        assert b.beta == pytest.approx(-a.beta) and b.eaf == pytest.approx(1 - a.eaf)
        assert b.neglog10_p == a.neglog10_p and b.effect_allele == "G"


def test_palindromic_variant_keeps_direction_in_order_and_withholds_it_reversed(offline_kwargs):
    offline_kwargs["readers"] = demo_readers(finngen_route="tabix")
    kept = lookup_variant("rs41279716", **offline_kwargs)
    k = _row(kept, "finngen_r12", "I9_CHD")
    assert k.palindromic and k.status == "ok" and k.beta is not None
    assert k.direction_withheld_reason is None
    offline_kwargs["readers"] = demo_readers(finngen_route="tabix")
    rev = lookup_variant("1_109270888_T_A", **offline_kwargs)
    for ds, oid in (("finngen_r12", "I9_CHD"), ("pan_ukbb_eur", LDL)):
        w = _row(rev, ds, oid)
        assert w.palindromic and w.status == "ok"
        assert w.beta is None and w.eaf is None and w.se is not None
        assert w.direction_withheld_reason == S.WITHHELD_PALINDROMIC_ORDER
    assert any("no direction shown" in n for n in rev.notes)


def test_a_different_variant_at_the_position_is_an_allele_mismatch(offline_kwargs):
    offline_kwargs["readers"] = demo_readers(finngen_route="tabix")
    r = lookup_variant("1_109274968_G_C", **offline_kwargs)
    assert _row(r, "pan_ukbb_eur", LDL).status == S.ST_ALLELE_MISMATCH
    assert _row(r, "finngen_r12", "I9_CHD").status == S.ST_ALLELE_MISMATCH


def test_finngen_api_route_labels_censored_values_not_untested(offline_kwargs):
    offline_kwargs["readers"] = demo_readers(finngen_route="api")
    offline_kwargs["outcomes"]["finngen_r12"] = FG_OUTCOMES + ["AB1_ACTINOMYCOSIS"]
    r = lookup_variant("rs12740374", **offline_kwargs)
    assert _row(r, "finngen_r12", "AB1_ACTINOMYCOSIS").status == S.ST_API_CENSORED
    assert _row(r, "finngen_r12", "I9_CHD").status == S.ST_OK


def test_hybrid_route_fills_censored_endpoints_from_the_files(offline_kwargs):
    api = DemoFinnGenApi()
    api.payload["results"] = [dict(x, beta=None, sebeta=None, pval=None, mlogp=None)
                              if x["phenocode"] == "E4_LIPOPROT" else x
                              for x in api.payload["results"]]
    offline_kwargs["readers"] = demo_readers(api=api)
    r = lookup_variant("rs12740374", **offline_kwargs)
    e4 = _row(r, "finngen_r12", "E4_LIPOPROT")
    assert e4.status == "ok" and e4.source_route == "tabix" and e4.beta == pytest.approx(-0.128558)
    assert _row(r, "finngen_r12", "I9_CHD").source_route == "api"
    assert r.provenance["finngen_r12"]["route"] == "api+tabix"


def test_hybrid_confirms_an_api_absence_against_a_file(offline_kwargs):
    r = lookup_variant("1_109274968_G_C", **offline_kwargs)  # API 404 both orders
    assert _row(r, "finngen_r12", "I9_CHD").status == S.ST_ALLELE_MISMATCH  # file disagrees
    api = DemoFinnGenApi()
    offline_kwargs["readers"] = demo_readers(api=api)
    r = lookup_variant("1_109284000_A_G", **offline_kwargs)  # in no file
    assert {x.status for x in r.dataset_rows("finngen_r12") if x.outcome_id in FG_OUTCOMES} == {S.ST_ABSENT}
    assert r.provenance["finngen_r12"]["n_files_read"] == 1


def test_an_outage_is_never_served_as_absent(offline_kwargs, monkeypatch, tmp_path):
    import snp_phewas_tabix as ct

    def down(*a, **k):
        raise SourceUnreachable("could not read: ConnectTimeout")

    monkeypatch.setattr(ct, "_http_get", down)
    offline_kwargs["readers"] = demo_readers(finngen_route="api", api=DemoFinnGenApi(outage=True),
                                             pan_index_base="https://example.invalid/",
                                             index_cache=str(tmp_path))
    r = lookup_variant("rs12740374", **offline_kwargs)
    fg = {x.status for x in r.dataset_rows("finngen_r12") if x.outcome_id in FG_OUTCOMES}
    assert fg == {S.ST_UNREACHABLE}
    pan = {x.status for x in r.dataset_rows("pan_ukbb_eur") if x.outcome_id in PAN_OUTCOMES}
    assert pan == {S.ST_UNREACHABLE}
    assert any("not evidence of no association" in n for n in r.notes)


def test_hybrid_reads_the_files_when_the_api_is_down(offline_kwargs):
    """Measured 2026-10-18: the FinnGen API answered HTTP 502 three times in a row during
    a live run while the endpoint files (another host) were readable."""
    offline_kwargs["readers"] = demo_readers(api=DemoFinnGenApi(outage=True))
    r = lookup_variant("rs12740374", **offline_kwargs)
    chd = _row(r, "finngen_r12", "I9_CHD")
    assert chd.status == S.ST_OK and chd.source_route == "tabix"
    assert chd.beta == pytest.approx(-0.0801145)
    prov = r.provenance["finngen_r12"]
    assert prov["api_outage_fallback"] and prov["n_files_read"] == 3
    assert any("web API could not be reached" in n for n in r.notes)


def test_an_api_and_file_outage_says_so_and_never_claims_the_files_were_read(offline_kwargs, tmp_path):
    from snp_phewas_catalog import get_dataset
    from snp_phewas_sources import FinnGenReader

    readers = demo_readers(api=DemoFinnGenApi(outage=True))
    readers["finngen_r12"] = FinnGenReader(get_dataset("finngen_r12"), route="hybrid",
                                           client=DemoFinnGenApi(outage=True),
                                           data_base=str(tmp_path) + "/", workers=1)
    offline_kwargs["readers"] = readers
    r = lookup_variant("rs12740374", **offline_kwargs)
    assert S.ST_ABSENT not in {x.status for x in r.dataset_rows("finngen_r12")
                               if x.outcome_id in FG_OUTCOMES}
    assert any("files could not be read either" in n for n in r.notes), r.notes
    assert not any("read from the files instead" in n for n in r.notes)


def test_build_mapping_failures_are_classified(offline_kwargs):
    def outage(*a, **k):
        raise SourceUnreachable("Ensembl assembly mapper: ReadTimeout")

    offline_kwargs["map_fn"] = outage
    r = lookup_variant("rs12740374", **offline_kwargs)
    assert {x.status for x in r.dataset_rows("pan_ukbb_eur") if x.outcome_id in PAN_OUTCOMES} == {S.ST_UNREACHABLE}
    offline_kwargs["map_fn"] = lambda *a, **k: []
    r = lookup_variant("rs12740374", **offline_kwargs)
    assert {x.status for x in r.dataset_rows("pan_ukbb_eur") if x.outcome_id in PAN_OUTCOMES} == {S.ST_BUILD_UNMAPPABLE}
    assert _row(r, "finngen_r12", "I9_CHD").status == "ok"


def test_grch37_input_is_mapped_to_grch38(offline_kwargs):
    r = lookup_variant("1:109817590:G:T", build="GRCh37", **offline_kwargs)
    assert r.variant.variant_id_grch38 == "1_109274968_G_T"


class _CountingReader:
    """Stands in for a dataset reader: records which outcomes were asked for and answers
    'absent' for each, so the scope tests need no network."""

    def __init__(self, finngen: bool = False):
        self.asked: list[str] = []
        self.finngen = finngen  # the FinnGen reader also returns results outside the catalog

    def read(self, outcomes, tgt, catalog_ids=None):
        from snp_phewas_sources import RawRecord

        self.asked = [o.outcome_id for o in outcomes]
        recs = {o.outcome_id: RawRecord(S.ST_ABSENT, route="tabix") for o in outcomes}
        return (recs, {"route": "tabix"}, []) if self.finngen else (recs, {"route": "tabix"})


def test_default_scope_reads_the_matched_outcomes_on_both_sides_and_says_so(offline_kwargs):
    """Two views only: the outcomes matched across the two biobanks (default) or every
    outcome in both (--all-outcomes)."""
    from snp_phewas_crosswalk import harmonized_finngen_endpoints, harmonized_pan_ukbb_outcomes

    pan, fg = _CountingReader(), _CountingReader(finngen=True)
    offline_kwargs.update(readers={"pan_ukbb_eur": pan, "finngen_r12": fg}, outcomes=None)
    r = lookup_variant("rs12740374", datasets="pan_ukbb_eur,finngen_r12", **offline_kwargs)
    assert len(pan.asked) == 589 and set(pan.asked) == set(harmonized_pan_ukbb_outcomes())
    assert len(fg.asked) == 622 and set(fg.asked) == set(harmonized_finngen_endpoints())
    assert r.status_counts("pan_ukbb_eur") == {"absent_in_dataset": 589, "not_queried": 7160 - 589}
    assert r.provenance["pan_ukbb_eur"]["scope"] == "harmonized"
    assert r.provenance["finngen_r12"]["scope"] == "harmonized"
    assert any("589 outcomes matched to the other biobank" in n and "--all-outcomes" in n
               for n in r.notes)
    assert any("622 outcomes matched to the other biobank" in n for n in r.notes)
    every, every_fg = _CountingReader(), _CountingReader(finngen=True)
    offline_kwargs.update(readers={"pan_ukbb_eur": every, "finngen_r12": every_fg})
    r = lookup_variant("rs12740374", datasets="pan_ukbb_eur,finngen_r12", scope="all",
                       **offline_kwargs)
    assert (len(every.asked), len(every_fg.asked)) == (7160, 2469)
    assert r.provenance["pan_ukbb_eur"]["scope"] == "all"


def test_pan_ukbb_alone_is_read_in_full_because_there_is_no_replication_view(offline_kwargs):
    """The matched-outcome cut exists for the cross-biobank comparison; with no FinnGen in the
    lookup its stated reason would be false, so Pan-UK Biobank alone reads every outcome."""
    pan = _CountingReader()
    offline_kwargs.update(readers={"pan_ukbb_eur": pan}, outcomes=None)
    r = lookup_variant("rs12740374", datasets="pan_ukbb_eur", **offline_kwargs)
    assert len(pan.asked) == 7160 and r.provenance["pan_ukbb_eur"]["scope"] == "all"
    assert not any("matched to the other biobank" in n for n in r.notes)
    with pytest.raises(ValueError, match="scope"):
        lookup_variant("rs12740374", scope="some", **offline_kwargs)


def test_an_explicit_outcome_list_overrides_the_default_scope(offline_kwargs):
    r = lookup_variant("rs12740374", **offline_kwargs)
    assert r.provenance["pan_ukbb_eur"]["scope"] == "requested"
    assert r.provenance["pan_ukbb_eur"]["n_queried"] == 3


# --- partition + source-id invariants (each constructed to fail) ------------------------

def _good_rows():
    r = lookup_variant("rs12740374", readers=demo_readers(), rsid_resolver=demo_rsid,
                       map_fn=demo_map, datasets="finngen_r12",
                       outcomes={"finngen_r12": FG_OUTCOMES}, harmonize=False)
    return r.rows, {o.outcome_id for o in load_outcome_catalog("finngen_r12")}


def test_partition_assert_fires_on_each_constructed_failure():
    rows, cat = _good_rows()
    S.assert_partition("finngen_r12", rows, cat)
    with pytest.raises(S.PartitionError, match="have no row"):
        S.assert_partition("finngen_r12", rows[1:], cat)
    with pytest.raises(S.PartitionError, match="duplicate"):
        S.assert_partition("finngen_r12", rows + rows[:1], cat)
    extra = S.OutcomeRow(**{**rows[0].__dict__, "outcome_id": "NOT_AN_ENDPOINT"})
    with pytest.raises(S.PartitionError, match="outside the catalog"):
        S.assert_partition("finngen_r12", rows[1:] + [extra], cat)
    bad = S.OutcomeRow(**{**rows[0].__dict__, "status": "vanished"})
    with pytest.raises(S.PartitionError, match="unknown statuses"):
        S.assert_partition("finngen_r12", [bad] + rows[1:], cat)


def test_source_id_assert_fires_when_the_native_id_is_dropped(tmp_path):
    rows, _ = _good_rows()
    ok = [r for r in rows if r.status == "ok"]
    S.write_rows_tsv(rows, str(tmp_path / "rows.tsv"))
    dropped = S.OutcomeRow(**{**ok[0].__dict__, "variant_id_native": None})
    with pytest.raises(S.SourceIdMissingError):
        S.write_rows_tsv([dropped], str(tmp_path / "bad.tsv"))


def test_a_reader_that_drops_a_queried_outcome_fails_loudly(offline_kwargs):
    class Lossy:
        def read(self, outcomes, tgt):
            return {}, {"route": "tabix"}

    offline_kwargs["readers"] = {**offline_kwargs["readers"], "pan_ukbb_eur": Lossy()}
    with pytest.raises(S.PartitionError, match="no record"):
        lookup_variant("rs12740374", **offline_kwargs)


# --- FinnGen API client ----------------------------------------------------------------

def test_pheweb_client_classifies_responses(monkeypatch):
    import snp_phewas_tabix as ct

    monkeypatch.setattr(ct, "RETRY_BACKOFF_S", 0.0)

    class Sess:
        def __init__(self, code=200, body=None, exc=None, codes=None):
            self.code, self.body, self.exc = code, body, exc
            self.codes, self.calls = list(codes or []), 0

        def get(self, url, timeout):
            self.calls += 1
            if self.exc:
                raise self.exc
            code = self.codes.pop(0) if self.codes else self.code

            class R:
                status_code = code

                def json(_s):
                    if self.body is None:
                        raise ValueError("no json")
                    return self.body
            return R()

    c = FinnGenPhewebClient(session=Sess(404))
    assert c.fetch_variant("1", 1, "A", "G") is None
    for s in (Sess(500), Sess(429), Sess(exc=requests.ReadTimeout()), Sess(200, None), Sess(200, {"x": 1})):
        with pytest.raises(SourceUnreachable):
            FinnGenPhewebClient(session=s).fetch_variant("1", 1, "A", "G")
    assert FinnGenPhewebClient(session=Sess(200, {"results": []})).fetch_variant("1", 1, "A", "G") == {"results": []}
    # One transient 503 is retried; a persistent one is an outage after 3 attempts.
    flaky = Sess(200, {"results": []}, codes=[503])
    assert FinnGenPhewebClient(session=flaky).fetch_variant("1", 1, "A", "G") == {"results": []}
    assert flaky.calls == 2
    down = Sess(503)
    with pytest.raises(SourceUnreachable):
        FinnGenPhewebClient(session=down).fetch_variant("1", 1, "A", "G")
    assert down.calls == 3
    gone = Sess(404)
    FinnGenPhewebClient(session=gone).fetch_variant("1", 1, "A", "G")
    assert gone.calls == 1  # a 404 is an answer, never retried


# --- build mapping on the opposite strand, and multi-allelic rsIDs ----------------------

def _minus_strand_map(from_build, chrom, start, end, to_build):
    return [{"mapped": {"coord_system": "chromosome", "seq_region_name": chrom,
                        "start": 5000, "end": 5000 + (end - start), "strand": -1}}]


def test_opposite_strand_mapping_reverse_complements_and_refuses_indels():
    from snp_phewas_variant import Placement, map_placement

    snv, _ = map_placement(Placement("GRCh38", "1", 100, "G", "T"), "GRCh37", map_fn=_minus_strand_map)
    assert (snv.ref, snv.alt, snv.strand_vs_grch38) == ("C", "A", -1)
    mnv, _ = map_placement(Placement("GRCh38", "1", 100, "AG", "CT"), "GRCh37", map_fn=_minus_strand_map)
    assert (mnv.ref, mnv.alt) == ("CT", "AG")  # reverse complement, not complement (TC/GA)
    indel, note = map_placement(Placement("GRCh38", "1", 100, "AG", "A"), "GRCh37",
                                map_fn=_minus_strand_map)
    assert indel is None and "insertion or deletion" in note


def test_a_multi_allelic_rsid_says_which_allele_it_answered_for():
    from snp_phewas_variant import resolve_variant

    v = resolve_variant("rs999", rsid_resolver=lambda q: ["1_100_G_T", "1_100_G_A", "2_5_C_T"])
    assert v.variant_id_grch38 == "1_100_G_T"
    assert "1_100_G_A" in v.resolution_note and "2_5_C_T" not in v.resolution_note
    single = resolve_variant("rs999", rsid_resolver=lambda q: "1_100_G_T")
    assert "other allele" not in single.resolution_note


def test_every_outcome_label_is_unique_within_a_dataset():
    """Two outcomes with one label read as a duplicate row (measured 2026-10-18 on the first
    live answer: two "LDL direct, adjusted by medication" rows, rank-normalized and raw)."""
    import collections

    for dsid in load_datasets():
        labels = collections.Counter(o.outcome_label for o in load_outcome_catalog(dsid))
        assert not [k for k, v in labels.items() if v > 1], dsid
    pan = {o.outcome_id: o.outcome_label for o in load_outcome_catalog("pan_ukbb_eur")}
    assert pan["categorical-20002-both_sexes-1473"] == "Non-cancer illness code, self-reported: high cholesterol"
    assert "rank-normalized" in pan["continuous-LDLC-both_sexes-medadj_irnt"]
    assert "raw units" in pan["continuous-LDLC-both_sexes-medadj_raw"]
    assert pan["biomarkers-30780-both_sexes-irnt"] == "LDL direct"


def test_phenome_wide_denominator_counts_what_was_asked_not_only_what_has_a_value():
    """An endpoint the API withholds above p = 0.05 was tested, and an outage was asked about:
    counting only rows with a value let either loosen the threshold (the --quick report said
    0.05 / 195 for a lookup over 2,469 endpoints)."""
    def row(status):
        return S.OutcomeRow(dataset_id="d", dataset_release="r", outcome_id=status,
                            outcome_label="o", outcome_category_native="", outcome_kind="binary",
                            n_cases=1, n_controls=1, n_total=2, variant_id_grch38="1_1_A_G",
                            variant_id_native=None, rs_id=None, effect_allele="G",
                            other_allele="A", status=status)

    counted = [S.ST_OK, S.ST_LOW_CONFIDENCE, S.ST_API_CENSORED, S.ST_UNREACHABLE,
               S.ST_FILE_MISSING, S.ST_RECORD_UNPARSABLE]
    excluded = [S.ST_ABSENT, S.ST_NOT_TESTED, S.ST_ALLELE_MISMATCH, S.ST_BUILD_UNMAPPABLE,
                S.ST_NOT_QUERIED, S.ST_NOT_IN_CATALOG]
    assert S.phenome_wide_n([row(s) for s in counted + excluded]) == len(counted)


def test_an_open_targets_reply_with_no_data_is_an_outage_not_no_variant(monkeypatch):
    import snp_phewas_ot_resolver as R

    class Resp:
        status_code = 200

        def json(self):
            return {"data": None}

    class Sess:
        def post(self, *a, **k):
            return Resp()

    with pytest.raises(R.OpenTargetsError, match="no data"):
        R.search_variant_ids("rs12740374", session=Sess())


def test_a_full_index_cache_disk_is_an_unreadable_outcome_not_a_crashed_lookup(tmp_path, monkeypatch):
    import errno

    import snp_phewas_sources as src

    def full(*a, **k):
        raise OSError(errno.ENOSPC, "No space left on device")

    monkeypatch.setattr(src, "ensure_index", full)
    rec = _read_fg(*_fg_file(tmp_path, b"1\t109274968\tG\tT\t-0.08\t0.0086\t19.7\t0.2\n"))
    assert rec.status == S.ST_UNREACHABLE and "No space left" in rec.detail, rec


def test_a_file_missing_a_contract_column_is_unparsable_not_a_crashed_lookup(tmp_path):
    data = tmp_path / "fg.gz"
    header = b"chrom\tpos\talt\tbeta\tsebeta\tmlogp\taf_alt\n"  # no ref column
    rec = b"1\t109274968\tT\t-0.08\t0.0086\t19.7\t0.2\n"
    data.write_bytes(_bgzf_block(header + rec) + _bgzf_block(b""))
    idx = tmp_path / "fg.gz.tbi"
    idx.write_bytes(_tbi(["1"], [len(header)]))
    r = _read_fg(str(data), str(idx))
    assert r.status == S.ST_RECORD_UNPARSABLE and "'ref'" in r.detail, r
