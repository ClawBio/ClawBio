"""Hypothesis properties for the two pieces of this skill that convert or parse:
re-expressing an effect on the target allele (`orient` + `apply`), and the tabix index
walker that reads one position out of a bgzipped file (`fetch_position`)."""

from __future__ import annotations

import gzip
import math
import struct
import tempfile
import zlib
from pathlib import Path

from hypothesis import example, given, settings
from hypothesis import strategies as st

import snp_phewas_schema as S
from snp_phewas_orient import PALINDROMIC_PAIRS, apply, orient
from snp_phewas_tabix import fetch_position

BASES = st.sampled_from("ACGT")
BETA = st.floats(-5, 5, allow_nan=False) | st.none()
EAF = st.floats(0, 1, allow_nan=False) | st.none()


def _pair():
    return st.tuples(BASES, BASES).filter(lambda p: p[0] != p[1])


@given(_pair(), BETA, EAF)
def test_same_allele_order_keeps_the_effect(target, beta, eaf):
    ref, alt = target
    o = orient(ref, alt, ref, alt)
    assert o.matches and not o.flip and o.withheld_reason is None
    assert apply(o, beta, eaf) == (beta, eaf)


@given(_pair(), BETA, EAF)
@example(("A", "T"), -0.1, 0.2)
@example(("G", "T"), -0.1, 0.2)
def test_swapped_order_negates_or_withholds_never_both(target, beta, eaf):
    """A swapped record is re-signed onto the target allele, except a strand-ambiguous
    pair, whose direction is withheld (never guessed)."""
    ref, alt = target
    o = orient(alt, ref, ref, alt)
    assert o.matches
    b, f = apply(o, beta, eaf)
    if (ref, alt) in PALINDROMIC_PAIRS:
        assert o.withheld_reason == S.WITHHELD_PALINDROMIC_ORDER and not o.flip
        assert (b, f) == (None, None)
    else:
        assert o.flip and o.withheld_reason is None
        assert b == (None if beta is None else -beta)
        assert f == (None if eaf is None else 1.0 - eaf)


@given(_pair(), BETA, EAF)
def test_flipping_twice_round_trips(target, beta, eaf):
    ref, alt = target
    o = orient(alt, ref, ref, alt)
    if o.withheld_reason:
        return
    b, f = apply(o, *apply(o, beta, eaf))
    assert b == beta
    assert f is None if eaf is None else math.isclose(f, eaf, abs_tol=1e-12)


@given(_pair(), _pair())
def test_a_record_with_other_alleles_never_matches(target, source):
    if set(source) == set(target):
        return
    assert not orient(source[0], source[1], target[0], target[1]).matches


@given(_pair(), st.sampled_from(["a", "c", "g", "t"]))
def test_orientation_ignores_case(target, _):
    ref, alt = target
    assert orient(ref.lower(), alt.lower(), ref, alt) == orient(ref, alt, ref, alt)


# --- the index walker --------------------------------------------------------------------

TILE = 1 << 14
HEADER = b"chr\tpos\tref\talt\n"


def _block(payload: bytes) -> bytes:
    comp = zlib.compressobj(6, zlib.DEFLATED, -15)
    data = comp.compress(payload) + comp.flush()
    head = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
    head += struct.pack("<H", 18 + len(data) + 8 - 1)
    return head + data + struct.pack("<II", zlib.crc32(payload) & 0xFFFFFFFF, len(payload))


def _tbi(linear: list[int]) -> bytes:
    nm = b"1\0"
    body = b"TBI\x01" + struct.pack("<8i", 1, 0, 1, 2, 2, ord("#"), 1, len(nm)) + nm
    body += struct.pack("<i", 0) + struct.pack("<i", len(linear))
    body += struct.pack(f"<{len(linear)}Q", *linear)
    return gzip.compress(body)


@st.composite
def bgzipped_files(draw):
    """Sorted records in random blocks, with empty blocks dropped between some of them, and
    a linear index that sometimes names the EMPTY block before a record (the Pan-UKBB
    index defect: the in-block offset only makes sense in the block after it)."""
    positions = sorted(draw(st.sets(st.integers(1, 6 * TILE), min_size=1, max_size=60)))
    lines = [f"1\t{p}\tA\tG\n".encode() for p in positions]
    cuts = sorted(draw(st.sets(st.integers(1, max(1, len(lines) - 1)), max_size=8)))
    chunks, prev = [], 0
    for c in cuts + [len(lines)]:
        if c > prev:
            chunks.append(lines[prev:c])
            prev = c
    blobs: list[bytes] = [_block(HEADER)]
    starts: dict[int, tuple[int, int]] = {}  # record index -> (coffset, uoffset)
    off = len(blobs[0])
    i = 0
    for chunk in chunks:
        empty_before = draw(st.booleans())
        empty_off = None
        if empty_before:
            e = _block(b"")
            blobs.append(e)
            empty_off = off
            off += len(e)
        u = 0
        for ln in chunk:
            defect = empty_off is not None and draw(st.booleans())
            starts[i] = ((empty_off if defect else off), u)
            u += len(ln)
            i += 1
        b = _block(b"".join(chunk))
        blobs.append(b)
        off += len(b)
    blobs.append(_block(b""))
    n_tiles = (positions[-1] - 1) // TILE + 1
    linear = []
    for t in range(n_tiles):
        first = next((k for k, p in enumerate(positions) if (p - 1) // TILE >= t), None)
        if first is None:
            linear.append(0)
        else:
            c, u = starts[first]
            linear.append((c << 16) | u)
    return positions, b"".join(blobs), _tbi(linear)


@settings(deadline=None, max_examples=60)
@given(bgzipped_files(), st.data())
def test_index_walker_agrees_with_a_full_scan(f, data):
    positions, blob, tbi = f
    q = data.draw(st.sampled_from(positions) | st.integers(1, 6 * TILE + 100))
    with tempfile.TemporaryDirectory() as d:
        bgz, idx = Path(d) / "x.tsv.bgz", Path(d) / "x.tsv.bgz.tbi"
        bgz.write_bytes(blob)
        idx.write_bytes(tbi)
        got = fetch_position(str(bgz), str(idx), "1", q)
    want = [["1", str(q), "A", "G"]] if q in positions else []
    assert got == want
