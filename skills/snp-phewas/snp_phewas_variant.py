"""Resolve the query variant to its GRCh38 identity and (when needed) its GRCh37 placement.

* rsID -> GRCh38 `chr_pos_ref_alt` through the Open Targets Platform variant search
  (`snp_phewas_ot_resolver`; the top hit wins and the answer says which variant id was
  used, and names the other alleles when the rsID is multi-allelic).
* A typed `chr_pos_ref_alt` (`_`, `:` or `-` separated, optional `chr`) passes through;
  with `build="GRCh37"` it is first mapped to GRCh38.
* GRCh38 -> GRCh37 for Pan-UKBB: the Ensembl REST assembly mapper
  (`/map/human/GRCh38/<region>/GRCh37`), which returns the mapped strand. A position with
  zero or several primary-chromosome placements, or one whose mapped length differs, is
  `build_unmappable` (a permanent boundary for this variant), while an Ensembl outage
  raises `SourceUnreachable` (an outage, never "unmappable").
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Callable

import requests

from snp_phewas_tabix import SourceUnreachable, get_with_retry
from snp_phewas_orient import complement


def reverse_complement(allele: str) -> str:
    """Opposite-strand spelling of a multi-base allele (read 5' to 3' on that strand)."""
    return complement(allele)[::-1]

ENSEMBL_REST = "https://rest.ensembl.org"
DEFAULT_TIMEOUT_S = 30.0

_VID_RE = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|MT)[_:\-]([0-9]+)[_:\-]([ACGT]+)[_:\-]([ACGT]+)$", re.I)
_RSID_RE = re.compile(r"^rs[0-9]+$", re.I)


class VariantResolutionError(ValueError):
    """The query names no variant we can resolve (plain-language message)."""


@dataclass(frozen=True)
class Placement:
    build: str
    chrom: str
    pos: int
    ref: str
    alt: str
    strand_vs_grch38: int = 1

    @property
    def variant_id(self) -> str:
        return f"{self.chrom}_{self.pos}_{self.ref}_{self.alt}"


@dataclass(frozen=True)
class ResolvedVariant:
    query: str
    grch38: Placement
    rs_id: str | None
    resolution_note: str

    @property
    def variant_id_grch38(self) -> str:
        return self.grch38.variant_id


def parse_variant_id(s: str) -> tuple[str, int, str, str] | None:
    m = _VID_RE.match((s or "").strip())
    if not m:
        return None
    chrom, pos, ref, alt = m.groups()
    return chrom.upper(), int(pos), ref.upper(), alt.upper()


MapFn = Callable[[str, str, int, int, str], list[dict[str, Any]]]


def ensembl_map(from_build: str, chrom: str, start: int, end: int, to_build: str, *,
                session: requests.Session | None = None,
                timeout_s: float = DEFAULT_TIMEOUT_S) -> list[dict[str, Any]]:
    """Raw Ensembl assembly mappings for a region (list of {original, mapped})."""
    url = f"{ENSEMBL_REST}/map/human/{from_build}/{chrom}:{start}..{end}:1/{to_build}"
    sess = session or requests.Session()
    try:
        r = get_with_retry(sess, url, headers={"Content-Type": "application/json"},
                           timeout=timeout_s)
    except requests.RequestException as e:
        raise SourceUnreachable(f"Ensembl assembly mapper: {type(e).__name__}") from e
    if r.status_code != 200:
        raise SourceUnreachable(f"Ensembl assembly mapper answered HTTP {r.status_code}")
    return list(r.json().get("mappings") or [])


@lru_cache(maxsize=4096)
def _cached_map(from_build: str, chrom: str, start: int, end: int, to_build: str) -> tuple:
    return tuple(_freeze(m) for m in ensembl_map(from_build, chrom, start, end, to_build))


def _freeze(m: dict) -> tuple:
    return tuple(sorted((k, tuple(sorted(v.items())) if isinstance(v, dict) else v)
                        for k, v in m.items()))


def _thaw(t: tuple) -> dict:
    return {k: dict(v) if isinstance(v, tuple) else v for k, v in t}


def _default_map(from_build: str, chrom: str, start: int, end: int, to_build: str) -> list[dict]:
    return [_thaw(t) for t in _cached_map(from_build, chrom, start, end, to_build)]


def map_placement(p: Placement, to_build: str, *, map_fn: MapFn | None = None
                  ) -> tuple[Placement | None, str]:
    """Map a placement to another build. Returns (placement or None, plain-language note)."""
    map_fn = map_fn or _default_map
    end = p.pos + len(p.ref) - 1
    maps = map_fn(p.build, p.chrom, p.pos, end, to_build)
    prim = [m for m in maps
            if m.get("mapped", {}).get("coord_system") == "chromosome"
            and str(m["mapped"].get("seq_region_name")) == p.chrom]
    if len(prim) != 1:
        return None, (f"{p.variant_id} ({p.build}) has {len(prim)} placements on "
                      f"chromosome {p.chrom} in {to_build}, so it cannot be looked up there")
    mp = prim[0]["mapped"]
    if int(mp["end"]) - int(mp["start"]) != end - p.pos:
        return None, (f"{p.variant_id} ({p.build}) maps to a region of a different length "
                      f"in {to_build}, so it cannot be looked up there")
    strand = int(mp.get("strand", 1))
    if strand == -1 and len(p.ref) != len(p.alt):
        # An indel on the opposite strand needs re-anchoring on a different flanking base,
        # which a coordinate mapping does not supply; looking up a guess would read as a
        # false allele mismatch.
        return None, (f"{p.variant_id} ({p.build}) is an insertion or deletion that maps to "
                      f"the opposite strand in {to_build}, so it cannot be placed there reliably")
    ref, alt = (p.ref, p.alt) if strand == 1 else (reverse_complement(p.ref),
                                                   reverse_complement(p.alt))
    out = Placement(build=to_build, chrom=p.chrom, pos=int(mp["start"]), ref=ref, alt=alt,
                    strand_vs_grch38=strand)
    note = f"{p.build} {p.chrom}:{p.pos} is {to_build} {out.chrom}:{out.pos}"
    if strand == -1:
        note += " on the opposite strand (alleles reverse-complemented)"
    return out, note


def resolve_variant(query: str, *, build: str = "GRCh38",
                    rsid_resolver: Callable[[str], str] | None = None,
                    map_fn: MapFn | None = None) -> ResolvedVariant:
    """Resolve an rsID or a typed variant id to its GRCh38 placement."""
    q = (query or "").strip()
    parsed = parse_variant_id(q)
    if parsed is not None:
        chrom, pos, ref, alt = parsed
        if build == "GRCh38":
            p = Placement("GRCh38", chrom, pos, ref, alt)
            return ResolvedVariant(q, p, None, f"variant {p.variant_id} (GRCh38) as given")
        if build != "GRCh37":
            raise VariantResolutionError(f"genome build {build!r} is not supported; use GRCh38 or GRCh37")
        p37 = Placement("GRCh37", chrom, pos, ref, alt)
        p38, note = map_placement(p37, "GRCh38", map_fn=map_fn)
        if p38 is None:
            raise VariantResolutionError(note)
        p38 = Placement("GRCh38", p38.chrom, p38.pos, p38.ref, p38.alt)
        return ResolvedVariant(q, p38, None, f"variant {p37.variant_id} given on GRCh37; {note}")
    if _RSID_RE.match(q):
        resolver = rsid_resolver or _ot_resolver()
        got = resolver(q.lower())
        vids = [got] if isinstance(got, str) else list(got)
        if not vids:
            raise VariantResolutionError(f"Open Targets has no variant for {q}")
        vid = vids[0]
        parsed = parse_variant_id(vid)
        if parsed is None:
            raise VariantResolutionError(f"{q} resolved to {vid!r}, which is not a variant id")
        chrom, pos, ref, alt = parsed
        p = Placement("GRCh38", chrom, pos, ref, alt)
        note = f"{q} is {p.variant_id} (GRCh38), resolved through Open Targets"
        others = [v for v in vids[1:] if (parse_variant_id(v) or (None, None))[:2] == (chrom, pos)]
        if others:
            note += (f"; {q} also names {len(others)} other allele(s) at this position "
                     f"({', '.join(others)}), and this answer is for {p.variant_id} only. "
                     f"Ask with one of those ids to look it up instead")
        return ResolvedVariant(q, p, q.lower(), note)
    raise VariantResolutionError(
        f"{query!r} is neither an rsID (rs123) nor a variant id (1_109274968_G_T)"
    )


def _ot_resolver() -> Callable[[str], list[str]]:
    from snp_phewas_ot_resolver import OpenTargetsError, search_variant_ids

    def _resolve(rsid: str) -> list[str]:
        try:
            return search_variant_ids(rsid)
        except OpenTargetsError as e:
            raise SourceUnreachable(f"Open Targets variant search: {e}") from e

    return _resolve
