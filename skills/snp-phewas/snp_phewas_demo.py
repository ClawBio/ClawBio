"""Offline demo data: rs12740374 over vendored source slices, through `readers=`.

The slices under `demo/` are +/-10 kb windows around rs12740374 cut from three Pan-UKBB
outcome files (GRCh37 1:109,807,590-109,827,590) and three FinnGen R12 endpoint files
(GRCh38 1:109,264,968-109,284,968), each keeping its file's own header line and re-indexed
with htslib, plus a trimmed copy of the FinnGen PheWeb variant API answer for
`1:109274968-G-T`. They are REAL source rows, not synthetic values (see `demo/README.md`
for provenance and the required acknowledgments). The window also holds the
strand-ambiguous A/T variant rs41279716, which the palindromic tests use.

`demo_readers()` wires the skill's own readers to these files, and `demo_rsid` /
`demo_map` stand in for the Open Targets and Ensembl calls, so a demo run touches no
network and writes nothing outside `--output`.
"""

from __future__ import annotations

import json
import os

from snp_phewas_catalog import get_dataset
from snp_phewas_sources import FinnGenReader, PanUkbbReader
from snp_phewas_tabix import SourceUnreachable

DEMO_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "demo")
PAN_DIR = os.path.join(DEMO_DIR, "pan_ukbb") + os.sep
FG_DIR = os.path.join(DEMO_DIR, "finngen") + os.sep
API_JSON = os.path.join(FG_DIR, "api_variant_1-109274968-G-T.json")

DEMO_VARIANT = "rs12740374"
PAN_OUTCOMES = ["biomarkers-30780-both_sexes-irnt", "icd10-E78-both_sexes", "phecode-411-both_sexes"]
FG_OUTCOMES = ["I9_CHD", "E4_LIPOPROT", "I9_IHD"]
DEMO_OUTCOMES = {"pan_ukbb_eur": PAN_OUTCOMES, "finngen_r12": FG_OUTCOMES}

# The demo window: GRCh37 = GRCh38 + 542,622 on chr1 near CELSR2 (Ensembl maps
# rs12740374 1:109274968 -> 1:109817590, strand 1; checked live 2026-10-17).
OFFSET_37 = 109_817_590 - 109_274_968

# rsID -> GRCh38 variant id as Open Targets resolves them (checked live 2026-10-17).
DEMO_RSIDS = {"rs12740374": "1_109274968_G_T", "rs41279716": "1_109270888_A_T"}


def demo_map(from_build, chrom, start, end, to_build):
    """Ensembl assembly-mapper answers for the demo window (a fixed offset, strand 1)."""
    d = OFFSET_37 if (from_build, to_build) == ("GRCh38", "GRCh37") else -OFFSET_37
    return [{"original": {"start": start, "end": end},
             "mapped": {"coord_system": "chromosome", "seq_region_name": chrom,
                        "start": start + d, "end": end + d, "strand": 1}}]


def demo_rsid(rsid: str) -> str:
    if rsid not in DEMO_RSIDS:
        raise KeyError(f"the offline demo knows only {', '.join(DEMO_RSIDS)}")
    return DEMO_RSIDS[rsid]


class DemoFinnGenApi:
    """Answers like the FinnGen PheWeb variant API from the vendored payload; 404 otherwise."""

    def __init__(self, *, outage: bool = False) -> None:
        with open(API_JSON, encoding="utf-8") as fh:
            self.payload = json.load(fh)
        self.outage = outage
        self.calls: list[str] = []

    def fetch_variant(self, chrom, pos, ref, alt):
        self.calls.append(f"{chrom}:{pos}-{ref}-{alt}")
        if self.outage:
            raise SourceUnreachable("FinnGen variant API: ConnectTimeout")
        if (str(chrom), int(pos), ref, alt) == ("1", 109274968, "G", "T"):
            return self.payload
        return None


def demo_readers(*, finngen_route: str = "hybrid", api: DemoFinnGenApi | None = None,
                 pan_index_base: str | None = None, index_cache: str | None = None):
    """The skill's own readers, pointed at the vendored slices."""
    return {
        "pan_ukbb_eur": PanUkbbReader(get_dataset("pan_ukbb_eur"), data_base=PAN_DIR,
                                      index_base=pan_index_base or PAN_DIR, workers=2,
                                      index_cache=index_cache),
        "finngen_r12": FinnGenReader(get_dataset("finngen_r12"), route=finngen_route,
                                     client=api or DemoFinnGenApi(), data_base=FG_DIR, workers=2),
    }


def demo_kwargs(**overrides):
    """Keyword arguments that make `lookup_variant` run fully offline on the demo data."""
    kw = dict(readers=demo_readers(), rsid_resolver=demo_rsid, map_fn=demo_map,
              outcomes={k: list(v) for k, v in DEMO_OUTCOMES.items()})
    kw.update(overrides)
    return kw
