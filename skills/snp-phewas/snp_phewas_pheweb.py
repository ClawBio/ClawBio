"""FinnGen PheWeb variant API: every endpoint's result for one variant in one request.

`GET https://r12.finngen.fi/api/variant/{chrom}:{pos}-{ref}-{alt}` returns a JSON object
whose `results` list carries one entry per endpoint (`phenocode`, `beta`, `sebeta`,
`pval`, `mlogp`, `maf`, `n_case`, `n_control`, ...). Measured 2026-10-17 for rs12740374
(`1:109274968-G-T`): one request, 2,470 results. The `gwas-lookup` skill calls the same
endpoint for its FinnGen leg; this client adds the censoring semantics and the outage
classification below.

**The API is NOT a full-null PheWAS.** It fills beta / SE / p only where p < 0.05 and
returns null for every other endpoint. Measured 2026-10-17 on four variants, each
returning 2,470 results: rs12740374 195 with values (max p 0.0499), rs7412
(`19:44908822-C-T`) 249 (max 0.0497), and two common variants with no I9_CHD signal,
`5:50145377-A-G` 154 (max 0.04997) and `5:50145508-G-A` 157 (max 0.0495). For rs12740374
the endpoint files DO hold the nulls (AB1_ACTINOMYCOSIS p 0.187, AB1_AMOEBIASIS p 0.239 by
tabix). So an API-only read must label the rest "above the API's reporting threshold",
never "not tested"; the default `hybrid` route reads those endpoints from the files.

A variant the API does not hold answers HTTP 404 (measured for `1:109274968-G-A`,
`1:109274968-T-G` and `1:1-A-G`), so a 404 is "not in FinnGen under this allele order",
while a timeout, connection error, 429 or 5xx is an outage, never an absence.
"""

from __future__ import annotations

from typing import Any

import requests

from snp_phewas_tabix import SourceUnreachable, get_with_retry

DEFAULT_API_BASE = "https://r12.finngen.fi/api/variant/"
DEFAULT_TIMEOUT_S = 60.0


class FinnGenPhewebClient:
    def __init__(self, api_base: str = DEFAULT_API_BASE, *,
                 timeout_s: float = DEFAULT_TIMEOUT_S,
                 session: requests.Session | None = None) -> None:
        self.api_base = api_base if api_base.endswith("/") else api_base + "/"
        self.timeout_s = timeout_s
        self.session = session or requests.Session()

    def variant_url(self, chrom: str, pos: int, ref: str, alt: str) -> str:
        return f"{self.api_base}{chrom}:{pos}-{ref}-{alt}"

    def fetch_variant(self, chrom: str, pos: int, ref: str, alt: str) -> dict[str, Any] | None:
        """The API payload, or None when the API holds no such variant (HTTP 404)."""
        url = self.variant_url(chrom, pos, ref, alt)
        try:
            r = get_with_retry(self.session, url, timeout=self.timeout_s)
        except requests.RequestException as e:
            raise SourceUnreachable(f"FinnGen variant API: {type(e).__name__}") from e
        if r.status_code == 404:
            return None
        if r.status_code != 200:
            raise SourceUnreachable(f"FinnGen variant API answered HTTP {r.status_code}")
        try:
            payload = r.json()
        except ValueError as e:
            raise SourceUnreachable("FinnGen variant API returned a non-JSON body") from e
        if not isinstance(payload, dict) or not isinstance(payload.get("results"), list):
            raise SourceUnreachable("FinnGen variant API returned an unexpected shape")
        return payload
