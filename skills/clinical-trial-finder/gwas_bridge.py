"""Bridge to GWAS Catalog -- resolves rsID to disease traits for trial search.

Queries the EBI GWAS Catalog REST API directly (free, no auth) to extract
genome-wide significant trait associations for a variant.  This is lighter
than calling the full gwas-lookup skill and avoids its import-path issues.

If gwas-lookup's result.json already exists (e.g. from a prior run), we
read that instead of re-querying.
"""

import json
import urllib.error
import urllib.request
from pathlib import Path

_GWAS_API = "https://www.ebi.ac.uk/gwas/rest/api"

# gwas-lookup output directory (if available from a prior run)
_GWAS_LOOKUP_DIR = Path(__file__).resolve().parent.parent / "gwas-lookup"


def resolve_rsid(rsid: str, max_traits: int = 5) -> dict:
    """Resolve an rsID to associated disease traits via the GWAS Catalog.

    Returns dict with keys: rsid, traits (list[str]), genes (list[str]).
    Raises ValueError if the API fails or returns no associations.
    """
    # Try reading cached gwas-lookup result first
    cached = _try_cached_result(rsid, max_traits)
    if cached:
        return cached

    # Query GWAS Catalog REST API directly
    return _query_gwas_catalog(rsid, max_traits)


def _try_cached_result(rsid: str, max_traits: int) -> dict | None:
    """Check if gwas-lookup has a cached result.json we can reuse."""
    for candidate in [
        _GWAS_LOOKUP_DIR / "data" / f"demo_{rsid}.json",
        Path(f"/tmp/gwas_lookup_{rsid}") / "result.json",
    ]:
        if candidate.exists():
            try:
                data = json.loads(candidate.read_text())
                merged = data.get("data", {}).get("merged", {})
                if merged:
                    return _extract_traits_and_genes(rsid, merged, max_traits)
            except (json.JSONDecodeError, KeyError):
                continue
    return None


def _query_gwas_catalog(rsid: str, max_traits: int) -> dict:
    """Query EBI GWAS Catalog API for variant-trait associations."""
    # projection=associationBySnp embeds efoTraits inline (default omits them)
    url = f"{_GWAS_API}/singleNucleotidePolymorphisms/{rsid}/associations?projection=associationBySnp"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})

    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())
    except urllib.error.HTTPError as e:
        if e.code == 404:
            raise ValueError(f"Variant {rsid} not found in GWAS Catalog") from None
        raise ValueError(f"GWAS Catalog API error: HTTP {e.code}") from None
    except (urllib.error.URLError, OSError) as e:
        raise ValueError(f"GWAS Catalog API unreachable: {e}") from None

    # Extract traits and genes from associations
    trait_pvals: dict[str, float] = {}
    genes: list[str] = []
    seen_genes: set[str] = set()

    for assoc in data.get("_embedded", {}).get("associations", []):
        pval = assoc.get("pvalue", 1.0)
        if isinstance(pval, str):
            try:
                pval = float(pval)
            except ValueError:
                continue

        # Extract trait names
        for trait in assoc.get("efoTraits", []):
            name = trait.get("trait", "").strip()
            if name and pval < 5e-8:  # genome-wide significance
                key = name.lower()
                if key not in trait_pvals or pval < trait_pvals[key]:
                    trait_pvals[key] = pval

        # Extract gene names from strongest risk allele -> gene mappings
        for locus in assoc.get("loci", []):
            for gene_entry in locus.get("authorReportedGenes", []):
                gene = gene_entry.get("geneName", "").strip()
                if gene and gene not in seen_genes and gene not in ("intergenic", "NR"):
                    seen_genes.add(gene)
                    genes.append(gene)

    if not trait_pvals:
        raise ValueError(
            f"No genome-wide significant associations found for {rsid}. "
            "Try --gene or --query instead."
        )

    # Prioritise disease traits over measurement/biomarker traits for trial relevance.
    # "measurement" and "level" traits rarely match CT.gov condition searches.
    _MEASUREMENT_WORDS = {"measurement", "level", "amount", "ratio", "concentration"}

    def _is_measurement(trait: str) -> bool:
        return any(w in trait for w in _MEASUREMENT_WORDS)

    sorted_traits = sorted(
        trait_pvals.items(),
        key=lambda x: (_is_measurement(x[0]), x[1]),  # diseases first, then by p-value
    )
    traits = [name.title() for name, _ in sorted_traits[:max_traits]]

    return {"rsid": rsid, "traits": traits, "genes": genes}


def _extract_traits_and_genes(rsid: str, merged: dict, max_traits: int) -> dict:
    """Extract traits and genes from gwas-lookup's merged results format."""
    trait_pvals: dict[str, float] = {}

    for assoc in merged.get("gwas_associations", []):
        trait = assoc.get("trait", "").strip()
        pval = assoc.get("pval", 1.0)
        if trait and pval < 5e-8:
            key = trait.lower()
            if key not in trait_pvals or pval < trait_pvals[key]:
                trait_pvals[key] = pval

    for source_hits in merged.get("phewas", {}).values():
        for hit in source_hits:
            trait = hit.get("phenostring", "").strip()
            pval = hit.get("pval", 1.0)
            if trait and pval < 5e-8:
                key = trait.lower()
                if key not in trait_pvals or pval < trait_pvals[key]:
                    trait_pvals[key] = pval

    sorted_traits = sorted(trait_pvals.items(), key=lambda x: x[1])
    traits = [k.title() for k, _ in sorted_traits[:max_traits]]

    genes: list[str] = []
    seen_genes: set[str] = set()
    for eqtl in merged.get("eqtl_associations", []):
        gene = eqtl.get("gene", "").strip()
        if gene and gene not in seen_genes:
            seen_genes.add(gene)
            genes.append(gene)

    if not traits and not genes:
        raise ValueError(f"No significant associations for {rsid}")

    return {"rsid": rsid, "traits": traits, "genes": genes}
