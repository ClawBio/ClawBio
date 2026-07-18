"""ClinGen PVS1 gene-mechanism gate: is loss-of-function an established disease mechanism for a gene?

Step 1 of the ClinGen PVS1 decision tree (Abou Tayoun et al. 2018, PMID 30192042). Naive ACMG
automation applies PVS1 (very strong, pathogenic) to *any* LoF consequence. That over-calls LoF in
genes where haploinsufficiency is not the disease mechanism. This gate withholds PVS1 unless the gene
has ClinGen-curated haploinsufficiency (HI=3) or is gnomAD LoF-constrained (LOEUF<0.35 or pLI>0.9),
and withholds it for recessive genes (HI=30) where a single heterozygous LoF is only a carrier.

Every decision is traceable to the ClinGen HI score + gnomAD constraint + source. Data ships with the
skill under data/ (public ClinGen + gnomAD reference tables); provenance in data/clingen/PROVENANCE.json.
No side effects at import.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

_DATA = Path(__file__).resolve().parent / "data"
_TSV = _DATA / "clingen" / "ClinGen_haploinsufficiency_GRCh38.tsv"
_PROV = _DATA / "clingen" / "PROVENANCE.json"
_CONSTRAINT_TSV = _DATA / "gnomad" / "gnomad_constraint_by_gene.tsv"

# gnomAD LoF-intolerance thresholds (standard clinical cut-offs): a gene is treated as
# haploinsufficient-by-constraint when LOEUF < 0.35 or pLI > 0.9.
_LOEUF_MAX = 0.35
_PLI_MIN = 0.9

_HI_CACHE: dict | None = None
_CONSTRAINT_CACHE: dict | None = None


def _source_url() -> str:
    try:
        return json.loads(_PROV.read_text()).get("source_url", "https://ftp.clinicalgenome.org/")
    except Exception:
        return "https://ftp.clinicalgenome.org/"


def _load() -> dict:
    """gene -> HI score string, parsed once. ClinGen files prefix the header with '#'."""
    global _HI_CACHE
    if _HI_CACHE is not None:
        return _HI_CACHE
    out: dict[str, str] = {}
    if _TSV.exists():
        lines = _TSV.read_text().splitlines()
        hdr = next((i for i, l in enumerate(lines)
                    if "Gene Symbol" in l and "Haploinsufficiency Score" in l), None)
        if hdr is not None:
            rows = csv.DictReader(lines[hdr:], delimiter="\t")
            for r in rows:
                g = (r.get("#Gene Symbol") or r.get("Gene Symbol") or "").strip()
                if g:
                    out[g] = (r.get("Haploinsufficiency Score") or "").strip()
    _HI_CACHE = out
    return out


def _load_constraint() -> dict:
    """gene -> (LOEUF, pLI), parsed once."""
    global _CONSTRAINT_CACHE
    if _CONSTRAINT_CACHE is not None:
        return _CONSTRAINT_CACHE
    out: dict[str, tuple] = {}
    if _CONSTRAINT_TSV.exists():
        for line in _CONSTRAINT_TSV.read_text().splitlines()[1:]:
            parts = line.split("\t")
            if len(parts) >= 3 and parts[0]:
                loeuf = float(parts[1]) if parts[1] else None
                pli = float(parts[2]) if parts[2] else None
                out[parts[0]] = (loeuf, pli)
    _CONSTRAINT_CACHE = out
    return out


def lof_constrained(gene: str) -> dict:
    """gnomAD LoF-intolerance signal for `gene`: constrained (LOEUF<0.35 or pLI>0.9) => a
    haploinsufficiency signal even when the gene is not ClinGen-dosage-curated."""
    loeuf, pli = _load_constraint().get((gene or "").strip(), (None, None))
    constrained = (loeuf is not None and loeuf < _LOEUF_MAX) or (pli is not None and pli > _PLI_MIN)
    return {"constrained": constrained, "loeuf": loeuf, "pli": pli}


def pvs1_applicability(gene: str) -> dict:
    """Whether PVS1 may be applied for LoF in `gene`, with a traceable basis.

    Applicable when LoF is an established or strongly-inferred disease mechanism:
      - ClinGen HI=3 (curated haploinsufficiency) -> applicable (AD).
      - else gnomAD LoF-constrained (LOEUF<0.35 or pLI>0.9) -> applicable (inferred).
    Withheld when:
      - ClinGen HI=30 (autosomal recessive) -> a single het LoF is a carrier.
      - otherwise no established/inferred haploinsufficiency.
    """
    g = (gene or "").strip()
    hi = _load().get(g)
    con = lof_constrained(g)
    src = _source_url()
    base = {"gene": gene, "hi_score": hi, "loeuf": con["loeuf"], "pli": con["pli"],
            "source": src, "inheritance": "unknown"}

    if hi == "30":
        return {**base, "applicable": False, "inheritance": "AR",
                "basis": f"ClinGen HI=30 ({gene} recessive); single het LoF is a carrier; "
                         f"PVS1 withheld without biallelic evidence [{src}]"}
    if hi == "3":
        return {**base, "applicable": True, "inheritance": "AD",
                "basis": f"ClinGen HI=3 (curated haploinsufficiency) for {gene} [{src}]"}
    if con["constrained"]:
        return {**base, "applicable": True, "inheritance": "AD_inferred",
                "basis": f"gnomAD LoF-constrained ({gene}: LOEUF={con['loeuf']}, pLI={con['pli']}); "
                         f"haploinsufficiency inferred from constraint (not ClinGen-curated)"}
    reason = {
        "40": "ClinGen HI=40 (dosage sensitivity unlikely)",
        "2": "ClinGen HI=2 (insufficient evidence)", "1": "ClinGen HI=1 (little evidence)",
        "0": "ClinGen HI=0 (no evidence)",
    }.get(hi, f"{gene} not ClinGen-dosage-curated")
    return {**base, "applicable": False,
            "basis": f"{reason} and not gnomAD LoF-constrained (LOEUF={con['loeuf']}, pLI={con['pli']}); "
                     f"LoF not an established disease mechanism, PVS1 withheld [{src}]"}
