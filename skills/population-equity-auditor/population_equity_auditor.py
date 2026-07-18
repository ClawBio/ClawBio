#!/usr/bin/env python3
"""Population Equity Auditor.

Audit a cohort's small-variant calls for *reference-panel referral inflation*: the tendency of
ACMG automation, run with a Euro-biased population reference (gnomAD), to mis-handle variants that
are common in an under-represented population but under-sampled in the reference.

For every variant the auditor classifies under a 2-axis configuration grid:

    frequency source   x   pathogenic-evidence safeguard
    ----------------       --------------------------------
    gnomAD (global)        naive     (PVS1 applied to any LoF; ACMG/AMP Richards 2015 combining)
    gnomAD (AMR panel)     hardened  (PVS1 gated by ClinGen gene-mechanism; the trustworthy config)
    population (cohort)

It then quantifies, ancestry-stratified and assay-scope-aware:
  - false *actionable* calls: variants called Pathogenic/Likely Pathogenic under naive + gnomAD-blind
    that are NOT actionable under the trustworthy config (hardened + population-aware);
  - which safeguard corrects each (the ClinGen PVS1 gate, population frequency, or both);
  - referral (VUS + actionable) per configuration, per-candidate and frequency-weighted per genome.

Safety invariant (never suppress a true pathogenic): population-aware frequency must NOT silently
downgrade a variant that is actionable from gene-mechanism-gated pathogenic evidence. Such variants
(a possible real founder pathogenic that is common in the cohort) are FLAGGED for expert review, not
auto-benigned.

Comparability: each cohort carries a scope descriptor (assay, caller, variant classes, build). Cross-
cohort comparisons are restricted to the shared variant subset (e.g. SNV-only) and the restriction is
reported. No side effects at import; the CLI/`run()` own all IO.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import dataclass, field
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
DISCLAIMER = ("ClawBio is a research and educational tool. It is not a medical device and does not "
              "provide clinical diagnoses. Consult a healthcare professional before making any "
              "medical decisions.")

VERSION = "0.1.0"
SKILL_NAME = "population-equity-auditor"


# ---------------------------------------------------------------------------
# Engine import: reuse the shipped clinical-variant-reporter ACMG engine (DRY).
# ---------------------------------------------------------------------------
def _load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    mod = module_from_spec(spec)
    # Register before exec: dataclasses with string annotations (from __future__ import annotations)
    # resolve fields via sys.modules[cls.__module__], which fails if the module is not registered.
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _acmg_engine():
    candidates = [
        SKILL_DIR.parent / "clinical-variant-reporter" / "acmg_engine.py",
        SKILL_DIR / "acmg_engine.py",  # vendored fallback
    ]
    for c in candidates:
        if c.exists():
            return _load_module("acmg_engine", c)
    raise FileNotFoundError(
        "acmg_engine.py not found; expected the clinical-variant-reporter skill alongside this one.")


ENG = _acmg_engine()
GATE = _load_module("pvs1_gate", SKILL_DIR / "pvs1_gate.py")

ACTIONABLE = frozenset({"Pathogenic", "Likely Pathogenic"})
BENIGN = frozenset({"Benign", "Likely Benign"})
VUS = "Uncertain Significance"
REFERRED = ACTIONABLE | {VUS}          # what triggers a diagnostic workup


# ---------------------------------------------------------------------------
# Cohort scope (comparability guard)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class CohortScope:
    cohort_id: str
    assay: str                          # WGS | WES | array
    caller: str
    variant_classes: frozenset          # e.g. {"SNV"} or {"SNV","indel"}
    build: str
    n_samples: int
    notes: str = ""

    def to_dict(self) -> dict:
        return {"cohort_id": self.cohort_id, "assay": self.assay, "caller": self.caller,
                "variant_classes": sorted(self.variant_classes), "build": self.build,
                "n_samples": self.n_samples, "notes": self.notes}

    def comparable_with(self, other: "CohortScope") -> dict:
        """The shared variant subset two cohorts may be compared on, plus explicit warnings for
        every dimension where they differ (variant classes, assay, build)."""
        shared = self.variant_classes & other.variant_classes
        warnings = []
        only_self = self.variant_classes - other.variant_classes
        only_other = other.variant_classes - self.variant_classes
        for missing, who in ((only_self, self.cohort_id), (only_other, other.cohort_id)):
            for vc in sorted(missing):
                present = self.cohort_id if who == other.cohort_id else other.cohort_id
                warnings.append(
                    f"{vc} present in {present} but absent in {who}; comparison restricted to "
                    f"{sorted(shared)}")
        if self.assay != other.assay:
            warnings.append(f"assay differs ({self.cohort_id}={self.assay}, "
                            f"{other.cohort_id}={other.assay}); referral rates are not directly "
                            f"comparable across assays")
        if self.build != other.build:
            warnings.append(f"genome build differs ({self.build} vs {other.build}); confirm liftover")
        return {"variant_classes": shared, "warnings": warnings}


# Cohort scope is built from each cohort's metadata sidecar (see validate_input); no cohort is hard-coded.


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------
def wilson_ci(k: float, n: float, z: float = 1.96) -> tuple:
    """Wilson 95% score interval for k successes in n (allele count in allele number)."""
    if n <= 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = (z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


# ---------------------------------------------------------------------------
# Classification under one configuration
# ---------------------------------------------------------------------------
def _evidence(rec: dict, af):
    """Build the engine's VariantEvidence with `af` in the frequency slot (drives BA1/BS1/PM2)."""
    cons = rec.get("consequence") or ""
    gc = rec.get("genomic_context") or {}
    return ENG.VariantEvidence(
        chrom=str(gc.get("chrom", rec.get("chrom", ""))).replace("chr", ""),
        pos=int(gc.get("pos", rec.get("pos", 0)) or 0),
        ref=gc.get("ref", rec.get("ref", "")), alt=gc.get("alt", rec.get("alt", "")),
        gene=rec.get("gene", gc.get("gene", "")),
        consequence=cons, gnomad_af=af, cadd_phred=rec.get("cadd_phred"),
        is_missense=(cons == "missense_variant"), is_synonymous=(cons == "synonymous_variant"),
        is_lof=(cons in ENG.LOF_CONSEQUENCES), clinvar_significance="", clinvar_review_stars=0)


def classify_config(rec: dict, af, hardened: bool, pm2_supporting: bool = False) -> dict:
    """Classify one variant with frequency `af`.

    hardened=True applies the ClinGen PVS1 gene-mechanism gate. pm2_supporting=True downgrades PM2
    from moderate (Richards 2015) to supporting (ClinGen SVI 2020, Pejaver et al.); "absent from
    gnomAD" is itself a reference-panel-bias signal for under-sampled ancestries, so the PM2 strength
    epoch is a first-class axis of the audit, not a fixed choice. The engine is not mutated: the PM2
    criterion's strength is adjusted on the local copy of the triggered criteria."""
    import dataclasses
    ev = _evidence(rec, af)
    crits = [c for c in ENG.evaluate_criteria(ev) if c.triggered]
    if pm2_supporting:
        crits = [dataclasses.replace(c, strength="supporting") if c.code == "PM2" else c
                 for c in crits]
    withheld = None
    if hardened and any(c.code == "PVS1" for c in crits):
        verdict = GATE.pvs1_applicability(rec.get("gene", ""))
        if not verdict["applicable"]:
            crits = [c for c in crits if c.code != "PVS1"]
            withheld = verdict["basis"]
    acmg_class = ENG.classify(crits)
    return {"acmg_class": acmg_class, "codes": [c.code for c in crits],
            "pvs1_withheld": withheld, "pm2_supporting": pm2_supporting}


def audit_variant(rec: dict) -> dict:
    """The configuration grid for one variant + the corrections/safety derived from it.

    Three reference-panel-bias mechanisms can each lift a common under-represented-population LoF SNV
    to *actionable* under naive gnomAD-blind automation, and three safeguards each address one:
      - PM2 strength (moderate 2015 -> supporting SVI 2020),
      - the ClinGen PVS1 gene-mechanism gate,
      - population-aware allele frequency (the epoch-invariant fix)."""
    g = rec.get("gnomad_af")
    amr = rec.get("gnomad_af_amr")
    pop = rec.get("cohort_af")
    configs = {
        # naive automation, PM2=moderate (Richards 2015) — the dangerous default of un-hardened tools
        "naive_blind": classify_config(rec, g, hardened=False),
        "naive_blind_amr": classify_config(rec, amr, hardened=False),
        # safeguard 1: PM2 downgraded to supporting (ClinGen SVI 2020), everything else naive+blind
        "naive_blind_pm2sup": classify_config(rec, g, hardened=False, pm2_supporting=True),
        # safeguard 2: ClinGen PVS1 gate, still gnomAD-blind, PM2=moderate
        "hardened_blind": classify_config(rec, g, hardened=True),
        # safeguard 3: population-aware frequency, still naive+moderate
        "naive_pop": classify_config(rec, pop, hardened=False),
        # trustworthy config: gate + population-aware (PM2=moderate, the hardest test)
        "hardened_pop": classify_config(rec, pop, hardened=True),
        # fully modern: gate + population-aware + PM2=supporting
        "hardened_pop_pm2sup": classify_config(rec, pop, hardened=True, pm2_supporting=True),
    }
    c = {k: v["acmg_class"] for k, v in configs.items()}
    naive_blind_actionable = c["naive_blind"] in ACTIONABLE
    trust_actionable = c["hardened_pop"] in ACTIONABLE

    false_actionable = naive_blind_actionable and not trust_actionable
    false_actionable_amr = (c["naive_blind_amr"] in ACTIONABLE) and not trust_actionable
    corrected_by_pm2_strength = naive_blind_actionable and (c["naive_blind_pm2sup"] not in ACTIONABLE)
    corrected_by_gate = naive_blind_actionable and (c["hardened_blind"] not in ACTIONABLE)
    corrected_by_frequency = naive_blind_actionable and (c["naive_pop"] not in ACTIONABLE)
    # SAFETY: gene-mechanism-gated evidence says actionable, population frequency pulls it below.
    # A possible real founder pathogenic masked by cohort frequency -> flag, never silently benign.
    population_masking_flag = (c["hardened_blind"] in ACTIONABLE) and not trust_actionable

    return {
        "gene": rec.get("gene"), "genomic_context": rec.get("genomic_context"),
        "consequence": rec.get("consequence"),
        "gnomad_af": g, "gnomad_af_amr": amr, "cohort_af": pop,
        "cohort_af_ci95": rec.get("cohort_af_ci95"),
        "inheritance_note": rec.get("inheritance_note"),
        "configs": configs,
        "false_actionable": false_actionable,
        "false_actionable_amr": false_actionable_amr,
        "corrected_by_pm2_strength": corrected_by_pm2_strength,
        "corrected_by_gate": corrected_by_gate,
        "corrected_by_frequency": corrected_by_frequency,
        "corrected_by_all_three": (corrected_by_pm2_strength and corrected_by_gate
                                   and corrected_by_frequency),
        "population_masking_flag": population_masking_flag,
        "referred": {k: (cls in REFERRED) for k, cls in c.items()},
    }


# ---------------------------------------------------------------------------
# Cohort-level audit
# ---------------------------------------------------------------------------
_CONFIG_KEYS = ("naive_blind", "naive_blind_amr", "naive_blind_pm2sup", "naive_pop",
                "hardened_blind", "hardened_pop", "hardened_pop_pm2sup")


def audit_cohort(records: list, scope: CohortScope) -> dict:
    audits = [audit_variant(r) for r in records]
    n = len(audits)

    def per_candidate_referral(key):
        return sum(1 for a in audits if a["referred"][key]) / n if n else None

    def fw_referral(key):
        num = den = 0.0
        for a in audits:
            w = 2.0 * (a["cohort_af"] or 0.0)
            den += w
            if a["referred"][key]:
                num += w
        return num / den if den else None

    masked = [a for a in audits if a["population_masking_flag"]]
    n_actionable_pm2sup = sum(
        1 for a in audits
        if a["configs"]["naive_blind_pm2sup"]["acmg_class"] in ACTIONABLE
        and a["configs"]["hardened_pop"]["acmg_class"] not in ACTIONABLE)
    report = {
        "skill": SKILL_NAME, "version": VERSION,
        "cohort_scope": scope.to_dict(),
        "n_variants": n, "variants_processed": n,
        # actionable over-call under naive automation, PM2=moderate (Richards 2015)
        "false_actionable_naive_blind": sum(1 for a in audits if a["false_actionable"]),
        "false_actionable_naive_blind_amr": sum(1 for a in audits if a["false_actionable_amr"]),
        # ... and the residual actionable over-call once PM2 is downgraded to supporting (SVI 2020):
        # if ~0, the *actionable* harm is specific to 2015-strength automation; the residual harm is
        # VUS/referral inflation, which population-aware frequency then clears.
        "false_actionable_pm2_supporting": n_actionable_pm2sup,
        # each safeguard's independent correction of the naive+moderate actionable over-calls
        "corrected_by_pm2_strength": sum(1 for a in audits if a["corrected_by_pm2_strength"]),
        "corrected_by_gate": sum(1 for a in audits if a["corrected_by_gate"]),
        "corrected_by_frequency": sum(1 for a in audits if a["corrected_by_frequency"]),
        "corrected_by_all_three": sum(1 for a in audits if a["corrected_by_all_three"]),
        "referral": {k: per_candidate_referral(k) for k in _CONFIG_KEYS},
        "referral_frequency_weighted": {k: fw_referral(k) for k in _CONFIG_KEYS},
        "safety": {
            "unsafe_frequency_downgrades": len(masked),
            "masked_variants": [{"gene": a["gene"], "genomic_context": a["genomic_context"],
                                 "cohort_af": a["cohort_af"]} for a in masked],
            # invariant holds because every masked variant is surfaced here, never auto-benigned
            "safety_invariant_holds": True,
        },
        "variants": [{
            "gene": a["gene"], "genomic_context": a["genomic_context"],
            "consequence": a["consequence"], "gnomad_af": a["gnomad_af"],
            "gnomad_af_amr": a["gnomad_af_amr"], "cohort_af": a["cohort_af"],
            "cohort_af_ci95": a["cohort_af_ci95"], "inheritance_note": a["inheritance_note"],
            "class_naive_blind": a["configs"]["naive_blind"]["acmg_class"],
            "class_naive_blind_pm2sup": a["configs"]["naive_blind_pm2sup"]["acmg_class"],
            "class_trustworthy": a["configs"]["hardened_pop"]["acmg_class"],
            "false_actionable": a["false_actionable"],
            "corrected_by": [s for s, on in (("pm2_strength", a["corrected_by_pm2_strength"]),
                                             ("gate", a["corrected_by_gate"]),
                                             ("frequency", a["corrected_by_frequency"])) if on],
            "population_masking_flag": a["population_masking_flag"],
        } for a in audits],
        "provenance": {
            "engine": "clinical-variant-reporter/acmg_engine (ACMG/AMP Richards 2015)",
            "pvs1_gate": "ClinGen HI + gnomAD constraint (Abou Tayoun 2018)",
            "clingen_source": GATE._source_url(),
            "reference_baselines": ["gnomAD_global", "gnomAD_AMR"],
        },
        "disclaimer": DISCLAIMER,
        "status": "ok",
    }
    report["content_sha256"] = hashlib.sha256(
        json.dumps(report, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return report


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------
# Optional per-gene inheritance/mechanism footnotes (public gene facts, cohort-agnostic). Extend freely.
INHERITANCE_NOTES = {
    "PIDD1": "recessive (biallelic LoF -> AR phenotype); single het is a carrier",
    "RNASEL": "cancer-susceptibility (HPC1), not high-penetrance Mendelian LoF",
}


def _f(x):
    x = (x or "").strip()
    if x in ("", "NULL", "NA", "."):
        return None
    try:
        return float(x)
    except ValueError:
        return None


def parse_cohort_tsv(path: Path, subpop_sizes: dict) -> list:
    """Parse an annotated per-subpopulation allele-count cohort TSV into auditor records. Cohort-agnostic.

    Columns: LOCATION REF ALT SYMBOL CADD_PHRED Consequence CLIN_SIG <subpopulation AC columns...>
             GNOMAD_AF GNOMAD_AF_AMR GNOMAD_AN GNOMAD_AN_AMR
    `subpop_sizes` maps each subpopulation column name to its sample count N (from the cohort metadata,
    not hard-coded). Subpopulation columns hold ALT allele counts; cohort AF = sum(AC) / (2 * N_total);
    per-subpopulation frequency = AC / (2 * N_subpop). Wilson CI on the cohort AF is attached."""
    total_an = 2 * sum(subpop_sizes.values())
    lines = [l for l in Path(path).read_text().splitlines() if l.strip()]
    header = lines[0].split("\t")
    idx = {name: i for i, name in enumerate(header)}
    subpops = [s for s in subpop_sizes if s in idx]
    records = []
    for ln in lines[1:]:
        f = ln.split("\t")
        chrom, rng = f[idx["LOCATION"]].split(":")
        pos = int(rng.split("-")[0])
        ac_total = sum(int(f[idx[s]] or 0) for s in subpops)
        cohort_af = ac_total / total_an if total_an else None
        strata = {}
        for s in subpops:
            ac = int(f[idx[s]] or 0)
            an = 2 * subpop_sizes[s]
            strata[s] = {"ac": ac, "an": an, "af": ac / an if an else None}
        records.append({
            "gene": f[idx["SYMBOL"]],
            "genomic_context": {"chrom": chrom, "pos": pos, "ref": f[idx["REF"]],
                                "alt": f[idx["ALT"]], "gene": f[idx["SYMBOL"]]},
            "consequence": f[idx["Consequence"]],
            "cadd_phred": _f(f[idx["CADD_PHRED"]]),
            "gnomad_af": _f(f[idx["GNOMAD_AF"]]) if "GNOMAD_AF" in idx else None,
            "gnomad_af_amr": _f(f[idx["GNOMAD_AF_AMR"]]) if "GNOMAD_AF_AMR" in idx else None,
            "cohort_af": cohort_af,
            "cohort_af_ci95": wilson_ci(ac_total, total_an),
            "inheritance_note": INHERITANCE_NOTES.get(f[idx["SYMBOL"]]),
            "strata": strata,
        })
    return records


def load_cohort_metadata(input_path: Path) -> dict | None:
    """Load the cohort metadata sidecar `<stem>.meta.json` (declares subpopulation_sizes + assay/caller/
    build). Keeping cohort structure out of the code is what makes the skill reusable by any hospital."""
    cand = Path(input_path).parent / (Path(input_path).stem + ".meta.json")
    if cand.exists():
        return json.loads(cand.read_text())
    return None


def validate_input(input_path: Path) -> dict:
    """Parse a cohort TSV plus its `<stem>.meta.json` sidecar into records + a CohortScope. The skill is
    cohort-agnostic: assay, caller, build, and subpopulation sizes all come from the metadata."""
    input_path = Path(input_path)
    if not input_path.exists():
        print(f"Error: input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)
    text_lines = [l for l in input_path.read_text().splitlines() if l.strip()]
    if not text_lines or "LOCATION" not in text_lines[0].upper():
        print("Error: expected an annotated cohort TSV with a LOCATION/SYMBOL header", file=sys.stderr)
        sys.exit(1)
    meta = load_cohort_metadata(input_path)
    if not meta or "subpopulation_sizes" not in meta:
        print(f"Error: missing metadata sidecar {input_path.stem}.meta.json with "
              "'subpopulation_sizes' (and assay/caller/build). See data/demo_cohort.meta.json.",
              file=sys.stderr)
        sys.exit(1)
    sizes = meta["subpopulation_sizes"]
    records = parse_cohort_tsv(input_path, sizes)
    scope = CohortScope(
        cohort_id=meta.get("cohort_id", input_path.stem),
        assay=meta.get("assay", "unspecified"), caller=meta.get("caller", "unspecified"),
        variant_classes=frozenset(meta.get("variant_classes", ["SNV"])),
        build=meta.get("build", "unspecified"), n_samples=sum(sizes.values()),
        notes=meta.get("notes", ""))
    return {"records": records, "scope": scope, "source": str(input_path)}


def run_analysis(data: dict) -> dict:
    result = audit_cohort(data["records"], data["scope"])
    result["source"] = data.get("source", "unknown")
    return result


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def _pct(x):
    return "n/a" if x is None else f"{x * 100:.1f}%"


def write_report(result: dict, output_dir: Path) -> None:
    output_dir = Path(output_dir)
    (output_dir / "tables").mkdir(parents=True, exist_ok=True)

    with open(output_dir / "result.json", "w") as fh:
        json.dump(result, fh, indent=2)

    scope = result["cohort_scope"]
    ref = result["referral"]
    fwr = result["referral_frequency_weighted"]
    md = [
        "# Population Equity Auditor Report",
        "",
        f"**Cohort**: {scope['cohort_id']} · {scope['assay']} · {scope['n_samples']} samples · "
        f"build {scope['build']}",
        f"**Caller**: {scope['caller']}",
        f"**Variant classes**: {', '.join(scope['variant_classes'])}"
        + (f"  \n**Scope note**: {scope['notes']}" if scope.get("notes") else ""),
        f"**Variants audited**: {result['n_variants']}",
        "",
        "## Reference-panel referral inflation",
        "",
        f"- False *actionable* (Likely Pathogenic) calls under **naive automation, PM2=moderate "
        f"(Richards 2015)**, gnomAD-global-blind: **{result['false_actionable_naive_blind']}** "
        f"(gnomAD-AMR baseline: {result['false_actionable_naive_blind_amr']})",
        f"- Residual actionable over-calls once **PM2 is downgraded to supporting (ClinGen SVI 2020)**: "
        f"**{result['false_actionable_pm2_supporting']}** "
        f"— if ~0, the *actionable* harm is specific to 2015-strength automation; the residual harm "
        f"is VUS/referral inflation, which population-aware frequency then clears.",
        "",
        "Each safeguard independently corrects the naive+moderate actionable over-calls:",
        "",
        f"- PM2 strength update (moderate → supporting): {result['corrected_by_pm2_strength']}",
        f"- ClinGen PVS1 gene-mechanism gate: {result['corrected_by_gate']}",
        f"- Population-aware allele frequency (epoch-invariant): {result['corrected_by_frequency']}",
        f"- All three: {result['corrected_by_all_three']}",
        "",
        "### Referral rate by configuration (per candidate variant / frequency-weighted per genome)",
        "",
        "| Configuration | Per-candidate | Freq-weighted |",
        "|---|---|---|",
        f"| naive + gnomAD-global, PM2=moderate (blind) | {_pct(ref['naive_blind'])} | "
        f"{_pct(fwr['naive_blind'])} |",
        f"| naive + gnomAD-AMR | {_pct(ref['naive_blind_amr'])} | {_pct(fwr['naive_blind_amr'])} |",
        f"| naive + gnomAD-global, PM2=supporting | {_pct(ref['naive_blind_pm2sup'])} | "
        f"{_pct(fwr['naive_blind_pm2sup'])} |",
        f"| naive + population-aware | {_pct(ref['naive_pop'])} | {_pct(fwr['naive_pop'])} |",
        f"| hardened + gnomAD-global | {_pct(ref['hardened_blind'])} | {_pct(fwr['hardened_blind'])} |",
        f"| **hardened + population-aware (trustworthy)** | **{_pct(ref['hardened_pop'])}** | "
        f"**{_pct(fwr['hardened_pop'])}** |",
        f"| hardened + population-aware + PM2=supporting | {_pct(ref['hardened_pop_pm2sup'])} | "
        f"{_pct(fwr['hardened_pop_pm2sup'])} |",
        "",
        "## Safety invariant",
        "",
        f"- Variants where population frequency would mask gene-mechanism-gated pathogenic evidence "
        f"(flagged for expert review, never auto-benigned): "
        f"**{result['safety']['unsafe_frequency_downgrades']}**",
        f"- Safety invariant holds: {result['safety']['safety_invariant_holds']}",
        "",
        "## Variants corrected by a safeguard (false-actionable under naive + gnomAD-blind, PM2=moderate)",
        "",
        "| Gene | Consequence | gnomAD AF | cohort AF | naive+blind | +PM2-supporting | trustworthy | "
        "corrected by | note |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for v in result["variants"]:
        if v["false_actionable"]:
            md.append(
                f"| {v['gene']} | {v['consequence']} | "
                f"{'absent' if v['gnomad_af'] is None else f'{v['gnomad_af']:.2e}'} | "
                f"{'n/a' if v['cohort_af'] is None else f'{v['cohort_af']:.3f}'} | "
                f"{v['class_naive_blind']} | {v['class_naive_blind_pm2sup']} | "
                f"{v['class_trustworthy']} | {', '.join(v['corrected_by'])} | "
                f"{v.get('inheritance_note') or ''} |")
    md += [
        "",
        f"**Provenance**: {result['provenance']['engine']}; PVS1 gate = "
        f"{result['provenance']['pvs1_gate']}; baselines = "
        f"{', '.join(result['provenance']['reference_baselines'])}.",
        f"**Content hash**: `{result['content_sha256']}`",
        "",
        f"*{DISCLAIMER}*",
        "",
    ]
    (output_dir / "report.md").write_text("\n".join(md))

    # per-variant table
    tsv = ["\t".join(["gene", "chrom", "pos", "consequence", "gnomad_af", "gnomad_af_amr",
                      "cohort_af", "class_naive_blind", "class_naive_blind_pm2sup",
                      "class_trustworthy", "false_actionable", "corrected_by",
                      "population_masking_flag", "inheritance_note"])]
    for v in result["variants"]:
        gc = v["genomic_context"] or {}
        tsv.append("\t".join(str(x) for x in [
            v["gene"], gc.get("chrom", ""), gc.get("pos", ""), v["consequence"],
            v["gnomad_af"], v["gnomad_af_amr"], v["cohort_af"], v["class_naive_blind"],
            v["class_naive_blind_pm2sup"], v["class_trustworthy"], v["false_actionable"],
            ";".join(v["corrected_by"]), v["population_masking_flag"],
            v.get("inheritance_note") or ""]))
    (output_dir / "tables" / "per_variant_audit.tsv").write_text("\n".join(tsv) + "\n")

    # reproducibility bundle: exact command + engine/data provenance + content hash
    repro = output_dir / "reproducibility"
    repro.mkdir(parents=True, exist_ok=True)
    prov = result.get("provenance", {})
    (repro / "commands.sh").write_text(
        "#!/usr/bin/env bash\n"
        "# Reproduce this population-equity-auditor report.\n"
        f"# source: {result.get('source', 'unknown')}\n"
        f"# engine: {prov.get('engine')}\n"
        f"# pvs1_gate: {prov.get('pvs1_gate')}  clingen_source: {prov.get('clingen_source')}\n"
        f"# reference_baselines: {', '.join(prov.get('reference_baselines', []))}\n"
        f"# result content_sha256: {result.get('content_sha256')}\n"
        "python3 population_equity_auditor.py --input "
        f"{result.get('source', '<cohort.tsv>')} --output <out_dir>\n")

    print(f"Report written to {output_dir / 'report.md'}")
    print(f"Results written to {output_dir / 'result.json'}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, dest="input_file", help="Annotated cohort TSV")
    p.add_argument("--output", type=Path, help="Output directory")
    p.add_argument("--demo", action="store_true", help="Run on the bundled Peru high-impact-LoF cohort")
    return p.parse_args()


def _demo_input() -> Path:
    return SKILL_DIR / "data" / "demo_cohort.tsv"


def main():
    args = parse_args()
    if args.demo:
        data = validate_input(_demo_input())
        result = run_analysis(data)
        write_report(result, args.output or Path("/tmp") / "population_equity_auditor" / "demo")
    elif args.input_file:
        data = validate_input(args.input_file)
        result = run_analysis(data)
        write_report(result, args.output or args.input_file.parent / "output")
    else:
        print("Error: provide --input <file> or --demo", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
