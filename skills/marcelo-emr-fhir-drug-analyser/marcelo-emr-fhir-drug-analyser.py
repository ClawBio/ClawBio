#!/usr/bin/env python3
"""
marcelo-emr-fhir-drug-analyser
Pulls clinical context from a HAPI FHIR server, intersects with PharmGx results,
and runs targeted ClinPGx lookups to generate a medication review report.
"""

import argparse
import json
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path

import requests

SKILL_DIR = Path(__file__).parent
PHARMGX_SCRIPT = SKILL_DIR.parent / "pharmgx-reporter" / "pharmgx_reporter.py"
CLINPGX_SCRIPT = SKILL_DIR.parent / "clinpgx" / "clinpgx.py"

DISCLAIMER = (
    "*ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions.*"
)

# ---------------------------------------------------------------------------
# FHIR connection
# ---------------------------------------------------------------------------

def load_fhir_config(path: Path) -> dict:
    """Parse patient_fhir.json and validate required fields."""
    raw = json.loads(path.read_text())
    fhir = raw.get("fhir", {})
    if not fhir.get("server_url") or not fhir.get("patient_id"):
        raise ValueError("fhir_connection JSON must contain fhir.server_url and fhir.patient_id")
    return fhir


def _build_session(auth: dict | None) -> requests.Session:
    session = requests.Session()
    session.headers["Accept"] = "application/fhir+json"
    if not auth or auth.get("type", "none") == "none":
        return session
    if auth["type"] == "bearer":
        session.headers["Authorization"] = f"Bearer {auth['token']}"
    elif auth["type"] == "basic":
        session.auth = (auth["username"], auth["password"])
    return session


def _paginate(session: requests.Session, url: str, params: dict) -> list[dict]:
    """Fetch all pages of a FHIR bundle and return the entry list."""
    entries = []
    resp = session.get(url, params=params, timeout=30)
    resp.raise_for_status()
    bundle = resp.json()
    entries.extend(bundle.get("entry", []))
    # follow pagination links
    while True:
        next_url = next(
            (lk["url"] for lk in bundle.get("link", []) if lk.get("relation") == "next"),
            None,
        )
        if not next_url:
            break
        resp = session.get(next_url, timeout=30)
        resp.raise_for_status()
        bundle = resp.json()
        entries.extend(bundle.get("entry", []))
    return [e["resource"] for e in entries if "resource" in e]


# ---------------------------------------------------------------------------
# FHIR resource retrieval
# ---------------------------------------------------------------------------

def fetch_fhir_resources(server_url: str, patient_id: str, auth: dict | None) -> dict:
    """Fetch all relevant FHIR resources for a patient."""
    base = server_url.rstrip("/")
    session = _build_session(auth)
    params_patient = {"patient": patient_id, "_count": "50"}

    print(f"  Connecting to {base} for patient {patient_id}...")

    resources: dict = {
        "patient": None,
        "medication_requests": [],
        "medication_statements": [],
        "conditions": [],
        "observations": [],
        "allergies": [],
        "diagnostic_reports": [],
        "query_timestamp": datetime.now().isoformat(),
        "server_url": server_url,
        "patient_id": patient_id,
    }

    # Patient
    r = session.get(f"{base}/Patient/{patient_id}", timeout=30)
    r.raise_for_status()
    resources["patient"] = r.json()
    print(f"  Patient found: {_patient_name(resources['patient'])}")

    # MedicationRequest
    resources["medication_requests"] = _paginate(
        session, f"{base}/MedicationRequest", params_patient
    )
    print(f"  MedicationRequests: {len(resources['medication_requests'])}")

    # MedicationStatement
    resources["medication_statements"] = _paginate(
        session, f"{base}/MedicationStatement", params_patient
    )
    print(f"  MedicationStatements: {len(resources['medication_statements'])}")

    # Condition
    resources["conditions"] = _paginate(
        session, f"{base}/Condition", params_patient
    )
    print(f"  Conditions: {len(resources['conditions'])}")

    # Observation — only clinically relevant codes
    obs_params = dict(params_patient)
    obs_params["_count"] = "100"
    resources["observations"] = _paginate(session, f"{base}/Observation", obs_params)
    print(f"  Observations: {len(resources['observations'])}")

    # AllergyIntolerance
    resources["allergies"] = _paginate(
        session, f"{base}/AllergyIntolerance", params_patient
    )
    print(f"  AllergyIntolerances: {len(resources['allergies'])}")

    # DiagnosticReport
    resources["diagnostic_reports"] = _paginate(
        session, f"{base}/DiagnosticReport", params_patient
    )
    print(f"  DiagnosticReports: {len(resources['diagnostic_reports'])}")

    return resources


def _patient_name(patient: dict) -> str:
    for n in patient.get("name", []):
        given = " ".join(n.get("given", []))
        family = n.get("family", "")
        full = f"{given} {family}".strip()
        if full:
            return full
    return patient.get("id", "Unknown")


# ---------------------------------------------------------------------------
# Medication normalization
# ---------------------------------------------------------------------------

def _extract_med_name(resource: dict) -> str:
    """Extract a normalized generic medication name from a FHIR resource."""
    # Try medicationCodeableConcept first
    cc = resource.get("medicationCodeableConcept", {})
    if cc:
        # Prefer text display over coded display
        if cc.get("text"):
            return cc["text"].strip().lower()
        for coding in cc.get("coding", []):
            if coding.get("display"):
                return coding["display"].strip().lower()
    # Try contained medication reference display
    if resource.get("medicationReference", {}).get("display"):
        return resource["medicationReference"]["display"].strip().lower()
    return "unknown"


def _med_status_class(resource: dict, resource_type: str) -> str:
    status = resource.get("status", "unknown")
    intent = resource.get("intent", "")  # MedicationRequest only
    if resource_type == "MedicationRequest":
        if intent in ("proposal", "plan") and status == "active":
            return "planned"
        if status in ("active", "on-hold"):
            return "current"
        if status in ("completed", "stopped", "cancelled"):
            return "past"
    if resource_type == "MedicationStatement":
        if status in ("active", "intended"):
            return "current"
        if status in ("completed", "stopped"):
            return "past"
    return "unknown"


def normalize_medications(
    med_requests: list[dict],
    med_statements: list[dict],
    conditions: list[dict],
) -> list[dict]:
    """Merge, deduplicate, and classify medications into current/past/planned."""
    seen: dict[str, dict] = {}  # name → entry

    def _add(resource: dict, rtype: str) -> None:
        name = _extract_med_name(resource)
        if name == "unknown":
            return
        status_class = _med_status_class(resource, rtype)
        if name in seen:
            # prefer more specific classification
            priority = {"current": 3, "planned": 2, "past": 1, "unknown": 0}
            if priority.get(status_class, 0) > priority.get(seen[name]["status"], 0):
                seen[name]["status"] = status_class
                seen[name]["source"].append(rtype)
        else:
            seen[name] = {
                "name": name,
                "status": status_class,
                "source": [rtype],
                "indication": None,
                "fhir_id": resource.get("id"),
            }

    for r in med_requests:
        _add(r, "MedicationRequest")
    for r in med_statements:
        _add(r, "MedicationStatement")

    # Link indications via condition text (simple keyword overlap)
    condition_names = [
        (c.get("code", {}).get("text", "") or "").lower()
        for c in conditions
        if c.get("clinicalStatus", {}).get("coding", [{}])[0].get("code") == "active"
    ]
    for med in seen.values():
        for cond in condition_names:
            if any(word in cond for word in med["name"].split() if len(word) > 3):
                med["indication"] = cond
                break

    return list(seen.values())


def check_allergies(medications: list[dict], allergies: list[dict]) -> list[str]:
    """Return medication names that match an allergy entry."""
    allergen_names = set()
    for a in allergies:
        for coding in a.get("code", {}).get("coding", []):
            if coding.get("display"):
                allergen_names.add(coding["display"].strip().lower())
        text = a.get("code", {}).get("text", "")
        if text:
            allergen_names.add(text.strip().lower())

    flagged = []
    for med in medications:
        name = med["name"]
        if any(name in allergen or allergen in name for allergen in allergen_names):
            flagged.append(name)
    return flagged


# ---------------------------------------------------------------------------
# PharmGx integration
# ---------------------------------------------------------------------------

def run_pharmgx(genotype_path: Path, output_dir: Path) -> Path:
    """Run pharmgx-reporter and return the path to result.json."""
    pharmgx_out = output_dir / "pharmgx_output"
    pharmgx_out.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, str(PHARMGX_SCRIPT),
        "--input", str(genotype_path),
        "--output", str(pharmgx_out),
        "--no-enrich",  # skip ClinPGx inside pharmgx; we do targeted calls ourselves
    ]
    print(f"  Running PharmGx: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"PharmGx failed:\n{result.stderr}")
    result_json = pharmgx_out / "result.json"
    if not result_json.exists():
        raise FileNotFoundError(f"PharmGx did not produce result.json at {result_json}")
    return result_json


def load_pharmgx_result(path: Path) -> dict:
    data = json.loads(path.read_text())
    # result.json wraps data under a "data" key
    return data.get("data", data)


# ---------------------------------------------------------------------------
# Medication-PGx intersection
# ---------------------------------------------------------------------------

def intersect_pgx_with_meds(medications: list[dict], pgx_data: dict) -> list[dict]:
    """
    Match FHIR medications against PharmGx drug_recommendations.
    Returns a list of intersection records with classification and priority.
    """
    drug_recs: dict = pgx_data.get("drug_recommendations", {})
    # Build a flat lookup: drug_name_lower → {classification, gene, phenotype, ...}
    pgx_lookup: dict[str, dict] = {}
    for cls in ("avoid", "caution", "standard", "indeterminate"):
        for entry in drug_recs.get(cls, []):
            drug_name = (entry.get("drug") or "").lower()
            brand_name = (entry.get("brand") or "").lower()
            for key in (drug_name, brand_name):
                if key:
                    pgx_lookup[key] = {**entry, "classification": cls}

    intersections: list[dict] = []
    for med in medications:
        name = med["name"]
        # Try exact match then partial match
        pgx_entry = pgx_lookup.get(name)
        if not pgx_entry:
            pgx_entry = next(
                (v for k, v in pgx_lookup.items() if k and (k in name or name in k)),
                None,
            )

        intersections.append({
            "name": med["name"],
            "status": med["status"],
            "indication": med["indication"],
            "pgx_matched": pgx_entry is not None,
            "classification": pgx_entry["classification"] if pgx_entry else "not_covered",
            "gene": pgx_entry.get("gene") if pgx_entry else None,
            "phenotype": pgx_entry.get("phenotype") if pgx_entry else None,
            "diplotype": pgx_entry.get("diplotype") if pgx_entry else None,
        })

    return intersections


# ---------------------------------------------------------------------------
# ClinPGx targeted expansion
# ---------------------------------------------------------------------------

def _clinpgx_targets(intersections: list[dict]) -> tuple[set[str], set[str]]:
    """Return (genes, drugs) that need ClinPGx evidence expansion."""
    genes: set[str] = set()
    drugs: set[str] = set()
    for ix in intersections:
        if ix["classification"] in ("avoid", "caution") and ix["pgx_matched"]:
            if ix["gene"]:
                genes.add(ix["gene"])
            drugs.add(ix["name"])
        # planned medications with any PGx signal
        if ix["status"] == "planned" and ix["pgx_matched"]:
            if ix["gene"]:
                genes.add(ix["gene"])
            drugs.add(ix["name"])
    return genes, drugs


def run_clinpgx(genes: set[str], drugs: set[str], output_dir: Path) -> Path | None:
    """Run clinpgx for targeted gene-drug pairs. Returns path to result.json or None."""
    if not genes and not drugs:
        return None
    clinpgx_out = output_dir / "clinpgx_output"
    clinpgx_out.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, str(CLINPGX_SCRIPT), "--output", str(clinpgx_out)]
    if genes:
        cmd += ["--genes", ",".join(sorted(genes))]
    if drugs:
        cmd += ["--drugs", ",".join(sorted(drugs))]
    print(f"  Running ClinPGx for genes={sorted(genes)} drugs={sorted(drugs)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  Warning: ClinPGx returned non-zero exit ({result.returncode}). Continuing without it.")
        return None
    result_json = clinpgx_out / "result.json"
    return result_json if result_json.exists() else None


def load_clinpgx_result(path: Path) -> dict:
    data = json.loads(path.read_text())
    return data.get("data", data)


# ---------------------------------------------------------------------------
# Insight building
# ---------------------------------------------------------------------------

PRIORITY_ORDER = {"avoid": 0, "caution": 1, "indeterminate": 2, "standard": 3, "not_covered": 4}
PRIORITY_LABEL = {
    "avoid": "🔴 High priority review",
    "caution": "🟡 Monitor closely",
    "standard": "🟢 Standard PGx profile",
    "indeterminate": "🟡 Insufficient PGx evidence",
    "not_covered": "⚪ No PGx data available",
}


def build_insights(
    intersections: list[dict],
    clinpgx_data: dict | None,
    allergy_flags: list[str],
) -> list[dict]:
    """Enrich each intersection with ClinPGx evidence and build the final insight list."""
    # Build clinpgx lookup by drug name
    cpgx_by_drug: dict[str, list] = {}
    if clinpgx_data:
        for dr in clinpgx_data.get("drug_results", []):
            dname = (dr.get("drug_name") or dr.get("name") or "").lower()
            if dname:
                cpgx_by_drug.setdefault(dname, []).append(dr)
        for gr in clinpgx_data.get("gene_results", []):
            for ann in gr.get("clinical_annotations", []):
                dname = (ann.get("drug") or "").lower()
                if dname:
                    cpgx_by_drug.setdefault(dname, []).append(ann)

    insights = []
    for ix in sorted(intersections, key=lambda x: PRIORITY_ORDER.get(x["classification"], 99)):
        cpgx_entries = cpgx_by_drug.get(ix["name"], [])
        cpgx_summary = _summarize_clinpgx(cpgx_entries) if cpgx_entries else None
        insights.append({
            **ix,
            "priority_label": PRIORITY_LABEL.get(ix["classification"], "⚪ Unknown"),
            "allergy_flag": ix["name"] in allergy_flags,
            "clinpgx_evidence": cpgx_summary,
            "clinician_note": _clinician_note(ix, cpgx_summary),
        })
    return insights


def _summarize_clinpgx(entries: list[dict]) -> str:
    parts = []
    for e in entries[:3]:  # cap at 3 entries to keep report concise
        level = e.get("evidence_level") or e.get("level_of_evidence") or ""
        guideline = e.get("guideline") or e.get("summary") or ""
        if level or guideline:
            parts.append(f"Evidence {level}: {guideline}".strip(": "))
    return "; ".join(parts) if parts else None


def _clinician_note(ix: dict, cpgx_summary: str | None) -> str:
    cls = ix["classification"]
    name = ix["name"]
    gene = ix["gene"] or "unknown gene"
    pheno = ix["phenotype"] or "unknown phenotype"
    if cls == "avoid":
        note = f"Current guidelines suggest avoiding {name} in {pheno} patients ({gene})."
    elif cls == "caution":
        note = f"Dose adjustment or enhanced monitoring may be warranted for {name} ({gene}, {pheno})."
    elif cls == "indeterminate":
        note = f"Insufficient PGx evidence to classify {name} for this patient's {gene} phenotype."
    elif cls == "standard":
        note = f"No actionable PGx concern for {name} based on {gene} phenotype ({pheno})."
    else:
        note = f"{name} is not covered in the current PGx database."
    if cpgx_summary:
        note += f" ClinPGx: {cpgx_summary}"
    if ix.get("allergy_flag"):
        note = f"⚠️ ALLERGY FLAG — {note}"
    return note


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def generate_report(
    patient_info: dict,
    medications: list[dict],
    insights: list[dict],
    allergy_flags: list[str],
    fhir_config: dict,
    output_dir: Path,
    context_only: bool = False,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    patient_name = _patient_name(patient_info)
    dob = patient_info.get("birthDate", "Unknown")
    gender = patient_info.get("gender", "Unknown")

    lines = [
        "# 🏥 FHIR EMR Pharmacogenomic Medication Review",
        "",
        f"**Generated**: {timestamp}  ",
        f"**FHIR Server**: {fhir_config['server_url']}  ",
        f"**Patient ID**: {fhir_config['patient_id']}  ",
        "",
        "---",
        "",
        "## Patient Context",
        "",
        f"| Field | Value |",
        f"|---|---|",
        f"| Name | {patient_name} |",
        f"| Date of Birth | {dob} |",
        f"| Gender | {gender} |",
        "",
    ]

    if context_only:
        lines += [
            "## PGx-Relevance Map (Context-Only Mode)",
            "",
            "> No genotype file was provided. The table below shows which current medications",
            "> are PGx-sensitive and which genes would be informative.",
            "",
            "| Medication | Status | PGx-Sensitive Genes to Check |",
            "|---|---|---|",
        ]
        pgx_sensitive_drugs = {
            "warfarin": "CYP2C9, VKORC1", "clopidogrel": "CYP2C19",
            "codeine": "CYP2D6", "tramadol": "CYP2D6",
            "simvastatin": "SLCO1B1", "atorvastatin": "SLCO1B1",
            "fluorouracil": "DPYD", "capecitabine": "DPYD",
            "mercaptopurine": "TPMT, NUDT15", "azathioprine": "TPMT, NUDT15",
            "irinotecan": "UGT1A1", "tacrolimus": "CYP3A5",
            "omeprazole": "CYP2C19", "esomeprazole": "CYP2C19",
            "sertraline": "CYP2D6", "amitriptyline": "CYP2D6, CYP2C19",
            "efavirenz": "CYP2B6",
        }
        for med in medications:
            if med["status"] in ("current", "planned"):
                genes = next(
                    (g for d, g in pgx_sensitive_drugs.items() if d in med["name"]),
                    "—",
                )
                lines.append(f"| {med['name']} | {med['status']} | {genes} |")
        lines += [
            "",
            "> Provide a 23andMe/AncestryDNA genotype file to complete full PGx analysis.",
            "",
        ]
    else:
        # Full analysis sections
        avoid = [i for i in insights if i["classification"] == "avoid" and i["status"] != "past"]
        caution = [i for i in insights if i["classification"] == "caution" and i["status"] != "past"]
        standard = [i for i in insights if i["classification"] == "standard"]
        not_covered = [i for i in insights if i["classification"] == "not_covered" and i["status"] in ("current", "planned")]
        planned_pgx = [i for i in insights if i["status"] == "planned" and i["pgx_matched"]]

        # High priority
        if avoid:
            lines += ["## 🔴 High Priority Review", ""]
            lines += _insight_table(avoid)
            lines.append("")

        # Monitor closely
        if caution:
            lines += ["## 🟡 Monitor Closely", ""]
            lines += _insight_table(caution)
            lines.append("")

        # Planned medications pre-emptive check
        if planned_pgx:
            lines += ["## 🔵 Planned Medications — Pre-emptive PGx Check", ""]
            lines += _insight_table(planned_pgx)
            lines.append("")

        # Standard
        if standard:
            lines += ["## 🟢 Standard PGx Profile", ""]
            lines.append("| Medication | Gene | Phenotype | Note |")
            lines.append("|---|---|---|---|")
            for i in standard:
                lines.append(
                    f"| {i['name']} | {i['gene'] or '—'} | {i['phenotype'] or '—'} | No actionable PGx concern |"
                )
            lines.append("")

        # No PGx data
        if not_covered:
            lines += ["## ⚪ Medications Without PGx Coverage", ""]
            lines.append("| Medication | Status |")
            lines.append("|---|---|")
            for i in not_covered:
                lines.append(f"| {i['name']} | {i['status']} |")
            lines.append("")

        # Allergy flags
        if allergy_flags:
            lines += [
                "## ⚠️ Allergy Conflicts",
                "",
                "The following medications appear in both the prescription list and the "
                "patient's AllergyIntolerance record:",
                "",
            ]
            for name in allergy_flags:
                lines.append(f"- **{name}** — review allergy record before any PGx guidance")
            lines.append("")

        # Clinician review note
        lines += [
            "## Notes for Clinician",
            "",
            "- Medications marked 🔴 have a pharmacogenomic signal suggesting potential harm or reduced efficacy.",
            "- Medications marked 🟡 may benefit from dose review or closer monitoring.",
            "- Candidate medication alternatives, if any, should be confirmed against the full clinical picture "
            "including allergies, renal/hepatic function, and co-morbidities before any prescribing decision.",
            "",
        ]

    # Limitations
    lines += [
        "## Limitations",
        "",
        "- PGx analysis covers only the genes and SNPs included in the pharmgx-reporter database (12 genes, ~31 SNPs).",
        "- Medication names are normalized from FHIR coded entries; brand-to-generic mapping may be incomplete.",
        "- ClinPGx evidence is retrieved only for flagged gene-drug pairs; other medications may have unpublished PGx associations.",
        "- FHIR data completeness depends on what the source EMR has recorded.",
        "",
        "---",
        "",
        DISCLAIMER,
    ]

    report_path = output_dir / "report.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")

    # result.json
    result = {
        "skill": "marcelo-emr-fhir-drug-analyser",
        "version": "0.1.0",
        "timestamp": datetime.now().isoformat(),
        "patient": {
            "name": patient_name,
            "dob": dob,
            "gender": gender,
            "fhir_id": patient_info.get("id"),
        },
        "medications": medications,
        "pgx_intersections": insights if not context_only else [],
        "allergy_flags": allergy_flags,
        "context_only": context_only,
        "metadata": {
            "server_url": fhir_config["server_url"],
            "patient_id": fhir_config["patient_id"],
        },
    }
    (output_dir / "result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    return report_path


def _insight_table(insights: list[dict]) -> list[str]:
    rows = [
        "| Medication | Status | Gene | Phenotype | Classification | Clinician Note |",
        "|---|---|---|---|---|---|",
    ]
    for i in insights:
        allergy = " ⚠️ ALLERGY" if i.get("allergy_flag") else ""
        rows.append(
            f"| {i['name']}{allergy} | {i['status']} | {i['gene'] or '—'} "
            f"| {i['phenotype'] or '—'} | {i['classification'].upper()} "
            f"| {i['clinician_note'][:120]}... |"
        )
    return rows


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def analyse(
    fhir_config: dict,
    genotype_path: Path | None,
    pharmgx_result_path: Path | None,
    context_only: bool,
    output_dir: Path,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Fetch FHIR resources
    print("\n[1/6] Fetching FHIR resources...")
    resources = fetch_fhir_resources(
        fhir_config["server_url"],
        fhir_config["patient_id"],
        fhir_config.get("auth"),
    )

    # 2. Normalize medications
    print("\n[2/6] Normalizing medications...")
    medications = normalize_medications(
        resources["medication_requests"],
        resources["medication_statements"],
        resources["conditions"],
    )
    allergy_flags = check_allergies(medications, resources["allergies"])
    current = [m for m in medications if m["status"] == "current"]
    planned = [m for m in medications if m["status"] == "planned"]
    print(f"  Current: {len(current)}  Planned: {len(planned)}  Allergy flags: {len(allergy_flags)}")

    # 3. Context-only mode
    if context_only or (not genotype_path and not pharmgx_result_path):
        if not context_only:
            print("\n  No genotype or PharmGx result provided — running in context-only mode.")
        print("\n[3/6] Skipping PGx analysis (context-only mode).")
        return generate_report(
            resources["patient"], medications, [], allergy_flags,
            fhir_config, output_dir, context_only=True,
        )

    # 4. Run or load PharmGx
    print("\n[3/6] Running PharmGx...")
    if pharmgx_result_path:
        print(f"  Loading pre-computed PharmGx result from {pharmgx_result_path}")
        pgx_data = load_pharmgx_result(pharmgx_result_path)
    else:
        result_json = run_pharmgx(genotype_path, output_dir)
        pgx_data = load_pharmgx_result(result_json)

    # 5. Intersect with FHIR medications
    print("\n[4/6] Intersecting PharmGx results with FHIR medication list...")
    intersections = intersect_pgx_with_meds(medications, pgx_data)
    flagged = [i for i in intersections if i["classification"] in ("avoid", "caution")]
    print(f"  Flagged (avoid/caution): {len(flagged)}")

    # 6. Targeted ClinPGx expansion
    print("\n[5/6] Running targeted ClinPGx expansion...")
    genes, drugs = _clinpgx_targets(intersections)
    clinpgx_result_path = run_clinpgx(genes, drugs, output_dir)
    clinpgx_data = load_clinpgx_result(clinpgx_result_path) if clinpgx_result_path else None

    # 7. Build insights and report
    print("\n[6/6] Building insights and generating report...")
    insights = build_insights(intersections, clinpgx_data, allergy_flags)
    return generate_report(
        resources["patient"], medications, insights, allergy_flags,
        fhir_config, output_dir,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="FHIR EMR Pharmacogenomic Drug Analyser — ClawBio skill"
    )
    parser.add_argument(
        "--fhir-connection", type=Path,
        help="JSON file with fhir.server_url and fhir.patient_id (e.g. demo_input/patient_fhir.json)",
    )
    parser.add_argument(
        "--genotype", type=Path,
        help="23andMe/AncestryDNA raw genotype file (.txt or .txt.gz)",
    )
    parser.add_argument(
        "--pharmgx-result", type=Path,
        help="Pre-computed PharmGx result.json (skips re-running pharmgx-reporter)",
    )
    parser.add_argument(
        "--context-only", action="store_true",
        help="Stop after FHIR retrieval; return PGx-relevance map without running PharmGx",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("/tmp/fhir_pgx_report"),
        help="Output directory (default: /tmp/fhir_pgx_report)",
    )
    parser.add_argument(
        "--demo", action="store_true",
        help="Run with bundled demo data (demo_input/patient_fhir.json + demo_input/genotype_data.txt)",
    )
    args = parser.parse_args()

    if args.demo:
        fhir_path = SKILL_DIR / "demo_input" / "patient_fhir.json"
        genotype_path = SKILL_DIR / "demo_input" / "genotype_data.txt"
        pharmgx_result_path = None
    elif args.fhir_connection:
        fhir_path = args.fhir_connection
        genotype_path = args.genotype
        pharmgx_result_path = args.pharmgx_result
    else:
        parser.error("Provide --fhir-connection <file> or use --demo")

    fhir_config = load_fhir_config(fhir_path)
    report = analyse(
        fhir_config=fhir_config,
        genotype_path=genotype_path,
        pharmgx_result_path=pharmgx_result_path,
        context_only=args.context_only,
        output_dir=args.output,
    )
    print(f"\nReport written to: {report}")


if __name__ == "__main__":
    main()
