"""Output writers for clinical-trial-finder.

Generates: report.md, summary.json, fhir_bundle.json, phase chart,
commands.sh, and checksums.sha256.  Each writer takes trial data and
an output directory, creates the file, and returns its Path.
"""

import argparse
import csv
import hashlib
import html as _html
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from constants import (
    CSV_COLUMNS,
    DEFAULT_PAGE_SIZE,
    DISCLAIMER,
    FHIR_PHASE,
    FHIR_PHASE_DISPLAY,
    FHIR_STATUS,
    FHIR_VALID_PHASES,
    FHIR_VALID_STATUSES,
    PHASE_DISPLAY,
    PHASE_ORDER,
    STATUS_COLOR,
    STATUS_LABEL,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def count_recruiting(trials: list[dict]) -> int:
    """Count trials with RECRUITING status."""
    return sum(1 for t in trials if t["status"] == "RECRUITING")


def _sha256(path: Path) -> str:
    """Return hex SHA-256 digest of a file.  Reads in 8KB chunks to handle
    large files without loading everything into memory."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------


def write_report(
    query_info: dict,
    trials: list[dict],
    output_dir: Path,
    gene_context: dict | None = None,
) -> Path:
    """Write a human-readable markdown report.

    Includes a summary table, per-trial details with NCT links, and the
    ClawBio safety disclaimer.  When --gene was used, the report header
    shows the gene context (symbol, diseases, OpenTargets score).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "report.md"

    recruiting = count_recruiting(trials)
    phases = sorted({t["phase"] for t in trials if t["phase"]})
    study_types = sorted({t["study_type"] for t in trials if t["study_type"]})
    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")

    lines = ["# Clinical Trial Finder Report", ""]

    # Header varies depending on whether this was a gene or text query
    if gene_context:
        lines += [
            f"**Gene**: `{gene_context['symbol']}` -- {gene_context['name']}  ",
            f"**Via**: OpenTargets Platform (association score >= {gene_context['min_score']})  ",
            f"**Associated diseases queried**: {', '.join(gene_context['diseases'])}  ",
        ]
    else:
        lines += [f"**Query**: `{query_info['query']}`  "]

    lines += [
        f"**Generated**: {timestamp}  ",
        "**Source**: ClinicalTrials.gov API v2  ",
        f"**Trials found**: {len(trials)} ({recruiting} recruiting)",
        "",
        "---",
        "",
        "## Summary",
        "",
        "| Metric | Value |",
        "|--------|-------|",
        f"| Total trials | {len(trials)} |",
        f"| Recruiting now | {recruiting} |",
        f"| Phases | {', '.join(phases) or 'N/A'} |",
        f"| Study types | {', '.join(study_types) or 'N/A'} |",
        "",
        "## Trials",
        "",
    ]

    # Per-trial detail sections with text status labels
    for t in trials:
        label = STATUS_LABEL.get(t["status"], "[UNKNOWN]")
        nct_url = f"https://clinicaltrials.gov/study/{t['nct_id']}"
        lines += [
            f"### {label} {t['title']}",
            "",
            f"- **NCT ID**: [{t['nct_id']}]({nct_url})",
            f"- **Status**: {t['status']}",
            f"- **Phase**: {t['phase'] or 'N/A'}",
            f"- **Type**: {t['study_type'] or 'N/A'}",
            f"- **Start**: {t['start_date'] or 'N/A'} | **Est. completion**: {t['completion_date'] or 'N/A'}",
        ]
        if t["conditions"]:
            lines.append(f"- **Conditions**: {', '.join(t['conditions'][:3])}")
        if t["interventions"]:
            lines.append(f"- **Interventions**: {', '.join(t['interventions'][:3])}")
        if t["summary"]:
            lines += ["", f"> {t['summary']}"]
        lines.append("")

    lines += [
        "---",
        "",
        "## Reproducibility",
        "",
        "- `commands.sh` -- exact command to reproduce",
        "- `checksums.sha256` -- SHA-256 of all outputs",
        "",
        "---",
        "",
        DISCLAIMER,
    ]
    report_path.write_text("\n".join(lines))
    return report_path


# ---------------------------------------------------------------------------
# JSON summary (consumed by profile-report and other skills)
# ---------------------------------------------------------------------------


def write_summary(query_info: dict, trials: list[dict], output_dir: Path) -> Path:
    """Write machine-readable JSON with all trial data.

    This is the primary interface for downstream skill chaining --
    profile-report reads summary.json to embed trials in the unified report.
    """
    payload = {
        "query": query_info["query"],
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "source": "clinicaltrials.gov/api/v2",
        "total": len(trials),
        "recruiting": count_recruiting(trials),
        "trials": trials,
    }
    path = output_dir / "summary.json"
    path.write_text(json.dumps(payload, indent=2))
    return path


# ---------------------------------------------------------------------------
# FHIR R4 Bundle
#
# Produces a valid FHIR R4 Bundle of ResearchStudy resources.
# Status and phase codes come from the published R4 value sets.
# Conditions use MeSH coding from CT.gov's derivedSection when available,
# falling back to free-text CodeableConcept otherwise.
# ---------------------------------------------------------------------------


def write_fhir_bundle(trials: list[dict], output_dir: Path) -> Path:
    """Write a FHIR R4 Bundle of ResearchStudy resources."""
    entries = [_trial_to_fhir(t) for t in trials]
    bundle = {
        "resourceType": "Bundle",
        "type": "searchset",
        # FHIR instant requires timezone -- we use UTC explicitly
        "timestamp": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S+00:00"),
        "total": len(entries),
        "entry": entries,
        "meta": {
            "tag": [{"system": "https://clawbio.ai", "code": "clinical-trial-finder"}]
        },
    }
    path = output_dir / "fhir_bundle.json"
    path.write_text(json.dumps(bundle, indent=2))
    return path


def _trial_to_fhir(t: dict) -> dict:
    """Map a normalised trial record to a FHIR R4 ResearchStudy Bundle entry.

    Each field maps to the corresponding FHIR R4 element:
      nct_id       -> identifier (official), id
      status       -> status (via FHIR_STATUS value set)
      phase        -> phase (CodeableConcept with system + code + display)
      study_type   -> category
      conditions   -> condition (MeSH-coded when available)
      interventions -> focus
      summary      -> description
      dates        -> period (start/end)
    """
    resource: dict = {
        "resourceType": "ResearchStudy",
        "id": t["nct_id"],
        "meta": {"profile": ["http://hl7.org/fhir/StructureDefinition/ResearchStudy"]},
        "identifier": [
            {
                "use": "official",
                "system": "https://clinicaltrials.gov",
                "value": t["nct_id"],
            }
        ],
        "title": t["title"],
        "status": FHIR_STATUS.get(t["status"], "in-review"),
        "phase": {
            "coding": [
                {
                    k: v
                    for k, v in {
                        "system": "http://terminology.hl7.org/CodeSystem/research-study-phase",
                        "code": FHIR_PHASE.get(t["phase"], "n-a"),
                        "display": FHIR_PHASE_DISPLAY.get(t["phase"], "N/A"),
                    }.items()
                    if v  # FHIR R4 forbids empty strings
                }
            ]
        },
    }

    # INTERVENTIONAL vs OBSERVATIONAL -- maps to ResearchStudy.category
    if t["study_type"]:
        resource["category"] = [
            {
                "coding": [
                    {
                        "system": "http://clinicaltrials.gov/study-type",
                        "code": t["study_type"],
                        "display": t["study_type"].replace("_", " ").title(),
                    }
                ]
            }
        ]

    # Prefer MeSH-coded conditions from derivedSection (authoritative).
    # Fall back to free-text if CT.gov hasn't indexed this trial yet.
    if t.get("condition_meshes"):
        resource["condition"] = [
            {
                "coding": [
                    {
                        "system": "http://id.nlm.nih.gov/mesh/",
                        "code": m["id"],
                        "display": m.get("term", ""),
                    }
                ],
                "text": m.get("term", m["id"]),
            }
            for m in t["condition_meshes"]
            if m.get("id")
        ]
    elif t["conditions"]:
        resource["condition"] = [{"text": c} for c in t["conditions"]]

    # Interventions (drug names, devices) -> ResearchStudy.focus
    if t["interventions"]:
        resource["focus"] = [{"text": i} for i in t["interventions"]]

    if t["summary"]:
        resource["description"] = t["summary"]

    # Period with start and optional end date
    if t["start_date"]:
        resource["period"] = {"start": t["start_date"]}
        if t["completion_date"]:
            resource["period"]["end"] = t["completion_date"]

    return {
        "fullUrl": f"https://clinicaltrials.gov/study/{t['nct_id']}",
        "resource": resource,
    }


# ---------------------------------------------------------------------------
# Phase distribution chart
# ---------------------------------------------------------------------------


def write_phase_chart(
    trials: list[dict], output_dir: Path, title: str = ""
) -> Path | None:
    """Write a stacked bar chart of trial counts by phase, coloured by status.

    Returns the Path to the PNG, or None if matplotlib is not installed.
    matplotlib is optional -- the skill works without it.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    # Aggregate: phase -> {status: count}
    phase_status: dict[str, Counter] = defaultdict(Counter)
    for t in trials:
        phase = PHASE_DISPLAY.get(t["phase"], t["phase"] or "Unknown")
        phase_status[phase][t["status"]] += 1

    # Use fixed phase order, then append any unexpected phases at the end
    phases = [p for p in PHASE_ORDER if p in phase_status]
    phases += [p for p in phase_status if p not in phases]
    if not phases:
        return None

    # Collect which statuses actually appear, then filter in defined order
    present_statuses = {s for p in phases for s in phase_status[p]}
    all_statuses = [s for s in STATUS_COLOR if s in present_statuses]

    # Width scales with phase count (1.4in per bar, min 7in for readability)
    fig, ax = plt.subplots(figsize=(max(7, len(phases) * 1.4), 5))
    # bottoms tracks cumulative height per bar for stacking
    bottoms = [0] * len(phases)

    # Each status becomes a layer in the stacked bar
    for status in all_statuses:
        counts = [phase_status[p].get(status, 0) for p in phases]
        ax.bar(
            phases,
            counts,
            bottom=bottoms,
            color=STATUS_COLOR[status],
            label=status.replace("_", " ").title(),
            edgecolor="white",
            linewidth=0.5,
        )
        bottoms = [b + c for b, c in zip(bottoms, counts)]

    # Total labels on top -- 2% proportional offset avoids overlap at any scale
    for i, total in enumerate(bottoms):
        if total:
            ax.text(
                i,
                total * 1.02,
                str(total),
                ha="center",
                va="bottom",
                fontsize=9,
                fontweight="bold",
            )

    ax.set_xlabel("Phase", fontsize=11)
    ax.set_ylabel("Number of Trials", fontsize=11)
    ax.set_title(
        title or "Clinical Trial Phase Distribution",
        fontsize=13,
        fontweight="bold",
        pad=14,
    )
    ax.set_ylim(0, max(bottoms, default=1) * 1.25 + 1)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.85)
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    path = fig_dir / "phase_distribution.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------


def write_csv(trials: list[dict], output_dir: Path) -> Path:
    """Write trials as a CSV table for import into Excel, R, or pandas.

    Columns match CSV_COLUMNS.  List fields (conditions, interventions) are
    pipe-delimited so they survive CSV parsing without quoting issues.
    """
    table_dir = output_dir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    path = table_dir / "trials.csv"

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for t in trials:
            row = dict(t)
            row["conditions"] = " | ".join(t.get("conditions", []))
            row["interventions"] = " | ".join(t.get("interventions", []))
            writer.writerow(row)
    return path


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------


def write_html(
    query_info: dict,
    trials: list[dict],
    output_dir: Path,
    gene_context: dict | None = None,
) -> Path:
    """Write a self-contained HTML report with styled trial cards.

    No external CSS/JS dependencies -- everything is inline so the file
    opens correctly from any file manager or browser.
    """
    recruiting = count_recruiting(trials)
    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")

    # Duplicated from STATUS_COLOR because HTML needs inline styles (no imports)
    _STATUS_CSS = {
        "RECRUITING": "#2ecc71",
        "COMPLETED": "#95a5a6",
        "TERMINATED": "#e74c3c",
        "WITHDRAWN": "#7f8c8d",
        "ACTIVE_NOT_RECRUITING": "#f1c40f",
        "NOT_YET_RECRUITING": "#3498db",
        "ENROLLING_BY_INVITATION": "#27ae60",
        "SUSPENDED": "#e67e22",
        "APPROVED_FOR_MARKETING": "#8e44ad",
        "UNKNOWN": "#bdc3c7",
    }

    _esc = _html.escape
    header = (
        f"<strong>Gene:</strong> {_esc(gene_context.get('symbol', '?'))} -- {_esc(gene_context.get('name', ''))}"
        if gene_context
        else f"<strong>Query:</strong> {_esc(query_info['query'])}"
    )

    cards = []
    for t in trials:
        color = _STATUS_CSS.get(t["status"], "#bdc3c7")
        nct_url = f"https://clinicaltrials.gov/study/{_esc(t['nct_id'])}"
        source_badge = (
            f' <span style="background:#e8daef;padding:2px 6px;border-radius:3px;font-size:0.75em">'
            f"{_esc(t.get('source', 'ctgov').upper())}</span>"
            if t.get("source")
            else ""
        )
        phase_label = _esc(PHASE_DISPLAY.get(t["phase"], t["phase"] or "N/A"))
        cards.append(f"""
        <div class="trial-card" data-status="{_esc(t["status"])}" data-phase="{
            phase_label
        }"
             style="border-left:4px solid {color};padding:12px 16px;margin:8px 0;
                    background:#fafafa;border-radius:0 4px 4px 0">
          <div style="font-weight:bold">
            <span style="background:{
            color
        };color:white;padding:2px 8px;border-radius:3px;
                         font-size:0.8em;margin-right:8px">{_esc(t["status"])}</span>
            <a href="{nct_url}" target="_blank">{_esc(t["title"])}</a>{source_badge}
          </div>
          <div style="margin-top:6px;color:#555;font-size:0.9em">
            {_esc(t["nct_id"])} | Phase: {_esc(t["phase"] or "N/A")} |
            {_esc(t["start_date"] or "?")} -- {_esc(t["completion_date"] or "?")}
          </div>
          {
            (
                '<div style="margin-top:4px;font-size:0.85em;color:#777">'
                + ", ".join(_esc(c) for c in t["conditions"][:3])
                + "</div>"
            )
            if t.get("conditions")
            else ""
        }
        </div>""")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Clinical Trial Finder Report</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
           max-width: 900px; margin: 40px auto; padding: 0 20px; color: #333; }}
    h1 {{ border-bottom: 2px solid #2ecc71; padding-bottom: 8px; }}
    table {{ border-collapse: collapse; width: 100%; margin: 16px 0; }}
    th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
    th {{ background: #f5f5f5; }}
    .disclaimer {{ background: #fff3cd; padding: 12px; border-radius: 4px;
                   margin-top: 24px; font-size: 0.9em; }}
    .filters {{ background: #f0f4f8; padding: 12px 16px; border-radius: 6px;
               margin: 16px 0; display: flex; gap: 12px; flex-wrap: wrap; align-items: center; }}
    .filters input, .filters select {{ padding: 6px 10px; border: 1px solid #ccc;
               border-radius: 4px; font-size: 0.9em; }}
    .filters input {{ flex: 1; min-width: 200px; }}
    .filter-count {{ font-size: 0.85em; color: #666; margin-left: auto; }}
    .trial-card {{ transition: opacity 0.15s; }}
    .trial-card.hidden {{ display: none; }}
  </style>
</head>
<body>
  <h1>Clinical Trial Finder Report</h1>
  <p>{header}<br>
     <strong>Generated:</strong> {timestamp}<br>
     <strong>Trials found:</strong> {len(trials)} ({recruiting} recruiting)</p>

  <table>
    <tr><th>Metric</th><th>Value</th></tr>
    <tr><td>Total trials</td><td>{len(trials)}</td></tr>
    <tr><td>Recruiting now</td><td>{recruiting}</td></tr>
  </table>

  <h2>Trials</h2>

  <div class="filters">
    <input type="text" id="searchBox" placeholder="Search trials..." oninput="filterTrials()">
    <select id="statusFilter" onchange="filterTrials()">
      <option value="">All statuses</option>
      {"".join(f'<option value="{_esc(s)}">{_esc(s)}</option>' for s in sorted({t["status"] for t in trials}))}
    </select>
    <select id="phaseFilter" onchange="filterTrials()">
      <option value="">All phases</option>
      {"".join(f'<option value="{_esc(p)}">{_esc(p)}</option>' for p in sorted({PHASE_DISPLAY.get(t["phase"], t["phase"] or "N/A") for t in trials}))}
    </select>
    <span class="filter-count" id="filterCount">{len(trials)} trials</span>
  </div>

  {"".join(cards)}

  <script>
  function filterTrials() {{
    const search = document.getElementById('searchBox').value.toLowerCase();
    const status = document.getElementById('statusFilter').value;
    const phase = document.getElementById('phaseFilter').value;
    let visible = 0;
    document.querySelectorAll('.trial-card').forEach(card => {{
      const text = card.textContent.toLowerCase();
      const matchSearch = !search || text.includes(search);
      const matchStatus = !status || card.dataset.status === status;
      const matchPhase = !phase || card.dataset.phase === phase;
      const show = matchSearch && matchStatus && matchPhase;
      card.classList.toggle('hidden', !show);
      if (show) visible++;
    }});
    document.getElementById('filterCount').textContent = visible + ' of {len(trials)} trials';
  }}
  </script>

  <div class="disclaimer">{DISCLAIMER.replace("*", "")}</div>
  <p style="color:#999;font-size:0.8em;margin-top:24px">
    Generated by ClawBio clinical-trial-finder v0.1.0</p>
</body>
</html>"""

    path = output_dir / "report.html"
    path.write_text(html)
    return path


# ---------------------------------------------------------------------------
# FHIR R4 validation
# ---------------------------------------------------------------------------


def validate_fhir_bundle(bundle: dict) -> list[str]:
    """Validate a FHIR R4 Bundle against basic structural rules.

    Returns a list of human-readable error strings.  Empty list = valid.
    This is NOT a full FHIR validator (use HAPI for that); it checks:
      - Bundle structure (resourceType, type, total, entries)
      - ResearchStudy required fields (id, status, title)
      - Status values against the R4 value set
      - Phase codes against the R4 CodeSystem
    """
    errors: list[str] = []

    if bundle.get("resourceType") != "Bundle":
        errors.append("Bundle.resourceType must be 'Bundle'")
    if bundle.get("type") != "searchset":
        errors.append("Bundle.type must be 'searchset'")
    if not bundle.get("timestamp"):
        errors.append("Bundle.timestamp is required")

    entries = bundle.get("entry", [])
    if bundle.get("total") != len(entries):
        errors.append(
            f"Bundle.total ({bundle.get('total')}) != entry count ({len(entries)})"
        )

    for i, entry in enumerate(entries):
        resource = entry.get("resource", {})
        prefix = f"entry[{i}]"

        if resource.get("resourceType") != "ResearchStudy":
            errors.append(f"{prefix}: resourceType must be 'ResearchStudy'")
        if not resource.get("id"):
            errors.append(f"{prefix}: id is required")
        if not resource.get("title"):
            errors.append(f"{prefix}: title is required")

        status = resource.get("status", "")
        if status not in FHIR_VALID_STATUSES:
            errors.append(f"{prefix}: invalid status '{status}'")

        phase = resource.get("phase", {})
        for coding in phase.get("coding", []):
            code = coding.get("code", "")
            if code and code not in FHIR_VALID_PHASES:
                errors.append(f"{prefix}: invalid phase code '{code}'")

        if not entry.get("fullUrl"):
            errors.append(f"{prefix}: fullUrl is required")

    return errors


# ---------------------------------------------------------------------------
# Reproducibility -- commands.sh and checksums.sha256
# ---------------------------------------------------------------------------


def write_commands(args: argparse.Namespace, output_dir: Path) -> Path:
    """Write commands.sh with the exact CLI invocation to reproduce the run.

    Reconstructs the command from parsed args (not sys.argv) so the output
    is always canonical, even if the user passed args in a different order.
    Only includes non-default values to keep the command minimal.
    """
    parts = ["python skills/clinical-trial-finder/clinical_trial_finder.py"]

    def _sq(val: str) -> str:
        """Shell-quote a value using single quotes (handles spaces and special chars)."""
        return "'" + str(val).replace("'", "'\\''") + "'"

    # Source mode -- mutually exclusive
    if args.demo:
        parts.append("--demo")
    elif args.input:
        parts.append(f"--input {_sq(args.input)}")
    elif args.query:
        parts.append(f"--query {_sq(args.query)}")
    elif args.gene:
        parts.append(f"--gene {_sq(args.gene)}")
    elif getattr(args, "rsid", None):
        parts.append(f"--rsid {_sq(args.rsid)}")

    # Optional filters and flags -- only if non-default
    if args.status:
        parts.append(f"--status {args.status}")
    if getattr(args, "country", None):
        parts.append(f"--country {_sq(args.country)}")
    if args.max_results != DEFAULT_PAGE_SIZE:
        parts.append(f"--max-results {args.max_results}")
    if args.fhir:
        parts.append("--fhir")
    if getattr(args, "euctr", False):
        parts.append("--euctr")
    if getattr(args, "ot_min_score", 0.6) != 0.6:
        parts.append(f"--ot-min-score {args.ot_min_score}")
    if getattr(args, "ot_max_diseases", 5) != 5:
        parts.append(f"--ot-max-diseases {args.ot_max_diseases}")

    parts.append(f"--output {_sq(output_dir)}")
    path = output_dir / "commands.sh"
    path.write_text(" \\\n  ".join(parts) + "\n")
    return path


def write_checksums(output_dir: Path) -> Path:
    """Write checksums.sha256 with SHA-256 digests of all generated outputs.

    Uses standard sha256sum format (hex digest, two spaces, relative path)
    so users can verify with: sha256sum -c checksums.sha256
    Skips files that don't exist (e.g. fhir_bundle.json when --fhir is off).
    """
    targets = [
        "report.md",
        "summary.json",
        "fhir_bundle.json",
        "report.html",
        "tables/trials.csv",
        "figures/phase_distribution.png",
        "commands.sh",
    ]
    lines = []
    for rel in targets:
        p = output_dir / rel
        try:
            lines.append(f"{_sha256(p)}  {rel}")
        except FileNotFoundError:
            pass  # Skip outputs that weren't generated (e.g. fhir_bundle.json without --fhir)
    path = output_dir / "checksums.sha256"
    path.write_text("\n".join(lines) + "\n")
    return path
