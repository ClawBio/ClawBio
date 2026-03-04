"""
ClawBio PharmGx — Rich-text HTML report generator.

Produces a single self-contained HTML file (no external dependencies) that
mirrors the markdown report with a clinical-style layout:
  - Persistent "DEMO — FOR RESEARCH USE ONLY" banner
  - Colour-coded drug-recommendation rows (red/amber/green)
  - Sticky table-of-contents sidebar with clickable section links
  - Collapsible variant-detail section
  - Print-friendly media query
"""

import hashlib
from datetime import datetime, timezone
from html import escape
from pathlib import Path


# ---------------------------------------------------------------------------
# Colour palette for drug status
# ---------------------------------------------------------------------------
_STATUS_COLOURS = {
    "avoid":         ("#b91c1c", "#fef2f2", "AVOID"),
    "caution":       ("#92400e", "#fffbeb", "CAUTION"),
    "standard":      ("#166534", "#f0fdf4", "OK"),
    "indeterminate": ("#6b7280", "#f9fafb", "INSUFFICIENT DATA"),
}

# ---------------------------------------------------------------------------
# CSS — embedded so the HTML is fully self-contained
# ---------------------------------------------------------------------------
_CSS = """\
:root {
  --clawbio-navy: #1e293b;
  --clawbio-teal: #0d9488;
  --banner-bg:    #dc2626;
  --sidebar-w:    220px;
}
*, *::before, *::after { box-sizing: border-box; }
html { scroll-behavior: smooth; }
body {
  margin: 0;
  font-family: "Inter", "Segoe UI", system-ui, -apple-system, sans-serif;
  color: #1e293b;
  background: #f8fafc;
  line-height: 1.6;
}

/* ── Banner ─────────────────────────────────────────────────────────────── */
.banner {
  position: sticky;
  top: 0;
  z-index: 100;
  background: var(--banner-bg);
  color: #fff;
  text-align: center;
  padding: 6px 12px;
  font-weight: 700;
  font-size: 14px;
  letter-spacing: 2px;
  text-transform: uppercase;
}

/* ── Sidebar nav ────────────────────────────────────────────────────────── */
.sidebar {
  position: fixed;
  top: 36px;  /* banner height */
  left: 0;
  width: var(--sidebar-w);
  height: calc(100vh - 36px);
  overflow-y: auto;
  background: var(--clawbio-navy);
  padding: 24px 0;
  z-index: 50;
}
.sidebar a {
  display: block;
  padding: 8px 20px;
  color: #94a3b8;
  text-decoration: none;
  font-size: 13px;
  border-left: 3px solid transparent;
  transition: all .15s;
}
.sidebar a:hover, .sidebar a.active {
  color: #fff;
  background: rgba(255,255,255,.06);
  border-left-color: var(--clawbio-teal);
}

/* ── Main content ───────────────────────────────────────────────────────── */
.main {
  margin-left: var(--sidebar-w);
  max-width: 960px;
  padding: 32px 40px 80px;
}
h1 {
  font-size: 28px;
  color: var(--clawbio-navy);
  margin-bottom: 4px;
}
.subtitle {
  color: #64748b;
  font-size: 14px;
  margin-bottom: 24px;
}

/* ── Meta grid ──────────────────────────────────────────────────────────── */
.meta-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
  gap: 12px;
  margin-bottom: 32px;
}
.meta-card {
  background: #fff;
  border: 1px solid #e2e8f0;
  border-radius: 8px;
  padding: 12px 16px;
}
.meta-label {
  font-size: 11px;
  text-transform: uppercase;
  color: #94a3b8;
  letter-spacing: .5px;
}
.meta-value {
  font-size: 16px;
  font-weight: 600;
  word-break: break-all;
}
.meta-value.mono { font-family: "Fira Code", "Consolas", monospace; font-size: 12px; }

/* ── Summary cards ──────────────────────────────────────────────────────── */
.summary-row {
  display: flex;
  gap: 16px;
  flex-wrap: wrap;
  margin-bottom: 32px;
}
.summary-card {
  flex: 1 1 140px;
  border-radius: 8px;
  padding: 16px 20px;
  text-align: center;
  min-width: 140px;
}
.summary-card .count { font-size: 32px; font-weight: 700; }
.summary-card .label { font-size: 12px; text-transform: uppercase; letter-spacing: .5px; }

.card-avoid   { background: #fef2f2; color: #b91c1c; border: 1px solid #fecaca; }
.card-caution { background: #fffbeb; color: #92400e; border: 1px solid #fde68a; }
.card-ok      { background: #f0fdf4; color: #166534; border: 1px solid #bbf7d0; }
.card-unknown { background: #f9fafb; color: #6b7280; border: 1px solid #e5e7eb; }

/* ── Section ────────────────────────────────────────────────────────────── */
section {
  background: #fff;
  border: 1px solid #e2e8f0;
  border-radius: 10px;
  padding: 24px 28px;
  margin-bottom: 24px;
}
section h2 {
  font-size: 18px;
  margin: 0 0 16px;
  padding-bottom: 8px;
  border-bottom: 2px solid #e2e8f0;
}

/* ── Tables ─────────────────────────────────────────────────────────────── */
table {
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;
}
thead th {
  background: var(--clawbio-navy);
  color: #fff;
  text-align: left;
  padding: 8px 12px;
  font-weight: 600;
  font-size: 12px;
  text-transform: uppercase;
  letter-spacing: .3px;
}
thead th:first-child { border-radius: 6px 0 0 0; }
thead th:last-child  { border-radius: 0 6px 0 0; }
tbody td {
  padding: 8px 12px;
  border-bottom: 1px solid #f1f5f9;
  vertical-align: top;
}
tbody tr:hover { background: #f8fafc; }

/* status badge */
.badge {
  display: inline-block;
  padding: 2px 8px;
  border-radius: 4px;
  font-size: 11px;
  font-weight: 700;
  letter-spacing: .3px;
}
.badge-avoid   { background: #fef2f2; color: #b91c1c; }
.badge-caution { background: #fffbeb; color: #92400e; }
.badge-ok      { background: #f0fdf4; color: #166534; }
.badge-unknown { background: #f9fafb; color: #6b7280; }

/* ── Alerts ─────────────────────────────────────────────────────────────── */
.alert-list { list-style: none; padding: 0; margin: 0; }
.alert-list li {
  padding: 10px 14px;
  border-left: 4px solid;
  margin-bottom: 8px;
  border-radius: 0 6px 6px 0;
  font-size: 14px;
}
.alert-avoid  { border-color: #dc2626; background: #fef2f2; }
.alert-caution{ border-color: #f59e0b; background: #fffbeb; }
.alert-list .drug-name { font-weight: 700; }
.alert-list .gene-tag {
  font-size: 11px;
  background: #e2e8f0;
  padding: 1px 6px;
  border-radius: 3px;
  margin-left: 6px;
}

/* ── Collapsible ────────────────────────────────────────────────────────── */
details summary {
  cursor: pointer;
  font-weight: 600;
  padding: 8px 0;
  user-select: none;
}
details summary:hover { color: var(--clawbio-teal); }

/* ── Disclaimer ─────────────────────────────────────────────────────────── */
.disclaimer {
  background: #fffbeb;
  border: 1px solid #fde68a;
  border-radius: 8px;
  padding: 16px 20px;
  font-size: 13px;
  color: #92400e;
}

/* ── Print ──────────────────────────────────────────────────────────────── */
@media print {
  .banner { position: relative; }
  .sidebar { display: none; }
  .main { margin-left: 0; }
  section { break-inside: avoid; }
}

/* ── Responsive ─────────────────────────────────────────────────────────── */
@media (max-width: 768px) {
  .sidebar { display: none; }
  .main { margin-left: 0; padding: 20px; }
}
"""


# ---------------------------------------------------------------------------
# HTML helpers
# ---------------------------------------------------------------------------

def _e(text):
    """HTML-escape."""
    return escape(str(text))


def _badge(classification):
    colour_key = classification
    label = _STATUS_COLOURS.get(colour_key, ("", "", classification.upper()))[2]
    css = f"badge-{colour_key}" if colour_key in _STATUS_COLOURS else "badge-unknown"
    return f'<span class="badge {css}">{_e(label)}</span>'


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_html_report(input_path, fmt, total_snps, pgx_snps, profiles,
                         drug_results, gene_defs, pgx_snp_defs):
    """Return a complete self-contained HTML string.

    Parameters match ``generate_report()`` in pharmgx_reporter.py, with
    two extras that let us avoid importing module-level constants:
      *gene_defs*   — the GENE_DEFS dict
      *pgx_snp_defs* — the PGX_SNPS dict (for total count)
    """
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    checksum = hashlib.sha256(Path(input_path).read_bytes()).hexdigest()
    fname = Path(input_path).name

    n_std = len(drug_results["standard"])
    n_cau = len(drug_results["caution"])
    n_avo = len(drug_results["avoid"])
    n_ind = len(drug_results.get("indeterminate", []))

    # ── sidebar links ───────────────────────────────────────────────────
    toc = [
        ("summary", "Drug Summary"),
        ("alerts", "Alerts"),
        ("genes", "Gene Profiles"),
        ("variants", "Variants"),
        ("drugs", "All Drugs"),
        ("methods", "Methods"),
        ("disclaimer", "Disclaimer"),
    ]
    sidebar_html = "\n".join(
        f'    <a href="#{sid}">{_e(label)}</a>' for sid, label in toc
    )

    # ── meta cards ──────────────────────────────────────────────────────
    meta_items = [
        ("Date", now, False),
        ("Input File", fname, True),
        ("Format", fmt, False),
        ("PGx SNPs", f"{len(pgx_snps)}/{len(pgx_snp_defs)}", False),
        ("Genes Profiled", str(len(profiles)), False),
        ("Drugs Assessed", str(n_std + n_cau + n_avo + n_ind), False),
        ("SHA-256", checksum[:16] + "...", True),
    ]
    meta_html = "\n".join(
        f'    <div class="meta-card"><div class="meta-label">{_e(lbl)}</div>'
        f'<div class="meta-value{"  mono" if mono else ""}">{_e(val)}</div></div>'
        for lbl, val, mono in meta_items
    )

    # ── summary cards ───────────────────────────────────────────────────
    summary_html = f"""\
  <div class="summary-row">
    <div class="summary-card card-avoid"><div class="count">{n_avo}</div><div class="label">Avoid</div></div>
    <div class="summary-card card-caution"><div class="count">{n_cau}</div><div class="label">Caution</div></div>
    <div class="summary-card card-ok"><div class="count">{n_std}</div><div class="label">Standard</div></div>
    {"" if n_ind == 0 else f'<div class="summary-card card-unknown"><div class="count">{n_ind}</div><div class="label">Insufficient Data</div></div>'}
  </div>"""

    # ── alert list ──────────────────────────────────────────────────────
    alert_items = []
    for d in drug_results["avoid"]:
        alert_items.append(
            f'    <li class="alert-avoid"><span class="drug-name">{_e(d["drug"])}</span>'
            f' ({_e(d["brand"])})<span class="gene-tag">{_e(d["gene"])}</span>'
            f' — {_e(d["recommendation"])}</li>'
        )
    for d in drug_results["caution"]:
        alert_items.append(
            f'    <li class="alert-caution"><span class="drug-name">{_e(d["drug"])}</span>'
            f' ({_e(d["brand"])})<span class="gene-tag">{_e(d["gene"])}</span>'
            f' — {_e(d["recommendation"])}</li>'
        )
    alerts_html = "\n".join(alert_items) if alert_items else "    <li>No actionable alerts.</li>"

    # ── gene profiles table ─────────────────────────────────────────────
    gene_rows = []
    for gene in gene_defs:
        if gene in profiles:
            p = profiles[gene]
            gene_rows.append(
                f'      <tr><td><strong>{_e(gene)}</strong></td>'
                f'<td>{_e(gene_defs[gene]["name"])}</td>'
                f'<td>{_e(p["diplotype"])}</td>'
                f'<td>{_e(p["phenotype"])}</td></tr>'
            )
    gene_table_html = "\n".join(gene_rows)

    # ── detected variants table ─────────────────────────────────────────
    var_rows = []
    for rsid, info in sorted(pgx_snps.items(), key=lambda x: x[1]["gene"]):
        var_rows.append(
            f'        <tr><td>{_e(rsid)}</td><td>{_e(info["gene"])}</td>'
            f'<td>{_e(info["allele"])}</td><td><code>{_e(info["genotype"])}</code></td>'
            f'<td>{_e(info["effect"])}</td></tr>'
        )
    variant_table_html = "\n".join(var_rows)

    # ── full drug table ─────────────────────────────────────────────────
    drug_rows = []
    for cat in ["avoid", "caution", "indeterminate", "standard"]:
        for d in sorted(drug_results.get(cat, []), key=lambda x: x["drug"]):
            drug_rows.append(
                f'      <tr><td><strong>{_e(d["drug"])}</strong></td>'
                f'<td>{_e(d["brand"])}</td><td>{_e(d["class"])}</td>'
                f'<td>{_e(d["gene"])}</td><td>{_badge(d["classification"])}</td>'
                f'<td>{_e(d["recommendation"])}</td></tr>'
            )
    drug_table_html = "\n".join(drug_rows)

    # ── assemble ────────────────────────────────────────────────────────
    html = f"""\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ClawBio PharmGx Report — {_e(fname)}</title>
  <style>
{_CSS}
  </style>
</head>
<body>

<!-- Banner -->
<div class="banner">DEMO — FOR RESEARCH USE ONLY</div>

<!-- Sidebar -->
<nav class="sidebar">
{sidebar_html}
</nav>

<!-- Main content -->
<div class="main">

  <h1>ClawBio PharmGx Report</h1>
  <div class="subtitle">Pharmacogenomic analysis &middot; {_e(now)}</div>

  <div class="meta-grid">
{meta_html}
  </div>

  <!-- Summary -->
  <section id="summary">
    <h2>Drug Response Summary</h2>
{summary_html}
  </section>

  <!-- Alerts -->
  <section id="alerts">
    <h2>Actionable Alerts</h2>
    <ul class="alert-list">
{alerts_html}
    </ul>
  </section>

  <!-- Gene Profiles -->
  <section id="genes">
    <h2>Gene Profiles</h2>
    <table>
      <thead><tr><th>Gene</th><th>Full Name</th><th>Diplotype</th><th>Phenotype</th></tr></thead>
      <tbody>
{gene_table_html}
      </tbody>
    </table>
  </section>

  <!-- Detected Variants -->
  <section id="variants">
    <h2>Detected Variants</h2>
    <details open>
      <summary>Show {len(pgx_snps)} variant(s)</summary>
      <table>
        <thead><tr><th>rsID</th><th>Gene</th><th>Star Allele</th><th>Genotype</th><th>Effect</th></tr></thead>
        <tbody>
{variant_table_html}
        </tbody>
      </table>
    </details>
  </section>

  <!-- Complete Drug Recommendations -->
  <section id="drugs">
    <h2>Complete Drug Recommendations</h2>
    <table>
      <thead><tr><th>Drug</th><th>Brand</th><th>Class</th><th>Gene</th><th>Status</th><th>Recommendation</th></tr></thead>
      <tbody>
{drug_table_html}
      </tbody>
    </table>
  </section>

  <!-- Methods -->
  <section id="methods">
    <h2>Methods</h2>
    <ul>
      <li><strong>Tool:</strong> ClawBio PharmGx Reporter v0.1.0</li>
      <li><strong>SNP panel:</strong> {len(pgx_snp_defs)} pharmacogenomic variants across {len(gene_defs)} genes</li>
      <li><strong>Star allele calling:</strong> Simplified DTC-compatible algorithm (single-SNP per allele)</li>
      <li><strong>Phenotype assignment:</strong> CPIC-based diplotype-to-phenotype mapping</li>
      <li><strong>Drug guidelines:</strong> 51 drugs from CPIC (cpicpgx.org), simplified for DTC context</li>
    </ul>
    <p style="margin-top:12px;font-size:13px;color:#64748b;">
      <strong>Input checksum (SHA-256):</strong> <code>{_e(checksum)}</code>
    </p>
    <p style="font-size:13px;color:#64748b;">
      <strong>Reproduce:</strong> <code>python pharmgx_reporter.py --input {_e(fname)} --output report</code>
    </p>
  </section>

  <!-- Disclaimer -->
  <section id="disclaimer">
    <h2>Disclaimer</h2>
    <div class="disclaimer">
      <p><strong>This report is for research and educational purposes only.</strong>
      It is NOT a diagnostic device and should NOT be used to make medication decisions
      without consulting a qualified healthcare professional.</p>
      <p>Pharmacogenomic recommendations are based on CPIC guidelines
      (<a href="https://cpicpgx.org" target="_blank">cpicpgx.org</a>).
      DTC genetic tests have limitations: they may not detect all relevant variants,
      and results should be confirmed by clinical-grade testing before clinical use.</p>
    </div>
  </section>

  <!-- References -->
  <section>
    <h2>References</h2>
    <ul style="font-size:13px;">
      <li>Corpas, M. (2026). ClawBio. <a href="https://github.com/ClawBio/ClawBio" target="_blank">github.com/ClawBio/ClawBio</a></li>
      <li>CPIC. Clinical Pharmacogenetics Implementation Consortium. <a href="https://cpicpgx.org" target="_blank">cpicpgx.org</a></li>
      <li>Caudle, K.E. et al. (2014). Standardizing terms for clinical pharmacogenetic test results. <em>Genet Med</em>, 16(9), 655-663.</li>
      <li>PharmGKB. <a href="https://www.pharmgkb.org" target="_blank">pharmgkb.org</a></li>
    </ul>
  </section>

</div><!-- .main -->

</body>
</html>
"""
    return html
