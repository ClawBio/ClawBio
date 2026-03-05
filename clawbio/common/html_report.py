"""Reusable HTML report builder for ClawBio skills.

Generates self-contained HTML reports with embedded CSS — no external
dependencies (no Jinja2). Reports open cleanly in Telegram's in-app browser,
regular browsers, and mobile devices.
"""

from __future__ import annotations

import html
import math
from datetime import datetime, timezone
from pathlib import Path

from clawbio.common.report import DISCLAIMER

# ---------------------------------------------------------------------------
# Embedded CSS
# ---------------------------------------------------------------------------

_CSS = """\
:root {
  /* Brand palette */
  --cb-green-900: #1b5e20;
  --cb-green-700: #2e7d32;
  --cb-green-500: #43a047;
  --cb-green-100: #e8f5e9;
  --cb-green-50:  #f1f8e9;

  /* Severity */
  --cb-red-700:   #c62828;
  --cb-red-100:   #ffebee;
  --cb-red-50:    #fff5f5;
  --cb-amber-700: #f57f17;
  --cb-amber-100: #fff8e1;
  --cb-grey-700:  #616161;
  --cb-grey-500:  #9e9e9e;
  --cb-grey-300:  #e0e0e0;
  --cb-grey-100:  #f5f5f5;
  --cb-grey-50:   #fafafa;

  /* Surfaces */
  --cb-bg:        #fafafa;
  --cb-surface:   #ffffff;
  --cb-text:      #212121;
  --cb-text-secondary: #616161;
  --cb-border:    #e0e0e0;

  /* Spacing */
  --cb-space-xs: 4px;
  --cb-space-sm: 8px;
  --cb-space-md: 16px;
  --cb-space-lg: 24px;
  --cb-space-xl: 32px;
  --cb-space-2xl: 48px;

  /* Radii */
  --cb-radius-sm: 6px;
  --cb-radius-md: 10px;
  --cb-radius-lg: 16px;

  /* Shadows */
  --cb-shadow-sm: 0 1px 3px rgba(0,0,0,0.08), 0 1px 2px rgba(0,0,0,0.06);
  --cb-shadow-md: 0 4px 6px rgba(0,0,0,0.07), 0 2px 4px rgba(0,0,0,0.06);

  /* Backward-compat aliases */
  --clawbio-green: var(--cb-green-700);
  --clawbio-amber: var(--cb-amber-700);
  --clawbio-red: var(--cb-red-700);
  --clawbio-grey: var(--cb-grey-700);
  --clawbio-bg: var(--cb-bg);
  --clawbio-card-bg: var(--cb-surface);
}
*, *::before, *::after { box-sizing: border-box; }
body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  line-height: 1.6;
  color: var(--cb-text);
  background: var(--cb-bg);
  margin: 0;
  padding: var(--cb-space-md);
  max-width: 960px;
  margin-left: auto;
  margin-right: auto;
}
h1 { color: var(--cb-green-700); border-bottom: 3px solid var(--cb-green-700); padding-bottom: 8px; }
h2 { color: #424242; margin-top: var(--cb-space-xl); }
h3 { color: var(--cb-grey-700); margin-top: var(--cb-space-lg); }

/* Branded header */
.report-header {
  background: linear-gradient(135deg, #1b5e20 0%, #2e7d32 50%, #388e3c 100%);
  color: white;
  padding: var(--cb-space-lg) var(--cb-space-xl);
  border-radius: var(--cb-radius-lg);
  margin: 0 0 var(--cb-space-xl) 0;
  position: relative;
  overflow: hidden;
}
.report-header .header-logo {
  position: absolute;
  right: -20px;
  top: 50%;
  transform: translateY(-50%);
  opacity: 0.10;
  width: 200px;
  height: 200px;
}
.report-header h1 {
  margin: 0; font-size: 1.8em; font-weight: 700; border: none; padding: 0;
  color: white; letter-spacing: -0.02em; position: relative;
}
.report-header .subtitle {
  margin: var(--cb-space-xs) 0 0 0; font-size: 0.95em; opacity: 0.9; font-weight: 400;
  position: relative;
}

/* Metadata block */
.metadata { background: var(--cb-green-100); border-radius: var(--cb-radius-sm); padding: 12px 16px; margin-bottom: var(--cb-space-lg); }
.metadata p { margin: 4px 0; font-size: 0.95em; }
.metadata strong { color: var(--cb-green-900); }

/* Tables */
table { width: 100%; border-collapse: collapse; margin: var(--cb-space-md) 0; font-size: 0.9em; }
th { background: var(--cb-green-100); color: var(--cb-green-900); text-align: left; padding: 10px 12px; border-bottom: 2px solid #a5d6a7; }
td { padding: 8px 12px; border-bottom: 1px solid var(--cb-border); }
tr:nth-child(even) { background: var(--cb-grey-100); }
tr:hover { background: var(--cb-green-100); }

/* Table wrapper for mobile scroll */
.table-wrap {
  overflow-x: auto; -webkit-overflow-scrolling: touch;
  margin: var(--cb-space-md) 0; border-radius: var(--cb-radius-md); border: 1px solid var(--cb-border);
}
.table-wrap table { margin: 0; border: none; }

/* Severity-banded table rows */
tr.row-avoid { background: var(--cb-red-50); }
tr.row-avoid:hover { background: var(--cb-red-100); }
tr.row-caution { background: #fffde7; }
tr.row-caution:hover { background: var(--cb-amber-100); }
tr.row-indeterminate { background: var(--cb-grey-50); }
tr.row-indeterminate:hover { background: var(--cb-grey-100); }
tr.row-standard { background: var(--cb-green-50); }
tr.row-standard:hover { background: var(--cb-green-100); }

/* Badges */
.badge { display: inline-flex; align-items: center; gap: 4px; padding: 3px 12px;
         border-radius: 20px; font-size: 0.78em; font-weight: 700;
         text-transform: uppercase; letter-spacing: 0.05em; white-space: nowrap; }
.badge-standard { background: #c8e6c9; color: #1b5e20; }
.badge-caution { background: #fff9c4; color: #e65100; }
.badge-avoid { background: #ffcdd2; color: #b71c1c; }
.badge-indeterminate { background: #e0e0e0; color: #424242; }

/* Alert boxes */
.alert-box { border-left: 4px solid; border-radius: var(--cb-radius-sm); padding: 12px 16px; margin: 12px 0; }
.alert-box-avoid { border-color: var(--cb-red-700); background: var(--cb-red-100); }
.alert-box-caution { border-color: var(--cb-amber-700); background: var(--cb-amber-100); }
.alert-box-info { border-color: var(--cb-grey-700); background: var(--cb-grey-100); }
.alert-box h4 { margin: 0 0 8px 0; }
.alert-box-avoid h4 { color: var(--cb-red-700); }
.alert-box-caution h4 { color: var(--cb-amber-700); }
.alert-box-info h4 { color: var(--cb-grey-700); }

/* Summary cards */
.summary-cards {
  display: grid; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
  gap: var(--cb-space-md); margin: var(--cb-space-md) 0 var(--cb-space-xl) 0;
}
.summary-card {
  background: var(--cb-surface); border-radius: var(--cb-radius-md);
  padding: var(--cb-space-lg) var(--cb-space-md); text-align: center;
  box-shadow: var(--cb-shadow-sm); border-top: 4px solid transparent;
  transition: box-shadow 0.2s ease, transform 0.2s ease;
}
.summary-card:hover { box-shadow: var(--cb-shadow-md); transform: translateY(-2px); }
.summary-card .count { font-size: 2.5em; font-weight: 800; display: block; line-height: 1.1; }
.summary-card .label {
  font-size: 0.8em; text-transform: uppercase; letter-spacing: 0.08em;
  color: var(--cb-text-secondary); margin-top: var(--cb-space-xs); display: block;
}
.summary-card.avoid { border-top-color: var(--cb-red-700); }
.summary-card.avoid .count { color: var(--cb-red-700); }
.summary-card.caution { border-top-color: var(--cb-amber-700); }
.summary-card.caution .count { color: var(--cb-amber-700); }
.summary-card.standard { border-top-color: var(--cb-green-700); }
.summary-card.standard .count { color: var(--cb-green-700); }
.summary-card.indeterminate { border-top-color: var(--cb-grey-500); }
.summary-card.indeterminate .count { color: var(--cb-grey-500); }

/* Executive summary */
.exec-summary {
  margin: var(--cb-space-lg) 0;
}
.exec-summary h3 { margin: 0 0 var(--cb-space-sm) 0; color: var(--cb-green-900); font-size: 1.1em; }
.exec-summary-grid { display: grid; grid-template-columns: 1fr 1fr; gap: var(--cb-space-md); }
.exec-stat {
  display: flex; align-items: flex-start; gap: var(--cb-space-md);
  background: var(--cb-surface); border: 1px solid var(--cb-border);
  border-left: 4px solid var(--cb-grey-300); border-radius: var(--cb-radius-sm);
  padding: var(--cb-space-md); box-shadow: var(--cb-shadow-sm);
}
.exec-stat.stat-avoid { border-left-color: var(--cb-red-700); }
.exec-stat.stat-caution { border-left-color: var(--cb-amber-700); }
.exec-stat.stat-ok { border-left-color: var(--cb-green-700); }
.exec-stat .stat-icon { font-size: 1.6em; flex-shrink: 0; line-height: 1; }
.exec-stat .stat-text { font-size: 0.9em; line-height: 1.5; color: var(--cb-text-secondary); }
.exec-stat .stat-text strong { display: block; font-size: 1.05em; color: var(--cb-text); margin-bottom: 2px; }

/* Donut chart */
.donut-chart-section {
  display: flex; align-items: center; gap: var(--cb-space-xl);
  margin: var(--cb-space-lg) 0; flex-wrap: wrap; justify-content: center;
}
.donut-chart-section svg { flex-shrink: 0; }
.donut-legend { display: flex; flex-direction: column; gap: var(--cb-space-sm); }
.donut-legend-item { display: flex; align-items: center; gap: var(--cb-space-sm); font-size: 0.9em; }
.donut-legend-swatch { width: 14px; height: 14px; border-radius: 3px; flex-shrink: 0; }

/* Progress bars */
.progress-bar-container {
  width: 100%; background: var(--cb-grey-300); border-radius: 20px;
  overflow: hidden; height: 12px; margin: var(--cb-space-xs) 0;
}
.progress-bar-fill { height: 100%; border-radius: 20px; }
.progress-bar-fill.fill-green { background: linear-gradient(90deg, #43a047, #2e7d32); }
.progress-bar-fill.fill-amber { background: linear-gradient(90deg, #ffa726, #f57f17); }
.progress-bar-fill.fill-red { background: linear-gradient(90deg, #ef5350, #c62828); }
.progress-bar-fill.fill-grey { background: linear-gradient(90deg, #bdbdbd, #757575); }

/* Disclaimer */
.disclaimer { background: #fff3e0; border: 1px solid #ffcc80; border-radius: var(--cb-radius-sm);
              padding: var(--cb-space-xs) var(--cb-space-md); margin: 0 0 var(--cb-space-sm) 0; font-size: 0.8em; color: #e65100; line-height: 1.4; }
.disclaimer-bottom { margin: var(--cb-space-xl) 0 0 0; }

/* Gene links (smaller text) */
.gene-links { font-size: 0.82em; line-height: 1.6; }

/* Collapsible details */
details { border: 1px solid var(--cb-border); border-radius: var(--cb-radius-md); margin: var(--cb-space-md) 0; overflow: hidden; }
details summary {
  cursor: pointer; padding: var(--cb-space-md); font-weight: 600; color: var(--cb-text);
  background: var(--cb-grey-50); border-radius: var(--cb-radius-md);
  list-style: none; display: flex; align-items: center; justify-content: space-between;
}
details summary::-webkit-details-marker { display: none; }
details summary::after {
  content: "\\25B6"; font-size: 0.7em; color: var(--cb-grey-500); transition: transform 0.2s ease;
}
details[open] summary::after { transform: rotate(90deg); }
details[open] summary { border-bottom: 1px solid var(--cb-border); border-radius: var(--cb-radius-md) var(--cb-radius-md) 0 0; }
details > :not(summary) { padding: 0 var(--cb-space-md); }

/* Branded footer */
.report-footer {
  margin-top: var(--cb-space-2xl); padding-top: var(--cb-space-lg);
  border-top: 2px solid var(--cb-green-100); text-align: center;
  color: var(--cb-text-secondary); font-size: 0.85em;
}
.report-footer .footer-brand { font-weight: 700; color: var(--cb-green-700); }

/* Responsive */
@media (max-width: 600px) {
  body { padding: var(--cb-space-sm); }
  table { font-size: 0.8em; }
  th, td { padding: 6px 8px; }
  .summary-cards { grid-template-columns: repeat(2, 1fr); }
  .summary-card { padding: var(--cb-space-md); }
  .summary-card .count { font-size: 1.8em; }
  .exec-summary-grid { grid-template-columns: 1fr; }
  .report-header { padding: var(--cb-space-md); }
  .report-header h1 { font-size: 1.4em; }
}

/* Evidence badges */
.badge-evidence-high { background: #c8e6c9; color: #1b5e20; }
.badge-evidence-moderate { background: #bbdefb; color: #0d47a1; }
.badge-evidence-low { background: #fff9c4; color: #e65100; }
.badge-evidence-minimal { background: #e0e0e0; color: #424242; }
.badge-evidence-na { background: #e0e0e0; color: #757575; }
.evidence-verified { color: #2e7d32; font-size: 0.82em; font-weight: 600; }
.evidence-unverified { color: #9e9e9e; font-size: 0.82em; font-style: italic; }
.evidence-source { color: #757575; font-size: 0.78em; display: block; margin-top: 2px; }
.rec-action { margin: 4px 0 0; font-size: 0.88em; line-height: 1.5; }
.evidence-recs { margin-top: 6px; padding-top: 6px; border-top: 1px solid var(--cb-grey-300); font-size: 0.78em; line-height: 1.5; }
.evidence-rec-source {
  display: inline-block; font-weight: 700; color: #fff; background: var(--cb-green-700);
  border-radius: 3px; padding: 0 5px; font-size: 0.85em; margin-right: 4px; vertical-align: middle;
}
.evidence-rec-text { color: var(--cb-text-secondary); }

/* Print */
@media print {
  body { background: white; max-width: 100%; padding: 0; }
  .report-header { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
  .summary-card { break-inside: avoid; box-shadow: none; border: 1px solid #ccc; }
  details { break-inside: avoid; }
  .disclaimer { break-inside: avoid; }
  a::after { content: " (" attr(href) ")"; font-size: 0.8em; color: #666; }
  tr { break-inside: avoid; }
}
"""

# Badge CSS class mapping
_BADGE_CLASS = {
    "standard": "badge-standard",
    "caution": "badge-caution",
    "avoid": "badge-avoid",
    "indeterminate": "badge-indeterminate",
}

# Badge display labels
_BADGE_LABEL = {
    "standard": "Standard",
    "caution": "Caution",
    "avoid": "Avoid",
    "indeterminate": "Insufficient",
}


class HtmlReportBuilder:
    """Builds a self-contained HTML report with embedded CSS."""

    def __init__(self, title: str, skill: str, extra_css: str = "") -> None:
        self._title = html.escape(title)
        self._skill = html.escape(skill)
        self._extra_css = extra_css
        self._sections: list[str] = []
        self._custom_header = False
        self._custom_footer = False

    # -- Building blocks ---------------------------------------------------

    # SVG claw logo for header background
    _CLAW_SVG = (
        '<svg class="header-logo" viewBox="0 0 200 200" xmlns="http://www.w3.org/2000/svg">'
        '<g fill="white" fill-rule="evenodd">'
        # Three claw marks
        '<path d="M60 30 C55 80, 40 130, 50 170 C55 175, 65 175, 68 168 '
        'C78 130, 70 80, 75 35 Z"/>'
        '<path d="M95 20 C90 75, 85 130, 90 175 C93 182, 103 182, 106 175 '
        'C111 130, 106 75, 101 20 Z"/>'
        '<path d="M132 30 C137 80, 152 130, 142 170 C139 175, 129 175, 126 168 '
        'C116 130, 124 80, 119 35 Z"/>'
        '</g></svg>'
    )

    def add_header_block(self, title: str, subtitle: str = "") -> "HtmlReportBuilder":
        """Add a branded gradient header with SVG claw logo background."""
        sub = f'<p class="subtitle">{html.escape(subtitle)}</p>' if subtitle else ""
        self._sections.append(
            f'<div class="report-header">'
            f"{self._CLAW_SVG}"
            f"<h1>{html.escape(title)}</h1>"
            f"{sub}"
            f"</div>"
        )
        self._custom_header = True
        return self

    def add_metadata(self, items: dict[str, str]) -> "HtmlReportBuilder":
        parts = []
        for key, val in items.items():
            parts.append(f"<p><strong>{html.escape(key)}:</strong> {html.escape(str(val))}</p>")
        self._sections.append(f'<div class="metadata">{"".join(parts)}</div>')
        return self

    def add_section(self, heading: str, level: int = 2) -> "HtmlReportBuilder":
        tag = f"h{min(max(level, 1), 6)}"
        self._sections.append(f"<{tag}>{html.escape(heading)}</{tag}>")
        return self

    def add_paragraph(self, text: str, css_class: str = "") -> "HtmlReportBuilder":
        cls = f' class="{html.escape(css_class)}"' if css_class else ""
        self._sections.append(f"<p{cls}>{html.escape(text)}</p>")
        return self

    def add_summary_cards(self, cards: list[tuple[str, int, str]]) -> "HtmlReportBuilder":
        parts = []
        for label, count, category in cards:
            cat_class = html.escape(category)
            parts.append(
                f'<div class="summary-card {cat_class}">'
                f'<span class="count">{int(count)}</span>'
                f'<span class="label">{html.escape(label)}</span>'
                f"</div>"
            )
        self._sections.append(f'<div class="summary-cards">{"".join(parts)}</div>')
        return self

    def add_alert_box(self, severity: str, title: str, body: str) -> "HtmlReportBuilder":
        sev = severity if severity in ("avoid", "caution", "info") else "info"
        self._sections.append(
            f'<div class="alert-box alert-box-{sev}">'
            f"<h4>{html.escape(title)}</h4>"
            f"<p>{html.escape(body)}</p>"
            f"</div>"
        )
        return self

    def add_executive_summary(
        self, items: list[tuple[str, str, str, str]]
    ) -> "HtmlReportBuilder":
        """Add a summary panel with icon/stat items.

        Each item is (icon, title, description, severity).
        severity: "avoid", "caution", "ok" — controls left border colour.
        For backward compat, 3-tuples (no severity) default to no colour.
        """
        stats = []
        for item in items:
            if len(item) == 4:
                icon, title, desc, severity = item
            else:
                icon, title, desc = item
                severity = ""
            sev_cls = f" stat-{html.escape(severity)}" if severity else ""
            stats.append(
                f'<div class="exec-stat{sev_cls}">'
                f'<span class="stat-icon">{icon}</span>'
                f'<div class="stat-text"><strong>{html.escape(title)}</strong>'
                f"<span>{html.escape(desc)}</span></div></div>"
            )
        self._sections.append(
            f'<div class="exec-summary">'
            f"<h3>Summary</h3>"
            f'<div class="exec-summary-grid">{"".join(stats)}</div></div>'
        )
        return self

    def add_donut_chart(
        self,
        segments: list[tuple[str, int, str]],
        size: int = 180,
    ) -> "HtmlReportBuilder":
        """Add an inline SVG donut chart with legend."""
        total = sum(s[1] for s in segments) or 1
        radius = 60
        circumference = 2 * math.pi * radius
        cx, cy = size // 2, size // 2

        circles = []
        offset = 0.0
        for _label, count, color in segments:
            if count == 0:
                continue
            dash = (count / total) * circumference
            gap = circumference - dash
            circles.append(
                f'<circle cx="{cx}" cy="{cy}" r="{radius}" fill="none" '
                f'stroke="{html.escape(color)}" stroke-width="24" '
                f'stroke-dasharray="{dash:.2f} {gap:.2f}" '
                f'stroke-dashoffset="{-offset:.2f}" '
                f'transform="rotate(-90 {cx} {cy})" />'
            )
            offset += dash

        svg = (
            f'<svg width="{size}" height="{size}" viewBox="0 0 {size} {size}" '
            f'xmlns="http://www.w3.org/2000/svg" role="img" '
            f'aria-label="Drug breakdown chart">'
            + "".join(circles)
            + f'<text x="{cx}" y="{cy}" text-anchor="middle" dy="0.35em" '
            f'font-size="28" font-weight="700" fill="currentColor">{total}</text>'
            f'<text x="{cx}" y="{cy + 18}" text-anchor="middle" '
            f'font-size="11" fill="#757575">drugs</text>'
            f"</svg>"
        )

        legend_items = []
        for label, count, color in segments:
            legend_items.append(
                f'<div class="donut-legend-item">'
                f'<span class="donut-legend-swatch" style="background:{html.escape(color)}"></span>'
                f"{html.escape(label)}: <strong>{count}</strong>"
                f"</div>"
            )

        self._sections.append(
            f'<div class="donut-chart-section">{svg}'
            f'<div class="donut-legend">{"".join(legend_items)}</div></div>'
        )
        return self

    def add_progress_bar(
        self, label: str, value: int, maximum: int, color: str = "green"
    ) -> "HtmlReportBuilder":
        """Add a labeled progress bar."""
        pct = min(100, (value / maximum * 100)) if maximum else 0
        self._sections.append(
            f'<div style="margin:8px 0;">'
            f'<div style="display:flex;justify-content:space-between;font-size:0.85em;">'
            f"<span>{html.escape(label)}</span>"
            f"<span>{value}/{maximum} ({pct:.0f}%)</span></div>"
            f'<div class="progress-bar-container">'
            f'<div class="progress-bar-fill fill-{html.escape(color)}" '
            f'style="width:{pct:.1f}%"></div></div></div>'
        )
        return self

    def add_table(
        self,
        headers: list[str],
        rows: list[list[str]],
        badge_col: int | None = None,
    ) -> "HtmlReportBuilder":
        parts = ["<table><thead><tr>"]
        for h in headers:
            parts.append(f"<th>{html.escape(h)}</th>")
        parts.append("</tr></thead><tbody>")
        for row in rows:
            parts.append("<tr>")
            for i, cell in enumerate(row):
                if i == badge_col:
                    badge_cls = _BADGE_CLASS.get(cell, "badge-indeterminate")
                    badge_lbl = _BADGE_LABEL.get(cell, html.escape(cell))
                    parts.append(f'<td><span class="badge {badge_cls}">{badge_lbl}</span></td>')
                else:
                    parts.append(f"<td>{html.escape(str(cell))}</td>")
            parts.append("</tr>")
        parts.append("</tbody></table>")
        self._sections.append("".join(parts))
        return self

    def add_table_wrapped(
        self,
        headers: list[str],
        rows: list[list[str]],
        badge_col: int | None = None,
        row_classes: list[str] | None = None,
    ) -> "HtmlReportBuilder":
        """Add a table wrapped for mobile scrolling with optional row classes."""
        parts = ['<div class="table-wrap"><table><thead><tr>']
        for h in headers:
            parts.append(f"<th>{html.escape(h)}</th>")
        parts.append("</tr></thead><tbody>")
        for idx, row in enumerate(rows):
            cls = ""
            if row_classes and idx < len(row_classes):
                cls = f' class="{html.escape(row_classes[idx])}"'
            parts.append(f"<tr{cls}>")
            for i, cell in enumerate(row):
                if i == badge_col:
                    badge_cls = _BADGE_CLASS.get(cell, "badge-indeterminate")
                    badge_lbl = _BADGE_LABEL.get(cell, html.escape(cell))
                    parts.append(f'<td><span class="badge {badge_cls}">{badge_lbl}</span></td>')
                else:
                    parts.append(f"<td>{html.escape(str(cell))}</td>")
            parts.append("</tr>")
        parts.append("</tbody></table></div>")
        self._sections.append("".join(parts))
        return self

    def add_details(self, summary_text: str, content_html: str) -> "HtmlReportBuilder":
        """Add a collapsible <details>/<summary> section (collapsed by default)."""
        self._sections.append(
            f"<details><summary>{html.escape(summary_text)}</summary>"
            f"{content_html}"
            f"</details>"
        )
        return self

    def add_raw_html(self, raw: str) -> "HtmlReportBuilder":
        self._sections.append(raw)
        return self

    def add_disclaimer(self) -> "HtmlReportBuilder":
        """Add disclaimer inline at the current position in the report."""
        self._sections.append(f'<div class="disclaimer"><strong>Disclaimer:</strong> {html.escape(DISCLAIMER)}</div>')
        self._has_disclaimer = True
        return self

    def add_footer_block(self, skill: str, version: str = "") -> "HtmlReportBuilder":
        """Add a branded footer with timestamp."""
        now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
        ver = f" v{html.escape(version)}" if version else ""
        self._sections.append(
            f'<div class="report-footer">'
            f'<p>Generated by <span class="footer-brand">ClawBio</span> '
            f"\u00b7 {html.escape(skill)}{ver} \u00b7 {now}</p>"
            f'<p style="font-size:0.8em;margin-top:4px;">Genetic data processed locally. '
            f"No data was transmitted to external servers.</p></div>"
        )
        self._custom_footer = True
        return self

    # -- Render ------------------------------------------------------------

    def render(self) -> str:
        now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
        body = "\n".join(self._sections)
        disclaimer_bottom = ""
        if getattr(self, "_has_disclaimer", False):
            disclaimer_bottom = (
                f'<div class="disclaimer disclaimer-bottom">'
                f'<strong>Disclaimer:</strong> {html.escape(DISCLAIMER)}</div>'
            )
        title_block = "" if self._custom_header else f"<h1>{self._title}</h1>"
        footer = "" if self._custom_footer else (
            f'<p style="color:#757575;font-size:0.9em;margin-top:32px;">'
            f"Generated by ClawBio &middot; {html.escape(self._skill)} &middot; {now}</p>"
        )
        return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{self._title}</title>
<style>
{_CSS}
{self._extra_css}
</style>
</head>
<body>
{title_block}
{body}
{disclaimer_bottom}
{footer}
</body>
</html>"""


def write_html_report(output_dir: str | Path, filename: str, content: str) -> Path:
    """Write an HTML string to *output_dir*/*filename* and return the path."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / filename
    path.write_text(content, encoding="utf-8")
    return path
