#!/usr/bin/env python3
"""ClawBio differential expression visualisation skill.

Port of ClawBio 1.x `diff-visualizer` for the ClawBio 2.0 `clawbio` CLI. The
analysis/plotting/reporting logic (bulk volcano/MA/bar/heatmap plots, scRNA
contrast/marker plots, AnnData-based enhancement, report.md/report.html/
result.json generation) is unchanged from the 1.x source.

Adapted for the CLI surface only:
  - Dropped the `from clawbio.common.*` imports and the `sys.path`
    PROJECT_ROOT hack (that module belongs to the 1.x repo layout).
    `sha256_file`, `DISCLAIMER`, `generate_report_header`,
    `generate_report_footer`, `write_result_json`, `HtmlReportBuilder`, and
    `write_html_report` below are local, dependency-free re-implementations
    of the same helpers with the same behavior — same pattern as
    skills/rnaseq_de/scripts/rnaseq_de.py and
    skills/claw_ancestry_pca/scripts/ancestry_pca.py in this repo.
  - `_examples_dir()` now resolves to `<skill_dir>/examples` (one level up
    from `scripts/`, matching this repo's skill layout convention) instead
    of `<script_dir>/examples`.
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import html
import json
import math
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Inlined from ClawBio 1.x `clawbio/common/checksums.py`, `report.py`,
# `html_report.py` (see module docstring for why these are local copies).
# ---------------------------------------------------------------------------


def sha256_file(filepath: str | Path) -> str:
    """Compute SHA-256 checksum of a file (full 64-char hex digest)."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)


def generate_report_header(
    title: str,
    skill_name: str,
    skill_version: str = "",
    input_files: list[Path] | None = None,
    extra_metadata: dict[str, str] | None = None,
) -> str:
    """Generate the standard markdown report header."""
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    checksums = []
    if input_files:
        for f in input_files:
            f = Path(f)
            if f.exists():
                checksums.append(f"- `{f.name}`: `{sha256_file(f)}`")
            else:
                checksums.append(f"- `{f.name}`: (file not found)")

    lines = [
        f"# {title}",
        "",
        f"**Date**: {now}",
        f"**Skill**: {skill_name}",
        f"**Version**: {skill_version}",
    ]
    if extra_metadata:
        for key, val in extra_metadata.items():
            lines.append(f"**{key}**: {val}")
    if checksums:
        lines.append("**Input files**:")
        lines.extend(checksums)
    lines.extend(["", "---", ""])

    return "\n".join(lines)


def generate_report_footer() -> str:
    """Generate the standard markdown report footer with disclaimer."""
    return f"""
---

## Disclaimer

*{DISCLAIMER}*
"""


def write_result_json(
    output_dir: str | Path,
    skill: str,
    version: str,
    summary: dict[str, Any],
    data: dict[str, Any],
    input_checksum: str = "",
    datasets: dict[str, str] | None = None,
    status: str | None = None,
    ok: bool | None = None,
) -> Path:
    """Write the standardized result.json envelope alongside report.md."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    envelope: dict[str, Any] = {}
    if status is not None:
        envelope["status"] = status
    if ok is not None:
        envelope["ok"] = ok
    envelope.update(
        {
            "skill": skill,
            "version": version,
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "input_checksum": f"sha256:{input_checksum}" if input_checksum else "",
            "datasets": datasets or {},
            "summary": summary,
            "data": data,
        }
    )

    result_path = output_dir / "result.json"
    result_path.write_text(json.dumps(envelope, indent=2, default=str))
    return result_path


_HTML_REPORT_CSS = """\
:root {
  --cb-green-900: #1b5e20;
  --cb-green-700: #2e7d32;
  --cb-green-500: #43a047;
  --cb-green-100: #e8f5e9;
  --cb-green-50:  #f1f8e9;
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
  --cb-bg:        #fafafa;
  --cb-surface:   #ffffff;
  --cb-text:      #212121;
  --cb-text-secondary: #616161;
  --cb-border:    #e0e0e0;
  --cb-space-xs: 4px;
  --cb-space-sm: 8px;
  --cb-space-md: 16px;
  --cb-space-lg: 24px;
  --cb-space-xl: 32px;
  --cb-space-2xl: 48px;
  --cb-radius-sm: 6px;
  --cb-radius-md: 10px;
  --cb-radius-lg: 16px;
  --cb-shadow-sm: 0 1px 3px rgba(0,0,0,0.08), 0 1px 2px rgba(0,0,0,0.06);
  --cb-shadow-md: 0 4px 6px rgba(0,0,0,0.07), 0 2px 4px rgba(0,0,0,0.06);
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
.metadata { background: var(--cb-green-100); border-radius: var(--cb-radius-sm); padding: 12px 16px; margin-bottom: var(--cb-space-lg); }
.metadata p { margin: 4px 0; font-size: 0.95em; }
.metadata strong { color: var(--cb-green-900); }
table { width: 100%; border-collapse: collapse; margin: var(--cb-space-md) 0; font-size: 0.9em; }
th { background: var(--cb-green-100); color: var(--cb-green-900); text-align: left; padding: 10px 12px; border-bottom: 2px solid #a5d6a7; }
td { padding: 8px 12px; border-bottom: 1px solid var(--cb-border); }
tr:nth-child(even) { background: var(--cb-grey-100); }
tr:hover { background: var(--cb-green-100); }
.table-wrap {
  overflow-x: auto; -webkit-overflow-scrolling: touch;
  margin: var(--cb-space-md) 0; border-radius: var(--cb-radius-md); border: 1px solid var(--cb-border);
}
.table-wrap table { margin: 0; border: none; }
tr.row-avoid { background: var(--cb-red-50); }
tr.row-avoid:hover { background: var(--cb-red-100); }
tr.row-caution { background: #fffde7; }
tr.row-caution:hover { background: var(--cb-amber-100); }
tr.row-indeterminate { background: var(--cb-grey-50); }
tr.row-indeterminate:hover { background: var(--cb-grey-100); }
tr.row-standard { background: var(--cb-green-50); }
tr.row-standard:hover { background: var(--cb-green-100); }
.badge { display: inline-flex; align-items: center; gap: 4px; padding: 3px 12px;
         border-radius: 20px; font-size: 0.78em; font-weight: 700;
         text-transform: uppercase; letter-spacing: 0.05em; white-space: nowrap; }
.badge-standard { background: #c8e6c9; color: #1b5e20; }
.badge-caution { background: #fff9c4; color: #e65100; }
.badge-avoid { background: #ffcdd2; color: #b71c1c; }
.badge-indeterminate { background: #e0e0e0; color: #424242; }
.alert-box { border-left: 4px solid; border-radius: var(--cb-radius-sm); padding: 12px 16px; margin: 12px 0; }
.alert-box-avoid { border-color: var(--cb-red-700); background: var(--cb-red-100); }
.alert-box-caution { border-color: var(--cb-amber-700); background: var(--cb-amber-100); }
.alert-box-info { border-color: var(--cb-grey-700); background: var(--cb-grey-100); }
.alert-box h4 { margin: 0 0 8px 0; }
.alert-box-avoid h4 { color: var(--cb-red-700); }
.alert-box-caution h4 { color: var(--cb-amber-700); }
.alert-box-info h4 { color: var(--cb-grey-700); }
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
.donut-chart-section {
  display: flex; align-items: center; gap: var(--cb-space-xl);
  margin: var(--cb-space-lg) 0; flex-wrap: wrap; justify-content: center;
}
.donut-chart-section svg { flex-shrink: 0; }
.donut-legend { display: flex; flex-direction: column; gap: var(--cb-space-sm); }
.donut-legend-item { display: flex; align-items: center; gap: var(--cb-space-sm); font-size: 0.9em; }
.donut-legend-swatch { width: 14px; height: 14px; border-radius: 3px; flex-shrink: 0; }
.progress-bar-container {
  width: 100%; background: var(--cb-grey-300); border-radius: 20px;
  overflow: hidden; height: 12px; margin: var(--cb-space-xs) 0;
}
.progress-bar-fill { height: 100%; border-radius: 20px; }
.progress-bar-fill.fill-green { background: linear-gradient(90deg, #43a047, #2e7d32); }
.progress-bar-fill.fill-amber { background: linear-gradient(90deg, #ffa726, #f57f17); }
.progress-bar-fill.fill-red { background: linear-gradient(90deg, #ef5350, #c62828); }
.progress-bar-fill.fill-grey { background: linear-gradient(90deg, #bdbdbd, #757575); }
.disclaimer { background: #fff3e0; border: 1px solid #ffcc80; border-radius: var(--cb-radius-sm);
              padding: var(--cb-space-xs) var(--cb-space-md); margin: 0 0 var(--cb-space-sm) 0; font-size: 0.8em; color: #e65100; line-height: 1.4; }
.disclaimer-bottom { margin: var(--cb-space-xl) 0 0 0; }
.gene-links { font-size: 0.82em; line-height: 1.6; }
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
.report-footer {
  margin-top: var(--cb-space-2xl); padding-top: var(--cb-space-lg);
  border-top: 2px solid var(--cb-green-100); text-align: center;
  color: var(--cb-text-secondary); font-size: 0.85em;
}
.report-footer .footer-brand { font-weight: 700; color: var(--cb-green-700); }
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

_BADGE_CLASS = {
    "standard": "badge-standard",
    "caution": "badge-caution",
    "avoid": "badge-avoid",
    "indeterminate": "badge-indeterminate",
}

_BADGE_LABEL = {
    "standard": "Standard",
    "caution": "Caution",
    "avoid": "Avoid",
    "indeterminate": "Indeterminate",
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

    _CLAW_SVG = (
        '<svg class="header-logo" viewBox="0 0 200 200" xmlns="http://www.w3.org/2000/svg">'
        '<g fill="white" fill-rule="evenodd">'
        '<path d="M60 30 C55 80, 40 130, 50 170 C55 175, 65 175, 68 168 '
        'C78 130, 70 80, 75 35 Z"/>'
        '<path d="M95 20 C90 75, 85 130, 90 175 C93 182, 103 182, 106 175 '
        'C111 130, 106 75, 101 20 Z"/>'
        '<path d="M132 30 C137 80, 152 130, 142 170 C139 175, 129 175, 126 168 '
        'C116 130, 124 80, 119 35 Z"/>'
        '</g></svg>'
    )

    def add_header_block(self, title: str, subtitle: str = "") -> "HtmlReportBuilder":
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
        self._sections.append(f'<div class="disclaimer"><strong>Disclaimer:</strong> {html.escape(DISCLAIMER)}</div>')
        self._has_disclaimer = True
        return self

    def add_footer_block(self, skill: str, version: str = "") -> "HtmlReportBuilder":
        now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
        ver = f" v{html.escape(version)}" if version else ""
        self._sections.append(
            f'<div class="report-footer">'
            f'<p>Generated by <span class="footer-brand">ClawBio</span> '
            f"· {html.escape(skill)}{ver} · {now}</p>"
            f'<p style="font-size:0.8em;margin-top:4px;">Genetic data processed locally. '
            f"No data was transmitted to external servers.</p></div>"
        )
        self._custom_footer = True
        return self

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
{_HTML_REPORT_CSS}
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


VERSION = "0.1.0"

try:
    plt.style.use("seaborn-v0_8-whitegrid")
except OSError:
    pass

PLOT_COLORS = {
    "up": "#D55E00",
    "down": "#0072B2",
    "neutral": "#B8C4D6",
    "teal": "#009E73",
    "gold": "#E3B505",
    "purple": "#7A5195",
    "ink": "#102A43",
    "muted": "#486581",
    "grid": "#D9E2EC",
    "panel": "#FCFDFE",
    "spine": "#BCCCDC",
}

plt.rcParams.update(
    {
        "figure.facecolor": "#FFFFFF",
        "axes.facecolor": PLOT_COLORS["panel"],
        "axes.edgecolor": PLOT_COLORS["spine"],
        "axes.labelcolor": PLOT_COLORS["ink"],
        "axes.titlecolor": PLOT_COLORS["ink"],
        "axes.titlesize": 16,
        "axes.titleweight": "semibold",
        "axes.labelsize": 12,
        "xtick.color": PLOT_COLORS["muted"],
        "ytick.color": PLOT_COLORS["muted"],
        "grid.color": PLOT_COLORS["grid"],
        "grid.linestyle": "-",
        "grid.linewidth": 0.8,
        "legend.frameon": True,
        "legend.edgecolor": PLOT_COLORS["spine"],
        "legend.facecolor": "#FFFFFF",
        "legend.framealpha": 0.94,
        "font.size": 11,
    }
)


@dataclass
class InputBundle:
    mode: str
    source_kind: str
    input_path: Path | None
    input_files: list[Path]
    bulk_df: pd.DataFrame | None = None
    scrna_contrast_df: pd.DataFrame | None = None
    scrna_within_cluster_contrast_df: pd.DataFrame | None = None
    scrna_markers_df: pd.DataFrame | None = None
    upstream_result: dict[str, Any] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


@dataclass
class RunArtifacts:
    notes: list[str] = field(default_factory=list)
    figure_paths: list[Path] = field(default_factory=list)
    table_paths: list[Path] = field(default_factory=list)
    n_input_rows: int = 0
    n_significant: int = 0
    top_table_preview: pd.DataFrame | None = None
    enhanced_inputs_used: list[str] = field(default_factory=list)
    groupby_used: str = ""
    mode_details: dict[str, Any] = field(default_factory=dict)


def _sep_for(path: Path) -> str:
    return "\t" if path.suffix.lower() == ".tsv" else ","


def _read_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=_sep_for(path))
    df.columns = [str(col).strip() for col in df.columns]
    return df


def _write_csv(path: Path, df: pd.DataFrame) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path


def _maybe_numeric(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in columns:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def _detect_table_kind(df: pd.DataFrame) -> str | None:
    cols = {str(col) for col in df.columns}
    if {"gene", "log2FoldChange"} <= cols and ("padj" in cols or "pvalue" in cols):
        return "bulk"
    if {"cluster", "comparison_id", "group1", "group2", "names", "scores"} <= cols:
        return "scrna_within_cluster_contrast"
    if {"cluster", "names", "scores"} <= cols:
        return "scrna_markers"
    if {"names", "scores"} <= cols:
        return "scrna_contrast"
    return None


def _validate_bulk_table(df: pd.DataFrame) -> pd.DataFrame:
    missing = {"gene", "log2FoldChange"} - set(df.columns)
    if missing:
        raise ValueError(f"Bulk DE table missing required columns: {', '.join(sorted(missing))}.")
    if "padj" not in df.columns and "pvalue" not in df.columns:
        raise ValueError("Bulk DE table requires either a 'padj' or 'pvalue' column.")
    out = _maybe_numeric(df, ["log2FoldChange", "padj", "pvalue", "baseMean"])
    out["gene"] = out["gene"].astype(str)
    out = out.dropna(subset=["log2FoldChange"]).reset_index(drop=True)
    if out.empty:
        raise ValueError("Bulk DE table has no finite log2FoldChange rows after validation.")
    return out


def _validate_scrna_contrast_table(df: pd.DataFrame) -> pd.DataFrame:
    missing = {"names", "scores"} - set(df.columns)
    if missing:
        raise ValueError(
            f"scRNA contrast table missing required columns: {', '.join(sorted(missing))}."
        )
    out = _maybe_numeric(df, ["scores", "logfoldchanges", "pvals_adj", "pvals"])
    out["names"] = out["names"].astype(str)
    out = out.dropna(subset=["scores"]).reset_index(drop=True)
    if out.empty:
        raise ValueError("scRNA contrast table has no finite score rows after validation.")
    return out


def _validate_scrna_within_cluster_contrast_table(df: pd.DataFrame) -> pd.DataFrame:
    missing = {"cluster", "comparison_id", "group1", "group2", "names", "scores"} - set(df.columns)
    if missing:
        raise ValueError(
            "scRNA within-cluster contrast table missing required columns: "
            f"{', '.join(sorted(missing))}."
        )
    out = _maybe_numeric(df, ["scores", "logfoldchanges", "pvals_adj", "pvals"])
    for col in ("cluster", "comparison_id", "group1", "group2", "groupby", "scope"):
        if col in out.columns:
            out[col] = out[col].astype(str)
    out["names"] = out["names"].astype(str)
    out = out.dropna(subset=["scores"]).reset_index(drop=True)
    if out.empty:
        raise ValueError("scRNA within-cluster contrast table has no finite score rows after validation.")
    return out


def _validate_scrna_markers_table(df: pd.DataFrame) -> pd.DataFrame:
    missing = {"cluster", "names", "scores"} - set(df.columns)
    if missing:
        raise ValueError(
            f"scRNA markers table missing required columns: {', '.join(sorted(missing))}."
        )
    out = _maybe_numeric(df, ["scores", "logfoldchanges", "pvals_adj", "pvals"])
    out["cluster"] = out["cluster"].astype(str)
    out["names"] = out["names"].astype(str)
    out = out.dropna(subset=["scores"]).reset_index(drop=True)
    if out.empty:
        raise ValueError("scRNA markers table has no finite score rows after validation.")
    return out


def _examples_dir() -> Path:
    # Script lives in <skill_dir>/scripts/; demo fixtures live in
    # <skill_dir>/examples/ (this repo's skill layout convention).
    return Path(__file__).resolve().parent.parent / "examples"


def _load_demo_bundle(mode: str) -> InputBundle:
    examples = _examples_dir()
    if mode == "scrna":
        contrast_path = examples / "demo_scrna_contrast.csv"
        markers_path = examples / "demo_scrna_markers.csv"
        return InputBundle(
            mode="scrna",
            source_kind="demo-scrna",
            input_path=contrast_path,
            input_files=[contrast_path, markers_path],
            scrna_contrast_df=_validate_scrna_contrast_table(_read_table(contrast_path)),
            scrna_markers_df=_validate_scrna_markers_table(_read_table(markers_path)),
            notes=["Using bundled scRNA demo tables for contrast and marker visualisation."],
        )

    bulk_path = examples / "demo_bulk_de_results.csv"
    return InputBundle(
        mode="bulk",
        source_kind="demo-bulk",
        input_path=bulk_path,
        input_files=[
            bulk_path,
            examples / "demo_bulk_counts.csv",
            examples / "demo_bulk_metadata.csv",
        ],
        bulk_df=_validate_bulk_table(_read_table(bulk_path)),
        notes=["Using bundled bulk RNA differential expression demo table."],
    )


def _load_upstream_result_json(path: Path) -> dict[str, Any]:
    result_path = path / "result.json"
    if not result_path.exists():
        return {}
    try:
        return json.loads(result_path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def resolve_input_bundle(args: argparse.Namespace) -> InputBundle:
    requested_mode = args.mode
    if args.demo:
        demo_mode = "scrna" if requested_mode == "scrna" else "bulk"
        return _load_demo_bundle(demo_mode)

    if not args.input:
        raise ValueError("diffviz requires --input unless --demo is used.")

    input_path = Path(args.input).expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input path not found: {input_path}")

    if input_path.is_dir():
        tables_dir = input_path / "tables"
        bulk_path = tables_dir / "de_results.csv"
        contrast_path = tables_dir / "contrastive_markers_full.csv"
        within_cluster_path = tables_dir / "within_cluster_contrastive_markers_full.csv"
        markers_path = tables_dir / "markers_top.csv"
        upstream_result = _load_upstream_result_json(input_path)

        has_bulk = bulk_path.exists()
        has_contrast = contrast_path.exists()
        has_within_cluster = within_cluster_path.exists()
        has_markers = markers_path.exists()

        if requested_mode == "bulk":
            if not has_bulk:
                raise ValueError(
                    "Bulk mode requested, but input directory does not contain tables/de_results.csv."
                )
            return InputBundle(
                mode="bulk",
                source_kind="rnaseq-output-dir",
                input_path=input_path,
                input_files=[bulk_path],
                bulk_df=_validate_bulk_table(_read_table(bulk_path)),
                upstream_result=upstream_result,
                notes=["Detected rnaseq-de output directory."],
            )

        if requested_mode == "scrna":
            if not has_contrast and not has_within_cluster and not has_markers:
                raise ValueError(
                    "scRNA mode requested, but input directory does not contain "
                    "tables/contrastive_markers_full.csv, "
                    "tables/within_cluster_contrastive_markers_full.csv, or tables/markers_top.csv."
                )
            return InputBundle(
                mode="scrna",
                source_kind="scrna-output-dir",
                input_path=input_path,
                input_files=[path for path in [contrast_path, within_cluster_path, markers_path] if path.exists()],
                scrna_contrast_df=_validate_scrna_contrast_table(_read_table(contrast_path))
                if has_contrast
                else None,
                scrna_within_cluster_contrast_df=_validate_scrna_within_cluster_contrast_table(
                    _read_table(within_cluster_path)
                )
                if has_within_cluster
                else None,
                scrna_markers_df=_validate_scrna_markers_table(_read_table(markers_path))
                if has_markers
                else None,
                upstream_result=upstream_result,
                notes=["Detected scrna-orchestrator output directory."],
            )

        if has_bulk and (has_contrast or has_within_cluster or has_markers):
            raise ValueError(
                "Input directory contains both bulk and scRNA result tables. "
                "Please re-run with --mode bulk or --mode scrna."
            )
        if has_bulk:
            return InputBundle(
                mode="bulk",
                source_kind="rnaseq-output-dir",
                input_path=input_path,
                input_files=[bulk_path],
                bulk_df=_validate_bulk_table(_read_table(bulk_path)),
                upstream_result=upstream_result,
                notes=["Detected rnaseq-de output directory."],
            )
        if has_contrast or has_within_cluster or has_markers:
            return InputBundle(
                mode="scrna",
                source_kind="scrna-output-dir",
                input_path=input_path,
                input_files=[path for path in [contrast_path, within_cluster_path, markers_path] if path.exists()],
                scrna_contrast_df=_validate_scrna_contrast_table(_read_table(contrast_path))
                if has_contrast
                else None,
                scrna_within_cluster_contrast_df=_validate_scrna_within_cluster_contrast_table(
                    _read_table(within_cluster_path)
                )
                if has_within_cluster
                else None,
                scrna_markers_df=_validate_scrna_markers_table(_read_table(markers_path))
                if has_markers
                else None,
                upstream_result=upstream_result,
                notes=["Detected scrna-orchestrator output directory."],
            )
        raise ValueError(
            "Could not detect supported upstream outputs in the directory. Expected one of "
            "tables/de_results.csv, tables/contrastive_markers_full.csv, "
            "tables/within_cluster_contrastive_markers_full.csv, or tables/markers_top.csv."
        )

    table_df = _read_table(input_path)
    detected_kind = _detect_table_kind(table_df)
    if requested_mode == "bulk":
        if detected_kind != "bulk":
            raise ValueError(
                "Bulk mode requested, but the input table is not a bulk DE result. "
                "Expected columns: gene, log2FoldChange, and padj or pvalue."
            )
        return InputBundle(
            mode="bulk",
            source_kind="bulk-table",
            input_path=input_path,
            input_files=[input_path],
            bulk_df=_validate_bulk_table(table_df),
        )
    if requested_mode == "scrna":
        if detected_kind == "scrna_within_cluster_contrast":
            return InputBundle(
                mode="scrna",
                source_kind="scrna-within-cluster-contrast-table",
                input_path=input_path,
                input_files=[input_path],
                scrna_within_cluster_contrast_df=_validate_scrna_within_cluster_contrast_table(table_df),
            )
        if detected_kind == "scrna_contrast":
            return InputBundle(
                mode="scrna",
                source_kind="scrna-contrast-table",
                input_path=input_path,
                input_files=[input_path],
                scrna_contrast_df=_validate_scrna_contrast_table(table_df),
            )
        if detected_kind == "scrna_markers":
            return InputBundle(
                mode="scrna",
                source_kind="scrna-markers-table",
                input_path=input_path,
                input_files=[input_path],
                scrna_markers_df=_validate_scrna_markers_table(table_df),
            )
        raise ValueError(
            "scRNA mode requested, but the input table is not a supported marker/contrast table. "
            "Expected columns: names + scores, with optional cluster for marker tables."
        )

    if detected_kind == "bulk":
        return InputBundle(
            mode="bulk",
            source_kind="bulk-table",
            input_path=input_path,
            input_files=[input_path],
            bulk_df=_validate_bulk_table(table_df),
        )
    if detected_kind == "scrna_within_cluster_contrast":
        return InputBundle(
            mode="scrna",
            source_kind="scrna-within-cluster-contrast-table",
            input_path=input_path,
            input_files=[input_path],
            scrna_within_cluster_contrast_df=_validate_scrna_within_cluster_contrast_table(table_df),
        )
    if detected_kind == "scrna_contrast":
        return InputBundle(
            mode="scrna",
            source_kind="scrna-contrast-table",
            input_path=input_path,
            input_files=[input_path],
            scrna_contrast_df=_validate_scrna_contrast_table(table_df),
        )
    if detected_kind == "scrna_markers":
        return InputBundle(
            mode="scrna",
            source_kind="scrna-markers-table",
            input_path=input_path,
            input_files=[input_path],
            scrna_markers_df=_validate_scrna_markers_table(table_df),
        )
    raise ValueError(
        "Could not auto-detect input table type. Supported tables are: "
        "bulk DE (gene/log2FoldChange/padj|pvalue), scRNA contrast (names/scores), "
        "or scRNA markers (cluster/names/scores)."
    )


def _pvalue_column(df: pd.DataFrame) -> str:
    if "padj" in df.columns:
        return "padj"
    if "pvalue" in df.columns:
        return "pvalue"
    if "pvals_adj" in df.columns:
        return "pvals_adj"
    if "pvals" in df.columns:
        return "pvals"
    return ""


def _sanitize_pvalues(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").replace([np.inf, -np.inf], np.nan)
    values = values.fillna(1.0).clip(lower=1e-300, upper=1.0)
    return values


def _style_axis(ax: Any, *, grid_axis: str = "both") -> None:
    ax.set_facecolor(PLOT_COLORS["panel"])
    if grid_axis != "none":
        ax.grid(axis=grid_axis, color=PLOT_COLORS["grid"], linestyle="-", linewidth=0.8, alpha=0.7)
    else:
        ax.grid(False)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(PLOT_COLORS["spine"])
        ax.spines[side].set_linewidth(1.0)


def _display_cap(
    series: pd.Series,
    *,
    quantile: float = 0.995,
    minimum: float = 4.0,
    fallback: float = 8.0,
) -> float:
    values = pd.to_numeric(series, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna().abs()
    if values.empty:
        return fallback
    cap = float(values.quantile(quantile))
    if not math.isfinite(cap) or cap <= 0:
        cap = fallback
    return max(minimum, cap)


def _save_current_figure(path: Path, dpi: int = 180) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close("all")
    return path


def _label_points(ax: Any, data: pd.DataFrame, x_col: str, y_col: str, label_col: str, n: int) -> None:
    if n <= 0 or data.empty:
        return
    top = data.head(n)
    for _, row in top.iterrows():
        x = row.get(x_col)
        y = row.get(y_col)
        label = str(row.get(label_col, ""))
        if pd.notna(x) and pd.notna(y):
            ax.text(float(x), float(y), label, fontsize=8)


def _bulk_display_frame(
    df: pd.DataFrame,
    *,
    min_basemean: float,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    meta: dict[str, Any] = {
        "rows_before": int(len(df)),
        "rows_after": int(len(df)),
        "applied": False,
        "min_basemean": float(min_basemean),
    }
    if "baseMean" not in df.columns or min_basemean <= 0:
        return df.copy(), meta

    base_mean = pd.to_numeric(df["baseMean"], errors="coerce")
    filtered = df.loc[base_mean > float(min_basemean)].copy()
    if filtered.empty:
        meta["note"] = (
            f"Requested baseMean > {min_basemean:g} display filter produced no rows; "
            "falling back to the unfiltered bulk table."
        )
        return df.copy(), meta

    meta["rows_after"] = int(len(filtered))
    meta["applied"] = True
    meta["rows_removed"] = int(len(df) - len(filtered))
    return filtered, meta


def plot_bulk_volcano(
    df: pd.DataFrame,
    outpath: Path,
    *,
    padj_threshold: float,
    lfc_threshold: float,
    label_top: int,
) -> Path:
    p_col = _pvalue_column(df)
    plot_df = df.copy()
    plot_df["log2FoldChange"] = pd.to_numeric(plot_df["log2FoldChange"], errors="coerce")
    plot_df["_neglog10"] = -np.log10(_sanitize_pvalues(plot_df[p_col]))
    plot_df["_significant"] = (
        (_sanitize_pvalues(plot_df[p_col]) <= padj_threshold)
        & (plot_df["log2FoldChange"].abs() >= lfc_threshold)
    )
    plot_df["_direction"] = np.where(
        plot_df["_significant"] & (plot_df["log2FoldChange"] >= lfc_threshold),
        "up",
        np.where(
            plot_df["_significant"] & (plot_df["log2FoldChange"] <= -lfc_threshold),
            "down",
            "neutral",
        ),
    )
    display_cap = _display_cap(
        plot_df["log2FoldChange"],
        minimum=max(4.0, lfc_threshold * 3.0),
    )
    plot_df["_display_lfc"] = plot_df["log2FoldChange"].clip(-display_cap, display_cap)
    plot_df["_label_score"] = plot_df["_neglog10"] * plot_df["log2FoldChange"].abs().replace(0, 0.05)
    label_df = plot_df.sort_values(["_significant", "_label_score"], ascending=[False, False])

    fig, ax = plt.subplots(figsize=(8.6, 6.2))
    ax.scatter(
        plot_df.loc[plot_df["_direction"] == "neutral", "_display_lfc"],
        plot_df.loc[plot_df["_direction"] == "neutral", "_neglog10"],
        s=18,
        color=PLOT_COLORS["neutral"],
        alpha=0.48,
        label="Not significant",
        rasterized=True,
    )
    ax.scatter(
        plot_df.loc[plot_df["_direction"] == "down", "_display_lfc"],
        plot_df.loc[plot_df["_direction"] == "down", "_neglog10"],
        s=24,
        color=PLOT_COLORS["down"],
        alpha=0.72,
        label="Downregulated",
        rasterized=True,
    )
    ax.scatter(
        plot_df.loc[plot_df["_direction"] == "up", "_display_lfc"],
        plot_df.loc[plot_df["_direction"] == "up", "_neglog10"],
        s=24,
        color=PLOT_COLORS["up"],
        alpha=0.74,
        label="Upregulated",
        rasterized=True,
    )
    ax.axvline(lfc_threshold, color=PLOT_COLORS["muted"], linestyle="--", linewidth=1.1)
    ax.axvline(-lfc_threshold, color=PLOT_COLORS["muted"], linestyle="--", linewidth=1.1)
    ax.axhline(-math.log10(max(padj_threshold, 1e-300)), color=PLOT_COLORS["muted"], linestyle="--", linewidth=1.1)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title("Bulk RNA Differential Expression Volcano")
    ax.set_xlim(-(display_cap * 1.08), display_cap * 1.08)
    _style_axis(ax, grid_axis="both")
    _label_points(ax, label_df, "_display_lfc", "_neglog10", "gene", label_top)
    ax.legend(loc="upper left", ncol=1)
    return _save_current_figure(outpath)


def plot_top_gene_bar(df: pd.DataFrame, outpath: Path, *, top_genes: int, title: str) -> Path:
    plot_df = df.head(top_genes).copy()
    if plot_df.empty:
        raise ValueError("No genes available for top-gene bar plot.")
    colors = [PLOT_COLORS["up"] if value >= 0 else PLOT_COLORS["down"] for value in plot_df.iloc[::-1]["log2FoldChange"]]
    fig, ax = plt.subplots(figsize=(8.2, max(4.4, 0.30 * len(plot_df) + 1.2)))
    ax.barh(
        plot_df.iloc[::-1]["gene"],
        plot_df.iloc[::-1]["log2FoldChange"],
        color=colors,
        edgecolor="#FFFFFF",
        linewidth=0.7,
    )
    ax.axvline(0, color=PLOT_COLORS["muted"], linewidth=1.1)
    ax.set_xlabel("log2 fold change")
    ax.set_title(title)
    _style_axis(ax, grid_axis="x")
    return _save_current_figure(outpath)


def plot_ma(
    df: pd.DataFrame,
    outpath: Path,
    *,
    padj_threshold: float,
    lfc_threshold: float,
) -> Path:
    if "baseMean" not in df.columns:
        raise ValueError("MA plot requires a baseMean column.")
    p_col = _pvalue_column(df)
    plot_df = df.dropna(subset=["baseMean"]).copy()
    if plot_df.empty:
        raise ValueError("MA plot requires at least one finite baseMean row.")
    significant = (
        (_sanitize_pvalues(plot_df[p_col]) <= padj_threshold)
        & (plot_df["log2FoldChange"].abs() >= lfc_threshold)
    )
    fig, ax = plt.subplots(figsize=(8.6, 6.2))
    ax.scatter(
        np.log10(plot_df.loc[~significant, "baseMean"] + 1.0),
        plot_df.loc[~significant, "log2FoldChange"],
        s=18,
        color=PLOT_COLORS["neutral"],
        alpha=0.45,
        rasterized=True,
    )
    sig_df = plot_df.loc[significant].copy()
    sig_colors = np.where(sig_df["log2FoldChange"] >= 0, PLOT_COLORS["up"], PLOT_COLORS["down"])
    ax.scatter(
        np.log10(sig_df["baseMean"] + 1.0),
        sig_df["log2FoldChange"],
        s=24,
        color=sig_colors,
        alpha=0.68,
        rasterized=True,
    )
    ax.axhline(lfc_threshold, linestyle="--", linewidth=1.1, color=PLOT_COLORS["muted"])
    ax.axhline(-lfc_threshold, linestyle="--", linewidth=1.1, color=PLOT_COLORS["muted"])
    ax.set_xlabel("log10(baseMean + 1)")
    ax.set_ylabel("log2 fold change")
    ax.set_title("MA Plot")
    _style_axis(ax, grid_axis="both")
    return _save_current_figure(outpath)


def _load_counts_matrix(path: Path) -> pd.DataFrame:
    counts = _read_table(path)
    if counts.shape[1] < 3:
        raise ValueError("Counts matrix must contain one gene column and at least two samples.")
    gene_col = counts.columns[0]
    counts = counts.set_index(gene_col)
    counts = counts.apply(pd.to_numeric, errors="coerce")
    if counts.isna().any().any():
        raise ValueError("Counts matrix contains non-numeric entries.")
    counts.index = counts.index.astype(str)
    return counts


def _load_bulk_metadata(path: Path) -> pd.DataFrame:
    metadata = _read_table(path)
    if "sample_id" not in metadata.columns:
        raise ValueError("Metadata must include a sample_id column.")
    metadata = metadata.copy()
    metadata["sample_id"] = metadata["sample_id"].astype(str)
    metadata = metadata.set_index("sample_id")
    return metadata


def _infer_bulk_group_column(metadata: pd.DataFrame, contrast: str) -> str:
    if contrast:
        factor = contrast.split(",")[0].strip()
        if factor and factor in metadata.columns:
            return factor
    for col in metadata.columns:
        series = metadata[col]
        if pd.api.types.is_numeric_dtype(series):
            continue
        nunique = series.astype(str).nunique(dropna=True)
        if 1 < nunique <= 20:
            return col
    return ""


def _prepare_bulk_matrix_inputs(
    counts_path: Path,
    metadata_path: Path,
    contrast: str,
) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    counts = _load_counts_matrix(counts_path)
    metadata = _load_bulk_metadata(metadata_path)
    missing_samples = [sample for sample in counts.columns if sample not in metadata.index]
    if missing_samples:
        raise ValueError(f"Metadata missing count-matrix samples: {', '.join(missing_samples[:5])}")
    metadata = metadata.loc[counts.columns].copy()
    group_col = _infer_bulk_group_column(metadata, contrast)
    return counts, metadata, group_col


def _evenly_spaced_items(items: list[str], take: int) -> list[str]:
    if take <= 0 or not items:
        return []
    if take >= len(items):
        return list(items)
    positions = np.linspace(0, len(items) - 1, num=take)
    return [items[int(round(pos))] for pos in positions]


def _select_heatmap_samples(
    order: list[str],
    metadata: pd.DataFrame,
    group_col: str,
    *,
    max_samples: int = 48,
) -> tuple[list[str], dict[str, Any]]:
    if len(order) <= max_samples:
        return list(order), {
            "original_sample_count": len(order),
            "selected_sample_count": len(order),
            "downsampled": False,
        }

    if not group_col:
        selected = _evenly_spaced_items(order, max_samples)
        return selected, {
            "original_sample_count": len(order),
            "selected_sample_count": len(selected),
            "downsampled": True,
        }

    groups = metadata.loc[order, group_col].fillna("Unknown").astype(str)
    group_order = list(dict.fromkeys(groups.tolist()))
    per_group_target = max(1, max_samples // max(1, len(group_order)))
    selected: list[str] = []
    leftovers: dict[str, list[str]] = {}

    for group in group_order:
        members = [sample for sample in order if groups.loc[sample] == group]
        chosen = _evenly_spaced_items(members, min(len(members), per_group_target))
        selected.extend(chosen)
        chosen_set = set(chosen)
        leftovers[group] = [sample for sample in members if sample not in chosen_set]

    remaining_slots = max_samples - len(selected)
    while remaining_slots > 0:
        available_groups = [group for group in group_order if leftovers[group]]
        if not available_groups:
            break
        available_groups.sort(key=lambda group: len(leftovers[group]), reverse=True)
        for group in available_groups:
            if remaining_slots <= 0:
                break
            selected.append(leftovers[group].pop(0))
            remaining_slots -= 1

    selected_set = set(selected)
    selected = [sample for sample in order if sample in selected_set]
    return selected, {
        "original_sample_count": len(order),
        "selected_sample_count": len(selected),
        "downsampled": len(selected) < len(order),
        "group_counts": groups.value_counts().sort_index().to_dict(),
    }


def _bulk_log_cpm_matrix(
    counts: pd.DataFrame,
    genes: list[str],
    sample_order: list[str],
) -> pd.DataFrame:
    available_genes = [gene for gene in genes if gene in counts.index]
    if not available_genes:
        raise ValueError("None of the selected top genes are present in the counts matrix.")
    matrix_counts = counts.loc[available_genes, sample_order]
    lib_sizes = matrix_counts.sum(axis=0)
    return np.log1p(matrix_counts.div(lib_sizes, axis=1) * 1_000_000.0)


def plot_bulk_heatmap(
    genes: list[str],
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    group_col: str,
    outpath: Path,
    *,
    max_samples: int = 48,
) -> tuple[Path, dict[str, Any]]:
    if group_col:
        order = metadata.sort_values(group_col).index.tolist()
    else:
        order = counts.columns.astype(str).tolist()

    selected_samples, selection_meta = _select_heatmap_samples(
        order,
        metadata,
        group_col,
        max_samples=max_samples,
    )
    log_cpm = _bulk_log_cpm_matrix(counts, genes, selected_samples)
    available_genes = log_cpm.index.astype(str).tolist()
    matrix = log_cpm.to_numpy(dtype=float)
    matrix = matrix - matrix.mean(axis=1, keepdims=True)
    std = matrix.std(axis=1, keepdims=True)
    std[std == 0] = 1.0
    z_matrix = matrix / std

    fig, ax = plt.subplots(
        figsize=(
            max(9, min(18, 0.28 * len(selected_samples) + 4)),
            max(6, min(14, 0.42 * len(available_genes) + 3)),
        )
    )
    im = ax.imshow(z_matrix, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    tick_stride = max(1, int(math.ceil(len(selected_samples) / 24)))
    tick_positions = np.arange(0, len(selected_samples), tick_stride)
    tick_labels: list[str] = []
    for pos in tick_positions:
        sample = selected_samples[int(pos)]
        short_sample = sample[:8] if len(sample) > 10 else sample
        if group_col and len(selected_samples) <= 24:
            tick_labels.append(f"{short_sample}\n{metadata.loc[sample, group_col]}")
        else:
            tick_labels.append(short_sample)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(available_genes)))
    ax.set_yticklabels(available_genes, fontsize=8)
    if selection_meta["downsampled"]:
        ax.set_title(
            "Top Gene Heatmap "
            f"({selection_meta['selected_sample_count']} representative samples, z-scored logCPM)"
        )
    else:
        ax.set_title("Top Gene Heatmap (z-scored logCPM)")
    _style_axis(ax, grid_axis="none")
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.set_ylabel("z-scored logCPM", rotation=90, va="bottom")
    return _save_current_figure(outpath), selection_meta


def plot_bulk_mean_comparison(
    genes: list[str],
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    group_col: str,
    outpath: Path,
) -> Path:
    if not group_col:
        raise ValueError("Mean comparison requires a grouping column in metadata.")
    groups = metadata[group_col].fillna("Unknown").astype(str)
    group_order = sorted(groups.unique().tolist())
    if len(group_order) != 2:
        raise ValueError("Mean comparison currently supports exactly 2 groups.")

    log_cpm = _bulk_log_cpm_matrix(counts, genes, counts.columns.astype(str).tolist())
    rows = []
    for gene in log_cpm.index.astype(str):
        row = {"gene": gene}
        for group in group_order:
            sample_ids = groups.index[groups == group].tolist()
            row[group] = float(log_cpm.loc[gene, sample_ids].mean())
        rows.append(row)
    compare_df = pd.DataFrame(rows)
    if compare_df.empty:
        raise ValueError("No genes available for mean comparison.")

    left_group, right_group = group_order
    plot_df = compare_df.head(min(15, len(compare_df))).copy().iloc[::-1].reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(9.2, max(5.2, 0.36 * len(plot_df) + 1.7)))
    for idx, row in plot_df.iterrows():
        ax.plot([row[left_group], row[right_group]], [idx, idx], color=PLOT_COLORS["neutral"], linewidth=2.0, alpha=0.9)
        ax.scatter(row[left_group], idx, color=PLOT_COLORS["down"], s=56, label=left_group if idx == 0 else None, zorder=3)
        ax.scatter(row[right_group], idx, color=PLOT_COLORS["up"], s=56, label=right_group if idx == 0 else None, zorder=3)
    ax.set_yticks(np.arange(len(plot_df)))
    ax.set_yticklabels(plot_df["gene"], fontsize=8)
    ax.set_xlabel("mean logCPM")
    ax.set_title(f"Mean Expression Comparison: {right_group} vs {left_group}")
    ax.legend(loc="lower right")
    _style_axis(ax, grid_axis="x")
    return _save_current_figure(outpath)


def plot_bulk_metric_heatmap(
    df: pd.DataFrame,
    outpath: Path,
    *,
    top_genes: int,
) -> Path:
    metric_df = df.head(min(top_genes, len(df))).copy()
    if metric_df.empty:
        raise ValueError("No genes available for gene metrics panel.")
    p_col = _pvalue_column(metric_df)
    metric_df["gene"] = metric_df["gene"].astype(str)
    metric_df["log2FC"] = pd.to_numeric(metric_df["log2FoldChange"], errors="coerce")
    metric_df["neglog10_p"] = -np.log10(_sanitize_pvalues(metric_df[p_col]))
    if "baseMean" in metric_df.columns:
        base_mean = pd.to_numeric(metric_df["baseMean"], errors="coerce").clip(lower=0)
        metric_df["log10_baseMean"] = np.log10(base_mean + 1.0)
    else:
        metric_df["log10_baseMean"] = np.nan

    plot_df = metric_df.iloc[::-1].reset_index(drop=True)
    genes = plot_df["gene"].tolist()
    y_pos = np.arange(len(plot_df))

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(13, max(6, 0.34 * len(plot_df) + 1.6)),
        sharey=True,
        gridspec_kw={"wspace": 0.08, "width_ratios": [1.15, 1.0, 1.0]},
    )

    fc_colors = [PLOT_COLORS["up"] if value >= 0 else PLOT_COLORS["down"] for value in plot_df["log2FC"]]
    axes[0].barh(y_pos, plot_df["log2FC"], color=fc_colors, alpha=0.9)
    axes[0].axvline(0, color=PLOT_COLORS["muted"], linewidth=1.1)
    axes[0].set_xlabel("log2FC")
    axes[0].set_title("Effect Size")
    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels(genes, fontsize=8)

    axes[1].barh(y_pos, plot_df["neglog10_p"], color=PLOT_COLORS["teal"], alpha=0.9)
    axes[1].set_xlabel("-log10(padj)")
    axes[1].set_title("Significance")
    axes[1].tick_params(axis="y", left=False, labelleft=False)

    axes[2].barh(y_pos, plot_df["log10_baseMean"].fillna(0.0), color=PLOT_COLORS["purple"], alpha=0.88)
    axes[2].set_xlabel("log10(baseMean+1)")
    axes[2].set_title("Mean Expression")
    axes[2].tick_params(axis="y", left=False, labelleft=False)

    for ax in axes:
        _style_axis(ax, grid_axis="x")

    fig.suptitle("Gene Metrics Panel", fontsize=16, y=0.995)
    return _save_current_figure(outpath)


def plot_bulk_effect_bubble(
    df: pd.DataFrame,
    outpath: Path,
    *,
    padj_threshold: float,
    lfc_threshold: float,
    label_top: int,
) -> Path:
    p_col = _pvalue_column(df)
    plot_df = df.copy()
    plot_df["_x"] = pd.to_numeric(plot_df["log2FoldChange"], errors="coerce")
    plot_df["_y"] = -np.log10(_sanitize_pvalues(plot_df[p_col]))
    if "baseMean" in plot_df.columns:
        size_metric = np.log10(pd.to_numeric(plot_df["baseMean"], errors="coerce").clip(lower=0) + 1.0)
    else:
        size_metric = plot_df["_x"].abs().clip(lower=0.2)
    plot_df["_size"] = size_metric.fillna(float(size_metric.median())).clip(lower=0.2) * 35
    plot_df["_significant"] = (
        (_sanitize_pvalues(plot_df[p_col]) <= padj_threshold)
        & (plot_df["_x"].abs() >= lfc_threshold)
    )
    plot_df["_label_score"] = plot_df["_y"] * plot_df["_x"].abs().replace(0, 0.05)
    label_df = plot_df.sort_values(["_significant", "_label_score"], ascending=[False, False])

    display_cap = _display_cap(plot_df["_x"], minimum=max(4.0, lfc_threshold * 3.0))
    plot_df["_display_x"] = plot_df["_x"].clip(-display_cap, display_cap)
    fig, ax = plt.subplots(figsize=(8.6, 6.2))
    colors = np.where(plot_df["_x"] >= 0, PLOT_COLORS["up"], PLOT_COLORS["down"])
    alphas = np.where(plot_df["_significant"], 0.72, 0.28)
    ax.scatter(
        plot_df["_display_x"],
        plot_df["_y"],
        s=plot_df["_size"],
        c=colors,
        alpha=alphas,
        edgecolors="white",
        linewidths=0.35,
        rasterized=True,
    )
    ax.axvline(lfc_threshold, color=PLOT_COLORS["muted"], linestyle="--", linewidth=1.1)
    ax.axvline(-lfc_threshold, color=PLOT_COLORS["muted"], linestyle="--", linewidth=1.1)
    ax.axhline(-math.log10(max(padj_threshold, 1e-300)), color=PLOT_COLORS["muted"], linestyle="--", linewidth=1.1)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title("Effect Bubble Plot")
    ax.set_xlim(-(display_cap * 1.08), display_cap * 1.08)
    _style_axis(ax, grid_axis="both")
    _label_points(ax, label_df.assign(_display_x=label_df["_x"].clip(-display_cap, display_cap)), "_display_x", "_y", "gene", label_top)
    return _save_current_figure(outpath)


def _sorted_significant_bulk(
    df: pd.DataFrame,
    *,
    padj_threshold: float,
    lfc_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    p_col = _pvalue_column(df)
    sig = df.loc[
        (_sanitize_pvalues(df[p_col]) <= padj_threshold)
        & (df["log2FoldChange"].abs() >= lfc_threshold)
    ].copy()
    sig = sig.sort_values([p_col, "log2FoldChange"], ascending=[True, False])
    top = sig.copy()
    if top.empty:
        top = df.copy()
        top["_rank_abs_lfc"] = top["log2FoldChange"].abs()
        top = top.sort_values([p_col, "_rank_abs_lfc"], ascending=[True, False])
        top = top.drop(columns=["_rank_abs_lfc"])
    return top, sig


def _sc_to_array(matrix: Any) -> np.ndarray:
    if hasattr(matrix, "toarray"):
        return np.asarray(matrix.toarray(), dtype=float)
    if hasattr(matrix, "A"):
        return np.asarray(matrix.A, dtype=float)
    return np.asarray(matrix, dtype=float)


def _load_adata(path: Path) -> Any:
    try:
        import anndata as ad
    except Exception as exc:  # pragma: no cover - import failure path
        raise RuntimeError("anndata is required for --adata-enhanced scRNA plots.") from exc
    return ad.read_h5ad(path)


def _infer_scrna_groupby(adata: Any) -> str:
    preferred = ["leiden", "cluster", "clusters", "cell_type", "annotation"]
    for col in preferred:
        if col in adata.obs.columns:
            return col
    for col in adata.obs.columns:
        series = adata.obs[col]
        if pd.api.types.is_categorical_dtype(series) or series.dtype == object:
            nunique = pd.Series(series.astype(str)).nunique(dropna=True)
            if 1 < nunique <= 30:
                return str(col)
    return ""


def _top_scrna_genes(
    contrast_df: pd.DataFrame | None,
    within_cluster_df: pd.DataFrame | None,
    markers_df: pd.DataFrame | None,
    top_genes: int,
) -> list[str]:
    genes: list[str] = []
    if contrast_df is not None and not contrast_df.empty:
        genes.extend(
            contrast_df.sort_values(
                ["pvals_adj", "scores"] if "pvals_adj" in contrast_df.columns else ["scores"],
                ascending=[True, False] if "pvals_adj" in contrast_df.columns else [False],
            )["names"].dropna().astype(str).tolist()
        )
    if within_cluster_df is not None and not within_cluster_df.empty:
        sort_cols = ["cluster", "comparison_id", "pvals_adj", "scores"] if "pvals_adj" in within_cluster_df.columns else ["cluster", "comparison_id", "scores"]
        ascending = [True, True, True, False] if "pvals_adj" in within_cluster_df.columns else [True, True, False]
        per_panel = (
            within_cluster_df.sort_values(sort_cols, ascending=ascending)
            .groupby(["cluster", "comparison_id"], as_index=False, group_keys=False)
            .head(min(max(top_genes // 2, 2), top_genes))
        )
        genes.extend(per_panel["names"].dropna().astype(str).tolist())
    if markers_df is not None and not markers_df.empty:
        per_cluster = (
            markers_df.sort_values(["cluster", "scores"], ascending=[True, False])
            .groupby("cluster", as_index=False, group_keys=False)
            .head(min(max(top_genes // 3, 2), top_genes))
        )
        genes.extend(per_cluster["names"].dropna().astype(str).tolist())
    ordered = list(dict.fromkeys(genes))
    return ordered[: max(top_genes, 4)]


def plot_scrna_contrast_volcano(
    df: pd.DataFrame,
    outpath: Path,
    *,
    top_genes: int,
    label_top: int,
    padj_threshold: float,
    lfc_threshold: float,
) -> tuple[Path, pd.DataFrame]:
    plot_df = df.copy()
    x_col = "logfoldchanges" if "logfoldchanges" in plot_df.columns and plot_df["logfoldchanges"].notna().any() else "scores"
    p_col = _pvalue_column(plot_df)
    if p_col:
        plot_df["_y"] = -np.log10(_sanitize_pvalues(plot_df[p_col]))
        plot_df["_significant"] = (
            (_sanitize_pvalues(plot_df[p_col]) <= padj_threshold)
            & (pd.to_numeric(plot_df[x_col], errors="coerce").abs() >= lfc_threshold)
        )
    else:
        plot_df["_y"] = plot_df["scores"].abs().clip(lower=0.01)
        plot_df["_significant"] = plot_df["scores"].abs() >= plot_df["scores"].abs().quantile(0.9)
    plot_df["_x"] = pd.to_numeric(plot_df[x_col], errors="coerce")
    plot_df = plot_df.dropna(subset=["_x", "_y"]).copy()
    plot_df["_label_score"] = plot_df["_y"] * plot_df["_x"].abs().replace(0, 0.05)
    top_df = plot_df.sort_values(["_significant", "_label_score"], ascending=[False, False]).head(top_genes).copy()

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(
        plot_df.loc[~plot_df["_significant"], "_x"],
        plot_df.loc[~plot_df["_significant"], "_y"],
        s=18,
        color="#90a4ae",
        alpha=0.75,
        label="Background",
    )
    ax.scatter(
        plot_df.loc[plot_df["_significant"], "_x"],
        plot_df.loc[plot_df["_significant"], "_y"],
        s=24,
        color="#ef6c00",
        alpha=0.85,
        label="Highlighted",
    )
    if x_col == "logfoldchanges":
        ax.axvline(lfc_threshold, linestyle="--", linewidth=1, color="#455a64")
        ax.axvline(-lfc_threshold, linestyle="--", linewidth=1, color="#455a64")
        ax.set_xlabel("log fold change")
    else:
        ax.set_xlabel("marker score")
    if p_col:
        ax.axhline(-math.log10(max(padj_threshold, 1e-300)), linestyle="--", linewidth=1, color="#455a64")
        ax.set_ylabel("-log10 p-value")
    else:
        ax.set_ylabel("absolute score")
    ax.set_title("scRNA Contrast Volcano")
    _label_points(ax, top_df, "_x", "_y", "names", label_top)
    ax.legend(loc="upper right")
    return _save_current_figure(outpath), top_df


def plot_scrna_top_markers_bar(
    df: pd.DataFrame,
    outpath: Path,
    *,
    top_genes: int,
    title: str,
) -> tuple[Path, pd.DataFrame]:
    plot_df = df.copy()
    sort_cols = ["pvals_adj", "scores"] if "pvals_adj" in plot_df.columns else ["scores"]
    ascending = [True, False] if "pvals_adj" in plot_df.columns else [False]
    top_df = plot_df.sort_values(sort_cols, ascending=ascending).head(top_genes).copy()
    values = top_df["logfoldchanges"] if "logfoldchanges" in top_df.columns and top_df["logfoldchanges"].notna().any() else top_df["scores"]
    colors = ["#6a1b9a" if value >= 0 else "#1565c0" for value in values.iloc[::-1]]
    fig, ax = plt.subplots(figsize=(8, max(4, 0.28 * len(top_df) + 1)))
    ax.barh(top_df.iloc[::-1]["names"], values.iloc[::-1], color=colors)
    ax.axvline(0, color="#455a64", linewidth=1)
    ax.set_title(title)
    ax.set_xlabel("log fold change" if values.name == "logfoldchanges" else "marker score")
    return _save_current_figure(outpath), top_df


def _count_significant_scrna_rows(df: pd.DataFrame, *, padj_threshold: float) -> int:
    p_col = _pvalue_column(df)
    if p_col:
        return int(len(df.loc[_sanitize_pvalues(df[p_col]) <= padj_threshold]))
    return int(len(df.loc[df["scores"].abs() >= df["scores"].abs().quantile(0.9)]))


def plot_within_cluster_marker_panels(
    df: pd.DataFrame,
    outpath: Path,
    *,
    top_genes: int,
) -> tuple[Path, pd.DataFrame]:
    sort_cols = ["cluster", "comparison_id", "pvals_adj", "scores"] if "pvals_adj" in df.columns else ["cluster", "comparison_id", "scores"]
    ascending = [True, True, True, False] if "pvals_adj" in df.columns else [True, True, False]
    grouped = (
        df.sort_values(sort_cols, ascending=ascending)
        .groupby(["cluster", "comparison_id"], as_index=False, group_keys=False)
        .head(max(1, top_genes))
        .reset_index(drop=True)
    )
    panel_keys = list(grouped[["cluster", "comparison_id"]].drop_duplicates().itertuples(index=False, name=None))
    fig, axes = plt.subplots(
        len(panel_keys),
        1,
        figsize=(8.5, max(3.4, 2.8 * len(panel_keys))),
        squeeze=False,
    )
    for ax, (cluster, comparison_id) in zip(axes.flatten(), panel_keys):
        subset = grouped.loc[
            (grouped["cluster"].astype(str) == str(cluster))
            & (grouped["comparison_id"].astype(str) == str(comparison_id))
        ].copy()
        value_col = "logfoldchanges" if "logfoldchanges" in subset.columns and subset["logfoldchanges"].notna().any() else "scores"
        values = pd.to_numeric(subset[value_col], errors="coerce").fillna(0.0)
        colors = [PLOT_COLORS["up"] if value >= 0 else PLOT_COLORS["down"] for value in values.iloc[::-1]]
        ax.barh(subset.iloc[::-1]["names"], values.iloc[::-1], color=colors)
        label = f"Cluster {cluster}: {comparison_id.replace('__vs__', ' vs ')}"
        ax.set_title(label)
        ax.axvline(0, color=PLOT_COLORS["spine"], linewidth=1)
        ax.set_xlabel("log fold change" if value_col == "logfoldchanges" else "marker score")
    fig.tight_layout()
    return _save_current_figure(outpath), grouped


def plot_marker_rank_bars(
    df: pd.DataFrame,
    outpath: Path,
    *,
    per_cluster: int,
) -> tuple[Path, pd.DataFrame]:
    grouped = (
        df.sort_values(["cluster", "scores"], ascending=[True, False])
        .groupby("cluster", as_index=False, group_keys=False)
        .head(per_cluster)
        .reset_index(drop=True)
    )
    clusters = grouped["cluster"].dropna().astype(str).unique().tolist()
    n_panels = len(clusters)
    fig, axes = plt.subplots(
        n_panels,
        1,
        figsize=(8, max(3.2, 2.5 * n_panels)),
        squeeze=False,
    )
    for ax, cluster in zip(axes.flatten(), clusters):
        subset = grouped.loc[grouped["cluster"].astype(str) == cluster].copy()
        ax.barh(subset.iloc[::-1]["names"], subset.iloc[::-1]["scores"], color="#00897b")
        ax.set_title(f"Cluster {cluster}")
        ax.set_xlabel("marker score")
    return _save_current_figure(outpath), grouped


def plot_manual_dotplot(adata: Any, genes: list[str], groupby: str, outpath: Path) -> Path:
    groups = adata.obs[groupby].astype(str)
    group_order = sorted(groups.unique().tolist())
    gene_indices = [int(np.where(adata.var_names == gene)[0][0]) for gene in genes]
    dots = []
    for group_i, group in enumerate(group_order):
        mask = (groups == group).to_numpy()
        group_matrix = _sc_to_array(adata[mask, gene_indices].X)
        pct = (group_matrix > 0).mean(axis=0)
        mean_expr = group_matrix.mean(axis=0)
        for gene_i, gene in enumerate(genes):
            dots.append((gene_i, group_i, float(pct[gene_i]), float(mean_expr[gene_i]), gene, group))
    dot_df = pd.DataFrame(dots, columns=["gene_i", "group_i", "pct", "mean_expr", "gene", "group"])
    fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(genes) + 2), max(4, 0.4 * len(group_order) + 2)))
    sc = ax.scatter(
        dot_df["gene_i"],
        dot_df["group_i"],
        s=dot_df["pct"].clip(lower=0.02) * 700,
        c=dot_df["mean_expr"],
        cmap="viridis",
        edgecolor="black",
        linewidth=0.3,
    )
    ax.set_xticks(np.arange(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(group_order)))
    ax.set_yticklabels(group_order, fontsize=8)
    ax.set_title(f"Marker Dotplot by {groupby}")
    plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04, label="mean expression")
    return _save_current_figure(outpath)


def plot_manual_heatmap(adata: Any, genes: list[str], groupby: str, outpath: Path) -> Path:
    groups = adata.obs[groupby].astype(str)
    group_order = sorted(groups.unique().tolist())
    gene_indices = [int(np.where(adata.var_names == gene)[0][0]) for gene in genes]
    rows = []
    for group in group_order:
        mask = (groups == group).to_numpy()
        group_matrix = _sc_to_array(adata[mask, gene_indices].X)
        rows.append(group_matrix.mean(axis=0))
    matrix = np.vstack(rows)
    if matrix.size == 0:
        raise ValueError("No values available for marker heatmap.")
    fig, ax = plt.subplots(figsize=(max(6, 0.5 * len(genes) + 2), max(4, 0.4 * len(group_order) + 2)))
    im = ax.imshow(matrix, aspect="auto", cmap="magma")
    ax.set_xticks(np.arange(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(group_order)))
    ax.set_yticklabels(group_order, fontsize=8)
    ax.set_title(f"Marker Heatmap by {groupby}")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="mean expression")
    return _save_current_figure(outpath)


def _ensure_umap(adata: Any) -> Any:
    if "X_umap" in getattr(adata, "obsm", {}) and np.asarray(adata.obsm["X_umap"]).shape[1] >= 2:
        return adata
    try:
        import scanpy as sc
    except Exception as exc:  # pragma: no cover - optional dependency path
        raise RuntimeError("scanpy is required to compute UMAP when X_umap is absent.") from exc
    adata = adata.copy()
    if "X_pca" not in getattr(adata, "obsm", {}):
        sc.pp.pca(adata)
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.umap(adata, random_state=0)
    return adata


def plot_umap_feature_panel(adata: Any, genes: list[str], outpath: Path) -> Path:
    adata = _ensure_umap(adata)
    coords = np.asarray(adata.obsm["X_umap"], dtype=float)
    gene_indices = [int(np.where(adata.var_names == gene)[0][0]) for gene in genes]
    expr = _sc_to_array(adata[:, gene_indices].X)
    n_panels = len(genes)
    n_cols = min(2, n_panels)
    n_rows = int(math.ceil(n_panels / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.0 * n_cols, 4.0 * n_rows), squeeze=False)
    axes_flat = axes.flatten()
    for ax, gene, gene_i in zip(axes_flat, genes, range(n_panels)):
        sc = ax.scatter(coords[:, 0], coords[:, 1], c=expr[:, gene_i], cmap="viridis", s=12, linewidth=0)
        ax.set_title(gene)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    for ax in axes_flat[n_panels:]:
        ax.axis("off")
    fig.suptitle("UMAP Feature Panel", fontsize=13)
    return _save_current_figure(outpath)


def render_bulk_outputs(
    bundle: InputBundle,
    args: argparse.Namespace,
    output_dir: Path,
    artifacts: RunArtifacts,
) -> None:
    figures_dir = output_dir / "figures"
    tables_dir = output_dir / "tables"
    bulk_df = bundle.bulk_df
    assert bulk_df is not None

    artifacts.n_input_rows = int(len(bulk_df))
    display_df, display_meta = _bulk_display_frame(
        bulk_df,
        min_basemean=float(args.min_basemean),
    )
    artifacts.mode_details["bulk_display_filter"] = display_meta
    upstream_summary = bundle.upstream_result.get("summary", {}) if bundle.upstream_result else {}
    if upstream_summary.get("lfc_shrinkage_applied"):
        coeff = str(upstream_summary.get("lfc_shrinkage_coeff", "")).strip()
        coeff_note = f" ({coeff})" if coeff else ""
        artifacts.notes.append(f"Upstream rnaseq-de applied PyDESeq2 log2 fold-change shrinkage{coeff_note}.")
    if display_meta.get("applied"):
        artifacts.notes.append(
            "Bulk display filter applied: "
            f"baseMean > {args.min_basemean:g} retained "
            f"{display_meta['rows_after']} of {display_meta['rows_before']} rows."
        )
    elif display_meta.get("note"):
        artifacts.notes.append(str(display_meta["note"]))

    top_df, sig_df = _sorted_significant_bulk(
        display_df,
        padj_threshold=args.padj_threshold,
        lfc_threshold=args.lfc_threshold,
    )
    artifacts.n_significant = int(len(sig_df))
    artifacts.mode_details["bulk_display_rows"] = int(len(display_df))
    artifacts.mode_details["bulk_display_significant_rows"] = int(len(sig_df))

    top_table = top_df.head(args.top_genes).copy()
    artifacts.top_table_preview = top_table
    artifacts.table_paths.append(_write_csv(tables_dir / "top_genes.csv", top_table))
    artifacts.table_paths.append(_write_csv(tables_dir / "significant_genes.csv", sig_df))

    artifacts.figure_paths.append(
        plot_bulk_volcano(
            display_df,
            figures_dir / "volcano.png",
            padj_threshold=args.padj_threshold,
            lfc_threshold=args.lfc_threshold,
            label_top=args.label_top,
        )
    )
    artifacts.figure_paths.append(
        plot_top_gene_bar(
            top_table.rename(columns={"gene": "gene", "log2FoldChange": "log2FoldChange"}),
            figures_dir / "top_genes_bar.png",
            top_genes=min(args.top_genes, len(top_table)),
            title="Top Bulk Differentially Expressed Genes",
        )
    )
    artifacts.figure_paths.append(
        plot_bulk_effect_bubble(
            display_df,
            figures_dir / "effect_bubble.png",
            padj_threshold=args.padj_threshold,
            lfc_threshold=args.lfc_threshold,
            label_top=args.label_top,
        )
    )
    artifacts.figure_paths.append(
        plot_bulk_metric_heatmap(
            top_df,
            figures_dir / "gene_metrics_heatmap.png",
            top_genes=min(args.top_genes, 20),
        )
    )
    if "baseMean" in bulk_df.columns:
        artifacts.figure_paths.append(
            plot_ma(
                display_df,
                figures_dir / "ma_plot.png",
                padj_threshold=args.padj_threshold,
                lfc_threshold=args.lfc_threshold,
            )
        )
    else:
        artifacts.notes.append("Skipped MA plot because the input did not contain a baseMean column.")

    counts_path = Path(args.counts).expanduser().resolve() if args.counts else None
    metadata_path = Path(args.metadata).expanduser().resolve() if args.metadata else None
    if bundle.source_kind == "demo-bulk":
        counts_path = counts_path or (_examples_dir() / "demo_bulk_counts.csv")
        metadata_path = metadata_path or (_examples_dir() / "demo_bulk_metadata.csv")
    if counts_path and metadata_path:
        if counts_path.exists() and metadata_path.exists():
            heatmap_genes = top_table["gene"].dropna().astype(str).tolist()[: min(args.top_genes, 20)]
            try:
                bulk_counts, bulk_metadata, group_col = _prepare_bulk_matrix_inputs(
                    counts_path,
                    metadata_path,
                    str(bundle.upstream_result.get("summary", {}).get("contrast", "")),
                )
                heatmap_path, heatmap_meta = plot_bulk_heatmap(
                    heatmap_genes,
                    bulk_counts,
                    bulk_metadata,
                    group_col,
                    figures_dir / "top_genes_heatmap.png",
                )
                artifacts.figure_paths.append(heatmap_path)
                artifacts.enhanced_inputs_used.extend(["counts", "metadata"])
                if group_col:
                    artifacts.mode_details["bulk_group_column"] = group_col
                    artifacts.mode_details["bulk_group_counts"] = heatmap_meta.get("group_counts", {})
                if heatmap_meta.get("downsampled"):
                    artifacts.notes.append(
                        "Bulk heatmap was downsampled to "
                        f"{heatmap_meta['selected_sample_count']} representative samples "
                        f"from {heatmap_meta['original_sample_count']} total samples for readability."
                    )
                if group_col and len(set(bulk_metadata[group_col].astype(str))) == 2:
                    artifacts.figure_paths.append(
                        plot_bulk_mean_comparison(
                            heatmap_genes,
                            bulk_counts,
                            bulk_metadata,
                            group_col,
                            figures_dir / "mean_expression_comparison.png",
                        )
                    )
            except Exception as exc:
                artifacts.notes.append(f"Skipped bulk heatmap: {exc}")
        else:
            artifacts.notes.append("Skipped bulk heatmap because --counts or --metadata path does not exist.")
    else:
        artifacts.notes.append(
            "Skipped bulk heatmap because both --counts and --metadata were not provided."
        )


def render_scrna_outputs(
    bundle: InputBundle,
    args: argparse.Namespace,
    output_dir: Path,
    artifacts: RunArtifacts,
) -> None:
    figures_dir = output_dir / "figures"
    tables_dir = output_dir / "tables"

    contrast_df = bundle.scrna_contrast_df
    within_cluster_df = bundle.scrna_within_cluster_contrast_df
    markers_df = bundle.scrna_markers_df
    total_rows = 0

    if contrast_df is not None:
        total_rows += int(len(contrast_df))
        volcano_path, contrast_top = plot_scrna_contrast_volcano(
            contrast_df,
            figures_dir / "contrast_volcano.png",
            top_genes=args.top_genes,
            label_top=args.label_top,
            padj_threshold=args.padj_threshold,
            lfc_threshold=args.lfc_threshold,
        )
        artifacts.figure_paths.append(volcano_path)
        bar_path, top_markers = plot_scrna_top_markers_bar(
            contrast_df,
            figures_dir / "top_markers_bar.png",
            top_genes=args.top_genes,
            title="Top scRNA Contrast Markers",
        )
        artifacts.figure_paths.append(bar_path)
        artifacts.table_paths.append(_write_csv(tables_dir / "top_markers.csv", top_markers))
        artifacts.n_significant += _count_significant_scrna_rows(
            contrast_df,
            padj_threshold=args.padj_threshold,
        )
        artifacts.top_table_preview = contrast_top.head(args.top_genes)

    if within_cluster_df is not None:
        total_rows += int(len(within_cluster_df))
        panel_path, top_by_panel = plot_within_cluster_marker_panels(
            within_cluster_df,
            figures_dir / "within_cluster_marker_panels.png",
            top_genes=max(1, min(args.top_genes, 6)),
        )
        artifacts.figure_paths.append(panel_path)
        artifacts.table_paths.append(_write_csv(tables_dir / "within_cluster_top_markers.csv", top_by_panel))
        artifacts.n_significant += _count_significant_scrna_rows(
            within_cluster_df,
            padj_threshold=args.padj_threshold,
        )
        artifacts.mode_details["within_cluster_comparisons"] = int(
            within_cluster_df[["cluster", "comparison_id"]].drop_duplicates().shape[0]
        )
        artifacts.mode_details["within_cluster_clusters"] = int(
            within_cluster_df["cluster"].astype(str).nunique()
        )
        if artifacts.top_table_preview is None:
            artifacts.top_table_preview = top_by_panel.head(args.top_genes)

    if markers_df is not None:
        total_rows += int(len(markers_df))
        marker_rank_path, top_by_cluster = plot_marker_rank_bars(
            markers_df,
            figures_dir / "marker_rank_bars.png",
            per_cluster=max(1, min(args.top_genes, 5)),
        )
        artifacts.figure_paths.append(marker_rank_path)
        artifacts.table_paths.append(_write_csv(tables_dir / "top_markers_by_cluster.csv", top_by_cluster))
        if artifacts.top_table_preview is None:
            artifacts.top_table_preview = top_by_cluster.head(args.top_genes)

    artifacts.n_input_rows = total_rows

    adata_path = Path(args.adata).expanduser().resolve() if args.adata else None
    if not adata_path:
        artifacts.notes.append("Skipped enhanced scRNA plots because --adata was not provided.")
        return
    if not adata_path.exists():
        artifacts.notes.append(f"Skipped enhanced scRNA plots because adata path does not exist: {adata_path}")
        return

    try:
        adata = _load_adata(adata_path)
    except Exception as exc:
        artifacts.notes.append(f"Skipped enhanced scRNA plots: {exc}")
        return

    selected_genes = _top_scrna_genes(contrast_df, within_cluster_df, markers_df, args.top_genes)
    available_genes = [gene for gene in selected_genes if gene in adata.var_names]
    if not available_genes:
        artifacts.notes.append("Skipped enhanced scRNA plots because selected genes were not present in AnnData.")
        return

    if within_cluster_df is not None and not within_cluster_df.empty:
        artifacts.notes.append(
            "Skipped scRNA dotplot and heatmap because within-cluster contrast inputs do not define a single grouping axis."
        )
    else:
        groupby = _infer_scrna_groupby(adata)
        if not groupby:
            artifacts.notes.append(
                "Skipped scRNA marker dotplot and heatmap because no suitable grouping column was found in AnnData."
            )
        else:
            genes_for_group_plots = available_genes[: min(10, len(available_genes))]
            try:
                artifacts.figure_paths.append(
                    plot_manual_dotplot(
                        adata,
                        genes_for_group_plots,
                        groupby,
                        figures_dir / "marker_dotplot.png",
                    )
                )
                artifacts.figure_paths.append(
                    plot_manual_heatmap(
                        adata,
                        genes_for_group_plots,
                        groupby,
                        figures_dir / "marker_heatmap.png",
                    )
                )
                artifacts.groupby_used = groupby
                artifacts.enhanced_inputs_used.append("adata")
            except Exception as exc:
                artifacts.notes.append(f"Skipped scRNA dotplot/heatmap: {exc}")

    try:
        artifacts.figure_paths.append(
            plot_umap_feature_panel(
                adata,
                available_genes[: min(4, len(available_genes))],
                figures_dir / "umap_feature_panel.png",
            )
        )
        if "adata" not in artifacts.enhanced_inputs_used:
            artifacts.enhanced_inputs_used.append("adata")
    except Exception as exc:
        artifacts.notes.append(f"Skipped UMAP feature panel: {exc}")


def _markdown_image_list(figure_paths: list[Path]) -> list[str]:
    labels = {
        "volcano.png": "Bulk Volcano Plot",
        "top_genes_bar.png": "Top Bulk Genes",
        "effect_bubble.png": "Effect Bubble Plot",
        "ma_plot.png": "Bulk MA Plot",
        "top_genes_heatmap.png": "Top Gene Heatmap",
        "gene_metrics_heatmap.png": "Gene Metrics Panel",
        "mean_expression_comparison.png": "Mean Expression Comparison",
        "contrast_volcano.png": "scRNA Contrast Volcano",
        "top_markers_bar.png": "Top scRNA Markers",
        "within_cluster_marker_panels.png": "Within-Cluster Marker Panels",
        "marker_rank_bars.png": "Marker Rank Bars",
        "marker_dotplot.png": "Marker Dotplot",
        "marker_heatmap.png": "Marker Heatmap",
        "umap_feature_panel.png": "UMAP Feature Panel",
    }
    lines: list[str] = []
    for figure_path in figure_paths:
        label = labels.get(figure_path.name, figure_path.name)
        lines.append(f"### {label}")
        lines.append("")
        lines.append(f"![{label}](figures/{figure_path.name})")
        lines.append("")
    return lines


def render_markdown_report(
    bundle: InputBundle,
    args: argparse.Namespace,
    output_dir: Path,
    artifacts: RunArtifacts,
) -> Path:
    metadata = {
        "Mode": bundle.mode,
        "Source kind": bundle.source_kind,
        "Input rows": str(artifacts.n_input_rows),
        "Significant rows": str(artifacts.n_significant),
        "Top genes requested": str(args.top_genes),
        "Label count": str(args.label_top),
        "Adjusted p-value threshold": str(args.padj_threshold),
        "Absolute log fold-change threshold": str(args.lfc_threshold),
    }
    bulk_display_filter = artifacts.mode_details.get("bulk_display_filter", {})
    if bulk_display_filter.get("applied"):
        metadata["Bulk display filter"] = f"baseMean > {args.min_basemean:g}"
        metadata["Bulk display rows"] = str(artifacts.mode_details.get("bulk_display_rows", artifacts.n_input_rows))
    if artifacts.groupby_used:
        metadata["scRNA grouping column"] = artifacts.groupby_used
    if artifacts.mode_details.get("within_cluster_comparisons"):
        metadata["Within-cluster comparisons"] = str(artifacts.mode_details["within_cluster_comparisons"])
        metadata["Within-cluster clusters"] = str(artifacts.mode_details.get("within_cluster_clusters", 0))
    header = generate_report_header(
        title="ClawBio Differential Visualisation Report",
        skill_name="diffviz",
        input_files=bundle.input_files,
        extra_metadata=metadata,
    )

    lines = ["## Summary", ""]
    lines.append(f"- Mode: **{bundle.mode}**")
    lines.append(f"- Source kind: **{bundle.source_kind}**")
    lines.append(f"- Input rows processed: **{artifacts.n_input_rows}**")
    lines.append(f"- Significant rows: **{artifacts.n_significant}**")
    lines.append(f"- Figures generated: **{len(artifacts.figure_paths)}**")
    if bulk_display_filter.get("applied"):
        lines.append(
            f"- Bulk display filter: **baseMean > {args.min_basemean:g}** "
            f"({bulk_display_filter['rows_after']} / {bulk_display_filter['rows_before']} rows retained)"
        )
    if artifacts.enhanced_inputs_used:
        lines.append(f"- Enhanced inputs used: **{', '.join(sorted(set(artifacts.enhanced_inputs_used)))}**")
    else:
        lines.append("- Enhanced inputs used: none")
    bulk_group_counts = artifacts.mode_details.get("bulk_group_counts", {})
    if bulk_group_counts:
        lines.append("")
        lines.append("### Cohort Composition")
        for group_name, group_count in bulk_group_counts.items():
            lines.append(f"- {group_name}: **{group_count}** samples")
    if artifacts.mode_details.get("within_cluster_comparisons"):
        lines.append(f"- Within-cluster comparisons: **{artifacts.mode_details['within_cluster_comparisons']}**")
        lines.append(f"- Clusters with within-cluster results: **{artifacts.mode_details.get('within_cluster_clusters', 0)}**")

    if artifacts.top_table_preview is not None and not artifacts.top_table_preview.empty:
        lines.extend(["", "## Top Features", ""])
        preview = artifacts.top_table_preview.head(min(10, len(artifacts.top_table_preview))).copy()
        preview_cols = [
            col
            for col in preview.columns
            if col
            in {
                "gene",
                "names",
                "cluster",
                "comparison_id",
                "group1",
                "group2",
                "groupby",
                "log2FoldChange",
                "logfoldchanges",
                "scores",
                "padj",
                "pvals_adj",
            }
        ]
        preview = preview[preview_cols]
        lines.append("| " + " | ".join(preview.columns) + " |")
        lines.append("|" + "|".join(["---"] * len(preview.columns)) + "|")
        for row in preview.itertuples(index=False):
            lines.append("| " + " | ".join(str(value) for value in row) + " |")

    lines.extend(["", "## Figures", ""])
    lines.extend(_markdown_image_list(artifacts.figure_paths))

    lines.extend(["## Tables", ""])
    for table_path in artifacts.table_paths:
        lines.append(f"- `tables/{table_path.name}`")

    if bundle.notes or artifacts.notes:
        lines.extend(["", "## Notes", ""])
        for note in bundle.notes + artifacts.notes:
            lines.append(f"- {note}")

    lines.extend(["", "## Reproducibility", ""])
    lines.append("- `reproducibility/commands.sh`")
    lines.append("- `reproducibility/environment.yml`")
    lines.append("- `reproducibility/checksums.sha256`")

    report_path = output_dir / "report.md"
    report_path.write_text(header + "\n".join(lines) + generate_report_footer(), encoding="utf-8")
    return report_path


def _img_to_data_uri(path: Path) -> str:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def _figure_gallery_html(figure_paths: list[Path]) -> str:
    if not figure_paths:
        return "<p>No figures generated.</p>"

    labels = {
        "volcano.png": "Bulk Volcano Plot",
        "top_genes_bar.png": "Top Bulk Genes",
        "effect_bubble.png": "Effect Bubble Plot",
        "ma_plot.png": "Bulk MA Plot",
        "top_genes_heatmap.png": "Top Gene Heatmap",
        "gene_metrics_heatmap.png": "Gene Metrics Panel",
        "mean_expression_comparison.png": "Mean Expression Comparison",
        "contrast_volcano.png": "scRNA Contrast Volcano",
        "top_markers_bar.png": "Top scRNA Markers",
        "marker_rank_bars.png": "Marker Rank Bars",
        "marker_dotplot.png": "Marker Dotplot",
        "marker_heatmap.png": "Marker Heatmap",
        "umap_feature_panel.png": "UMAP Feature Panel",
    }

    first_path = figure_paths[0]
    first_uri = _img_to_data_uri(first_path)
    first_title = labels.get(first_path.name, first_path.stem.replace("_", " ").title())
    cards: list[str] = []
    for idx, path in enumerate(figure_paths):
        uri = _img_to_data_uri(path)
        title = labels.get(path.name, path.stem.replace("_", " ").title())
        is_active = " is-active" if idx == 0 else ""
        cards.append(
            f'<button type="button" class="gallery-thumb{is_active}" '
            f'onclick="window.clawbioSelectFigure(this)" '
            f'data-title="{title}" data-src="{uri}">'
            f'<img src="{uri}" alt="{title}" />'
            f"<span>{title}</span>"
            f"</button>"
        )

    viewer = (
        '<div class="gallery-shell">'
        '<div class="gallery-sidebar">'
        '<div class="gallery-thumb-grid">'
        + "".join(cards)
        + "</div></div>"
        '<div class="gallery-stage">'
        f'<h3 id="gallery-active-title">{first_title}</h3>'
        f'<img id="gallery-active-image" src="{first_uri}" alt="{first_title}" />'
        "</div>"
        "</div>"
        "<script>"
        "window.clawbioSelectFigure = function(button) {"
        "  const title = button.getAttribute('data-title') || '';"
        "  const src = button.getAttribute('data-src') || '';"
        "  const titleEl = document.getElementById('gallery-active-title');"
        "  const imageEl = document.getElementById('gallery-active-image');"
        "  if (titleEl) titleEl.textContent = title;"
        "  if (imageEl) { imageEl.src = src; imageEl.alt = title; }"
        "  document.querySelectorAll('.gallery-thumb').forEach((el) => el.classList.remove('is-active'));"
        "  button.classList.add('is-active');"
        "};"
        "</script>"
    )
    return viewer


def render_html_report(
    bundle: InputBundle,
    args: argparse.Namespace,
    output_dir: Path,
    artifacts: RunArtifacts,
) -> Path:
    extra_css = """
:root {
  --cb-green-900: #0f3d3e;
  --cb-green-700: #0f766e;
  --cb-green-500: #14b8a6;
  --cb-green-100: #dff7f3;
  --cb-green-50: #f3fbfa;
  --cb-amber-700: #c2410c;
  --cb-bg: #f4f7fb;
  --cb-surface: #ffffff;
  --cb-text: #102a43;
  --cb-text-secondary: #486581;
  --cb-border: #d9e2ec;
  --cb-accent: #2563eb;
}
body {
  max-width: 1180px;
  background: radial-gradient(circle at top left, #f8fbff 0%, #f4f7fb 48%, #eef3f9 100%);
}
.report-header {
  background: linear-gradient(135deg, #0f766e 0%, #0f4c81 52%, #1d4ed8 100%);
}
.gallery-shell {
  display: grid;
  grid-template-columns: minmax(220px, 300px) minmax(0, 1fr);
  gap: 18px;
  margin-top: 12px;
  align-items: start;
}
.gallery-sidebar {
  max-height: 860px;
  overflow: auto;
  padding-right: 4px;
}
.gallery-thumb-grid {
  display: grid;
  gap: 12px;
}
.gallery-thumb {
  appearance: none;
  background: white;
  border: 1px solid var(--cb-border);
  border-radius: 12px;
  padding: 10px;
  box-shadow: var(--cb-shadow-sm);
  cursor: pointer;
  text-align: left;
  transition: transform 120ms ease, box-shadow 120ms ease, border-color 120ms ease;
}
.gallery-thumb:hover {
  transform: translateY(-1px);
  box-shadow: var(--cb-shadow-md);
}
.gallery-thumb.is-active {
  border-color: var(--cb-accent);
  box-shadow: 0 0 0 2px rgba(37, 99, 235, 0.14);
}
.gallery-thumb img {
  width: 100%;
  border-radius: 8px;
  border: 1px solid var(--cb-border);
  display: block;
  aspect-ratio: 16 / 10;
  object-fit: cover;
  background: #f8fafc;
}
.gallery-thumb span {
  display: block;
  margin-top: 8px;
  font-size: 0.95rem;
  font-weight: 600;
  color: var(--cb-text);
}
.gallery-stage {
  background: white;
  border: 1px solid var(--cb-border);
  border-radius: 14px;
  box-shadow: var(--cb-shadow-sm);
  padding: 16px;
}
.gallery-stage h3 {
  margin-top: 0;
  margin-bottom: 12px;
  font-size: 1.1rem;
}
.gallery-stage img {
  width: 100%;
  border-radius: 10px;
  border: 1px solid var(--cb-border);
  background: #f8fafc;
}
@media (max-width: 980px) {
  .gallery-shell {
    grid-template-columns: 1fr;
  }
  .gallery-sidebar {
    max-height: none;
  }
}
"""
    builder = HtmlReportBuilder(
        "ClawBio Differential Visualisation Report",
        "diffviz",
        extra_css=extra_css,
    )
    builder.add_disclaimer()
    builder.add_header_block(
        "Differential Expression Visualisation",
        "Bulk RNA-seq and single-cell figure/report packaging",
    )
    builder.add_metadata(
        {
            "Mode": bundle.mode,
            "Source kind": bundle.source_kind,
            "Input rows": str(artifacts.n_input_rows),
            "Significant rows": str(artifacts.n_significant),
            "Figures generated": str(len(artifacts.figure_paths)),
            "Enhanced inputs": ", ".join(sorted(set(artifacts.enhanced_inputs_used))) or "none",
        }
    )
    builder.add_summary_cards(
        [
            ("Input rows", artifacts.n_input_rows, "standard"),
            ("Significant", artifacts.n_significant, "caution" if artifacts.n_significant else "indeterminate"),
            ("Figures", len(artifacts.figure_paths), "standard"),
        ]
    )
    bulk_display_filter = artifacts.mode_details.get("bulk_display_filter", {})
    bulk_group_counts = artifacts.mode_details.get("bulk_group_counts", {})
    if bulk_group_counts:
        builder.add_section("Cohort Composition")
        builder.add_table(
            ["Group", "Samples"],
            [[str(group), str(count)] for group, count in bulk_group_counts.items()],
        )
    builder.add_section("Thresholds")
    builder.add_table(
        ["Parameter", "Value"],
        [
            ["Top genes", str(args.top_genes)],
            ["Label count", str(args.label_top)],
            ["Adjusted p-value threshold", str(args.padj_threshold)],
            ["Absolute log fold-change threshold", str(args.lfc_threshold)],
            ["Bulk display filter", f"baseMean > {args.min_basemean:g}" if bulk_display_filter.get("applied") else "not applied"],
        ],
    )
    if artifacts.top_table_preview is not None and not artifacts.top_table_preview.empty:
        builder.add_section("Top Features")
        preview = artifacts.top_table_preview.head(min(10, len(artifacts.top_table_preview))).copy()
        preview_cols = [
            col
            for col in preview.columns
            if col
            in {
                "gene",
                "names",
                "cluster",
                "comparison_id",
                "group1",
                "group2",
                "groupby",
                "log2FoldChange",
                "logfoldchanges",
                "scores",
                "padj",
                "pvals_adj",
            }
        ]
        preview = preview[preview_cols]
        builder.add_table_wrapped(preview.columns.tolist(), preview.astype(str).values.tolist())
    builder.add_section("Figure Gallery")
    builder.add_raw_html(_figure_gallery_html(artifacts.figure_paths))
    if bundle.notes or artifacts.notes:
        builder.add_section("Notes")
        for note in bundle.notes + artifacts.notes:
            builder.add_paragraph(note)
    builder.add_footer_block("diffviz", VERSION)
    html_content = builder.render()
    return write_html_report(output_dir, "report.html", html_content)


def _bundle_input_checksum(paths: list[Path]) -> str:
    digest = hashlib.sha256()
    for path in sorted((p for p in paths if p.exists()), key=lambda item: str(item)):
        digest.update(path.name.encode("utf-8"))
        digest.update(sha256_file(path).encode("utf-8"))
    return digest.hexdigest()


def build_result_json(
    bundle: InputBundle,
    args: argparse.Namespace,
    output_dir: Path,
    artifacts: RunArtifacts,
) -> Path:
    data = {
        "mode": bundle.mode,
        "source_kind": bundle.source_kind,
        "input_files": [path.name for path in bundle.input_files],
        "figures": [path.name for path in artifacts.figure_paths],
        "tables": [path.name for path in artifacts.table_paths],
        "thresholds": {
            "top_genes": args.top_genes,
            "label_top": args.label_top,
            "padj_threshold": args.padj_threshold,
            "lfc_threshold": args.lfc_threshold,
        },
        "notes": bundle.notes + artifacts.notes,
        "enhanced_inputs_used": sorted(set(artifacts.enhanced_inputs_used)),
        "groupby_used": artifacts.groupby_used,
        "disclaimer": DISCLAIMER,
    }
    if artifacts.top_table_preview is not None:
        data["top_preview"] = artifacts.top_table_preview.head(10).to_dict(orient="records")
    if bundle.upstream_result:
        data["upstream_result_meta"] = {
            "skill": bundle.upstream_result.get("skill", ""),
            "summary": bundle.upstream_result.get("summary", {}),
        }
    return write_result_json(
        output_dir=output_dir,
        skill="diffviz",
        version=VERSION,
        summary={
            "mode": bundle.mode,
            "source_kind": bundle.source_kind,
            "n_input_rows": artifacts.n_input_rows,
            "n_significant": artifacts.n_significant,
            "bulk_display_rows": artifacts.mode_details.get("bulk_display_rows"),
            "bulk_display_filter": artifacts.mode_details.get("bulk_display_filter", {}),
            "within_cluster_comparisons": artifacts.mode_details.get("within_cluster_comparisons"),
            "within_cluster_clusters": artifacts.mode_details.get("within_cluster_clusters"),
            "figures_generated": [path.name for path in artifacts.figure_paths],
            "enhanced_inputs_used": sorted(set(artifacts.enhanced_inputs_used)),
        },
        data=data,
        input_checksum=_bundle_input_checksum(bundle.input_files),
    )


def build_repro_command(args: argparse.Namespace) -> str:
    parts = ["python3", "skills/diff_visualizer/scripts/diff_visualizer.py"]
    if args.demo:
        parts.append("--demo")
    else:
        parts.extend(["--input", args.input])
    parts.extend(["--output", args.output])
    if args.mode != "auto":
        parts.extend(["--mode", args.mode])
    if args.counts:
        parts.extend(["--counts", args.counts])
    if args.metadata:
        parts.extend(["--metadata", args.metadata])
    if args.adata:
        parts.extend(["--adata", args.adata])
    if args.top_genes != 25:
        parts.extend(["--top-genes", str(args.top_genes)])
    if args.label_top != 10:
        parts.extend(["--label-top", str(args.label_top)])
    if args.padj_threshold != 0.05:
        parts.extend(["--padj-threshold", str(args.padj_threshold)])
    if args.lfc_threshold != 1.0:
        parts.extend(["--lfc-threshold", str(args.lfc_threshold)])
    if args.min_basemean != 10.0:
        parts.extend(["--min-basemean", str(args.min_basemean)])
    return " ".join(parts)


def write_reproducibility(
    bundle: InputBundle,
    args: argparse.Namespace,
    output_dir: Path,
    artifacts: RunArtifacts,
    report_path: Path,
    html_path: Path,
    result_path: Path,
) -> None:
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)
    commands = f"""#!/usr/bin/env bash
set -euo pipefail

{build_repro_command(args)}
"""
    (repro_dir / "commands.sh").write_text(commands, encoding="utf-8")
    environment = """name: clawbio-diffviz
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.11
  - pandas=2.2
  - numpy=1.26
  - matplotlib=3.8
  - pip
  - pip:
      - anndata>=0.10
      - scanpy>=1.10
"""
    (repro_dir / "environment.yml").write_text(environment, encoding="utf-8")
    checksum_targets = [*bundle.input_files, report_path, html_path, result_path, *artifacts.figure_paths, *artifacts.table_paths]
    lines = []
    for path in checksum_targets:
        if not path.exists():
            continue
        rel = path.relative_to(output_dir) if path.is_relative_to(output_dir) else path.name
        lines.append(f"{sha256_file(path)}  {rel}")
    (repro_dir / "checksums.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> dict[str, Any]:
    output_dir = Path(args.output).expanduser().resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(
            f"Output directory '{output_dir}' is not empty. Choose a new --output path to avoid overwriting reports."
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)

    bundle = resolve_input_bundle(args)
    artifacts = RunArtifacts()

    if bundle.mode == "bulk":
        render_bulk_outputs(bundle, args, output_dir, artifacts)
    else:
        render_scrna_outputs(bundle, args, output_dir, artifacts)

    report_path = render_markdown_report(bundle, args, output_dir, artifacts)
    html_path = render_html_report(bundle, args, output_dir, artifacts)
    result_path = build_result_json(bundle, args, output_dir, artifacts)
    write_reproducibility(bundle, args, output_dir, artifacts, report_path, html_path, result_path)

    return {
        "output_dir": str(output_dir),
        "report_path": str(report_path),
        "html_path": str(html_path),
        "result_path": str(result_path),
        "figures": [path.name for path in artifacts.figure_paths],
        "tables": [path.name for path in artifacts.table_paths],
        "mode": bundle.mode,
        "source_kind": bundle.source_kind,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="ClawBio differential expression visualiser for bulk RNA-seq and scRNA outputs."
    )
    parser.add_argument("--input", help="Input directory or result table")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--mode", choices=["auto", "bulk", "scrna"], default="auto")
    parser.add_argument("--counts", help="Optional bulk counts matrix for heatmap enrichment")
    parser.add_argument("--metadata", help="Optional bulk sample metadata for heatmap enrichment")
    parser.add_argument("--adata", help="Optional AnnData .h5ad for scRNA enhancement plots")
    parser.add_argument("--top-genes", type=int, default=25, help="Top genes/markers to prioritise")
    parser.add_argument("--label-top", type=int, default=10, help="Top points to label on volcano plots")
    parser.add_argument("--padj-threshold", type=float, default=0.05, help="Adjusted p-value threshold")
    parser.add_argument("--lfc-threshold", type=float, default=1.0, help="Absolute log fold-change threshold")
    parser.add_argument("--min-basemean", type=float, default=10.0, help="Minimum baseMean retained in bulk display plots/tables")
    parser.add_argument("--demo", action="store_true", help="Run with bundled demo inputs")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    try:
        result = run(args)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

    print("ClawBio diffviz completed")
    print(f"  Output: {result['output_dir']}")
    print(f"  Mode:   {result['mode']}")
    print(f"  Source: {result['source_kind']}")
    print(f"  Figures: {', '.join(result['figures'])}")


if __name__ == "__main__":
    main()
