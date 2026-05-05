#!/usr/bin/env python3
"""Regenerate the dynamic sections of benchmarks.html from a clawbio_bench
aggregate_report.json snapshot.

Usage:
    python scripts/generate_leaderboard.py \
        --aggregate bench_runs/latest/aggregate_report.json \
        --html benchmarks.html

The HTML template uses HTML comment markers to bound each replaceable
section: BENCH_META, BENCH_SUMMARY, BENCH_TABLE, BENCH_FOOTER, BENCH_JSONLD.
Static content (hero, principles, CTA, design tokens) is left untouched.

The script is idempotent and safe to run repeatedly. It exits 0 on
success and prints a diff summary to stdout. If the HTML is unchanged
the workflow downstream can skip the PR step.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

BENCH_TO_FOLDER = {
    "bio-orchestrator": "bio-orchestrator",
    "equity-scorer": "equity-scorer",
    "nutrigx-advisor": "nutrigx-advisor",
    "pharmgx-reporter": "pharmgx-reporter",
    "claw-metagenomics": "claw-metagenomics",
    "clawbio-finemapping": "fine-mapping",
    "clinical-variant-reporter": "clinical-variant-reporter",
    "cvr-variant-identity": "clinical-variant-reporter",
    "cvr-acmg-correctness": "clinical-variant-reporter",
    "gwas-prs": "gwas-prs",
}

DISPLAY_NAME = {
    "clawbio-finemapping": "fine-mapping",
}

ORIGINAL_AUDIT_DATE = "2026-04-05"
ORIGINAL_AUDIT_COMMIT = "1481fb4"
ORIGINAL_AUDIT_PASS = 80
ORIGINAL_AUDIT_TOTAL = 140


def status_pill(rate: float, harness_errors: int, evaluated: int) -> tuple[str, str]:
    if evaluated == 0 and harness_errors > 0:
        return "status-p0", "Infra"
    if rate >= 75.0:
        return "status-clear", "Clear"
    if rate >= 50.0:
        return "status-watch", "Watch"
    if rate >= 25.0:
        return "status-p1", "P1"
    return "status-p0", "P0"


def rate_color(rate: float, evaluated: int, harness_errors: int) -> str:
    if evaluated == 0 and harness_errors > 0:
        return "var(--text-muted)"
    if rate >= 75.0:
        return "var(--accent)"
    if rate >= 50.0:
        return "var(--warn)"
    return "var(--danger)"


def fail_findings(harness: dict) -> list[str]:
    fail_cats = set(harness.get("fail_categories", []))
    counts = harness.get("categories", {})
    findings = [(c, n) for c, n in counts.items() if c in fail_cats and n > 0]
    findings.sort(key=lambda kv: -kv[1])
    if not findings:
        return ["none"]
    return [f"{cat} ({n})" if n > 1 else cat for cat, n in findings[:3]]


def render_table_rows(harnesses: dict) -> str:
    items = sorted(
        harnesses.items(),
        key=lambda kv: (
            -(kv[1].get("pass_rate", 0.0) if kv[1].get("evaluated", 0) > 0 else -1),
            kv[0],
        ),
    )
    rows = []
    for name, h in items:
        evaluated = h.get("evaluated", 0)
        passed = h.get("pass_count", 0)
        rate = h.get("pass_rate", 0.0)
        errs = h.get("harness_errors", 0)
        folder = BENCH_TO_FOLDER.get(name, name)
        display = DISPLAY_NAME.get(name, name)

        pill_cls, pill_label = status_pill(rate, errs, evaluated)
        bar_color = rate_color(rate, evaluated, errs)
        bar_width = 0.0 if (evaluated == 0 and errs > 0) else rate
        rate_text = "infra" if (evaluated == 0 and errs > 0) else f"{rate:.1f}%"
        total_str = f"{passed} / {evaluated}" if evaluated else f"0 / 0 ({errs} infra)"

        findings = fail_findings(h) if evaluated else [f"harness_error ({errs})"]
        finding_html = "".join(
            f'<span class="finding-tag">{f}</span>' for f in findings
        )

        rows.append(
            f"""            <tr>
              <td><a class="skill-link" href="https://github.com/ClawBio/ClawBio/tree/main/skills/{folder}" target="_blank">{display}</a></td>
              <td>{total_str}</td>
              <td><div class="rate-cell" style="color:{bar_color};">{rate_text}</div><span class="rate-bar-wrap"><span class="rate-bar" style="width:{bar_width}%;background:{bar_color};"></span></span></td>
              <td>{finding_html}</td>
              <td><span class="pill-status {pill_cls}">{pill_label}</span></td>
            </tr>"""
        )
    return "\n".join(rows)


def render_meta(report: dict) -> str:
    suite = report.get("benchmark_suite_version", "?")
    date = report.get("date", "?")
    commit = report.get("clawbio_commit", "?")[:7]
    return f"""        <div class="meta-card">
          <div class="meta-label">Last Run</div>
          <div class="meta-value">{date}</div>
        </div>
        <div class="meta-card">
          <div class="meta-label">Bench</div>
          <div class="meta-value"><a href="https://github.com/biostochastics/clawbio_bench" target="_blank">clawbio_bench v{suite}</a></div>
        </div>
        <div class="meta-card">
          <div class="meta-label">Bench Author</div>
          <div class="meta-value">Biostochastics LLC</div>
        </div>
        <div class="meta-card">
          <div class="meta-label">ClawBio Commit</div>
          <div class="meta-value"><a href="https://github.com/ClawBio/ClawBio/commit/{commit}" target="_blank"><code>{commit}</code></a></div>
        </div>"""


def render_summary(report: dict) -> str:
    total_p = sum(h.get("pass_count", 0) for h in report["harnesses"].values())
    total_t = sum(h.get("evaluated", 0) for h in report["harnesses"].values())
    n_skills = len(report["harnesses"])
    rate = (100 * total_p / total_t) if total_t else 0.0
    return f"""        <div>
          <div class="summary-headline"><span class="pass">{total_p}</span> / {total_t} tests passing <span style="color:var(--text-muted);font-weight:600;">({rate:.1f}%)</span></div>
          <div class="summary-rest">{n_skills} skills audited across 3 dimensions: safety, correctness, honesty. Up from {ORIGINAL_AUDIT_PASS} / {ORIGINAL_AUDIT_TOTAL} ({100*ORIGINAL_AUDIT_PASS/ORIGINAL_AUDIT_TOTAL:.1f}%) at the original {ORIGINAL_AUDIT_DATE} audit.</div>
        </div>"""


def render_footer(report: dict) -> str:
    date = report.get("date", "?")
    commit = report.get("clawbio_commit", "?")[:7]
    return f"""    <p>Last benchmark run: {date} against commit <a href="https://github.com/ClawBio/ClawBio/commit/{commit}" target="_blank"><code>{commit}</code></a>. Original audit baseline {ORIGINAL_AUDIT_DATE} at commit <a href="https://github.com/ClawBio/ClawBio/commit/{ORIGINAL_AUDIT_COMMIT}" target="_blank"><code>{ORIGINAL_AUDIT_COMMIT}</code></a>: {ORIGINAL_AUDIT_PASS} / {ORIGINAL_AUDIT_TOTAL} ({100*ORIGINAL_AUDIT_PASS/ORIGINAL_AUDIT_TOTAL:.1f}%).</p>
    <p style="margin-top:0.5rem;">ClawBio is open source under the <a href="https://github.com/ClawBio/ClawBio/blob/main/LICENSE" target="_blank">MIT License</a>.</p>"""


def render_jsonld(report: dict) -> str:
    return json.dumps(
        {
            "@context": "https://schema.org",
            "@type": "Dataset",
            "name": "ClawBio Benchmark Leaderboard",
            "description": "Public scientific-correctness leaderboard for ClawBio bioinformatics skills, audited by an independent third-party benchmark suite.",
            "url": "https://clawbio.ai/benchmarks.html",
            "creator": {"@type": "Organization", "name": "ClawBio", "url": "https://clawbio.ai"},
            "license": "https://opensource.org/licenses/MIT",
            "isAccessibleForFree": True,
            "dateModified": report.get("date", ""),
            "keywords": "bioinformatics, benchmark, scientific-correctness, AI-agents, genomics, equity-scorer, fine-mapping, pharmgx, nutrigx, clawbio_bench",
        },
        indent=2,
    )


def replace_section(html: str, marker: str, payload: str) -> str:
    pattern = re.compile(
        rf"(<!--\s*{marker}_START\s*-->)(.*?)(<!--\s*{marker}_END\s*-->)",
        re.DOTALL,
    )
    if not pattern.search(html):
        raise SystemExit(f"marker {marker} not found in HTML")
    return pattern.sub(rf"\1\n{payload}\n        \3", html)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--aggregate", required=True, type=Path)
    ap.add_argument("--html", required=True, type=Path)
    args = ap.parse_args()

    report = json.loads(args.aggregate.read_text())
    html = args.html.read_text()
    original = html

    html = replace_section(html, "BENCH_META", render_meta(report))
    html = replace_section(html, "BENCH_SUMMARY", render_summary(report))
    html = replace_section(html, "BENCH_TABLE", render_table_rows(report["harnesses"]))
    html = replace_section(html, "BENCH_FOOTER", render_footer(report))
    html = replace_section(html, "BENCH_JSONLD", render_jsonld(report))

    if html == original:
        print("benchmarks.html unchanged")
        return 0

    args.html.write_text(html)
    total_p = sum(h.get("pass_count", 0) for h in report["harnesses"].values())
    total_t = sum(h.get("evaluated", 0) for h in report["harnesses"].values())
    print(f"benchmarks.html regenerated: {total_p}/{total_t} ({100*total_p/total_t:.1f}%)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
