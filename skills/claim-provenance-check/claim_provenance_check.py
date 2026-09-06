#!/usr/bin/env python3
"""
claim_provenance_check.py — ClawBio Claim Provenance Check Skill
==================================================================
Bind every citation tag in a claim to a retrieved evidence row, and report which evidence
classes were reached for the subject versus NO COVERAGE. Refuses claims with an invented or
missing citation instead of scoring them — no partial credit.

Author:  Sippar (Nuru-AI)
Version: 0.1.0

Usage:
    python claim_provenance_check.py --input <case.json> --output <output_dir>
    python claim_provenance_check.py --demo --output /tmp/claim-provenance-check_demo
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

from provenance_core import (
    ClaimState,
    CoverageReport,
    EvidenceRow,
    SourceOutcome,
    state_for_claim,
)

# No external dependencies required — standard library only.

SAFETY_DISCLAIMER = (
    "This report is a mechanical citation check, not a scientific, clinical, or factual "
    "judgement. SUPPORTED / CONTESTED describe what the cited evidence rows say, not whether "
    "the underlying evidence is correct. NO COVERAGE means the claim could not be bound to any "
    "retrieved evidence row and must be treated as a refusal, never as a weaker form of "
    "support. ClawBio is a research and educational tool, not a medical device, and does not "
    "provide clinical diagnoses. Consult a healthcare professional before making any medical "
    "decisions."
)

REQUIRED_CASE_KEYS = ("subject", "rows", "outcomes", "claims")


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def load_case(input_path: Path) -> dict:
    """Parse and shape-check a case JSON file."""
    try:
        case = json.loads(input_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"{input_path} is not valid JSON: {exc}") from exc

    missing = [key for key in REQUIRED_CASE_KEYS if key not in case]
    if missing:
        raise ValueError(f"{input_path} is missing required key(s): {', '.join(missing)}")
    return case


def build_rows(raw_rows: list[dict]) -> list[EvidenceRow]:
    return [EvidenceRow(**row) for row in raw_rows]


def build_outcomes(raw_outcomes: list[dict]) -> list[SourceOutcome]:
    return [SourceOutcome(**outcome) for outcome in raw_outcomes]


def evaluate_claims(claims: list[str], rows: list[EvidenceRow]) -> list[dict]:
    """Bind every claim and return its verdict, in the order given."""
    verdicts = []
    for claim in claims:
        state, issues = state_for_claim(claim, rows)
        verdicts.append(
            {
                "claim": claim,
                "state": state.value,
                "issues": [{"code": i.code, "detail": i.detail} for i in issues],
            }
        )
    return verdicts


def render_report(subject: str, coverage: CoverageReport, verdicts: list[dict]) -> str:
    shipped = sum(1 for v in verdicts if v["state"] != ClaimState.NO_COVERAGE.value)
    refused = len(verdicts) - shipped

    lines = [
        f"# Claim Provenance Report — {subject}",
        "",
        f"> {SAFETY_DISCLAIMER}",
        "",
        "## Coverage",
        "",
        "```",
        coverage.render(),
        "```",
        "",
        "## Claims",
        "",
    ]
    for v in verdicts:
        lines.append(f"- **[{v['state']}]** {v['claim']}")
        for issue in v["issues"]:
            lines.append(f"    - refused: `{issue['code']}` — {issue['detail']}")
    lines += [
        "",
        f"**{shipped} shipped, {refused} refused.** A refusal is an outcome, not an error.",
        "",
    ]
    return "\n".join(lines)


def run(input_path: Path, output_dir: Path) -> dict:
    """
    Core skill logic: bind every claim in the case file to its retrieved evidence, and render
    the coverage map plus per-claim verdicts.

    Args:
        input_path:  Path to a case JSON file (subject, rows, outcomes, claims).
        output_dir:  Directory where report.md and result.json will be written.

    Returns:
        dict: Machine-readable results (also written to result.json).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    case = load_case(input_path)
    rows = build_rows(case["rows"])
    outcomes = build_outcomes(case["outcomes"])
    coverage = CoverageReport(case["subject"], outcomes)
    verdicts = evaluate_claims(case["claims"], rows)

    report_md = render_report(case["subject"], coverage, verdicts)
    (output_dir / "report.md").write_text(report_md, encoding="utf-8")

    results: dict = {
        "skill": "claim-provenance-check",
        "input": str(input_path),
        "generated_at": datetime.now().isoformat(),
        "coverage": coverage.as_dict(),
        "claims": verdicts,
    }
    (output_dir / "result.json").write_text(
        json.dumps(results, indent=2, default=str), encoding="utf-8"
    )

    return results


# ---------------------------------------------------------------------------
# Demo mode
# ---------------------------------------------------------------------------

def run_demo(output_dir: Path) -> dict:
    """Run on the bundled synthetic demo_input.txt to verify the skill works end-to-end."""
    output_dir.mkdir(parents=True, exist_ok=True)
    demo_input = Path(__file__).with_name("demo_input.txt")
    results = run(demo_input, output_dir)

    print(f"  Demo complete. Output: {output_dir}")
    print(f"  Files: {', '.join(sorted(f.name for f in output_dir.rglob('*') if f.is_file()))}")
    return results


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Bind every citation in a claim to a retrieved evidence row and report which "
            "evidence classes were reached, refusing invented citations with no partial credit."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", dest="input_path", help="Path to a case JSON file")
    parser.add_argument(
        "--output",
        dest="output_dir",
        default="/tmp/claim-provenance-check_output",
        help="Output directory (default: /tmp/claim-provenance-check_output)",
    )
    parser.add_argument("--demo", action="store_true", help="Run on the bundled synthetic demo_input.txt")
    args = parser.parse_args()

    out = Path(args.output_dir)

    if args.demo:
        run_demo(out)
    elif args.input_path:
        run(Path(args.input_path), out)
        print(f"  Done. Output: {out}")
        report = out / "report.md"
        if report.exists():
            print(f"  Report: {report}")
    else:
        parser.error("Provide --input <file> or --demo")


if __name__ == "__main__":
    main()
