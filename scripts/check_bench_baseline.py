#!/usr/bin/env python3
"""Fail CI when clawbio-bench regresses against a committed baseline.

Provenance: #342. The `scientific-audit` job ran the bench as

    uv run clawbio-bench --smoke --repo . || true

so the exit code was discarded and every later step was `if: always()`. The job
reported SUCCESS on every PR for months while the bench itself returned
`pass: false` with seven blocking harnesses, and `nutrigx-advisor` sat at 0/10
in the shipped catalog the whole time (#237). A verification step that reports
success for something nobody checked converts an unknown into a false
assurance, and then nobody looks again.

Why this gate does not simply assert `overall.pass`
---------------------------------------------------
Seven harnesses are below 100% today. Requiring a clean bench would block every
unrelated PR on pre-existing debt, and the pressure to re-add `|| true` would be
immediate and would win. So the gate fails on *regression* from a recorded
baseline: new breakage goes red, existing debt stays visible in `bench_baseline.json`
and is tracked in #106.

Raising a baseline number is a deliberate, reviewable commit. Lowering one
should be treated as a bug report.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Bench rates are floats rendered to one decimal place. Anything smaller than
# this is representation noise, not a regression.
TOLERANCE = 0.05


def load_baseline(path: Path) -> dict:
    """Load the committed baseline. A malformed or absent baseline is fatal:
    silently treating it as empty would make the gate vacuous, which is the
    #342 failure mode again."""
    data = json.loads(Path(path).read_text())
    if not isinstance(data, dict) or not isinstance(data.get("harnesses"), dict):
        raise ValueError(f"{path}: expected an object with a 'harnesses' mapping")
    return data


def compare(
    report: dict,
    baseline: dict,
    baseline_errors: dict | None = None,
) -> list[str]:
    """Return a list of human-readable regressions. Empty means clean.

    Three distinct regressions are detected, and the second matters as much as
    the first: deleting a failing harness must not turn the build green.
    """
    failures: list[str] = []
    harnesses = report.get("harnesses", {})
    expected = baseline.get("harnesses", {})

    for name, want in expected.items():
        entry = harnesses.get(name)
        if entry is None:
            failures.append(
                f"{name}: missing from this run, baseline expects {want}% "
                "(a harness that disappears is a regression, not a pass)"
            )
            continue

        got = entry.get("pass_rate")
        if got is None:
            failures.append(f"{name}: no pass_rate reported, baseline expects {want}%")
            continue

        if got < want - TOLERANCE:
            failures.append(
                f"{name}: pass rate fell from {want}% to {got}% "
                f"(-{round(want - got, 2)} points)"
            )

    if baseline_errors is not None:
        for name, want_errors in baseline_errors.items():
            entry = harnesses.get(name)
            if entry is None:
                continue
            got_errors = entry.get("harness_errors", 0)
            if got_errors > want_errors:
                failures.append(
                    f"{name}: harness errors rose from {want_errors} to {got_errors} "
                    "(an errored harness proves nothing, so this is not a pass)"
                )

    return failures


def _render(report: dict, baseline: dict) -> str:
    """Per-harness table, printed on pass as well as on failure so the rates are
    visible in the job log without downloading an artifact."""
    expected = baseline.get("harnesses", {})
    lines = ["| Harness | Baseline | This run | Delta |", "|---|---|---|---|"]
    for name, entry in sorted(report.get("harnesses", {}).items()):
        got = entry.get("pass_rate")
        want = expected.get(name)
        if want is None:
            lines.append(f"| {name} | (new) | {got}% | - |")
            continue
        delta = round((got or 0) - want, 2)
        marker = "" if delta >= -TOLERANCE else " **REGRESSION**"
        sign = "+" if delta > 0 else ""
        lines.append(f"| {name} | {want}% | {got}%{marker} | {sign}{delta} |")
    overall = report.get("overall", {})
    if overall:
        lines.append("")
        lines.append(
            f"Bench overall: {overall.get('total_pass')}/{overall.get('total_cases')} "
            f"= {overall.get('total_pass_rate')}% "
            f"(bench's own verdict: pass={overall.get('pass')})"
        )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", help="path to the bench aggregate_report.json")
    parser.add_argument(
        "--baseline",
        default=str(Path(__file__).resolve().parents[1] / "bench_baseline.json"),
    )
    parser.add_argument(
        "--summary",
        help="append the rendered table to this file (e.g. $GITHUB_STEP_SUMMARY)",
    )
    args = parser.parse_args(argv)

    report_path = Path(args.report)
    if not report_path.is_file():
        # The bench producing no report at all is the #342 failure mode wearing
        # a different hat: absence of evidence read as evidence of absence.
        print(
            f"ERROR: no bench report at {report_path}. The bench did not run to "
            "completion; this is a failure, not a pass.",
            file=sys.stderr,
        )
        return 2

    try:
        report = json.loads(report_path.read_text())
        baseline = load_baseline(Path(args.baseline))
    except (json.JSONDecodeError, ValueError, OSError) as exc:
        print(f"ERROR: could not read bench report or baseline: {exc}", file=sys.stderr)
        return 2

    table = _render(report, baseline)
    print(table)
    if args.summary:
        with open(args.summary, "a", encoding="utf-8") as handle:
            handle.write("\n## Scientific Correctness Audit\n\n" + table + "\n")

    baseline_errors = baseline.get("harness_errors")
    failures = compare(report, baseline, baseline_errors=baseline_errors)
    if failures:
        print("\nBench regressed against the committed baseline:\n", file=sys.stderr)
        for line in failures:
            print(f"  - {line}", file=sys.stderr)
        print(
            "\nIf a drop is intentional and understood, update bench_baseline.json "
            "in the same PR and say why in the description.",
            file=sys.stderr,
        )
        return 1

    print("\nNo regression against bench_baseline.json.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
