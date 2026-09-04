#!/usr/bin/env python3
"""completeness-aware-caller — three-state gene calls and ANI gating for
incomplete genomes.

Turns (gene search results, a completeness estimate, optional ANI
comparisons) into calls that refuse to overreach: a gene not found in a
genome that is only 70% complete is NOT absent — it is `cannot_conclude`.
Likewise, ranking two references by ANI is refused when the margin between
them is within the drift that incompleteness alone can induce.

Domain decisions (see SKILL.md for citations):
- Absence may be asserted only when completeness >= 0.95, so that under a
  random-loss model the probability of having missed a truly present gene
  stays at or below 5%.
- ANI drift is modelled as DRIFT_COEFF * (1 - completeness), calibrated on
  the STM815 degradation benchmark (0.41 ANI points of drift at 50%
  retention). A ranking must exceed SAFETY_FACTOR times that drift.
- ANI values inside the 94-96% species-boundary zone with completeness
  below 0.9 are flagged as uncertain species assignments.
"""

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

MIN_COMPLETENESS_FOR_ABSENT = 0.95
DRIFT_COEFF = 0.82          # ANI points per unit incompleteness (benchmark)
SAFETY_FACTOR = 2.0
SPECIES_BOUNDARY = (94.0, 96.0)
BOUNDARY_MIN_COMPLETENESS = 0.90

DISCLAIMER = (
    "*ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions.*"
)


def parse_completeness(value):
    """Accept a fraction (0-1] or a percentage (1-100]; return a fraction."""
    c = float(value)
    if 1.0 < c <= 100.0:
        c /= 100.0
    if not 0.0 < c <= 1.0:
        raise ValueError(
            f"completeness {value} outside (0, 1] as fraction or (1, 100] as percent"
        )
    return c


def call_gene(found, completeness, min_completeness=MIN_COMPLETENESS_FOR_ABSENT):
    """Three-state call for one gene given a completeness estimate."""
    if found:
        return {
            "status": "present",
            "confidence": 1.0,
            "message": "Detected in the assembly; detection is positive evidence.",
        }
    miss_probability = 1.0 - completeness
    if completeness >= min_completeness:
        return {
            "status": "absent",
            "confidence": completeness,
            "message": (
                f"Not found at {completeness:.1%} completeness; residual chance "
                f"of a missed gene is ~{miss_probability:.1%}."
            ),
        }
    return {
        "status": "cannot_conclude",
        "confidence": completeness,
        "message": (
            f"CANNOT CONCLUDE absence: the assembly is only "
            f"{completeness * 100:.0f}% complete, so a truly present gene "
            f"would be missed with ~{miss_probability:.0%} probability. "
            f"Re-sequence or close the assembly before asserting loss of "
            f"function."
        ),
    }


def ani_margin_gate(ani_a, ani_b, completeness,
                    drift_coeff=DRIFT_COEFF, safety=SAFETY_FACTOR):
    """Decide whether ANI to reference A vs B may be ranked at this
    completeness. Margin must exceed safety * modelled drift."""
    margin = abs(ani_a - ani_b)
    drift = drift_coeff * (1.0 - completeness)
    lo, hi = SPECIES_BOUNDARY
    boundary_uncertain = (
        completeness < BOUNDARY_MIN_COMPLETENESS
        and any(lo <= ani <= hi for ani in (ani_a, ani_b))
    )
    if margin >= safety * drift:
        decision = "rank"
        message = (
            f"Margin {margin:.2f} ANI points exceeds {safety:.0f}x modelled "
            f"drift ({drift:.2f}); ranking is supported."
        )
    else:
        decision = "cannot_conclude"
        message = (
            f"CANNOT CONCLUDE ranking: margin {margin:.2f} ANI points is "
            f"within {safety:.0f}x the drift ({drift:.2f}) that "
            f"{(1 - completeness):.0%} incompleteness alone can induce. "
            f"Do not report which reference is closer."
        )
    return {
        "decision": decision,
        "margin": margin,
        "drift": drift,
        "species_boundary_uncertain": boundary_uncertain,
        "message": message,
    }


# ------------------------------------------------------------------ reporting
def write_outputs(outdir, completeness, gene_calls, ani_gate, argv):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    result = {
        "completeness": completeness,
        "gene_calls": gene_calls,
        "ani_gate": ani_gate,
        "thresholds": {
            "min_completeness_for_absent": MIN_COMPLETENESS_FOR_ABSENT,
            "drift_coeff": DRIFT_COEFF,
            "safety_factor": SAFETY_FACTOR,
        },
        "generated": datetime.now(timezone.utc).isoformat(),
    }
    (outdir / "result.json").write_text(json.dumps(result, indent=2))

    lines = [
        "# Completeness-Aware Caller Report",
        "",
        f"**Assembly completeness**: {completeness:.1%}",
        "",
        "## Gene calls",
        "",
        "| Gene | Status | Confidence | Rationale |",
        "|------|--------|-----------|-----------|",
    ]
    for g in gene_calls:
        lines.append(
            f"| {g['gene']} | {g['status'].upper().replace('_', ' ')} | "
            f"{g['confidence']:.2f} | {g['message']} |"
        )
    if ani_gate is not None:
        lines += [
            "",
            "## ANI ranking gate",
            "",
            f"**Decision**: {ani_gate['decision'].upper().replace('_', ' ')}",
            "",
            ani_gate["message"],
        ]
        if ani_gate["species_boundary_uncertain"]:
            lines.append(
                "\n⚠️ At least one ANI value falls in the 94-96% species-"
                "boundary zone while completeness is below 90%: species "
                "assignment is uncertain."
            )
    lines += ["", "---", "", DISCLAIMER, ""]
    (outdir / "report.md").write_text("\n".join(lines))

    (outdir / "commands.sh").write_text(
        "#!/bin/sh\n# Replay command\n" + " ".join(argv) + "\n"
    )


# ------------------------------------------------------------------ demo
# Mirrors the real frag70 level of the STM815 benchmark. nodB/nodC/nodS/nodU
# lose every copy and vanish from the assembly; nifH survives only because
# STM815 carries two copies of it on pBPHY02 — copy number, not importance,
# decides which genes a fragmented assembly can still show you.
DEMO_GENES = [
    {"gene": "nodC", "found": False},
    {"gene": "nodB", "found": False},
    {"gene": "nifH", "found": True},
    {"gene": "nifD", "found": True},
]
DEMO_COMPLETENESS = 0.70
DEMO_ANI = (81.45, 80.75)  # STM815@70% vs LB400 and vs J2315 (benchmark)


def run(args, argv):
    if args.demo:
        completeness = DEMO_COMPLETENESS
        gene_specs = DEMO_GENES
        ani_pair = DEMO_ANI
    else:
        if args.busco_json:
            busco = json.loads(Path(args.busco_json).read_text())
            completeness = parse_completeness(busco["C"])
        elif args.completeness:
            completeness = parse_completeness(args.completeness)
        else:
            raise SystemExit("Provide --completeness, --busco-json, or --demo")
        gene_specs = (
            json.loads(Path(args.genes).read_text()) if args.genes else []
        )
        ani_pair = (
            (args.ani_a, args.ani_b)
            if args.ani_a is not None and args.ani_b is not None
            else None
        )

    gene_calls = [
        {"gene": spec["gene"], **call_gene(spec["found"], completeness)}
        for spec in gene_specs
    ]
    ani_gate = (
        ani_margin_gate(ani_pair[0], ani_pair[1], completeness)
        if ani_pair else None
    )
    write_outputs(args.output, completeness, gene_calls, ani_gate, argv)
    print(f"Report written to: {Path(args.output) / 'report.md'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--completeness", help="Fraction (0-1] or percent (1-100]")
    ap.add_argument("--busco-json",
                    help="Path to a busco-assessor result.json (reads C%%)")
    ap.add_argument("--genes",
                    help='JSON file: [{"gene": "nifH", "found": false}, ...]')
    ap.add_argument("--ani-a", type=float,
                    help="ANI to reference A (e.g. from FastANI)")
    ap.add_argument("--ani-b", type=float, help="ANI to reference B")
    ap.add_argument("--output", required=True, help="Output directory")
    ap.add_argument("--demo", action="store_true",
                    help="Run the bundled STM815@70%% benchmark scenario")
    args = ap.parse_args()
    run(args, sys.argv)


if __name__ == "__main__":
    main()
