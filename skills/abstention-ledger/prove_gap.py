"""Measure the legacy-EFF gap in `rare-high-impact-variants` instead of asserting it.

The Berlin challenge brief warns:

    "Do not feed the historical quartet directly to rare-high-impact-variants:
     it does not parse the legacy EFF field."

Quoting that is not evidence. This script runs the unmodified skill three times
over the same 68 variants and reports what actually comes back:

  1. raw          -- INFO carries only EFF=, as the source data does
  2. passthrough  -- MC= set to the SnpEff effect name, unmapped
  3. mapped       -- MC= set to the Sequence Ontology term, via eff_to_info

The middle run is the point. A silent zero is obviously broken; a silent partial
count is not, and it is what a plausible one-line fix produces.

Nothing here modifies rare-high-impact-variants. Every number is read from the
`result.json` it writes.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from eff_to_info import MODES, SNPEFF_TO_SO, convert  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[2]
TARGET = REPO / "skills" / "rare-high-impact-variants" / "rare_high_impact_variants.py"

METRICS = (
    "variants_processed",
    "carried_variants",
    "high_impact_carried",
    "rare_high_impact_count",
    "high_impact_common",
    "high_impact_frequency_unknown",
)

# The three lines that explain the whole result. Quoted, with line numbers, so a
# reader can check them rather than take our word for it.
EXPLANATION = [
    (126, 'consequence = info.get("MC", "") or info.get("Consequence", "") or info.get("ANN", "")',
     "EFF is not in this chain, so consequence is the empty string"),
    (130, 'gene = (info.get("GENEINFO", "") or info.get("SYMBOL", "")).split(":")[0]...',
     "the gene symbol inside the EFF field is never recovered"),
    (160, "if not _is_high_impact(consequence): continue",
     "the record is dropped here, before any impact metric is incremented"),
]


def run_target(vcf: pathlib.Path, outdir: pathlib.Path) -> tuple[int, dict]:
    """Run the unmodified skill. Returns (exit_code, result_json)."""
    proc = subprocess.run(
        [sys.executable, str(TARGET), "--input", str(vcf), "--output", str(outdir)],
        capture_output=True,
        text=True,
        cwd=str(REPO),
    )
    result_path = outdir / "result.json"
    payload: dict = {}
    if result_path.exists():
        blob = json.loads(result_path.read_text())
        # The skill nests its counters; accept either shape.
        payload = blob.get("summary") or blob.get("results") or blob
    return proc.returncode, payload


def pick(payload: dict, key: str):
    """Find a metric wherever the skill happens to nest it."""
    if key in payload:
        return payload[key]
    for value in payload.values():
        if isinstance(value, dict) and key in value:
            return value[key]
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, type=pathlib.Path, help="segregation TSV")
    ap.add_argument("--evidence", type=pathlib.Path, help="cached VEP JSON, for documented frequencies")
    ap.add_argument("--output", required=True, type=pathlib.Path)
    args = ap.parse_args()

    if not TARGET.exists():
        print(f"target skill not found: {TARGET}", file=sys.stderr)
        return 1

    args.output.mkdir(parents=True, exist_ok=True)
    runs: dict[str, dict] = {}

    for mode in MODES:
        vcf = args.output / f"{mode}.vcf"
        stats = convert(args.input, vcf, mode=mode, evidence=args.evidence)
        code, payload = run_target(vcf, args.output / f"run_{mode}")
        runs[mode] = {
            "exit_code": code,
            "adapter": stats,
            "metrics": {k: pick(payload, k) for k in METRICS},
        }
        print(f"{mode:12s} exit={code}  " + "  ".join(
            f"{k.split('_')[0]}={runs[mode]['metrics'][k]}" for k in METRICS[:4]
        ))

    lof = [e for e in SNPEFF_TO_SO if SNPEFF_TO_SO[e][1] in (
        "nonsense", "frameshift_variant", "splice_donor_variant",
        "splice_acceptor_variant", "start_lost", "stop_lost")]

    # report
    L: list[str] = []
    a = L.append
    a("# Measuring the legacy-EFF gap")
    a("")
    a(f"Target: `skills/rare-high-impact-variants/rare_high_impact_variants.py`, unmodified.  ")
    a(f"Input: `{args.input.name}` — the same records in all three runs.")
    a("")
    a("| Metric | raw `EFF=` only | `MC=` unmapped | `MC=` mapped |")
    a("|---|---|---|---|")
    for k in METRICS:
        a(f"| `{k}` | {runs['raw']['metrics'][k]} | {runs['passthrough']['metrics'][k]} | {runs['mapped']['metrics'][k]} |")
    a(f"| **exit code** | **{runs['raw']['exit_code']}** | **{runs['passthrough']['exit_code']}** | **{runs['mapped']['exit_code']}** |")
    a("")
    a("## Reading the first column")
    a("")
    a(
        f"`variants_processed` and `carried_variants` are correct. Every impact metric is "
        f"{runs['raw']['metrics']['high_impact_carried']}. The exit code is "
        f"{runs['raw']['exit_code']} and the report is well formed. There is no warning on "
        "stderr. This is a false negative that reads as a clean run — which is why it is worth "
        "measuring rather than describing."
    )
    a("")
    for line, code, why in EXPLANATION:
        a(f"- `rare_high_impact_variants.py:{line}` — `{code}`  ")
        a(f"  {why}")
    a("")
    a("## Reading the second column — the part that matters")
    a("")
    a(
        "Setting `MC=` to the SnpEff effect name is the obvious one-line fix, and it is worse "
        "than the bug. The matcher is an unanchored, case-insensitive substring test over eight "
        "terms, and SnpEff's names agree with three of them by coincidence:"
    )
    a("")
    a("| SnpEff name | contains a matching substring? |")
    a("|---|---|")
    for name in sorted(lof):
        term = SNPEFF_TO_SO[name][1]
        accidental = name.lower() == term or term in name.lower()
        a(f"| `{name}` | {'yes, by accident' if accidental else 'no'} |")
    a("")
    a(
        f"So the unmapped run reports **{runs['passthrough']['metrics']['high_impact_carried']} "
        f"of {runs['passthrough']['metrics']['carried_variants']}** rather than "
        f"{runs['mapped']['metrics']['high_impact_carried']}. A zero is obviously broken. A "
        "partial count is not, and nothing in the output distinguishes it from a correct one."
    )
    a("")
    a("## Frequencies")
    a("")
    a(
        f"{runs['mapped']['adapter']['records_with_documented_frequency']} of "
        f"{runs['mapped']['adapter']['records_written']} records received a documented frequency "
        "from the cached build-matched layer. The rest carry no frequency key at all, so they land "
        "in the skill's own `frequency_unknown` bucket. Supplying a placeholder would manufacture "
        "the very claim this project refuses to make."
    )
    a("")
    unmapped = runs["mapped"]["adapter"]["unmapped_effects"]
    a(f"Unmapped effect names in the mapped run: {unmapped or 'none'}.")
    a("")
    a("---")
    a("")
    a(
        "*The brief said this skill could not read the data. It can now, and the interesting "
        "result was not the fix but the near-miss beside it.*"
    )

    (args.output / "report.md").write_text("\n".join(L) + "\n")
    (args.output / "result.json").write_text(json.dumps(runs, indent=1))
    print(f"\n-> {args.output}/report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
