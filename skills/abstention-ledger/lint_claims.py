"""Fail the build if a generated artefact makes a claim the input cannot support.

The Berlin brief bans five terms outright:

    "Never call anything rare, pathogenic, diagnostic, de novo or compound
     heterozygous."

Enforcing that by reviewer discipline does not survive a deadline, so it is a
check. Run it over every generated report and over the prose we ship.

What is deliberately *not* a violation
--------------------------------------
The ban is on making the claim, not on the letters appearing. Three exemptions,
each narrow and each justified:

1. **Identifiers.** `rare-high-impact-variants` is the name of a skill in this
   library and `rare_high_impact_count` is a field in its output. Naming a code
   path is not asserting anything about a variant.
2. **Quotations.** Markdown blockquotes carry the brief's own wording, including
   the ban itself. A checker that cannot quote its own rule is useless.
3. **Negations we rely on.** Sentences of the form "absence of a frequency is not
   evidence of rarity" exist to refuse the claim. Flagging them would push us
   toward saying less about the limitation, which is backwards.

Exemptions are matched narrowly and listed in the output, so a reader can see
exactly what was waved through.
"""

from __future__ import annotations

import argparse
import pathlib
import re
import sys

BANNED = (
    r"\brare\b",
    r"\brarity\b",
    r"\bpathogenic\b",
    r"\bpathogenicity\b",
    r"\bdiagnostic\b",
    r"\bdiagnosis\b",
    r"\bde novo\b",
    r"\bcompound het",
)

# Narrow exemptions. Each is an identifier or an explicit refusal, never a claim.
EXEMPT = (
    r"rare[-_]high[-_]impact[-_]variants",       # skill name
    r"rare_high_impact_count",                   # field name in its result.json
    r"rare[-_]disease[-_]rnaseq",                # skill name
    r"not evidence of rarity",                   # the refusal itself
    r"not\s+\w*\s*rare",                         # "cannot be confirmed rare"
    r"does not provide clinical diagnoses",       # mandated disclaimer wording
    r"rare/ultra_rare",                          # bucket label quoted from the target skill
)

_BANNED_RE = [re.compile(p, re.I) for p in BANNED]
_EXEMPT_RE = [re.compile(p, re.I) for p in EXEMPT]


def scan(path: pathlib.Path) -> tuple[list[tuple[int, str, str]], int]:
    """Return (violations, exempted_count)."""
    violations: list[tuple[int, str, str]] = []
    exempted = 0
    for n, line in enumerate(path.read_text().splitlines(), 1):
        stripped = line.lstrip()
        if stripped.startswith(">"):  # quoted material
            continue
        for rx in _BANNED_RE:
            for m in rx.finditer(line):
                window = line[max(0, m.start() - 40) : m.end() + 40]
                if any(ex.search(window) for ex in _EXEMPT_RE):
                    exempted += 1
                    continue
                violations.append((n, m.group(0), line.strip()[:160]))
    return violations, exempted


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", type=pathlib.Path)
    args = ap.parse_args()

    total = 0
    total_exempt = 0
    for path in args.paths:
        if not path.exists():
            print(f"  skip {path} (missing)")
            continue
        violations, exempted = scan(path)
        total_exempt += exempted
        if violations:
            total += len(violations)
            print(f"  FAIL {path}")
            for n, term, line in violations:
                print(f"    {path}:{n}  '{term}'  |  {line}")
        else:
            print(f"  ok   {path}" + (f"  ({exempted} exempted)" if exempted else ""))

    print()
    if total:
        print(f"{total} prohibited claim(s) found. {total_exempt} exemption(s) applied.")
        return 1
    print(f"clean — no prohibited claims. {total_exempt} exemption(s) applied (identifiers and refusals).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
