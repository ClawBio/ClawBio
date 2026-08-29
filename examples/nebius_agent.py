#!/usr/bin/env python3
"""Drive ClawBio skills from a model hosted on Nebius Token Factory.

ClawBio skills are deterministic Python: they do not need a language model and
they consume no inference credit. This script is the missing half. It gives a
Token Factory model a set of ClawBio skills as callable tools, lets it decide
which to run, and prints the tool calls so the reasoning is auditable.

    export NEBIUS_API_KEY=...            # tokenfactory.nebius.com/project/api-keys
    uv run python examples/nebius_agent.py \
        --ask "Which variants are genuinely rare and high-impact, and which can you not tell?"

Use --dry-run to exercise the tool wiring with no API call and no spend.

Only skills in SKILLS below can be run, and only with the flags listed there.
The model chooses among them; it cannot compose an arbitrary command line.

When run without --dry-run, the report text returned by a ClawBio skill is sent
to Nebius as tool output in the model conversation. The whitelist currently pins
every skill to --demo, so repo-bundled demo data is used and patient input paths
cannot reach the subprocesses.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_BASE_URL = "https://api.tokenfactory.nebius.com/v1/"
# Verified present in the Token Factory catalogue on 2026-08-12. Note that the
# model ID in Nebius's own function-calling docs
# (meta-llama/Meta-Llama-3.1-8B-Instruct-fast) 404s: the docs are stale, so list
# the catalogue rather than trusting them. Use --list-models to see what is live.
DEFAULT_MODEL = "Qwen/Qwen3-235B-A22B-Instruct-2507"

# Whitelist. Every entry was run end to end on 2026-08-12.
# Keep these args demo-only: reports are sent to Nebius, so do not add real
# patient input paths here. `_require_demo_only_whitelist` enforces that at
# import time so a later edit cannot silently add `--input`.
SKILLS: dict[str, dict] = {
    "rare_high_impact_variants": {
        "script": "skills/rare-high-impact-variants/rare_high_impact_variants.py",
        "args": ["--demo"],
        "report": "report.md",
        "description": (
            "Count rare, high-impact loss-of-function variants carried in a VCF. "
            "Reports separately those with NO population frequency data, which cannot "
            "be confirmed rare. Absence of a frequency is not evidence of rarity."
        ),
    },
    "vcf_annotator": {
        "script": "skills/vcf-annotator/vcf_annotator.py",
        "args": ["--demo"],
        "report": "report.md",
        "description": (
            "Annotate VCF variants with Ensembl VEP, ClinVar and gnomAD, ranked by "
            "predicted impact (HIGH/MODERATE/LOW/MODIFIER)."
        ),
    },
    "clinical_variant_reporter": {
        "script": "skills/clinical-variant-reporter/clinical_variant_reporter.py",
        "args": ["--demo"],
        "report": "report.md",
        "description": (
            "Classify germline variants under the ACMG/AMP 2015 28-criteria framework "
            "and return the triggered evidence codes."
        ),
    },
    "gwas_prs": {
        "script": "skills/gwas-prs/gwas_prs.py",
        "args": ["--demo"],
        "report": "prs_report.md",
        "description": (
            "Compute polygenic risk scores from PGS Catalog scores. Note the reference "
            "population stated in the report when interpreting any percentile."
        ),
    },
}


def _require_demo_only_whitelist(skills: dict[str, dict]) -> None:
    """Refuse to start if any whitelist entry can point at patient input."""
    bad: list[str] = []
    for name, spec in skills.items():
        args = list(spec.get("args") or [])
        if "--demo" not in args or "--input" in args:
            bad.append(name)
    if bad:
        raise SystemExit(
            "SKILLS whitelist must pin every entry to --demo and must not "
            f"include --input (reports are sent to Nebius). Offending: {', '.join(bad)}"
        )


_require_demo_only_whitelist(SKILLS)

SYSTEM_PROMPT = """You are a genomics analyst working through ClawBio skills.

Rules you must follow:
1. Do not answer from memory. Call a skill and base every claim on its output.
2. Quote the specific numbers the skill reported. Do not round or embellish them.
3. State plainly what the evidence does NOT support. If variants lack population
   frequency data, say they cannot be confirmed rare rather than assuming they are.
4. If no available skill can answer the question, say so instead of guessing.

Being explicit about what you cannot conclude is worth more than a confident answer."""


def run_skill(name: str) -> str:
    """Run one whitelisted skill and return its report text."""
    _require_demo_only_whitelist(SKILLS)
    spec = SKILLS.get(name)
    if spec is None:
        return f"ERROR: unknown skill {name!r}. Available: {', '.join(SKILLS)}"

    script = REPO_ROOT / spec["script"]
    if not script.exists():
        return f"ERROR: {spec['script']} not found in this checkout"

    with tempfile.TemporaryDirectory() as td:
        cmd = [sys.executable, str(script), *spec["args"], "--output", td]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if proc.returncode != 0:
            tail = proc.stderr.strip().splitlines()[-3:]
            return f"ERROR: {name} exited {proc.returncode}\n" + "\n".join(tail)

        report = Path(td) / spec["report"]
        if report.exists():
            return report.read_text()[:6000]
        return proc.stdout[-4000:] or f"{name} produced no report at {spec['report']}"


def tool_schemas() -> list[dict]:
    return [{
        "type": "function",
        "function": {
            "name": "run_clawbio_skill",
            "description": (
                "Run one ClawBio bioinformatics skill on its bundled demo data and "
                "return the report it produces. Available skills:\n" +
                "\n".join(f"- {n}: {s['description']}" for n, s in SKILLS.items())
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "skill": {"type": "string", "enum": list(SKILLS)},
                },
                "required": ["skill"],
            },
        },
    }]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ask", default=(
        "Which variants in this genome are genuinely rare and high-impact, "
        "and which ones can you not tell either way?"))
    ap.add_argument("--model", default=DEFAULT_MODEL)
    ap.add_argument("--base-url", default=os.environ.get("NEBIUS_BASE_URL", DEFAULT_BASE_URL))
    ap.add_argument("--max-turns", type=int, default=6)
    ap.add_argument("--max-tokens", type=int, default=1500,
                    help="Per-response cap. Reasoning models need headroom: GLM-5.2 "
                         "returns an empty answer if it spends the budget thinking")
    ap.add_argument("--dry-run", action="store_true",
                    help="Exercise the tool wiring with no API call and no spend")
    ap.add_argument("--list-models", action="store_true",
                    help="Print the live Token Factory catalogue and exit")
    args = ap.parse_args()

    if args.list_models:
        from openai import OpenAI
        key = os.environ.get("NEBIUS_API_KEY")
        if not key:
            sys.exit("NEBIUS_API_KEY is not set")
        client = OpenAI(base_url=args.base_url, api_key=key)
        for mid in sorted(m.id for m in client.models.list().data):
            print(mid)
        return

    if args.dry_run:
        print("DRY RUN: no API call, no credit spent.\n")
        print(f"Base URL : {args.base_url}")
        print(f"Model    : {args.model}")
        print(f"Skills   : {', '.join(SKILLS)}\n")
        schema = tool_schemas()[0]["function"]
        print(f"Tool exposed to the model: {schema['name']}")
        print(f"Enum offered            : {schema['parameters']['properties']['skill']['enum']}\n")
        probe = "rare_high_impact_variants"
        print(f"Calling {probe} locally to prove dispatch works:\n")
        out = run_skill(probe)
        print("\n".join(out.splitlines()[:14]))
        print("\nDry run complete. Set NEBIUS_API_KEY and drop --dry-run to go live.")
        return

    api_key = os.environ.get("NEBIUS_API_KEY")
    if not api_key:
        sys.exit("NEBIUS_API_KEY is not set. Get one at "
                 "https://tokenfactory.nebius.com/project/api-keys")

    try:
        from openai import OpenAI
    except ImportError:
        sys.exit("pip install openai")

    client = OpenAI(base_url=args.base_url, api_key=api_key)
    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": args.ask},
    ]

    print(f"Model    : {args.model}")
    print(f"Question : {args.ask}\n")

    calls = 0
    for turn in range(args.max_turns):
        response = client.chat.completions.create(
            model=args.model, messages=messages,
            tools=tool_schemas(), tool_choice="auto", temperature=0.2,
            max_tokens=args.max_tokens,
        )
        choice = response.choices[0]
        msg = choice.message

        if choice.finish_reason == "length" and not msg.content and not msg.tool_calls:
            sys.exit(
                f"{args.model} used all {args.max_tokens} tokens without producing an "
                f"answer. Reasoning models do this. Raise --max-tokens or pick another "
                f"model (--list-models)."
            )
        messages.append(msg.model_dump(exclude_none=True))

        if not msg.tool_calls:
            print("--- ANSWER ---")
            print(msg.content)
            usage = getattr(response, "usage", None)
            if usage:
                print(f"\nTokens: {usage.prompt_tokens} in, "
                      f"{usage.completion_tokens} out, over {calls} skill call(s).")
            return

        for call in msg.tool_calls:
            skill = json.loads(call.function.arguments).get("skill", "")
            calls += 1
            print(f"[tool call {calls}] run_clawbio_skill(skill={skill!r})")
            result = run_skill(skill)
            first = next((l for l in result.splitlines() if l.strip()), "")
            print(f"             -> {len(result)} chars, starting {first[:60]!r}\n")
            messages.append({
                "role": "tool", "tool_call_id": call.id, "content": result,
            })

    print(f"Stopped after {args.max_turns} turns without a final answer.")


if __name__ == "__main__":
    main()
