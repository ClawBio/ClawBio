#!/usr/bin/env python3
"""ClawBio Skill Scaffold - generate benchmark-ready skill directories from a one-line spec.

Usage:
    python scaffold_skill.py hla-typing "HLA allele typing from WGS/WES VCF data"
    python scaffold_skill.py --list-existing
"""

import argparse
import json
import sys
import textwrap
from pathlib import Path

CLAWBIO_ROOT = Path(__file__).resolve().parent
SKILLS_DIR = CLAWBIO_ROOT / "skills"
TEMPLATE_PATH = CLAWBIO_ROOT / "templates" / "SKILL-TEMPLATE.md"

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Scaffold a new ClawBio skill directory with all required files."
    )
    parser.add_argument("name", nargs="?", help="Skill name (lowercase-hyphenated, e.g. hla-typing)")
    parser.add_argument("description", nargs="?", help="One-line skill description")
    parser.add_argument("--list-existing", action="store_true", help="List existing skills and exit")
    parser.add_argument("--force", action="store_true", help="Overwrite existing skill directory")
    parser.add_argument("--bench-dir", type=Path, default=None,
                        help="Path to clawbio_bench test_cases/ dir (optional)")
    return parser.parse_args()


def to_python_name(name: str) -> str:
    """hla-typing -> hla_typing"""
    return name.replace("-", "_")


def to_title(name: str) -> str:
    """hla-typing -> Hla Typing"""
    return name.replace("-", " ").title()


def to_prefix(name: str) -> str:
    """hla-typing -> ht, pharmgx-reporter -> pr, equity-scorer -> es"""
    parts = name.split("-")
    if len(parts) == 1:
        return parts[0][:2]
    return "".join(p[0] for p in parts[:3])


def generate_skill_md(name: str, description: str) -> str:
    """Generate a conformant SKILL.md from the description."""
    py_name = to_python_name(name)
    title = to_title(name)

    # Derive trigger keywords from description
    words = [w.lower() for w in description.split() if len(w) > 3]
    keywords = words[:3] if len(words) >= 3 else words + ["genomics", "bioinformatics", "analysis"]
    keywords = keywords[:5]

    return textwrap.dedent(f"""\
---
name: {name}
description: >-
  {description}
version: 0.1.0
author: Manuel Corpas
domain: genomics
license: MIT

inputs:
  - name: input_file
    type: file
    format: [vcf, csv, tsv, txt]
    description: Primary input data file
    required: true

outputs:
  - name: report
    type: file
    format: md
    description: Analysis report
  - name: result
    type: file
    format: json
    description: Machine-readable results

dependencies:
  python: ">=3.11"
  packages:
    - pandas>=2.0

tags: [{', '.join(keywords)}]

demo_data:
  - path: demo_input.txt
    description: Synthetic test data

endpoints:
  cli: python skills/{name}/{py_name}.py --input {{input_file}} --output {{output_dir}}

metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    homepage: https://github.com/ClawBio/ClawBio
    os: [macos, linux]
    install:
      - kind: pip
        package: pandas
        bins: []
    trigger_keywords:
      - {keywords[0]}
      - {keywords[1] if len(keywords) > 1 else name}
      - {keywords[2] if len(keywords) > 2 else 'analyze ' + name}
---

# {title}

You are **{title}**, a specialised ClawBio agent for genomics. Your role is to {description.lower().rstrip('.')}.

## Trigger

**Fire this skill when the user says any of:**
- "{description.lower()}"
- "run {name}"
- "{' '.join(keywords[:2])}"
- "analyze {keywords[0]}"

**Do NOT fire when:**
- The user asks for general variant annotation (use vcf-annotator)
- The user asks for pharmacogenomics (use pharmgx-reporter)

**Design notes:** The trigger must be loud, not subtle. Models skip subdued
descriptions. Use exact phrases, domain-specific terms, and multiple synonyms.

## Why This Exists

- **Without it**: Users must manually {description.lower().rstrip('.')} using command-line tools and custom scripts
- **With it**: Automated analysis in seconds with a structured, reproducible report
- **Why ClawBio**: Grounded in real databases and algorithms, not LLM guessing

## Core Capabilities

1. **Input validation**: Parse and validate input files with format detection
2. **Analysis**: {description}
3. **Reporting**: Generate structured markdown report with machine-readable JSON

## Scope

**One skill, one task.** This skill does {description.lower().rstrip('.')} and nothing else.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| VCF | `.vcf` | CHROM, POS, REF, ALT, GT | `demo_input.txt` |
| TSV | `.tsv` | variant columns | `sample.tsv` |

## Workflow

When the user asks for {title.lower()}:

1. **Validate**: Check input format and required fields
2. **Parse**: Extract relevant variants and annotations
3. **Analyze**: Apply {title.lower()} algorithm
4. **Generate**: Write result.json with structured findings
5. **Report**: Write report.md with findings, tables, and disclaimer

**Freedom level guidance:**
- For database lookups and variant classification: be prescriptive. Every step must be exact.
- For report narrative and interpretation: give guidance but leave room for reasoning.

## CLI Reference

```bash
# Standard usage
python skills/{name}/{py_name}.py \\
  --input <input_file> --output <report_dir>

# Demo mode (synthetic data, no user files needed)
python skills/{name}/{py_name}.py --demo --output /tmp/{py_name}_demo

# Via ClawBio runner
python clawbio.py run {name} --input <file> --output <dir>
python clawbio.py run {name} --demo
```

## Demo

To verify the skill works:

```bash
python clawbio.py run {name} --demo
```

Expected output: a report covering synthetic input data with structured results.

## Algorithm / Methodology

1. **Parse input**: Read VCF/TSV and extract relevant loci
2. **Lookup**: Query reference databases for annotations
3. **Score**: Apply scoring algorithm to classify findings
4. **Report**: Generate structured output

**Key thresholds / parameters**:
- TODO: define thresholds with citations

## Example Queries

- "{description.lower()}"
- "run {name} on my VCF"
- "analyze my sample with {name}"

## Example Output

```markdown
# {title} Report

**Input**: demo_input.txt (5 variants)
**Date**: 2026-04-06

| Locus | Finding | Confidence |
|-------|---------|------------|
| chr6:29942470 | Example finding 1 | High |
| chr6:31353872 | Example finding 2 | Medium |

## Summary
Analysis completed on 5 variants. 2 findings reported.

*{DISCLAIMER}*
```

## Output Structure

```
output_directory/
├── report.md              # Primary markdown report
├── result.json            # Machine-readable results
├── tables/
│   └── results.csv        # Tabular data (optional)
└── reproducibility/
    ├── commands.sh         # Exact commands to reproduce
    ├── environment.yml     # Environment snapshot
    └── checksums.sha256    # SHA-256 of outputs, relative to output_directory
```

Verify the run with
`cd <output_directory> && sha256sum -c reproducibility/checksums.sha256`.

## Dependencies

**Required**:
- `pandas` >= 2.0; data manipulation

**Optional**:
- `biopython`; sequence handling (graceful degradation without it)

## Gotchas

- **Gotcha 1**: The model tends to infer results from gene names alone. Instead, always require actual genotype data from the input file. Why: inferred results are unreliable and clinically dangerous.
- **Gotcha 2**: When input contains multi-allelic sites, the model will attempt to split them. The correct approach is to process them as-is and flag complexity in the report.
- **Gotcha 3**: Empty or malformed VCF lines cause silent failures. Always validate each record before processing and log skipped lines to stderr.

## Safety

- **Local-first**: No data upload without explicit consent
- **Disclaimer**: Every report includes: *"{DISCLAIMER}"*
- **Audit trail**: Log all operations to reproducibility bundle
- **No hallucinated science**: All parameters trace to cited databases

## Agent Boundary

The agent (LLM) dispatches and explains. The skill (Python) executes.
The agent must NOT override thresholds or invent associations.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- User mentions {keywords[0]} or {name}
- Input file contains relevant loci

**Chaining partners**: this skill connects with:
- `pharmgx-reporter`: downstream pharmacogenomic implications
- `profile-report`: feeds into unified patient profile

## Maintenance

- **Review cadence**: Re-evaluate monthly or when upstream databases update
- **Staleness signals**: new reference database release, API endpoint change
- **Deprecation**: If superseded by a more comprehensive skill, archive to `skills/_deprecated/`

## Citations

- TODO: Add relevant database and paper citations
""")


def generate_python_script(name: str, description: str) -> str:
    """Generate the main Python script skeleton."""
    py_name = to_python_name(name)
    title = to_title(name)

    return textwrap.dedent(f'''\
#!/usr/bin/env python3
"""{title} - {description}."""

import argparse
import json
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = SKILL_DIR.parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

{f'DISCLAIMER = ("{DISCLAIMER}")'}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, dest="input_file", help="Input file path")
    parser.add_argument("--output", type=Path, help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run with synthetic demo data")
    return parser.parse_args()


def validate_input(input_path: Path) -> dict:
    """Validate and parse the input file. Returns parsed data dict."""
    if not input_path.exists():
        print(f"Error: input file not found: {{input_path}}", file=sys.stderr)
        sys.exit(1)
    lines = []
    with open(input_path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                lines.append(line)
    return {{"lines": lines, "source": str(input_path)}}


def run_analysis(data: dict) -> dict:
    """Core analysis logic. Returns result dict."""
    # TODO: implement core {name} logic
    return {{
        "skill": "{name}",
        "version": "0.1.0",
        "source": data.get("source", "unknown"),
        "variants_processed": len(data.get("lines", [])),
        "findings": [],
        "status": "skeleton"
    }}


def write_reproducibility_bundle(output_dir: Path, invocation: str) -> None:
    """Write output_dir/reproducibility/ via the shared clawbio.common layer.

    Do not hand-roll these writers: the shared helpers keep line endings, the
    checksum format and the bundle layout identical across every ClawBio skill.
    """
    write_commands_sh(
        output_dir,
        "# {title} - reproducibility\\n"
        "set -euo pipefail\\n"
        "\\n"
        "# 1. Recreate the environment\\n"
        "conda env create -f environment.yml\\n"
        "conda activate {name}\\n"
        "\\n"
        "# 2. Re-run the analysis\\n"
        f"python {py_name}.py {{invocation}}\\n"
        "\\n"
        "# 3. Verify output checksums (labels are relative to the output directory)\\n"
        '( cd "$(dirname "$0")/.." && sha256sum -c reproducibility/checksums.sha256 )',
    )
    write_environment_yml(
        output_dir,
        env_name="{name}",
        python_version="3.11",
        conda_deps=["pandas>=2.0"],
        pip_deps=[],
    )
    # anchor=output_dir keeps labels relative, so `cd <output_dir> &&
    # sha256sum -c reproducibility/checksums.sha256` resolves every entry.
    write_checksums(
        [output_dir / "report.md", output_dir / "result.json"],
        output_dir,
        anchor=output_dir,
    )


def write_report(result: dict, output_dir: Path, invocation: str = "--demo") -> None:
    """Write report.md, result.json and the reproducibility bundle to output_dir."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # result.json
    with open(output_dir / "result.json", "w") as f:
        json.dump(result, f, indent=2)

    # report.md
    n = result.get("variants_processed", 0)
    findings = result.get("findings", [])
    report = [
        "# {title} Report",
        "",
        f"**Input**: {{result.get('source', 'unknown')}}",
        f"**Variants processed**: {{n}}",
        f"**Findings**: {{len(findings)}}",
        "",
        "## Results",
        "",
        "| Locus | Finding | Confidence |",
        "|-------|---------|------------|",
    ]
    for f_ in findings:
        report.append(f"| {{f_.get('locus', '-')}} | {{f_.get('finding', '-')}} | {{f_.get('confidence', '-')}} |")
    if not findings:
        report.append("| - | No findings (skeleton implementation) | - |")
    report.extend([
        "",
        "## Summary",
        "",
        f"Analysis completed on {{n}} variants. {{len(findings)}} findings reported.",
        "",
        f"*{{DISCLAIMER}}*",
        "",
    ])
    with open(output_dir / "report.md", "w") as f:
        f.write("\\n".join(report))

    # Written last: checksums must cover the finished report and result.
    write_reproducibility_bundle(output_dir, invocation)

    print(f"Report written to {{output_dir / 'report.md'}}")
    print(f"Results written to {{output_dir / 'result.json'}}")
    print(f"Reproducibility bundle written to {{output_dir / 'reproducibility'}}")


def run_demo(output_dir: Path) -> None:
    """Run with built-in synthetic demo data."""
    demo_input = SKILL_DIR / "demo_input.txt"
    if not demo_input.exists():
        print("Error: demo data not found", file=sys.stderr)
        sys.exit(1)
    data = validate_input(demo_input)
    result = run_analysis(data)
    write_report(result, output_dir, invocation=f"--demo --output {{output_dir}}")


def main():
    args = parse_args()
    if args.demo:
        output = args.output or Path("/tmp") / "{py_name}" / "demo"
        run_demo(output)
    elif args.input_file:
        data = validate_input(args.input_file)
        result = run_analysis(data)
        output = args.output or args.input_file.parent / "output"
        write_report(
            result, output,
            invocation=f"--input {{args.input_file}} --output {{output}}",
        )
    else:
        print("Error: provide --input <file> or --demo", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
''')


def generate_api(name: str) -> str:
    """Generate api.py wrapper."""
    py_name = to_python_name(name)
    return textwrap.dedent(f'''\
"""API entry point for skill registry and orchestrator integration."""

from pathlib import Path
from importlib.util import spec_from_file_location, module_from_spec

_mod_path = Path(__file__).parent / "{py_name}.py"
_spec = spec_from_file_location("{py_name}", _mod_path)
_mod = module_from_spec(_spec)
_spec.loader.exec_module(_mod)


def run(input_path: str, output_dir: str = "/tmp/{name}") -> dict:
    """Run the skill programmatically. Returns result dict."""
    data = _mod.validate_input(Path(input_path))
    result = _mod.run_analysis(data)
    _mod.write_report(result, Path(output_dir))
    return result
''')


def generate_demo_data(name: str, description: str) -> str:
    """Generate synthetic demo input data."""
    return textwrap.dedent(f"""\
# Synthetic demo data for {name}
# Format: TSV with header
# This is NOT real patient data
# Description: {description}
#
# chrom\tpos\tref\talt\tgenotype
chr6\t29942470\tA\tG\t0/1
chr6\t31353872\tC\tT\t1/1
chr6\t32609227\tG\tA\t0/1
chr6\t32631934\tT\tC\t0/0
chr6\t32659879\tA\tC\t0/1
""")


# Appended verbatim to every generated test file. No per-skill substitution:
# it reads the skill's own SKILL.md generically. Kept as a plain (non-f) string
# so its braces and regex escapes pass through untouched.
_OUTPUT_CONTRACT_BLOCK = '''

def _parse_output_contract(skill_md):
    """Extract files promised in the SKILL.md '## Output Structure' tree.

    Returns output-relative file paths. Directory lines, and any entry whose
    inline comment contains 'optional', are skipped. Returns [] when there is no
    parseable section, so skills without the section are simply not gated.
    """
    if not skill_md.exists():
        return []
    text = skill_md.read_text()
    m = re.search(r"##\\s*Output Structure\\s*\\n+```[^\\n]*\\n(.*?)\\n```", text, re.S)
    if not m:
        return []
    files = []
    parents = {}
    for raw in m.group(1).splitlines():
        if not raw.strip():
            continue
        parts = re.split(r"\\s+#", raw, maxsplit=1)
        entry, comment = parts[0], (parts[1] if len(parts) > 1 else "")
        mm = re.match(r"^([\\s│├└─]*)(.*)$", entry)
        prefix, name = mm.group(1), mm.group(2).strip()
        if not name:
            continue
        depth = len(prefix) // 4
        if depth == 0:
            continue  # the root output_directory/ line
        if name.endswith("/"):
            parents[depth] = name.rstrip("/")
            for d in [k for k in parents if k > depth]:
                del parents[d]
            continue
        if "optional" in comment.lower():
            continue
        rel = "/".join(parents[d] for d in sorted(parents) if d < depth)
        files.append(rel + "/" + name if rel else name)
    return files


class TestOutputContract:
    """Every artifact promised in SKILL.md '## Output Structure' must be produced.

    Guards against doc/code drift: a skill documenting an output it never writes.
    Mark conditional artifacts with '(optional)' in the tree comment to exempt them.
    """

    def test_documented_outputs_are_produced(self, tmp_path):
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        if not promised:
            pytest.skip("No parseable '## Output Structure' section in SKILL.md")
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        assert result.returncode == 0, f"demo run failed: {result.stderr}"
        missing = [p for p in promised if not (tmp_path / p).exists()]
        assert not missing, (
            "SKILL.md Output Structure promises artifacts the skill did not "
            "produce: " + ", ".join(missing) + ". Write them, mark them "
            "'(optional)' in the SKILL.md tree, or remove them from the "
            "documented Output Structure."
        )
'''


def generate_tests(name: str) -> str:
    """Generate test suite."""
    py_name = to_python_name(name)
    header = textwrap.dedent(f'''\
"""Tests for {name}. Red/green TDD: these should fail until implementation is complete."""

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "{py_name}.py"
DEMO_INPUT = SKILL_DIR / "demo_input.txt"


class TestCLI:
    """CLI interface tests."""

    def test_no_args_exits_nonzero(self):
        result = subprocess.run(
            [sys.executable, str(SCRIPT)],
            capture_output=True, text=True
        )
        assert result.returncode != 0

    def test_demo_mode_produces_output(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        assert result.returncode == 0, f"stderr: {{result.stderr}}"
        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "result.json").exists()

    def test_input_mode_produces_output(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT),
             "--input", str(DEMO_INPUT),
             "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        assert result.returncode == 0, f"stderr: {{result.stderr}}"
        assert (tmp_path / "report.md").exists()

    def test_missing_input_exits_nonzero(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT),
             "--input", str(tmp_path / "nonexistent.txt"),
             "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        assert result.returncode != 0


class TestOutputFormat:
    """Output format validation."""

    def test_result_json_is_valid(self, tmp_path):
        subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        result = json.loads((tmp_path / "result.json").read_text())
        assert isinstance(result, dict)
        assert "skill" in result
        assert result["skill"] == "{name}"

    def test_report_contains_disclaimer(self, tmp_path):
        subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        report = (tmp_path / "report.md").read_text()
        assert "not a medical device" in report.lower()

    def test_result_has_variants_count(self, tmp_path):
        subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        result = json.loads((tmp_path / "result.json").read_text())
        assert "variants_processed" in result
        assert result["variants_processed"] > 0


class TestDemoData:
    """Demo data integrity."""

    def test_demo_input_exists(self):
        assert DEMO_INPUT.exists(), f"Demo data missing: {{DEMO_INPUT}}"

    def test_demo_input_has_content(self):
        content = DEMO_INPUT.read_text()
        lines = [l for l in content.splitlines() if l.strip() and not l.startswith("#")]
        assert len(lines) > 0, "Demo input has no data lines"
''')
    return header + _OUTPUT_CONTRACT_BLOCK


def generate_bench_test_cases(name: str, description: str) -> dict[str, dict[str, str]]:
    """Generate clawbio_bench-compatible test case directories.

    Returns: {dir_name: {filename: content}}
    """
    prefix = to_prefix(name)
    cases = {}

    # Case 1: happy path
    cases[f"{prefix}_01_basic_correct"] = {
        "ground_truth.txt": textwrap.dedent(f"""\
# BENCHMARK: {name} v0.1.0
# PAYLOAD: input.txt
# GROUND_TRUTH_COUNT: 5
# COUNT_TOLERANCE: 0
# FINDING_CATEGORY: correct
# FINDING: {prefix.upper()}-01 basic correct output
# HAZARD_METRIC: Tool produces correct variant count from clean input
# DERIVATION: 5 non-comment, non-empty lines in input file
# CITATION: Synthetic test data
"""),
        "input.txt": textwrap.dedent(f"""\
# Synthetic test input for {name}
chr6\t29942470\tA\tG\t0/1
chr6\t31353872\tC\tT\t1/1
chr6\t32609227\tG\tA\t0/1
chr6\t32631934\tT\tC\t0/0
chr6\t32659879\tA\tC\t0/1
"""),
    }

    # Case 2: edge case - empty input
    cases[f"{prefix}_02_empty_input"] = {
        "ground_truth.txt": textwrap.dedent(f"""\
# BENCHMARK: {name} v0.1.0
# PAYLOAD: input.txt
# GROUND_TRUTH_COUNT: 0
# FINDING_CATEGORY: correct_with_warnings
# FINDING: {prefix.upper()}-02 empty input handled gracefully
# HAZARD_METRIC: Tool handles empty input without crashing
# DERIVATION: Input file contains only comments and blank lines, 0 data lines
# CITATION: Synthetic test data
"""),
        "input.txt": textwrap.dedent(f"""\
# Synthetic empty test input for {name}
# No data lines follow
"""),
    }

    # Case 3: malformed input
    cases[f"{prefix}_03_malformed_input"] = {
        "ground_truth.txt": textwrap.dedent(f"""\
# BENCHMARK: {name} v0.1.0
# PAYLOAD: input.txt
# EXPECTED_EXIT_CODE: 1
# FINDING_CATEGORY: crash
# FINDING: {prefix.upper()}-03 malformed input causes crash
# HAZARD_METRIC: Tool should exit non-zero on malformed input with clear error message
# DERIVATION: Input contains lines with wrong column count; tool should reject
# CITATION: Synthetic test data
"""),
        "input.txt": textwrap.dedent("""\
# Malformed input: wrong number of columns, mixed delimiters
chr6\t29942470
chr6 31353872 C T 1/1
\t\t\t\t
chr6\t32609227\tG\t\t0/1
"""),
    }

    return cases


def scaffold(name: str, description: str, force: bool = False,
             bench_dir: Path | None = None) -> None:
    """Generate all skill files."""
    skill_dir = SKILLS_DIR / name
    py_name = to_python_name(name)

    if skill_dir.exists() and not force:
        print(f"Error: {skill_dir} already exists. Use --force to overwrite.", file=sys.stderr)
        sys.exit(1)

    # Create directories
    skill_dir.mkdir(parents=True, exist_ok=True)
    (skill_dir / "tests").mkdir(exist_ok=True)
    (skill_dir / "examples").mkdir(exist_ok=True)

    # Write files
    files = {
        "SKILL.md": generate_skill_md(name, description),
        f"{py_name}.py": generate_python_script(name, description),
        "api.py": generate_api(name),
        "demo_input.txt": generate_demo_data(name, description),
        f"tests/test_{py_name}.py": generate_tests(name),
    }

    for filename, content in files.items():
        path = skill_dir / filename
        path.write_text(content)
        print(f"  Created: skills/{name}/{filename}")

    # Benchmark test cases
    bench_cases = generate_bench_test_cases(name, description)
    if bench_dir:
        bench_skill_dir = bench_dir / name
    else:
        bench_skill_dir = skill_dir / "bench_test_cases"

    bench_skill_dir.mkdir(parents=True, exist_ok=True)
    for case_name, case_files in bench_cases.items():
        case_dir = bench_skill_dir / case_name
        case_dir.mkdir(parents=True, exist_ok=True)
        for fname, fcontent in case_files.items():
            (case_dir / fname).write_text(fcontent)
        print(f"  Created: bench test case {case_name}/")

    # Line count check
    skill_md_lines = len((skill_dir / "SKILL.md").read_text().splitlines())

    # Print conformance checklist
    print(f"\n{'='*60}")
    print(f"SKILL.md Conformance Checklist for: {name}")
    print(f"{'='*60}")

    skill_md = (skill_dir / "SKILL.md").read_text()
    checks = [
        (f"YAML: name matches folder ({name})", f"name: {name}" in skill_md),
        ("YAML: version semver", "version: 0.1.0" in skill_md),
        ("YAML: author present", "author: Manuel Corpas" in skill_md),
        ("YAML: description one line", "description:" in skill_md),
        ("YAML: inputs with format + required", "required: true" in skill_md),
        ("YAML: outputs with format", "format: md" in skill_md and "format: json" in skill_md),
        ("YAML: trigger_keywords >= 3", skill_md.count("      - ") >= 3),
        ("Section: ## Trigger", "## Trigger" in skill_md),
        ("Section: ## Scope", "## Scope" in skill_md),
        ("Section: ## Workflow (numbered)", "1. **" in skill_md),
        ("Section: ## Example Output", "## Example Output" in skill_md),
        ("Section: ## Gotchas >= 3", skill_md.count("**Gotcha") >= 3),
        ("Section: ## Safety + disclaimer", "not a medical device" in skill_md),
        ("Section: ## Agent Boundary", "## Agent Boundary" in skill_md),
        ("File: demo data exists", (skill_dir / "demo_input.txt").exists()),
        ("File: tests/ with test file", (skill_dir / "tests" / f"test_{py_name}.py").exists()),
        ("File: reproducibility bundle via clawbio.common",
         "clawbio.common.reproducibility" in (skill_dir / f"{py_name}.py").read_text()
         and "reproducibility/" in skill_md),
        (f"Line count: {skill_md_lines} < 500", skill_md_lines < 500),
    ]

    all_pass = True
    for label, passed in checks:
        status = "PASS" if passed else "FAIL"
        if not passed:
            all_pass = False
        print(f"  [{status}] {label}")

    print(f"\n{'='*60}")
    print("Benchmark Readiness")
    print(f"{'='*60}")

    bench_checks = [
        (f">= 3 test cases ({len(bench_cases)})", len(bench_cases) >= 3),
        ("PAYLOAD field in ground_truth.txt", True),
        ("FINDING_CATEGORY defined", True),
        ("DERIVATION present", True),
        ("Pass/fail categories defined", True),
    ]
    for label, passed in bench_checks:
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {label}")

    # Print registry and routing entries (not auto-written)
    print(f"\n{'='*60}")
    print("Manual integration (copy-paste into the right files):")
    print(f"{'='*60}")

    registry_entry = {
        name: {
            "description": description,
            "gate_traits": {
                "analytical_reasoning": 0.80,
                "pattern_recognition": 0.80,
            },
            "config_traits": {},
            "module_path": f"skills/{name}/api.py",
            "function": "run",
            "input_type": "vcf_path",
            "cost_tier": "medium",
            "returns": "Analysis report with findings",
        }
    }

    print(f"\n--- skill_registry.json entry ---")
    print(json.dumps(registry_entry, indent=2))

    print(f"\n--- CLAUDE.md routing table row ---")
    print(f"| {description}, {name} | `skills/{name}/` | Run `{py_name}.py` |")

    if bench_dir:
        print(f"\nBenchmark test cases written to: {bench_skill_dir}")
    else:
        print(f"\nBenchmark test cases at: skills/{name}/bench_test_cases/")
        print(f"  Move to clawbio_bench/src/clawbio_bench/test_cases/{name}/ when ready.")

    print(f"\nDone. Skill scaffolded at: skills/{name}/")
    if all_pass:
        print("All 17 conformance checks PASSED.")
    else:
        print("WARNING: Some conformance checks failed. Review above.")


def list_existing():
    """List existing skills."""
    if not SKILLS_DIR.exists():
        print("No skills directory found.")
        return
    skills = sorted(d.name for d in SKILLS_DIR.iterdir()
                    if d.is_dir() and (d / "SKILL.md").exists())
    print(f"Existing skills ({len(skills)}):")
    for s in skills:
        print(f"  {s}")


def main():
    args = parse_args()

    if args.list_existing:
        list_existing()
        return

    if not args.name or not args.description:
        print("Error: provide skill name and description.", file=sys.stderr)
        print('Usage: python scaffold_skill.py hla-typing "HLA allele typing from WGS/WES VCF data"',
              file=sys.stderr)
        sys.exit(1)

    # Validate name
    if not all(c.isalnum() or c == "-" for c in args.name):
        print("Error: skill name must be lowercase-hyphenated (e.g. hla-typing)", file=sys.stderr)
        sys.exit(1)

    print(f"Scaffolding skill: {args.name}")
    print(f"Description: {args.description}")
    print()

    scaffold(args.name, args.description, force=args.force, bench_dir=args.bench_dir)


if __name__ == "__main__":
    main()
