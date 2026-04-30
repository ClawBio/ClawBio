#!/usr/bin/env python3
"""
skill_builder.py — ClawBio Skill Builder
=========================================
Scaffold a new ClawBio skill from a spec file (JSON/YAML) or interactively.

Usage:
    # Spec-driven (JSON — no extra deps)
    python skill_builder.py --input spec.json

    # Spec-driven (YAML — requires pyyaml)
    python skill_builder.py --input spec.yaml

    # Interactive mode
    python skill_builder.py --interactive

    # Demo (scaffolds hello-bioinformatics)
    python skill_builder.py --demo --output /tmp/skill_builder_demo

    # Validate an existing SKILL.md
    python skill_builder.py --validate-only --input skills/my-skill/SKILL.md

    # Dry run — print without writing
    python skill_builder.py --input spec.json --dry-run

    # Via ClawBio runner
    python clawbio.py run skill-builder --demo
    python clawbio.py run skill-builder --input spec.json

    # Agent / LLM pipeline — JSON manifest to stdout, no interactive prompts
    python skill_builder.py --agent --input spec.json
    python skill_builder.py --agent --input spec.json --dry-run
    echo '{"name":"foo",...}' | python skill_builder.py --agent --input -
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import textwrap
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

VERSION = "0.1.0"
SCRIPT_DIR = Path(__file__).resolve().parent
DEMO_SPEC_PATH = SCRIPT_DIR / "demo_spec.json"

# Sentinel file that marks the ClawBio repo root
REPO_SENTINEL = "clawbio.py"

# Sections required in every SKILL.md body (from CONTRIBUTING.md checklist)
REQUIRED_SECTIONS = [
    ("YAML frontmatter",             r"^---"),
    ("trigger_keywords",             r"trigger_keywords"),
    ("Why This Exists",              r"##\s+Why This Exists"),
    ("Core Capabilities",            r"##\s+Core Capabilities"),
    ("Input Formats",                r"##\s+Input Formats"),
    ("Workflow",                     r"##\s+Workflow"),
    ("CLI Reference",                r"##\s+CLI Reference"),
    ("Demo section",                 r"##\s+Demo"),
    ("Output Structure",             r"##\s+Output Structure"),
    ("Dependencies",                 r"##\s+Dependencies"),
    ("Safety",                       r"##\s+Safety"),
    ("Integration with Bio Orchestrator", r"##\s+Integration with Bio Orchestrator"),
    ("Citations",                    r"##\s+Citations"),
]

# ---------------------------------------------------------------------------
# ANSI colours
# ---------------------------------------------------------------------------

def _use_color() -> bool:
    return sys.stdout.isatty() and os.environ.get("NO_COLOR") is None


_C = _use_color()
BOLD  = "\033[1m"  if _C else ""
DIM   = "\033[2m"  if _C else ""
GREEN = "\033[32m" if _C else ""
RED   = "\033[31m" if _C else ""
CYAN  = "\033[36m" if _C else ""
RESET = "\033[0m"  if _C else ""

# ---------------------------------------------------------------------------
# Draft management (auto-save between interactive sessions)
# ---------------------------------------------------------------------------

DRAFT_SPEC_PATH = Path.home() / ".clawbio_skill_builder_draft.json"


class _Jump(Exception):
    """Raised inside _ask()/_ask_list() when the user wants to navigate."""
    def __init__(self, target: int) -> None:
        self.target = target   # 0-based field index to jump to


def _save_draft(spec: dict) -> None:
    try:
        DRAFT_SPEC_PATH.write_text(json.dumps(spec, indent=2), encoding="utf-8")
    except OSError:
        pass  # non-fatal


def _load_draft() -> dict | None:
    if not DRAFT_SPEC_PATH.exists():
        return None
    try:
        return json.loads(DRAFT_SPEC_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _clear_draft() -> None:
    try:
        DRAFT_SPEC_PATH.unlink(missing_ok=True)
    except OSError:
        pass


# ---------------------------------------------------------------------------
# Spec loading
# ---------------------------------------------------------------------------

def _try_load_yaml(path: Path) -> dict:
    """Attempt to parse a YAML file; falls back gracefully if pyyaml absent."""
    try:
        import yaml  # type: ignore
        with path.open() as fh:
            return yaml.safe_load(fh) or {}
    except ImportError:
        print(
            f"{RED}pyyaml is not installed. Install it with `pip install pyyaml` "
            f"or convert your spec to JSON.{RESET}",
            file=sys.stderr,
        )
        sys.exit(1)


def load_spec(path: Path) -> dict:
    """Load and return a skill spec from a JSON or YAML file."""
    suffix = path.suffix.lower()
    if suffix == ".json":
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            print(f"{RED}Invalid JSON in {path}: {exc}{RESET}", file=sys.stderr)
            sys.exit(1)
    elif suffix in (".yaml", ".yml"):
        return _try_load_yaml(path)
    else:
        # Sniff content
        text = path.read_text(encoding="utf-8").lstrip()
        if text.startswith("{"):
            try:
                return json.loads(text)
            except json.JSONDecodeError:
                pass
        return _try_load_yaml(path)


# ---------------------------------------------------------------------------
# Module-level interactive helpers
# ---------------------------------------------------------------------------

# Ordered field definitions: (key, label, kind, default, help_text)
# kind: "str" | "str_required" | "list" | "deps" | "partners"
_FIELDS: list[tuple[str, str, str, str, str]] = [
    (
        "name",
        "Skill name (lowercase-hyphen)",
        "str_required",
        "",
        "Unique identifier for the skill.\n"
        "  Affects: directory name (skills/<name>/), Python script (<name>.py),\n"
        "           test file (test_<name>.py), and the key in catalog.json.",
    ),
    (
        "description",
        "One-line description",
        "str_required",
        "",
        "One-sentence summary of what this skill does.\n"
        "  Affects: SKILL.md YAML frontmatter 'description:', the opening paragraph\n"
        "           of SKILL.md, and the catalog.json description field.",
    ),
    (
        "author",
        "Your name or handle",
        "str_required",
        "",
        "Credit for the skill author.\n"
        "  Affects: SKILL.md YAML frontmatter 'author:' and the Python script header.",
    ),
    (
        "domain",
        "Scientific domain (e.g. genomics, metagenomics)",
        "str",
        "bioinformatics",
        "The scientific domain this skill belongs to.\n"
        "  Affects: SKILL.md preamble — 'specialised ClawBio agent for <domain>'.",
    ),
    (
        "cli_alias",
        "CLI alias for clawbio.py run <alias>",
        "str",
        "",
        "Short command name users type to invoke the skill.\n"
        "  Affects: the key in the clawbio.py SKILLS dict and catalog.json\n"
        "           demo_command field (python clawbio.py run <alias>).\n"
        "  Defaults to the first hyphen-segment of the skill name.",
    ),
    (
        "capabilities",
        "Core capabilities",
        "list",
        "",
        "What this skill can do — one capability per line.\n"
        "  Affects: SKILL.md '## Core Capabilities' section (one numbered bullet\n"
        "           per entry).",
    ),
    (
        "trigger_keywords",
        "Trigger keywords (phrases that route here)",
        "list",
        "",
        "Phrases the ClawBio router matches to dispatch to this skill.\n"
        "  Affects: SKILL.md YAML frontmatter 'trigger_keywords:' and the\n"
        "           catalog.json trigger_keywords field.",
    ),
    (
        "tags",
        "Tags (optional, e.g. alignment, QC)",
        "list",
        "",
        "Short discoverability labels for the skill.\n"
        "  Affects: SKILL.md YAML frontmatter 'tags:' and the catalog.json\n"
        "           tags field.",
    ),
    (
        "_deps_raw",
        "Required pip dependencies (comma-separated, e.g. pysam>=0.22, biopython)",
        "str",
        "",
        "Python packages this skill needs at runtime.\n"
        "  Affects: SKILL.md YAML 'metadata.openclaw.install:', the\n"
        "           '## Dependencies' section, and import hints in the Python script.",
    ),
    (
        "chaining_partners",
        "Chaining partners (comma-separated skill names, e.g. fastqc, multiqc)",
        "str",
        "",
        "Other ClawBio skills this one is commonly used alongside.\n"
        "  Affects: SKILL.md '## Integration with Bio Orchestrator' section\n"
        "           and the catalog.json chaining_partners field.",
    ),
]
_NUM_FIELDS = len(_FIELDS)


def _find_resume_cursor(spec: dict) -> int:
    """Return the 0-based index of the first unfilled required-ish field."""
    for i, (key, _label, kind, _default, _help) in enumerate(_FIELDS):
        if kind == "str_required" and not spec.get(key, ""):
            return i
        if kind == "list" and not spec.get(key):
            return i
    return _NUM_FIELDS  # all filled; go straight to review


def _ask(prompt: str, default: str = "", required: bool = False,
         field_num: int = 0) -> str:
    """Prompt for a single string value; support !back / !N navigation."""
    hint = f"  {DIM}(Field {field_num}/{_NUM_FIELDS}  — !back to go back, !N to jump){RESET}"
    while True:
        print(hint)
        display = f"{prompt} [{default}]: " if default else f"{prompt}: "
        value = input(display).strip()

        # Navigation commands
        if value.startswith("!"):
            cmd = value[1:].lower()
            if cmd == "back":
                raise _Jump(max(0, field_num - 2))   # field_num is 1-based, index is 0-based
            try:
                target = int(cmd)
                if 1 <= target <= _NUM_FIELDS:
                    raise _Jump(target - 1)
                print(f"  {RED}Jump target must be 1–{_NUM_FIELDS}.{RESET}")
            except ValueError:
                print(f"  {RED}Unknown command '{value}'. Use !back or !N (1–{_NUM_FIELDS}).{RESET}")
            continue

        if not value:
            value = default
        if required and not value:
            print(f"  {RED}This field is required.{RESET}")
            continue
        return value


def _ask_list(prompt: str, field_num: int = 0) -> list[str]:
    """Prompt for a list of strings; support !back / !N navigation."""
    hint = f"  {DIM}(Field {field_num}/{_NUM_FIELDS}  — !back to go back, !N to jump){RESET}"
    print(f"{hint}")
    print(f"{prompt} (one per line, blank line to finish):")
    items: list[str] = []
    while True:
        line = input("  > ").strip()
        if not line:
            break
        if line.startswith("!"):
            cmd = line[1:].lower()
            if items:
                print(f"  {RED}Partial input discarded — previous value kept.{RESET}")
            if cmd == "back":
                raise _Jump(max(0, field_num - 2))
            try:
                target = int(cmd)
                if 1 <= target <= _NUM_FIELDS:
                    raise _Jump(target - 1)
                print(f"  {RED}Jump target must be 1–{_NUM_FIELDS}.{RESET}")
            except ValueError:
                print(f"  {RED}Unknown command '{line}'. Use !back or !N (1–{_NUM_FIELDS}).{RESET}")
            continue
        items.append(line)
    return items


def _fmt_list(items: list[str]) -> str:
    if not items:
        return f"{DIM}(none){RESET}"
    return "\n".join(f"      • {item}" for item in items)


def _print_review(s: dict) -> None:
    deps_req = s.get("dependencies", {}).get("required", [])
    print(f"\n{CYAN}{'─' * 56}{RESET}")
    print(f"{BOLD}  Review your spec{RESET}")
    print(f"{CYAN}{'─' * 56}{RESET}")
    print(f"  {BOLD}[1] name{RESET}             {s.get('name','')}")
    print(f"  {BOLD}[2] description{RESET}      {s.get('description','')}")
    print(f"  {BOLD}[3] author{RESET}           {s.get('author','')}")
    print(f"  {BOLD}[4] domain{RESET}           {s.get('domain','')}")
    print(f"  {BOLD}[5] cli_alias{RESET}        {s.get('cli_alias','')}")
    print(f"  {BOLD}[6] capabilities{RESET}")
    print(_fmt_list(s.get("capabilities", [])))
    print(f"  {BOLD}[7] trigger_keywords{RESET}")
    print(_fmt_list(s.get("trigger_keywords", [])))
    print(f"  {BOLD}[8] tags{RESET}")
    print(_fmt_list(s.get("tags", [])))
    print(f"  {BOLD}[9] dependencies{RESET}     {', '.join(deps_req) if deps_req else '(none)'}")
    print(f"  {BOLD}[10] chaining_partners{RESET} {', '.join(s.get('chaining_partners', [])) or '(none)'}")
    print(f"{CYAN}{'─' * 56}{RESET}")


def _prompt_field(s: dict, idx: int) -> dict:
    """Re-prompt field at 0-based index *idx* and update spec in-place."""
    key, label, kind, default, help_text = _FIELDS[idx]
    field_num = idx + 1

    if help_text:
        print(f"\n  {DIM}{help_text}{RESET}\n")

    if kind in ("str", "str_required"):
        current = s.get(key, default or "")
        # A stored list (e.g. chaining_partners already set) must be joined back to a
        # comma-separated string so that _ask() receives a str default, not a list.
        if isinstance(current, list):
            current = ", ".join(current)
        new_val = _ask(label, default=current, required=(kind == "str_required"),
                       field_num=field_num)
        if key == "_deps_raw":
            s["dependencies"] = {
                "required": [d.strip() for d in new_val.split(",") if d.strip()],
                "optional": s.get("dependencies", {}).get("optional", []),
            }
        elif key == "chaining_partners":
            s["chaining_partners"] = [p.strip() for p in new_val.split(",") if p.strip()]
        else:
            s[key] = new_val
    elif kind == "list":
        new_val = _ask_list(label, field_num=field_num)
        s[key] = new_val
    elif kind == "deps":
        current = ", ".join(s.get("dependencies", {}).get("required", []))
        raw = _ask(label, default=current, field_num=field_num)
        s["dependencies"] = {
            "required": [d.strip() for d in raw.split(",") if d.strip()],
            "optional": s.get("dependencies", {}).get("optional", []),
        }
    elif kind == "partners":
        current = ", ".join(s.get("chaining_partners", []))
        raw = _ask(label, default=current, field_num=field_num)
        s["chaining_partners"] = [p.strip() for p in raw.split(",") if p.strip()]

    # Keep cli_alias default in sync if name changed. Match the convention
    # used elsewhere in the code (_clawbio_entry, generate_catalog_entry):
    # the alias is the first hyphen-segment of the skill name.
    if key == "name" and not s.get("cli_alias"):
        s["cli_alias"] = s["name"].split("-")[0]

    return s


def _run_review_loop(spec: dict) -> dict:
    """
    Show review, let user confirm / edit individual fields / abort.
    Returns the final confirmed spec.
    Jumping to a field from within the review edit re-opens that field and
    returns to review afterward (nested jumps are suppressed).
    """
    while True:
        _print_review(spec)
        print(f"\n  {BOLD}Proceed?{RESET}  Y = yes   n = abort   1-10 = edit a field")
        choice = input("  > ").strip().lower()

        if choice in ("", "y", "yes"):
            break
        elif choice in ("n", "no", "abort"):
            print(f"\n{RED}Aborted — no files written.{RESET}")
            _clear_draft()
            sys.exit(0)
        elif choice in {str(i) for i in range(1, _NUM_FIELDS + 1)}:
            idx = int(choice) - 1
            try:
                spec = _prompt_field(spec, idx)
            except _Jump as j:
                # User typed !N inside an edit — honour the jump
                try:
                    spec = _prompt_field(spec, j.target)
                except _Jump:
                    pass  # suppress nested jumps; return to review
        else:
            print(f"  {DIM}Enter Y, n, or a field number (1-{_NUM_FIELDS}).{RESET}")

    print()
    return spec


# ---------------------------------------------------------------------------
# Interactive spec collection
# ---------------------------------------------------------------------------

def collect_spec_interactive() -> dict:
    """Collect skill spec interactively via stdin prompts, with draft resume and review/edit loop."""
    print(f"\n{CYAN}{BOLD}ClawBio Skill Builder — Interactive Mode{RESET}")
    print(f"{DIM}Press Enter to accept defaults shown in [brackets].{RESET}")
    print(f"{DIM}Type !back to go back one field, !N to jump to field N.{RESET}")
    print(f"{DIM}Progress is auto-saved after each field to: {DRAFT_SPEC_PATH}{RESET}\n")

    # ----------------------------------------------------------------
    # Draft resume
    # ----------------------------------------------------------------
    spec: dict = {}
    start_cursor = 0
    draft = _load_draft()
    if draft:
        print(f"{CYAN}A saved draft was found: {DRAFT_SPEC_PATH}{RESET}")
        resume_choice = input("  Resume? [Y/n]: ").strip().lower()
        if resume_choice in ("", "y", "yes"):
            spec = draft
            start_cursor = _find_resume_cursor(spec)
            print(f"  {DIM}Resuming from field {start_cursor + 1}.{RESET}\n")
        else:
            _clear_draft()
            print(f"  {DIM}Starting fresh.{RESET}\n")

    # ----------------------------------------------------------------
    # Cursor-based collection loop
    # ----------------------------------------------------------------
    cursor = start_cursor
    while cursor < _NUM_FIELDS:
        try:
            spec = _prompt_field(spec, cursor)
            _save_draft(spec)
            cursor += 1
        except _Jump as j:
            cursor = j.target

    # ----------------------------------------------------------------
    # Review / edit / confirm loop
    # ----------------------------------------------------------------
    spec = _run_review_loop(spec)
    return spec


# ---------------------------------------------------------------------------
# Spec validation
# ---------------------------------------------------------------------------

_VALID_NAME_RE = re.compile(r"^[a-z][a-z0-9-]*$")


def validate_spec(spec: dict) -> list[str]:
    """Return a list of error strings (empty = valid)."""
    errors: list[str] = []

    if not spec.get("name"):
        errors.append("'name' is required")
    elif not _VALID_NAME_RE.match(spec["name"]):
        errors.append(
            f"'name' must be lowercase letters, digits, and hyphens — got '{spec['name']}'"
        )

    if not spec.get("description"):
        errors.append("'description' is required")

    if not spec.get("author"):
        errors.append("'author' is required")

    return errors


# ---------------------------------------------------------------------------
# SKILL.md generation
# ---------------------------------------------------------------------------

def _yaml_list(items: list[str], indent: int = 6) -> str:
    """Render a YAML block list with given indent."""
    pad = " " * indent
    return "\n".join(f"{pad}- {item}" for item in items) if items else f"{' ' * indent}[]"


def _md_input_table(formats: list) -> str:
    """Render the Input Formats markdown table.

    Accepts either a list of dicts (full metadata) or a list of strings
    (bare format names / extensions).
    """
    if not formats:
        formats = [
            {
                "format": "TODO: format name",
                "extension": ".ext",
                "required_fields": "field1, field2",
                "example": "demo_input.ext",
            }
        ]
    # Normalise: strings → minimal dicts
    normalised = []
    for f in formats:
        if isinstance(f, str):
            ext = f if f.startswith(".") else f"*.{f}"
            normalised.append({
                "format": f.upper(),
                "extension": ext,
                "required_fields": "—",
                "example": f"demo_input.{f.lstrip('.')}",
            })
        else:
            normalised.append(f)
    header = "| Format | Extension | Required Fields | Example |\n|--------|-----------|-----------------|---------|"
    rows = "\n".join(
        f"| {f.get('format','TODO')} | `{f.get('extension','.ext')}` "
        f"| {f.get('required_fields','TODO')} | `{f.get('example','demo')}` |"
        for f in normalised
    )
    return f"{header}\n{rows}"


def _md_chaining(partners: list[str]) -> str:
    if not partners:
        return "- *(none yet — add partners as you build)*"
    return "\n".join(f"- `{p}`: TODO — describe how they connect" for p in partners)


def _md_deps(deps: dict) -> str:
    req = deps.get("required", [])
    opt = deps.get("optional", [])
    lines = ["**Required** (in `requirements.txt` or skill-level install):"]
    if req:
        lines += [f"- `{d}` — TODO: purpose" for d in req]
    else:
        # Only claim stdlib-only when *nothing* is declared. If optional deps
        # are present the skill is plainly not stdlib-only, so the phrase is
        # misleading — say "none required" instead.
        lines.append(
            "- Python 3.11+ standard library only — no external packages"
            if not opt
            else "- None required (pure stdlib)"
        )
    if opt:
        lines += ["", "**Optional**:"]
        lines += [
            f"- `{d}` — TODO: purpose (graceful degradation without it)" for d in opt
        ]
    return "\n".join(lines)


def _yaml_inline_list(items: list[str]) -> str:
    if not items:
        return "[]"
    return "[" + ", ".join(items) + "]"


def generate_skill_md(spec: dict) -> str:
    """Generate a complete SKILL.md from a skill spec dict."""
    name         = spec["name"]
    description  = spec["description"]
    author       = spec.get("author", "TODO")
    version      = spec.get("version", "0.1.0")
    domain       = spec.get("domain", "bioinformatics")
    emoji        = spec.get("emoji", "🦖")
    tags         = spec.get("tags", [])
    trigger_kws  = spec.get("trigger_keywords", [f"TODO: keyword for {name}"])
    capabilities = spec.get("capabilities", ["TODO: describe capability 1", "TODO: describe capability 2"])
    formats      = spec.get("input_formats", [])
    deps         = spec.get("dependencies", {"required": [], "optional": []})
    partners     = spec.get("chaining_partners", [])
    license_     = spec.get("license", "MIT")
    cli_alias    = spec.get("cli_alias", name.split("-")[0])
    script_name  = name.replace("-", "_")
    title        = " ".join(w.capitalize() for w in name.split("-"))

    caps_md = "\n".join(
        f"{i+1}. **Capability {i+1}**: {cap}" for i, cap in enumerate(capabilities)
    )

    install_section: str
    req_pkgs = deps.get("required", [])
    if req_pkgs:
        install_lines = "\n".join(
            f"      - kind: pip\n        package: {p.split()[0]}\n        bins: []"
            for p in req_pkgs
        )
        install_section = f"    install:\n{install_lines}"
    else:
        install_section = "    install: []"

    return f"""\
---
name: {name}
description: >-
  {description}
version: {version}
author: {author}
license: {license_}
tags: {_yaml_inline_list(tags)}
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "{emoji}"
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
{install_section}
    trigger_keywords:
{_yaml_list(trigger_kws, indent=8)}
---

# {emoji} {title}

You are **{title}**, a specialised ClawBio agent for {domain}. Your role is to {description.rstrip('.')}.

## Why This Exists

- **Without it**: Users must TODO: describe the manual process this replaces.
- **With it**: TODO: describe the automated outcome in seconds/minutes.
- **Why ClawBio**: TODO: explain what makes this better than an LLM guessing — grounded in real databases/algorithms.

## Core Capabilities

{caps_md}

## Input Formats

{_md_input_table(formats)}

## Workflow

When the user asks for TODO: describe the task type:

1. **Validate**: Check input format and required fields
2. **Process**: TODO: core computation — be specific about algorithm/database used
3. **Generate**: TODO: output generation — what gets written where
4. **Report**: Write `report.md` with findings and reproducibility bundle

## CLI Reference

```bash
# Standard usage
python skills/{name}/{script_name}.py \\
  --input <input_file> --output <report_dir>

# Demo mode (synthetic data, no user files needed)
python skills/{name}/{script_name}.py --demo --output /tmp/{name}_demo

# Via ClawBio runner
python clawbio.py run {cli_alias} --input <file> --output <dir>
python clawbio.py run {cli_alias} --demo
```

## Demo

```bash
python clawbio.py run {cli_alias} --demo
```

Expected output: TODO: brief description of what the demo produces.

## Algorithm / Methodology

Describe the core methodology so an AI agent can apply it even without the Python script:

1. **Step**: TODO: detail
2. **Step**: TODO: detail
3. **Step**: TODO: detail

**Key thresholds / parameters**:
- TODO: Parameter 1: value (source: database/paper)

## Example Queries

- "TODO: example query 1 that would route here"
- "TODO: example query 2"

## Output Structure

```
output_directory/
├── report.md              # Primary markdown report
├── result.json            # Machine-readable results
├── figures/
│   └── plot.png           # Visualisation(s)
├── tables/
│   └── results.csv        # Tabular data
└── reproducibility/
    ├── commands.sh        # Exact commands to reproduce
    └── environment.yml    # Conda/pip environment snapshot
```

## Dependencies

{_md_deps(deps)}

## Safety

- **Local-first**: No data upload without explicit consent
- **Disclaimer**: Every report includes the ClawBio medical disclaimer
- **Audit trail**: Log all operations to reproducibility bundle
- **No hallucinated science**: All parameters trace to cited databases

## Integration with Bio Orchestrator

**Trigger conditions** — the orchestrator routes here when:
{chr(10).join(f"- {kw}" for kw in trigger_kws[:3])}

**Chaining partners** — this skill connects with:
{_md_chaining(partners)}

## Citations

- TODO: [Database/Paper 1](URL) — what it provides
- TODO: [Database/Paper 2](URL) — what it provides
"""


# ---------------------------------------------------------------------------
# Python skeleton generation
# ---------------------------------------------------------------------------

def generate_skill_py(spec: dict) -> str:
    """Generate a Python skeleton for the skill."""
    name        = spec["name"]
    description = spec["description"]
    author      = spec.get("author", "TODO")
    script_name = name.replace("-", "_")
    title       = " ".join(w.capitalize() for w in name.split("-"))
    deps        = spec.get("dependencies", {"required": [], "optional": []})
    req_pkgs    = deps.get("required", [])
    import_hint = (
        "\n".join(f"# import {p.split()[0].split('=')[0].split('>')[0].strip()}" for p in req_pkgs)
        if req_pkgs else "# No external dependencies required"
    )

    return textwrap.dedent(f'''\
        #!/usr/bin/env python3
        """
        {script_name}.py — ClawBio {title} Skill
        {"=" * (len(script_name) + len(title) + 22)}
        {description}

        Author:  {author}
        Version: {spec.get("version", "0.1.0")}

        Usage:
            python {script_name}.py --input <input_file> --output <output_dir>
            python {script_name}.py --demo --output /tmp/{name}_demo
        """

        from __future__ import annotations

        import argparse
        import json
        import sys
        from datetime import datetime
        from pathlib import Path

        {import_hint}

        # ---------------------------------------------------------------------------
        # Core logic
        # ---------------------------------------------------------------------------

        def run(input_path: Path, output_dir: Path) -> dict:
            """
            Core skill logic.

            Args:
                input_path:  Path to the input file.
                output_dir:  Directory where all output files will be written.

            Returns:
                dict: Machine-readable results (also written to result.json).
            """
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / "figures").mkdir(exist_ok=True)
            (output_dir / "tables").mkdir(exist_ok=True)
            repro_dir = output_dir / "reproducibility"
            repro_dir.mkdir(exist_ok=True)

            # ------------------------------------------------------------------
            # TODO: Implement your skill logic here.
            #
            # Typical pattern:
            #   data = parse_input(input_path)
            #   results = analyse(data)
            #   write_figures(results, output_dir / "figures")
            #   write_tables(results, output_dir / "tables")
            # ------------------------------------------------------------------

            results: dict = {{
                "skill": "{name}",
                "input": str(input_path),
                "generated_at": datetime.now().isoformat(),
                # TODO: add your result fields here
            }}

            # --- Write report.md ---
            report_lines = [
                f"# {title} Report",
                "",
                f"**Input**: `{{input_path.name}}`",
                f"**Generated**: {{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}}",
                "",
                "## Summary",
                "",
                "TODO: Add your analysis summary here.",
                "",
                "---",
                "",
                "> **Disclaimer**: This output is for research and educational use only.",
                "> It is not a medical device and should not replace professional advice.",
            ]
            (output_dir / "report.md").write_text("\\n".join(report_lines), encoding="utf-8")

            # --- Write result.json ---
            (output_dir / "result.json").write_text(
                json.dumps(results, indent=2, default=str), encoding="utf-8"
            )

            # --- Write reproducibility bundle ---
            cmd = (
                f"python skills/{name}/{script_name}.py "
                f"--input {{input_path}} --output {{output_dir}}"
            )
            (repro_dir / "commands.sh").write_text(
                f"#!/bin/bash\\n# Reproduced: {{datetime.now().isoformat()}}\\n{{cmd}}\\n",
                encoding="utf-8",
            )
            (repro_dir / "environment.yml").write_text(
                "name: {name}\\nchannels:\\n  - conda-forge\\ndependencies:\\n  - python=3.11\\n",
                encoding="utf-8",
            )

            return results


        # ---------------------------------------------------------------------------
        # Demo mode
        # ---------------------------------------------------------------------------

        def run_demo(output_dir: Path) -> None:
            """Run with synthetic demo data to verify the skill works end-to-end."""
            output_dir.mkdir(parents=True, exist_ok=True)

            # ------------------------------------------------------------------
            # TODO: Create synthetic demo input that exercises your core logic.
            # The demo should run in < 10 seconds with no external downloads.
            # ------------------------------------------------------------------

            demo_input = output_dir / "demo_input.txt"
            demo_input.write_text(
                "# Demo input for {title}\\n"
                "# TODO: replace with realistic synthetic data\\n"
                "sample_id\\tvalue\\n"
                "SAMPLE_01\\t42\\n"
                "SAMPLE_02\\t17\\n",
                encoding="utf-8",
            )

            results = run(demo_input, output_dir)

            print(f"  Demo complete. Output: {{output_dir}}")
            print(f"  Files: {{', '.join(f.name for f in output_dir.rglob('*') if f.is_file())}}")
            return results


        # ---------------------------------------------------------------------------
        # CLI entry point
        # ---------------------------------------------------------------------------

        def main() -> None:
            parser = argparse.ArgumentParser(
                description="{description}",
                formatter_class=argparse.RawDescriptionHelpFormatter,
            )
            parser.add_argument("--input",  dest="input_path", help="Path to input file")
            parser.add_argument(
                "--output",
                dest="output_dir",
                default="/tmp/{name}_output",
                help="Output directory (default: /tmp/{name}_output)",
            )
            parser.add_argument("--demo", action="store_true", help="Run with synthetic demo data")
            args = parser.parse_args()

            out = Path(args.output_dir)

            if args.demo:
                run_demo(out)
            elif args.input_path:
                results = run(Path(args.input_path), out)
                print(f"  Done. Output: {{out}}")
                report = out / "report.md"
                if report.exists():
                    print(f"  Report: {{report}}")
            else:
                parser.error("Provide --input <file> or --demo")


        if __name__ == "__main__":
            main()
    ''')


# ---------------------------------------------------------------------------
# Test skeleton generation
# ---------------------------------------------------------------------------

def generate_test_py(spec: dict) -> str:
    """Generate a pytest skeleton for the skill."""
    name        = spec["name"]
    script_name = name.replace("-", "_")
    title       = " ".join(w.capitalize() for w in name.split("-"))

    return textwrap.dedent(f'''\
        """
        Tests for the {title} skill.

        Run:
            pytest skills/{name}/tests/ -v
        or via ClawBio runner:
            python -m pytest skills/{name}/tests/ -v
        """

        import json
        from pathlib import Path

        import pytest
        import sys

        # Make the skill importable regardless of working directory
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

        from {script_name} import run, run_demo


        # ---------------------------------------------------------------------------
        # Fixtures
        # ---------------------------------------------------------------------------

        @pytest.fixture
        def tmp_output(tmp_path: Path) -> Path:
            """Provide a temporary output directory per test."""
            return tmp_path / "output"


        # ---------------------------------------------------------------------------
        # Tests
        # ---------------------------------------------------------------------------

        def test_demo_runs(tmp_output: Path) -> None:
            """Demo mode should complete without raising an exception."""
            run_demo(tmp_output)


        def test_report_generated(tmp_output: Path) -> None:
            """Demo should produce a non-empty report.md."""
            run_demo(tmp_output)
            report = tmp_output / "report.md"
            assert report.exists(), "report.md was not created"
            assert report.stat().st_size > 0, "report.md is empty"


        def test_result_json_valid(tmp_output: Path) -> None:
            """result.json should be present and parse as valid JSON."""
            run_demo(tmp_output)
            result_file = tmp_output / "result.json"
            assert result_file.exists(), "result.json was not created"
            data = json.loads(result_file.read_text())
            assert isinstance(data, dict), "result.json should be a JSON object"
            assert data.get("skill") == "{name}", "result.json should include skill name"


        def test_reproducibility_bundle(tmp_output: Path) -> None:
            """Reproducibility bundle (commands.sh + environment.yml) should be present."""
            run_demo(tmp_output)
            repro = tmp_output / "reproducibility"
            assert (repro / "commands.sh").exists(), "commands.sh missing"
            assert (repro / "environment.yml").exists(), "environment.yml missing"


        # ---------------------------------------------------------------------------
        # TODO: Add domain-specific tests below
        # ---------------------------------------------------------------------------

        # def test_core_logic(tmp_output: Path) -> None:
        #     """TODO: Test the actual analysis logic with known input/output."""
        #     demo_input = tmp_output / "test_input.txt"
        #     demo_input.parent.mkdir(parents=True, exist_ok=True)
        #     demo_input.write_text("...")
        #     results = run(demo_input, tmp_output)
        #     assert results["..."] == expected_value
    ''')


# ---------------------------------------------------------------------------
# Catalog.json helpers
# ---------------------------------------------------------------------------

def generate_catalog_entry(spec: dict) -> dict:
    """Build a catalog.json entry dict for the skill."""
    name        = spec["name"]
    cli_alias   = spec.get("cli_alias") or name.split("-")[0]
    script_name = name.replace("-", "_")  # noqa: F841 — kept for future schema
    deps        = spec.get("dependencies") or {}
    return {
        "name": name,
        "cli_alias": cli_alias,
        "description": spec["description"],
        "version": spec.get("version", "0.1.0"),
        "status": "mvp",
        "has_script": True,
        "has_tests": True,
        "has_demo": True,
        "demo_command": f"python clawbio.py run {cli_alias} --demo",
        "dependencies": list(deps.get("required", []) or []),
        "tags": spec.get("tags", []),
        "trigger_keywords": spec.get("trigger_keywords", []),
        "chaining_partners": spec.get("chaining_partners", []),
    }


def update_catalog_json(catalog_path: Path, entry: dict) -> bool:
    """
    Append the new skill entry to skills/catalog.json.
    Returns True on success, False if skill already exists.

    Uses ensure_ascii=False so non-ASCII characters already present in other
    skills' descriptions (em-dashes etc.) are not re-encoded to \\uXXXX, which
    would produce a noisy diff across every unrelated entry.
    """
    try:
        catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        print(f"  {RED}Could not read catalog.json: {exc}{RESET}", file=sys.stderr)
        return False

    existing_names = {s["name"] for s in catalog.get("skills", [])}
    if entry["name"] in existing_names:
        print(f"  {DIM}catalog.json: '{entry['name']}' already exists — skipping.{RESET}")
        return False

    catalog.setdefault("skills", []).append(entry)
    catalog["skill_count"] = len(catalog["skills"])
    catalog_path.write_text(
        json.dumps(catalog, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"  {GREEN}Updated catalog.json (+1 skill){RESET}")
    return True


def regenerate_catalog(repo_root: Path) -> bool:
    """
    Run ``scripts/generate_catalog.py`` to rebuild catalog.json from scratch.

    The catalog is derived from each SKILL.md's frontmatter plus hand-maintained
    tables inside the generator (aliases, trigger keywords, chaining). That
    script is the repo's single source of truth, so calling it here keeps
    skill-builder's output in lockstep with the rest of the tooling.
    """
    import subprocess

    script = repo_root / "scripts" / "generate_catalog.py"
    if not script.exists():
        return False
    try:
        subprocess.run(
            [sys.executable, str(script)],
            cwd=repo_root,
            check=True,
            capture_output=True,
        )
    except (subprocess.CalledProcessError, OSError) as exc:
        stderr = getattr(exc, "stderr", b"")
        stderr_text = stderr.decode("utf-8", errors="replace") if stderr else ""
        print(
            f"  {RED}generate_catalog.py failed: {exc}{RESET}\n{stderr_text}",
            file=sys.stderr,
        )
        return False
    print(f"  {GREEN}Regenerated catalog.json via scripts/generate_catalog.py{RESET}")
    return True


# ---------------------------------------------------------------------------
# clawbio.py patching
# ---------------------------------------------------------------------------

def _clawbio_entry(spec: dict) -> str:
    """Render the Python dict literal for the SKILLS entry."""
    name        = spec["name"]
    cli_alias   = spec.get("cli_alias") or name.split("-")[0]
    description = spec["description"]
    script_name = name.replace("-", "_")

    return textwrap.dedent(f'''\
        "{cli_alias}": {{
            "script": SKILLS_DIR / "{name}" / "{script_name}.py",
            "demo_args": ["--demo"],
            "description": "{description}",
            "allowed_extra_flags": set(),
            "accepts_genotypes": False,
        }},
    ''')


# Pattern used to insert new entries *inside* the SKILLS dict in clawbio.py.
# Matches the closing brace of SKILLS followed by the comment that demarcates
# the next top-level block. We capture the brace so the replacement can put
# the new entry immediately before it.
#
#     …last entry…
#         },
#     }                                         <-- this }
#                                               <-- blank line
#     # Skills that run in the full-profile pipeline   <-- anchor comment
#
# An older version of this function anchored only on the comment and inserted
# *after* the closing brace, producing invalid Python (entry outside the dict).
_CLAWBIO_DICT_CLOSE_RE = re.compile(
    r"\n(?P<close>\})\n(?P<blank>\s*)(?P<anchor>#\s*Skills that run in the full-profile pipeline)",
)


def _alias_already_present(source: str, cli_alias: str) -> bool:
    """True iff cli_alias appears as a top-level SKILLS dict key.

    Plain substring search gives false positives when the alias happens to
    appear inside a description or another string literal, so match the
    dict-key syntax anchored to the 4-space indent used in SKILLS.
    """
    key_re = re.compile(rf'^\s{{4}}"{re.escape(cli_alias)}":\s*\{{', re.MULTILINE)
    return bool(key_re.search(source))


def patch_clawbio_py(clawbio_path: Path, spec: dict) -> bool:
    """
    Insert a new entry into the SKILLS dict in clawbio.py.
    Returns True on success, False if alias already exists or layout unrecognised.
    """
    try:
        source = clawbio_path.read_text(encoding="utf-8")
    except OSError as exc:
        print(f"  {RED}Could not read clawbio.py: {exc}{RESET}", file=sys.stderr)
        return False

    cli_alias = spec.get("cli_alias") or spec["name"].split("-")[0]

    if _alias_already_present(source, cli_alias):
        print(f"  {DIM}clawbio.py: alias '{cli_alias}' already present — skipping.{RESET}")
        return False

    m = _CLAWBIO_DICT_CLOSE_RE.search(source)
    if m is None:
        print(
            f"  {RED}clawbio.py: could not locate SKILLS dict close marker — "
            f"skipping auto-patch.{RESET}",
            file=sys.stderr,
        )
        return False

    new_entry = _clawbio_entry(spec)
    # Indent the entry to match the SKILLS dict (4 spaces).
    indented = textwrap.indent(new_entry, "    ")
    # Insert the indented entry immediately before the closing }.
    insert_at = m.start("close")
    patched = source[:insert_at] + indented + source[insert_at:]
    clawbio_path.write_text(patched, encoding="utf-8")
    print(f"  {GREEN}Patched clawbio.py (+1 SKILLS entry for '{cli_alias}'){RESET}")
    return True


# ---------------------------------------------------------------------------
# Repo root detection
# ---------------------------------------------------------------------------

def find_repo_root() -> Path | None:
    """
    Walk up from the script location looking for clawbio.py (repo root sentinel).
    Returns the repo root Path, or None if not found.
    """
    candidate = SCRIPT_DIR
    for _ in range(6):  # max 6 levels up
        if (candidate / REPO_SENTINEL).exists():
            return candidate
        candidate = candidate.parent
    return None


def _is_inside(child: Path, parent: Path) -> bool:
    """True if ``child`` is the same path as ``parent`` or nested inside it.

    Both paths are resolved so relative vs. absolute, or trailing slashes,
    don't trip the comparison. Used to decide whether a scaffold should
    touch the repo-wide registry files.
    """
    try:
        child_r  = child.resolve()
        parent_r = parent.resolve()
    except OSError:
        return False
    return child_r == parent_r or parent_r in child_r.parents


def _run_skill_lint(repo_root: Path, skill_name: str) -> None:
    """
    Run the repo's own ``scripts/lint_skills.py`` and echo any row that
    mentions the freshly scaffolded skill.

    This is advisory — it never fails the scaffold, it just surfaces issues
    the same CI job will flag later (wrong ``os:`` value, missing emoji etc.).
    """
    import subprocess

    lint = repo_root / "scripts" / "lint_skills.py"
    if not lint.exists():
        return
    try:
        result = subprocess.run(
            [sys.executable, str(lint)],
            cwd=repo_root,
            capture_output=True,
            text=True,
        )
    except OSError:
        return

    rows = [
        line for line in result.stdout.splitlines()
        if line.startswith("|") and skill_name in line
    ]
    if not rows:
        return
    # Each row looks like: | name | Status | OS | Emoji | Issues |
    for row in rows:
        if " FAIL " in row:
            print(f"  {RED}lint: {row.strip()}{RESET}")
        elif " Warn " in row:
            print(f"  {DIM}lint: {row.strip()}{RESET}")


# ---------------------------------------------------------------------------
# SKILL.md validation
# ---------------------------------------------------------------------------

def validate_skill_md_content(content: str) -> list[tuple[str, bool, str]]:
    """
    Validate the *text* of a SKILL.md against the CONTRIBUTING checklist.

    Returns a list of (section_name, passed, detail) tuples.
    Accepts a plain string so callers that already hold the content in memory
    (e.g. agent dry-run mode) do not need to touch the filesystem.
    """
    results: list[tuple[str, bool, str]] = []
    for section_name, pattern in REQUIRED_SECTIONS:
        found = bool(re.search(pattern, content, re.MULTILINE))
        results.append((section_name, found, "present" if found else "MISSING"))
    return results


def validate_skill_md(md_path: Path) -> list[tuple[str, bool, str]]:
    """
    Validate a SKILL.md against the CONTRIBUTING checklist.

    Returns a list of (section_name, passed, detail) tuples.
    """
    try:
        content = md_path.read_text(encoding="utf-8")
    except OSError as exc:
        return [("File readable", False, str(exc))]
    return validate_skill_md_content(content)


def print_validation_report(md_path: Path, results: list[tuple[str, bool, str]]) -> int:
    """Print the validation report and return the number of failures."""
    passed = sum(1 for _, ok, _ in results if ok)
    total  = len(results)
    print(f"\n{BOLD}SKILL.md Validation — {md_path.name}{RESET}")
    print(f"{'─' * 50}")
    for name, ok, detail in results:
        icon  = f"{GREEN}PASS{RESET}" if ok else f"{RED}FAIL{RESET}"
        print(f"  [{icon}] {name}")
    print(f"{'─' * 50}")
    score = f"{passed}/{total}"
    if passed == total:
        print(f"{GREEN}{BOLD}Score: {score} — all checks passed!{RESET}")
    else:
        failures = total - passed
        print(f"{RED}{BOLD}Score: {score} — fix {failures} issue(s) before submitting PR{RESET}")
    print()
    return total - passed


# ---------------------------------------------------------------------------
# Scaffold orchestration
# ---------------------------------------------------------------------------

def scaffold_skill(
    spec: dict,
    output_base: Path,
    dry_run: bool = False,
    repo_root: Path | None = None,
    ask_overwrite: bool = False,
) -> dict:
    """
    Create all skill files and (if repo_root found) update registry.
    Returns a manifest dict summarising what was created.

    ask_overwrite: if True, prompt the user before overwriting existing files
                   instead of silently skipping them.
    """
    name        = spec["name"]
    script_name = name.replace("-", "_")
    skill_dir   = output_base / name

    files_written: list[str] = []
    files_skipped: list[str] = []

    def write_file(path: Path, content: str, label: str, force: bool = False) -> None:
        """Write content to path, respecting dry_run / overwrite policy."""
        rel = str(path.relative_to(output_base) if output_base in path.parents else path)
        if path.exists() and not force:
            if ask_overwrite:
                ans = input(f"  Overwrite {rel}? [y/N]: ").strip().lower()
                if ans not in ("y", "yes"):
                    print(f"  {DIM}kept:          {rel}{RESET}")
                    files_skipped.append(rel)
                    return
                # Fall through to write (updated)
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(content, encoding="utf-8")
                print(f"  {CYAN}updated:{RESET}       {rel}")
                files_written.append(rel)
                return
            else:
                print(f"  {DIM}skip (exists):  {rel}{RESET}")
                files_skipped.append(rel)
                return
        if dry_run:
            print(f"\n{CYAN}{'─' * 60}{RESET}")
            print(f"{BOLD}{label}: {rel}{RESET}")
            print(f"{CYAN}{'─' * 60}{RESET}")
            print(content[:2000])
            if len(content) > 2000:
                print(f"{DIM}... ({len(content) - 2000} more chars){RESET}")
            files_written.append(f"[dry-run] {rel}")
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")
            print(f"  {GREEN}created:{RESET}       {rel}")
            files_written.append(rel)

    print(f"\n{CYAN}{BOLD}🦖 Scaffolding skill: {name}{RESET}")
    if dry_run:
        print(f"{DIM}Dry run — no files will be written.{RESET}\n")

    # 1. SKILL.md
    write_file(skill_dir / "SKILL.md", generate_skill_md(spec), "SKILL.md")

    # 2. Python script
    write_file(skill_dir / f"{script_name}.py", generate_skill_py(spec), "Python script")

    # 3. Test skeleton
    write_file(
        skill_dir / "tests" / f"test_{script_name}.py",
        generate_test_py(spec),
        "Test skeleton",
    )

    # 4. Example spec copy (JSON, for reproducibility)
    # In --edit mode (ask_overwrite=True) always refresh the spec — no prompt.
    write_file(
        skill_dir / "examples" / "example_spec.json",
        json.dumps(spec, indent=2),
        "Example spec",
        force=ask_overwrite,
    )

    # 5. Validate the generated SKILL.md (only in real mode)
    if not dry_run:
        skill_md_path = skill_dir / "SKILL.md"
        if skill_md_path.exists():
            results = validate_skill_md(skill_md_path)
            failures = print_validation_report(skill_md_path, results)
            if failures:
                print(f"  {DIM}(Fill in the TODO: placeholders to resolve these.){RESET}\n")

    # 6. Update catalog.json + clawbio.py if:
    #      - we found the repo root,
    #      - we're writing a real scaffold (not a dry run),
    #      - AND the output is landing inside that repo's skills/ directory.
    #
    # The last guard matters: running `--demo --output /tmp/foo` from inside a
    # clone used to patch the real repo's registry to point at /tmp files,
    # leaving the registry pointing at nonexistent skills. Only register when
    # the files we just wrote are actually discoverable from the repo.
    catalog_updated = False
    clawbio_updated = False
    writes_into_repo = (
        repo_root is not None
        and not dry_run
        and _is_inside(output_base, repo_root / "skills")
    )
    if writes_into_repo:
        catalog_path  = repo_root / "skills" / "catalog.json"
        clawbio_path  = repo_root / "clawbio.py"

        if clawbio_path.exists():
            clawbio_updated = patch_clawbio_py(clawbio_path, spec)
        else:
            print(f"  {DIM}clawbio.py not found at {clawbio_path} — skipping.{RESET}")

        # Prefer the canonical generator over hand-appending: catalog.json's
        # own header says "generated_by: scripts/generate_catalog.py", and the
        # generator consults tables inside that script that we don't have here.
        # Fall back to hand-append only if the generator is unavailable.
        if (repo_root / "scripts" / "generate_catalog.py").exists():
            catalog_updated = regenerate_catalog(repo_root)
        elif catalog_path.exists():
            catalog_updated = update_catalog_json(
                catalog_path, generate_catalog_entry(spec)
            )
        else:
            print(f"  {DIM}catalog.json not found at {catalog_path} — skipping.{RESET}")

        # Surface OpenClaw discoverability lint for the freshly written skill
        # so CI surprises don't happen at PR time.
        _run_skill_lint(repo_root, spec["name"])
    elif repo_root and not dry_run:
        print(
            f"  {DIM}Output is outside repo skills/ directory — "
            f"leaving catalog.json and clawbio.py untouched.{RESET}"
        )

    manifest = {
        "skill": name,
        "output_dir": str(skill_dir),
        "files_written": files_written,
        "files_skipped": files_skipped,
        "catalog_updated": catalog_updated,
        "clawbio_updated": clawbio_updated,
        "dry_run": dry_run,
        "generated_at": datetime.now().isoformat(),
    }

    return manifest


# ---------------------------------------------------------------------------
# Report writing
# ---------------------------------------------------------------------------

def write_skill_builder_report(
    manifest: dict,
    output_dir: Path,
    spec: dict,
    repo_root: Path | None,
) -> None:
    """Write the skill-builder's own report.md + result.json to output_dir."""
    name      = spec["name"]
    cli_alias = spec.get("cli_alias") or name.split("-")[0]
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    # report.md — numbered fresh after deciding which conditional steps apply
    steps: list[str] = []
    if not manifest.get("catalog_updated"):
        steps.append(
            "Rebuild `skills/catalog.json`: "
            "`python scripts/generate_catalog.py` (from the repo root)"
        )
    if not manifest.get("clawbio_updated"):
        steps.append(
            "Add the skill to the `SKILLS` dict in `clawbio.py` — see "
            "`examples/example_spec.json` for the alias and script path"
        )
    steps.append(
        f"Fill in the `TODO:` placeholders in `skills/{name}/SKILL.md` and "
        f"`skills/{name}/{name.replace('-','_')}.py`"
    )
    steps.append(f"Run tests: `python -m pytest skills/{name}/tests/ -v`")
    steps.append(
        f"Lint the skill: `python scripts/lint_skills.py` "
        "(CI will reject any new skill whose frontmatter fails this)"
    )
    steps.append(
        f"Submit PR: `git checkout -b add-{name} && git add skills/{name}/ "
        f"&& git commit -m 'Add {name} skill' && git push`"
    )
    next_steps = [f"{i}. {s}" for i, s in enumerate(steps, start=1)]

    report_lines = [
        f"# Skill Builder Report",
        "",
        f"**Skill scaffolded**: `{name}`",
        f"**Generated**: {manifest['generated_at']}",
        f"**Output directory**: `{manifest['output_dir']}`",
        "",
        "## Files Created",
        "",
    ]
    for f in manifest["files_written"]:
        report_lines.append(f"- `{f}`")
    if manifest["files_skipped"]:
        report_lines += ["", "## Files Skipped (already exist)", ""]
        for f in manifest["files_skipped"]:
            report_lines.append(f"- `{f}`")

    report_lines += [
        "",
        "## Registry Updates",
        "",
        f"- catalog.json updated: {'yes' if manifest['catalog_updated'] else 'no'}",
        f"- clawbio.py updated:   {'yes' if manifest['clawbio_updated'] else 'no'}",
        "",
        "## Next Steps",
        "",
    ]
    report_lines += next_steps
    report_lines += [
        "",
        "---",
        "",
        "> This report was generated by the **ClawBio Skill Builder**.",
        "> Skill Builder is itself a ClawBio skill — built with itself.",
    ]

    (output_dir / "report.md").write_text("\n".join(report_lines), encoding="utf-8")
    (output_dir / "result.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    cmd_input = spec.get("_source_path", "spec.json")
    (repro_dir / "commands.sh").write_text(
        f"#!/bin/bash\n# Reproduced: {manifest['generated_at']}\n"
        f"python skills/skill-builder/skill_builder.py --input {cmd_input}\n",
        encoding="utf-8",
    )


# ---------------------------------------------------------------------------
# Demo mode
# ---------------------------------------------------------------------------

def run_demo(output_dir: Path, dry_run: bool = False) -> None:
    """Scaffold the hello-bioinformatics demo skill."""
    if not DEMO_SPEC_PATH.exists():
        print(
            f"{RED}Demo spec not found at {DEMO_SPEC_PATH}. "
            f"Ensure demo_spec.json is present alongside skill_builder.py.{RESET}",
            file=sys.stderr,
        )
        sys.exit(1)

    spec = load_spec(DEMO_SPEC_PATH)
    spec["_source_path"] = str(DEMO_SPEC_PATH)

    # Demo always writes into a subdirectory so it's self-contained
    skill_output_base = output_dir
    repo_root = find_repo_root()

    manifest = scaffold_skill(spec, skill_output_base, dry_run=dry_run, repo_root=repo_root)

    if not dry_run:
        write_skill_builder_report(manifest, output_dir, spec, repo_root)
        print(f"\n  {GREEN}{BOLD}Demo complete!{RESET}")
        print(f"  Skill directory: {manifest['output_dir']}")
        print(f"  Report:          {output_dir / 'report.md'}")


# ---------------------------------------------------------------------------
# Agent mode
# ---------------------------------------------------------------------------

def run_agent(
    spec: dict,
    output_base: Path,
    dry_run: bool = False,
    repo_root: Path | None = None,
) -> None:
    """
    Machine-readable execution path for LLM / CI callers.

    Scaffolds the skill with all human-readable progress output suppressed,
    then emits a single JSON object to stdout:

    Success::

        {
            "status": "ok",
            "skill": "my-skill",
            "output_dir": "/abs/path/skills/my-skill",
            "files_written": ["my-skill/SKILL.md", ...],
            "files_skipped": [],
            "catalog_updated": false,
            "clawbio_updated": false,
            "validation": {"passed": 13, "failed": 0, "failures": []},
            "dry_run": false,
            "generated_at": "2026-03-19T12:00:00",
            // only present when dry_run=true:
            "generated_content": {
                "SKILL.md": "...",
                "my_skill.py": "...",
                "tests/test_my_skill.py": "...",
                "examples/example_spec.json": "..."
            }
        }

    Error (spec validation failed — exits non-zero)::

        {"status": "error", "errors": ["name: must be lowercase-hyphen"]}
    """
    import io as _io

    script_name = spec["name"].replace("-", "_")

    # ------------------------------------------------------------------
    # Scaffold with stdout silenced (all decorative print() calls go to /dev/null)
    # ------------------------------------------------------------------
    _saved_stdout = sys.stdout
    sys.stdout = _io.StringIO()
    try:
        manifest = scaffold_skill(spec, output_base, dry_run=dry_run, repo_root=repo_root)
    finally:
        sys.stdout = _saved_stdout

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------
    if dry_run:
        # No files written — validate the generated content string directly.
        skill_md_content = generate_skill_md(spec)
        val_results = validate_skill_md_content(skill_md_content)
    else:
        skill_md_path = output_base / spec["name"] / "SKILL.md"
        val_results = (
            validate_skill_md(skill_md_path)
            if skill_md_path.exists()
            else [("File readable", False, "SKILL.md not found after scaffold")]
        )

    passed   = sum(1 for _, ok, _ in val_results if ok)
    failures = [name for name, ok, _ in val_results if not ok]
    manifest["validation"] = {"passed": passed, "failed": len(failures), "failures": failures}
    manifest["status"] = "ok"

    # ------------------------------------------------------------------
    # In dry-run mode, include the generated file contents in the JSON so
    # the calling agent can inspect them without hitting the filesystem.
    # ------------------------------------------------------------------
    if dry_run:
        manifest["generated_content"] = {
            "SKILL.md":                        generate_skill_md(spec),
            f"{script_name}.py":               generate_skill_py(spec),
            f"tests/test_{script_name}.py":    generate_test_py(spec),
            "examples/example_spec.json":      json.dumps(spec, indent=2),
        }

    print(json.dumps(manifest, indent=2))


# ---------------------------------------------------------------------------
# Edit mode
# ---------------------------------------------------------------------------

def edit_skill(skill_dir: Path) -> None:
    """
    Re-scaffold an existing skill directory after editing its spec.

    Reads ``<skill_dir>/examples/example_spec.json``, drops directly into
    the review/edit loop so the user can adjust any field, then calls
    scaffold_skill() with ask_overwrite=True so existing files can be
    selectively updated.  The example_spec.json is always refreshed.
    """
    spec_path = skill_dir / "examples" / "example_spec.json"
    if not spec_path.exists():
        print(
            f"{RED}No example_spec.json found in {skill_dir / 'examples'}.{RESET}\n"
            f"  Create one or use --input to scaffold from scratch.",
            file=sys.stderr,
        )
        sys.exit(1)

    spec = load_spec(spec_path)
    skill_name = spec.get("name", skill_dir.name)

    print(f"\n{CYAN}{BOLD}ClawBio Skill Builder — Edit Mode: {skill_name}{RESET}")
    print(f"{DIM}Review or edit any field, then confirm to re-scaffold.{RESET}\n")

    spec = _run_review_loop(spec)

    # Resolve the output base (parent of skill_dir)
    output_base = skill_dir.parent
    repo_root   = find_repo_root()

    manifest = scaffold_skill(
        spec,
        output_base,
        dry_run=False,
        repo_root=repo_root,
        ask_overwrite=True,
    )

    print(f"\n  {GREEN}{BOLD}Edit complete!{RESET}")
    print(f"  Skill directory: {manifest['output_dir']}")



# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="ClawBio Skill Builder — scaffold a new ClawBio skill from a spec",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Examples:
              python skill_builder.py --input spec.json
              python skill_builder.py --interactive
              python skill_builder.py --demo --output /tmp/demo
              python skill_builder.py --validate-only --input skills/my-skill/SKILL.md
              python skill_builder.py --input spec.json --dry-run
              python skill_builder.py --edit skills/my-skill/

              # Agent / LLM pipeline — JSON manifest to stdout, no interactive prompts:
              python skill_builder.py --agent --input spec.json
              python skill_builder.py --agent --input spec.json --dry-run
              echo '{"name":"foo",...}' | python skill_builder.py --agent --input -

            Notes:
              Interactive mode (--interactive) auto-saves your progress after every
              field to:
                ~/.clawbio_skill_builder_draft.json
              You will be offered the option to resume this draft the next time you
              run --interactive. To discard it and start fresh, delete the file:
                rm ~/.clawbio_skill_builder_draft.json
        """),
    )

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--input",         dest="input_path", metavar="FILE",
                      help="Path to skill spec (JSON/YAML), SKILL.md (--validate-only), "
                           "or '-' to read JSON from stdin")
    mode.add_argument("--interactive",   action="store_true",
                      help="Collect skill spec interactively")
    mode.add_argument("--demo",          action="store_true",
                      help="Scaffold the hello-bioinformatics demo skill")
    mode.add_argument("--edit",          dest="edit_dir", metavar="DIR",
                      help="Re-scaffold an existing skill directory (reads examples/example_spec.json)")

    parser.add_argument("--output",      dest="output_dir",  metavar="DIR",
                        help="Output directory (default: ./generated_skills/ or skills/ in repo)")
    parser.add_argument("--validate-only", action="store_true",
                        help="Validate a SKILL.md against the CONTRIBUTING checklist and exit")
    parser.add_argument("--dry-run",     action="store_true",
                        help="Print generated content without writing any files")
    parser.add_argument("--agent",       action="store_true",
                        help="Machine-readable mode: suppress all progress output and emit a "
                             "single JSON manifest to stdout. Compatible with --input FILE or "
                             "--input - (spec JSON read from stdin). Exits non-zero on error.")

    args = parser.parse_args()

    repo_root   = find_repo_root()
    default_out = (repo_root / "skills") if repo_root else Path("generated_skills")

    # ------------------------------------------------------------------
    # Validate-only mode (reads a SKILL.md, does not scaffold)
    # ------------------------------------------------------------------
    if args.validate_only:
        if not args.input_path:
            parser.error("--validate-only requires --input <path/to/SKILL.md>")
        md_path = Path(args.input_path)
        if not md_path.exists():
            print(f"{RED}File not found: {md_path}{RESET}", file=sys.stderr)
            sys.exit(1)
        results  = validate_skill_md(md_path)
        failures = print_validation_report(md_path, results)
        sys.exit(0 if failures == 0 else 1)

    # ------------------------------------------------------------------
    # Edit mode
    # ------------------------------------------------------------------
    if args.edit_dir:
        skill_dir = Path(args.edit_dir).resolve()
        if not skill_dir.is_dir():
            print(f"{RED}Directory not found: {skill_dir}{RESET}", file=sys.stderr)
            sys.exit(1)
        edit_skill(skill_dir)
        return

    # ------------------------------------------------------------------
    # Demo mode
    # ------------------------------------------------------------------
    if args.demo:
        out = Path(args.output_dir) if args.output_dir else Path("/tmp/skill_builder_demo")
        run_demo(out, dry_run=args.dry_run)
        return

    # ------------------------------------------------------------------
    # Spec-driven or interactive
    # ------------------------------------------------------------------
    if args.interactive:
        spec = collect_spec_interactive()
    else:
        if args.input_path == "-":
            # Read spec JSON from stdin (agent / pipeline use)
            try:
                spec = json.load(sys.stdin)
            except json.JSONDecodeError as exc:
                if args.agent:
                    print(json.dumps({"status": "error", "errors": [f"stdin JSON parse error: {exc}"]}))
                else:
                    print(f"{RED}JSON parse error reading from stdin: {exc}{RESET}", file=sys.stderr)
                sys.exit(1)
        else:
            input_path = Path(args.input_path)
            if not input_path.exists():
                if args.agent:
                    print(json.dumps({"status": "error", "errors": [f"spec file not found: {input_path}"]}))
                else:
                    print(f"{RED}Spec file not found: {input_path}{RESET}", file=sys.stderr)
                sys.exit(1)
            spec = load_spec(input_path)
            spec["_source_path"] = str(input_path)

    # Validate spec fields
    errors = validate_spec(spec)
    if errors:
        if args.agent:
            print(json.dumps({"status": "error", "errors": errors}))
        else:
            print(f"{RED}{BOLD}Spec validation errors:{RESET}", file=sys.stderr)
            for err in errors:
                print(f"  • {err}", file=sys.stderr)
        sys.exit(1)

    # Resolve output directory
    if args.output_dir:
        output_base = Path(args.output_dir)
    else:
        output_base = default_out

    # ------------------------------------------------------------------
    # Agent mode — emit JSON manifest and return
    # ------------------------------------------------------------------
    if args.agent:
        run_agent(spec, output_base, dry_run=args.dry_run, repo_root=repo_root)
        return

    # Scaffold (human mode)
    manifest = scaffold_skill(spec, output_base, dry_run=args.dry_run, repo_root=repo_root)

    if not args.dry_run:
        _clear_draft()  # remove any lingering interactive draft

        # Write skill-builder's own output report into a *hidden* subdirectory
        # so it's skipped by scripts/lint_skills.py and scripts/generate_catalog.py
        # (both of which skip directories starting with a dot). Older leading-
        # underscore names were treated as real skills and failed lint for
        # missing a SKILL.md.
        sb_report_dir = output_base / f".skill_builder_report_{spec['name']}"
        write_skill_builder_report(manifest, sb_report_dir, spec, repo_root)

        print(f"\n  {GREEN}{BOLD}Scaffold complete!{RESET}")
        print(f"  Skill directory: {manifest['output_dir']}")
        if (sb_report_dir / "report.md").exists():
            print(f"  Builder report:  {sb_report_dir / 'report.md'}")
        if repo_root:
            print(f"  Repo root found: {repo_root}")
        else:
            print(f"\n  {DIM}Repo root not found — catalog.json and clawbio.py were not updated.")
            print(f"  To update them manually, run this script from inside the ClawBio repo.{RESET}")


if __name__ == "__main__":
    main()
