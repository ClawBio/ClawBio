#!/usr/bin/env python3
"""
generate_skills_summary.py — Build SKILLS_SUMMARY.md from the skill registry
============================================================================
Renders a domain-categorised inventory of every ClawBio skill, combining the
live CLI registry with the curated taxonomy in ``scripts/skill_domains.json``.

Usage:
    python scripts/generate_skills_summary.py           # write SKILLS_SUMMARY.md
    python scripts/generate_skills_summary.py --check   # verify it is current
    python scripts/generate_skills_summary.py --output - # print to stdout

Exit codes:
    0  success (or, under --check, the file is up to date)
    1  the domain taxonomy does not cover the registry exactly
    2  under --check, SKILLS_SUMMARY.md is missing or stale

Why the registry and not skills/catalog.json
--------------------------------------------
``catalog.json`` is a per-folder index and deliberately drops entries via
``generate_catalog.EXCLUDED_FOLDERS``, so it currently describes 97 folders
while ``clawbio.py list`` reports 100 skills. It also keys on folder names,
whereas the runnable surface is keyed on CLI aliases (``drugphoto`` and
``drug-photo`` are one folder but two distinct list entries). Deriving the
tiers from the same functions the CLI uses keeps this report consistent with
what a user actually sees. The catalog is still read, but only to report drift
between the two views.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

CLAWBIO_DIR = Path(__file__).resolve().parents[1]
SKILLS_DIR = CLAWBIO_DIR / "skills"
CATALOG_PATH = SKILLS_DIR / "catalog.json"
DOMAINS_PATH = Path(__file__).resolve().parent / "skill_domains.json"
OUTPUT_PATH = CLAWBIO_DIR / "SKILLS_SUMMARY.md"

sys.path.insert(0, str(CLAWBIO_DIR))

# Tier labels used throughout the report.
TIER_REGISTERED = "R"
TIER_AGENT = "A"


# --------------------------------------------------------------------------- #
# Registry introspection
# --------------------------------------------------------------------------- #


def _frontmatter_description(skill_md: Path) -> str:
    """Pull the `description` value out of a SKILL.md YAML frontmatter block.

    Deliberately tolerant: descriptions in this repo are variously quoted,
    folded with `>-`, and wrapped across lines, and a hard parse failure here
    would block the whole report over cosmetic YAML.
    """
    try:
        text = skill_md.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""

    match = re.match(r"^---\s*\n(.*?)\n---", text, re.DOTALL)
    block = match.group(1) if match else text[:2000]

    found = re.search(r"^description:\s*(.+?)(?=\n\w+:|\Z)", block, re.DOTALL | re.MULTILINE)
    if not found:
        return ""

    value = " ".join(found.group(1).split())
    # Strip folded-scalar markers and surrounding quotes.
    value = re.sub(r"^[>|][-+]?\s*", "", value)
    return value.strip().strip("\"'").strip()


def _has_python_script(skill_dir: Path) -> bool:
    """Mirror the CLI's `[has script]` / `[spec only]` indicator."""
    for path in skill_dir.rglob("*.py"):
        if path.name == "__init__.py" or "__pycache__" in path.parts:
            continue
        relative = path.relative_to(skill_dir).parts
        if relative and relative[0] == "tests":
            continue
        return True
    return False


def collect_skills() -> dict[str, dict]:
    """Return every skill keyed by its list-facing name.

    Registered skills are keyed by CLI alias, agent-readable skills by folder
    name. `assert_flat_namespace` guards the assumption that these cannot
    collide.
    """
    from clawbio.cli import SKILLS, _skill_md_only_directories

    skills: dict[str, dict] = {}

    for alias, info in SKILLS.items():
        script: Path = info["script"]
        skills[alias] = {
            "name": alias,
            "tier": TIER_REGISTERED,
            "description": info.get("description", "").strip(),
            "script_name": script.name,
            "script_exists": script.exists(),
            "folder": script.parent.name,
            "has_script": True,
            "spec_only": False,
        }

    for skill_dir in _skill_md_only_directories():
        has_script = _has_python_script(skill_dir)
        skills[skill_dir.name] = {
            "name": skill_dir.name,
            "tier": TIER_AGENT,
            "description": _frontmatter_description(skill_dir / "SKILL.md"),
            "script_name": None,
            "script_exists": None,
            "folder": skill_dir.name,
            "has_script": has_script,
            "spec_only": not has_script,
        }

    return skills


def assert_flat_namespace(skills: dict[str, dict]) -> None:
    """Fail if a CLI alias ever collides with an agent-readable folder name.

    The taxonomy file keys on a single flat namespace. If the two ever
    overlap, a domain entry would silently describe the wrong skill, so this
    is a hard error rather than a warning.
    """
    from clawbio.cli import SKILLS, _skill_md_only_directories

    collisions = sorted(set(SKILLS) & {p.name for p in _skill_md_only_directories()})
    if collisions:
        raise SystemExit(
            "CLI aliases collide with agent-readable folder names: "
            + ", ".join(collisions)
            + "\nThe domain taxonomy assumes a flat namespace; disambiguate before regenerating."
        )


# --------------------------------------------------------------------------- #
# Taxonomy
# --------------------------------------------------------------------------- #


def load_taxonomy() -> dict:
    return json.loads(DOMAINS_PATH.read_text(encoding="utf-8"))


def validate_coverage(taxonomy: dict, skills: dict[str, dict]) -> None:
    """Every skill in exactly one domain, and no phantom entries.

    Exits non-zero with an actionable diff. A new skill should break this
    rather than quietly vanish from the report.
    """
    assigned: list[str] = []
    for domain in taxonomy["domains"]:
        assigned.extend(domain["members"])

    duplicates = sorted({name for name in assigned if assigned.count(name) > 1})
    unclassified = sorted(set(skills) - set(assigned))
    unknown = sorted(set(assigned) - set(skills))

    problems: list[str] = []
    if unclassified:
        problems.append(
            f"{len(unclassified)} skill(s) missing from {DOMAINS_PATH.name}: "
            + ", ".join(unclassified)
        )
    if unknown:
        problems.append(
            f"{len(unknown)} entry/entries in {DOMAINS_PATH.name} match no skill: "
            + ", ".join(unknown)
        )
    if duplicates:
        problems.append("assigned to more than one domain: " + ", ".join(duplicates))

    if problems:
        for problem in problems:
            print(f"error: {problem}", file=sys.stderr)
        raise SystemExit(1)


def catalog_drift(skills: dict[str, dict]) -> tuple[int, list[str]]:
    """Compare folder coverage against skills/catalog.json.

    Returns the catalog's skill count and the folders this report covers that
    the catalog omits. Reported rather than reconciled: the exclusions are
    intentional on the catalog side.
    """
    if not CATALOG_PATH.exists():
        return 0, []
    try:
        catalog = json.loads(CATALOG_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return 0, []

    catalogued = {entry.get("name") for entry in catalog.get("skills", [])}
    folders = {skill["folder"] for skill in skills.values()}
    return int(catalog.get("skill_count", len(catalogued))), sorted(folders - catalogued)


# --------------------------------------------------------------------------- #
# Rendering
# --------------------------------------------------------------------------- #


def _format_member(skill: dict) -> str:
    description = skill["description"]
    warning = ""
    if skill["tier"] == TIER_REGISTERED and not skill["script_exists"]:
        warning = f" **⚠️ script missing: `{skill['script_name']}`**"
    if description:
        return f"`{skill['name']}` — {description}{warning}"
    return f"`{skill['name']}`{warning}"


def render(taxonomy: dict, skills: dict[str, dict]) -> str:
    total = len(skills)
    registered = [s for s in skills.values() if s["tier"] == TIER_REGISTERED]
    agent = [s for s in skills.values() if s["tier"] == TIER_AGENT]
    domains = taxonomy["domains"]

    lines: list[str] = []
    add = lines.append

    add("# ClawBio Skill Library — Domain Report")
    add("")
    add(
        f"**{total} skills across {len(domains)} domains.** "
        f"{len(registered)} registered, {len(agent)} agent-readable."
    )
    add("")
    add(
        "Generated by `scripts/generate_skills_summary.py` from the live skill "
        "registry in `clawbio/cli.py` and the domain taxonomy in "
        "`scripts/skill_domains.json`. Do not edit by hand — rerun the script."
    )
    add("")
    add(
        "**Tier key:** **[R]** registered, runnable via `clawbio run <skill>` · "
        "**[A]** agent-readable, SKILL.md methodology only"
    )
    add("")
    add("---")
    add("")

    for index, domain in enumerate(domains, start=1):
        members = [skills[name] for name in domain["members"]]
        n_reg = sum(1 for m in members if m["tier"] == TIER_REGISTERED)
        n_agt = len(members) - n_reg

        add(f"## {index}. {domain['title']} — {len(members)} ({n_reg}R / {n_agt}A)")
        add("")
        if domain.get("note"):
            add(domain["note"])
            add("")

        for tier, label in ((TIER_REGISTERED, "[R]"), (TIER_AGENT, "[A]")):
            tier_members = [m for m in members if m["tier"] == tier]
            if not tier_members:
                continue
            add(f"**{label}**")
            add("")
            for member in tier_members:
                add(f"- {_format_member(member)}")
            add("")

    add("---")
    add("")
    add("# Cross-cutting findings")
    add("")

    # Broken registrations — computed, not asserted.
    broken = sorted(s["name"] for s in registered if not s["script_exists"])
    add("## Registered skills with a missing script")
    add("")
    if broken:
        for name in broken:
            skill = skills[name]
            add(f"- `{name}` — `{skill['script_name']}` not found in `skills/{skill['folder']}/`")
        add("")
        add(
            f"{len(registered) - len(broken)} of {len(registered)} registered "
            "scripts resolve on disk."
        )
    else:
        add("None. Every registered skill resolves to a script on disk.")
    add("")

    # The registration gap.
    with_script = sorted(s["name"] for s in agent if s["has_script"])
    spec_only = sorted(s["name"] for s in agent if s["spec_only"])
    add("## The registration gap")
    add("")
    add(
        f"{len(with_script)} of {len(agent)} agent-readable skills ship a Python "
        "script but have no entry in the `SKILLS` dict in `clawbio/cli.py`, so "
        "they cannot be reached through `clawbio run`."
    )
    add("")
    if total:
        add(f"That is {len(with_script) / total:.0%} of the library's implemented functionality.")
        add("")
    if spec_only:
        add(
            f"Genuinely spec-only ({len(spec_only)}): "
            + ", ".join(f"`{name}`" for name in spec_only)
            + "."
        )
        add("")

    # Coverage asymmetry — computed from the taxonomy.
    all_registered_domains = [
        d for d in domains
        if all(skills[m]["tier"] == TIER_REGISTERED for m in d["members"])
    ]
    all_agent_domains = [
        d for d in domains
        if all(skills[m]["tier"] == TIER_AGENT for m in d["members"])
    ]
    add("## Coverage asymmetry")
    add("")
    if all_registered_domains:
        add(
            "Fully runnable domains: "
            + ", ".join(d["title"] for d in all_registered_domains)
            + "."
        )
        add("")
    if all_agent_domains:
        add(
            "Domains with no CLI path at all: "
            + ", ".join(d["title"] for d in all_agent_domains)
            + "."
        )
        add("")

    # Curated observations that cannot be derived from metadata.
    findings = taxonomy.get("curated_findings", {})
    pairs = findings.get("duplicate_pairs", [])
    if pairs:
        add("## Duplicate pairs worth consolidating")
        add("")
        for pair in pairs:
            add(f"- {pair}")
        add("")

    notes = findings.get("notes", [])
    if notes:
        add("## Other observations")
        add("")
        for note in notes:
            add(f"- {note}")
        add("")

    # Catalog drift.
    catalog_count, missing_from_catalog = catalog_drift(skills)
    if catalog_count:
        add("## Drift against `skills/catalog.json`")
        add("")
        add(
            f"The catalog indexes {catalog_count} skill folders; this report covers "
            f"{total} list entries. The two counts differ by design — the catalog is "
            "keyed on folders and applies `EXCLUDED_FOLDERS`, while the CLI is keyed "
            "on aliases."
        )
        if missing_from_catalog:
            add("")
            add(
                "Folders present here but absent from the catalog: "
                + ", ".join(f"`{name}`" for name in missing_from_catalog)
                + "."
            )
        add("")

    # Domain totals.
    add("---")
    add("")
    add("## Domain totals")
    add("")
    ranked = sorted(domains, key=lambda d: (-len(d["members"]), d["title"]))
    for domain in ranked:
        members = [skills[name] for name in domain["members"]]
        n_reg = sum(1 for m in members if m["tier"] == TIER_REGISTERED)
        add(f"- {domain['title']} — {len(members)} ({n_reg}R / {len(members) - n_reg}A)")
    add("")
    add(f"**Total: {total} skills ({len(registered)} registered, {len(agent)} agent-readable)**")
    add("")

    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify SKILLS_SUMMARY.md matches the registry without writing",
    )
    parser.add_argument(
        "--output",
        default=str(OUTPUT_PATH),
        help="output path, or - for stdout (default: SKILLS_SUMMARY.md)",
    )
    args = parser.parse_args()

    skills = collect_skills()
    assert_flat_namespace(skills)
    taxonomy = load_taxonomy()
    validate_coverage(taxonomy, skills)

    content = render(taxonomy, skills)

    if args.check:
        if not OUTPUT_PATH.exists():
            print(f"error: {OUTPUT_PATH.name} does not exist; run without --check", file=sys.stderr)
            raise SystemExit(2)
        if OUTPUT_PATH.read_text(encoding="utf-8") != content:
            print(
                f"error: {OUTPUT_PATH.name} is out of date; "
                "run python scripts/generate_skills_summary.py",
                file=sys.stderr,
            )
            raise SystemExit(2)
        print(f"{OUTPUT_PATH.name} is up to date — {len(skills)} skills.")
        return

    if args.output == "-":
        sys.stdout.write(content)
        return

    destination = Path(args.output)
    destination.write_text(content, encoding="utf-8")
    registered = sum(1 for s in skills.values() if s["tier"] == TIER_REGISTERED)
    try:
        shown = destination.relative_to(CLAWBIO_DIR)
    except ValueError:
        shown = destination
    print(
        f"Wrote {shown} — {len(skills)} skills "
        f"({registered} registered, {len(skills) - registered} agent-readable) "
        f"across {len(taxonomy['domains'])} domains"
    )


if __name__ == "__main__":
    main()
