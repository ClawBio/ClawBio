"""docs/data-handling.md must name every skill that can send data off the machine.

A prose promise that "networked skills are individually labelled" governs
nothing. This test scans every skill for outbound-call code, or for a declared
remote endpoint, or for use of the hosted Genomic Intelligence client, and
fails if any such skill is missing from the data-handling page. Add a row to
the page when you add a network call; the page is the contract institutions
read.
"""

from __future__ import annotations

import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
SKILLS = ROOT / "skills"
DOC = ROOT / "docs" / "data-handling.md"

CALL = re.compile(
    r"requests\.(get|post|put|delete|request|Session)\(|urllib\.request\.urlopen\(|urlopen\(|"
    r"httpx\.|aiohttp|\bsubprocess\.[a-z_]+\([^)]*\b(curl|wget|nextflow)\b|download\.file\(|"
    r"httr::|bioblend|BioBlend|AsyncOpenAI\(|OpenAI\(|from clawbio\.gi|import clawbio\.gi|"
    r"session\.(get|post)\(|voyageai|google\.cloud\.bigquery|huggingface_hub|hf_hub_download|"
    r"from_pretrained\(|(import|from) labstep|synapseclient|[\"']uvx[\"']|[\"']nextflow[\"']\s*,\s*[\"']run[\"']"
)
SOURCE_SUFFIXES = (".py", ".R", ".r", ".sh")


def skills_that_reach_the_network() -> set[str]:
    found: set[str] = set()
    for skill_dir in sorted(p for p in SKILLS.iterdir() if p.is_dir()):
        name = skill_dir.name
        for path in skill_dir.rglob("*"):
            if not path.is_file() or path.suffix not in SOURCE_SUFFIXES:
                continue
            if "tests" in path.parts or "node_modules" in path.parts:
                continue
            text = path.read_text(errors="ignore")
            if CALL.search(text):
                found.add(name)
                break
        skill_md = skill_dir / "SKILL.md"
        if skill_md.exists():
            match = re.match(r"^---\n(.*?)\n---", skill_md.read_text(errors="ignore"), re.S)
            try:
                front = yaml.safe_load(match.group(1)) if match else {}
            except yaml.YAMLError:
                front = {}
            endpoints = ((front or {}).get("metadata") or {}).get("endpoints")
            if isinstance(endpoints, (list, dict)):
                values = endpoints if isinstance(endpoints, list) else endpoints.values()
                if any(isinstance(v, str) and v.startswith("http") for v in values):
                    found.add(name)
    return found


def documented_skills() -> set[str]:
    text = DOC.read_text()
    return set(re.findall(r"`([a-z0-9][a-z0-9-]*)`", text))


def test_data_handling_page_exists():
    assert DOC.exists(), "docs/data-handling.md is missing"


def test_every_networked_skill_is_documented():
    networked = skills_that_reach_the_network()
    assert networked, "scanner found nothing, which means the scanner is broken"
    missing = sorted(networked - documented_skills())
    assert not missing, (
        "skills that reach the network but are absent from docs/data-handling.md: "
        + ", ".join(missing)
    )


def documented_rows() -> set[str]:
    """Skill names in the first cell of every table row (several rows list more than one)."""
    names: set[str] = set()
    for line in DOC.read_text().splitlines():
        if line.startswith("| `"):
            first_cell = line.split("|")[1]
            names.update(re.findall(r"`([a-z0-9][a-z0-9-]*)`", first_cell))
    return names


def test_documented_skills_exist():
    """A row for a skill that no longer exists is a stale claim."""
    existing = {p.name for p in SKILLS.iterdir() if p.is_dir()}
    stale = sorted(documented_rows() - existing)
    assert not stale, f"docs/data-handling.md has rows for skills that do not exist: {stale}"


def test_rows_cover_every_scanned_skill_or_the_checked_list():
    """Every scanned skill must have a table row, except those the page explicitly
    lists as checked and clean, which must still be named."""
    networked = skills_that_reach_the_network()
    rows = documented_rows()
    text = DOC.read_text()
    checked_section = text.split("## Checked and found to make no outbound call", 1)[-1]
    checked = set(re.findall(r"`([a-z0-9][a-z0-9-]*)`", checked_section.split("## ", 1)[0]))
    unaccounted = sorted(networked - rows - checked)
    assert not unaccounted, f"scanned skills with neither a row nor a checked entry: {unaccounted}"


def test_scanner_still_fires_on_a_known_networked_skill():
    """If the regex rots, the page silently stops being enforced."""
    networked = skills_that_reach_the_network()
    for known in ("vcf-annotator", "clinical-variant-reporter", "gi-annotation", "pathway-enricher"):
        assert known in networked, f"scanner no longer detects {known}"
