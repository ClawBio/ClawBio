#!/usr/bin/env python3
"""Lint SKILL.md frontmatter for OpenClaw discoverability.

Checks every skills/*/SKILL.md for the metadata.openclaw block and validates
that the os field uses Node.js process.platform values (darwin, linux, win32)
rather than informal names (macos, windows).

Exit code 0  = all skills pass
Exit code 1  = at least one skill has a warning or error

Outputs a Markdown summary table to stdout (suitable for GITHUB_STEP_SUMMARY).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

VALID_OS_VALUES = {"darwin", "linux", "win32", "freebsd", "openbsd", "sunos", "aix"}
COMMON_MISTAKES = {
    "macos": "darwin",
    "mac": "darwin",
    "osx": "darwin",
    "windows": "win32",
    "win": "win32",
    "ubuntu": "linux",
    "debian": "linux",
    "centos": "linux",
    "rhel": "linux",
    "fedora": "linux",
}

REQUIRED_FIELDS = ["os"]
RECOMMENDED_FIELDS = ["emoji", "homepage"]

SKILLS_DIR = Path(__file__).resolve().parent.parent / "skills"


def parse_yaml_simple(yaml_str: str) -> dict:
    """Minimal YAML parser — enough for flat/nested keys we care about.

    We avoid importing pyyaml so this script runs with zero dependencies.
    """
    try:
        import yaml  # noqa: F811

        return yaml.safe_load(yaml_str) or {}
    except ImportError:
        pass

    result: dict = {}

    m = re.search(r"^name:\s*(.+)", yaml_str, re.MULTILINE)
    if m:
        result["name"] = m.group(1).strip().strip('"').strip("'")

    m = re.search(r"^description:\s*(.+)", yaml_str, re.MULTILINE)
    if m:
        result["description"] = m.group(1).strip().strip('"').strip("'")

    metadata_match = re.search(r"metadata:\s*\n((?:[ \t]+.+\n)*)", yaml_str)
    metadata_block = metadata_match.group(1) if metadata_match else ""

    def _extract_block(parent: str, key: str) -> str:
        match = re.search(rf"{key}:\s*\n((?:[ \t]+.+\n)*)", parent)
        return match.group(1) if match else ""

    if metadata_block:
        metadata: dict = {}
        for field in ("version", "author", "domain"):
            field_match = re.search(rf"^\s*{field}:\s*(.+)", metadata_block, re.MULTILINE)
            if field_match:
                metadata[field] = field_match.group(1).strip().strip('"').strip("'")

        openclaw_block = _extract_block(metadata_block, "openclaw")
        if openclaw_block:
            openclaw: dict = {}

            inline_os_match = re.search(r"os:\s*\[([^\]]*)\]", openclaw_block)
            if inline_os_match:
                openclaw["os"] = [
                    v.strip().strip('"').strip("'")
                    for v in inline_os_match.group(1).split(",")
                    if v.strip()
                ]
            else:
                block_os_match = re.search(r"os:\s*\n((?:\s*-\s*.+\n)+)", openclaw_block)
                if block_os_match:
                    openclaw["os"] = [
                        line.strip().lstrip("- ").strip('"').strip("'")
                        for line in block_os_match.group(1).splitlines()
                        if line.strip()
                    ]

            em_match = re.search(r'emoji:\s*["\']?([^"\'\n]+)', openclaw_block)
            if em_match:
                openclaw["emoji"] = em_match.group(1).strip()

            hp_match = re.search(r"homepage:\s*(\S+)", openclaw_block)
            if hp_match:
                openclaw["homepage"] = hp_match.group(1).strip()

            al_match = re.search(r"always:\s*(\S+)", openclaw_block)
            if al_match:
                openclaw["always"] = al_match.group(1).strip().lower() == "true"

            kw_section = re.search(
                r"trigger_keywords:\s*\n((?:\s+-\s+.+\n)*)", openclaw_block
            )
            if kw_section:
                openclaw["trigger_keywords"] = [
                    line.strip().lstrip("- ").strip('"').strip("'")
                    for line in kw_section.group(1).strip().splitlines()
                    if line.strip()
                ]

            metadata["openclaw"] = openclaw

        result["metadata"] = metadata

    if "version" not in result:
        m = re.search(r"^version:\s*(.+)", yaml_str, re.MULTILINE)
        if m:
            result["version"] = m.group(1).strip().strip('"').strip("'")
    if "author" not in result:
        m = re.search(r"^author:\s*(.+)", yaml_str, re.MULTILINE)
        if m:
            result["author"] = m.group(1).strip().strip('"').strip("'")
    if "domain" not in result:
        m = re.search(r"^domain:\s*(.+)", yaml_str, re.MULTILINE)
        if m:
            result["domain"] = m.group(1).strip().strip('"').strip("'")

    return result


def normalize_skill_metadata(raw: dict) -> dict:
    metadata = raw.get("metadata") if isinstance(raw.get("metadata"), dict) else {}
    openclaw = metadata.get("openclaw") if isinstance(metadata.get("openclaw"), dict) else {}
    return {
        "version": raw.get("version", metadata.get("version", "0.1.0")),
        "author": raw.get("author", metadata.get("author", "")),
        "domain": raw.get("domain", metadata.get("domain", "")),
        "openclaw": openclaw,
    }


def lint_skill(skill_dir: Path) -> dict:
    """Lint a single skill. Returns a result dict."""
    skill_md = skill_dir / "SKILL.md"
    name = skill_dir.name
    result = {
        "name": name,
        "has_skill_md": skill_md.exists(),
        "errors": [],
        "warnings": [],
        "os": [],
        "emoji": "",
        "discoverable": False,
    }

    if not skill_md.exists():
        result["errors"].append("No SKILL.md found")
        return result

    text = skill_md.read_text(encoding="utf-8")
    fm = re.match(r"^---\n(.*?\n)---", text, re.DOTALL)
    fm = fm.group(1) if fm else None
    if fm is None:
        result["errors"].append("No YAML frontmatter (--- markers)")
        return result

    data = normalize_skill_metadata(parse_yaml_simple(fm))
    openclaw = data.get("openclaw") or {}

    if not openclaw:
        result["errors"].append("No metadata.openclaw block")
        return result

    # Check os field
    os_values = openclaw.get("os", [])
    if not os_values:
        result["errors"].append("Missing `os` field — OpenClaw cannot gate by platform")
    else:
        for val in os_values:
            if val in COMMON_MISTAKES:
                result["errors"].append(
                    f'os: "{val}" should be "{COMMON_MISTAKES[val]}" '
                    f"(Node.js process.platform)"
                )
            elif val not in VALID_OS_VALUES:
                result["warnings"].append(f'os: "{val}" is not a known platform')
        result["os"] = os_values

    # Check recommended fields
    if not openclaw.get("emoji"):
        result["warnings"].append("Missing `emoji`")
    else:
        result["emoji"] = openclaw["emoji"]

    if not openclaw.get("homepage"):
        result["warnings"].append("Missing `homepage`")

    # Discoverable = has openclaw block + valid os field with no errors
    result["discoverable"] = len(result["errors"]) == 0

    return result


def main() -> int:
    if not SKILLS_DIR.is_dir():
        print(f"ERROR: skills directory not found: {SKILLS_DIR}", file=sys.stderr)
        return 1

    skill_dirs = sorted(
        d for d in SKILLS_DIR.iterdir()
        if d.is_dir() and not d.name.startswith(".")
    )

    # Skip non-skill directories (e.g. catalog.json is a file, not a dir)
    results = [lint_skill(d) for d in skill_dirs]

    # Print Markdown table
    total = len(results)
    discovered = sum(1 for r in results if r["discoverable"])
    errored = sum(1 for r in results if r["errors"])
    warned = sum(1 for r in results if r["warnings"] and not r["errors"])

    print("## OpenClaw Skill Discovery Report")
    print()
    print(f"**{discovered}/{total}** skills discoverable | "
          f"**{errored}** with errors | **{warned}** with warnings")
    print()
    print("| Skill | Status | OS | Emoji | Issues |")
    print("|-------|--------|----|-------|--------|")

    for r in results:
        if r["discoverable"] and not r["warnings"]:
            status = "Pass"
        elif r["discoverable"]:
            status = "Warn"
        else:
            status = "FAIL"

        os_str = ", ".join(r["os"]) if r["os"] else "-"
        emoji = r["emoji"] or "-"
        issues = "; ".join(r["errors"] + r["warnings"]) or "-"

        print(f"| {r['name']} | {status} | {os_str} | {emoji} | {issues} |")

    print()

    if errored:
        print(f"**{errored} skill(s) will NOT be discovered by OpenClaw.** "
              "Fix the errors above.")

    return 1 if errored else 0


if __name__ == "__main__":
    sys.exit(main())
