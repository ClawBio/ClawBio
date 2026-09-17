"""Hatchling build hook: bundle a curated subset of skills/ into the wheel.

The repository keeps all 84 skills (and their heavy demo data: gzipped genomes,
expected-output images, per-skill data dirs) at the repository root. Shipping all
of that in the pip wheel would push it past 35 MB, so this hook force-includes a
lean subset under the ``clawbio/`` namespace:

  * every skill's logic (.md / .py / .yaml / .sh) so agents can read and call it,
  * small data files (< 256 KB) that ship cheaply,
  * the full demo payload for the headline skills only, so the README quick-start
    (``clawbio run pharmgx --demo``) works offline straight after install,
  * the whole examples/ directory (small, used by several --demo flows).

Files are read from the repository root and mapped to ``clawbio/skills/...`` and
``clawbio/examples/...`` in the wheel only; nothing moves on disk, so the Claude
Code plugin marketplace and every public ``skills/...`` URL stay intact.
"""

from __future__ import annotations

from pathlib import Path

from hatchling.builders.hooks.plugin.interface import BuildHookInterface

# Skills whose full demo payload ships so their --demo runs offline out of the box.
HEADLINE_SKILLS = {"pharmgx-reporter", "drug-photo", "gwas-lookup", "just-prs-mcp"}

# Code licences that may ship inside the MIT wheel.
WHEEL_LICENCES = {"MIT", "Apache-2.0"}

# Skills whose code licence cannot be redistributed inside the MIT wheel. They
# stay in the repository as folders labelled by licence and are excluded here.
# Kept in step with FOLDER_ONLY in tests/test_licence_policy.py by a test.
#
# This list is a declaration, not the gate: _excluded_skill_folders() below
# re-reads every SKILL.md at build time, so a folder that never reaches
# skills/catalog.json (generate_catalog.py has its own EXCLUDED_FOLDERS) still
# cannot slip into the wheel.
WHEEL_EXCLUDED_SKILLS = {"fastreer", "wes-clinical-report-en", "wes-clinical-report-es"}

# Logic/source files included for every skill regardless of size.
LOGIC_SUFFIXES = {".md", ".py", ".yaml", ".yml", ".sh", ".cff", ".toml", ".cfg"}
# Data files included for non-headline skills only when small.
SMALL_DATA_SUFFIXES = {".json", ".csv", ".tsv", ".txt"}
SMALL_DATA_MAX_BYTES = 256 * 1024

SKIP_DIR_NAMES = {
    "__pycache__",
    ".pytest_cache",
    ".ipynb_checkpoints",
    ".vscode",
    ".idea",
    ".codex",
    ".tessl-plugin",
}
SKIP_FILE_NAMES = {".DS_Store", ".gitignore", ".gitkeep"}
SKIP_SUFFIXES = {".pyc", ".pyo"}


def _declared_licence(skill_md: Path) -> str:
    """Read `license:` from a SKILL.md YAML frontmatter block.

    Deliberately a local copy of scripts/generate_catalog.py's frontmatter
    reader: scripts/ is not in the sdist (see tool.hatch.build.targets.sdist),
    so importing it here would break any wheel built from an unpacked sdist.
    """
    lines = skill_md.read_text(encoding="utf-8", errors="ignore").splitlines()
    if not lines or lines[0].strip() != "---":
        return ""
    for line in lines[1:]:
        if line.strip() == "---":
            break
        # lstrip so a licence nested under `metadata:` is found too; the
        # startswith still cannot match `model_license:`, which begins `model_`.
        stripped = line.lstrip()
        if stripped.startswith("license:"):
            return stripped.split(":", 1)[1].strip().strip("\"'")
    return ""


def _excluded_skill_folders(skills_dir: Path) -> set[str]:
    """Every folder under skills/ that must not be published.

    Two reasons, both read off disk because that is what the build walks:

    - No SKILL.md: it is not a skill. Untracked scratch work in a maintainer's
      tree, or a RETIRED/ holding pen. Publishing whatever happens to sit in
      skills/ is how local data reaches PyPI.
    - A code licence that is not MIT or Apache-2.0, blank included: the wheel is
      MIT, and a folder that will not say what it is cannot be shown to be
      redistributable.
    """
    excluded = set(WHEEL_EXCLUDED_SKILLS)
    if not skills_dir.is_dir():
        return excluded
    for folder in skills_dir.iterdir():
        if not folder.is_dir():
            continue
        skill_md = folder / "SKILL.md"
        if not skill_md.is_file() or _declared_licence(skill_md) not in WHEEL_LICENCES:
            excluded.add(folder.name)
    return excluded


def _keep_for_non_headline(path: Path) -> bool:
    """Include skill logic and only small data files for non-headline skills."""
    suffix = path.suffix.lower()
    if suffix in LOGIC_SUFFIXES:
        return True
    if suffix in SMALL_DATA_SUFFIXES:
        try:
            return path.stat().st_size <= SMALL_DATA_MAX_BYTES
        except OSError:
            return False
    return False


class CustomBuildHook(BuildHookInterface):
    PLUGIN_NAME = "custom"

    def initialize(self, version, build_data):  # noqa: D401 - hatchling interface
        root = Path(self.root)
        force_include = build_data.setdefault("force_include", {})
        excluded_skills = _excluded_skill_folders(root / "skills")

        for base in ("skills", "examples"):
            base_dir = root / base
            if not base_dir.is_dir():
                continue
            for path in base_dir.rglob("*"):
                if not path.is_file():
                    continue
                if any(part in SKIP_DIR_NAMES for part in path.parts):
                    continue
                if path.name in SKIP_FILE_NAMES or path.suffix.lower() in SKIP_SUFFIXES:
                    continue

                rel = path.relative_to(root)

                if base == "skills":
                    skill = rel.parts[1] if len(rel.parts) > 1 else ""
                    if skill in excluded_skills:
                        continue
                    if skill not in HEADLINE_SKILLS and not _keep_for_non_headline(path):
                        continue

                force_include[str(path)] = f"clawbio/{rel.as_posix()}"
