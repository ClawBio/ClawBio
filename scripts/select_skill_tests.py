"""Pick the skills whose tests a pull request must run.

Reads changed paths (one per line) on stdin and prints skill folder names.
A skill is selected when the PR changes:
  - a file inside the skill's folder, or inside another skill it reaches by path;
  - a clawbio module it imports, directly or through other clawbio modules;
  - a file that shapes every test environment (ENVIRONMENT_FILES): all skills.
Explicit imports only: clawbio/__init__.py imports the runner, so counting the
package imports Python runs implicitly would select every skill for any clawbio
change. Import-time breakage there fails tests/ and clawbio/tests, which run on
every PR. Every skill's tests still run on the scheduled sweep.
"""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

ENVIRONMENT_FILES = {"uv.lock", "pyproject.toml", "conftest.py", "pytest.ini"}


def _module_name(root: Path, path: Path) -> str:
    parts = list(path.relative_to(root).with_suffix("").parts)
    if parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def _imports(path: Path, module: str, known: set[str]) -> set[str]:
    """clawbio modules a file imports, plus `clawbio.cli` if it runs clawbio.py."""
    text = path.read_text(encoding="utf-8", errors="replace")
    found: set[str] = set()
    try:
        tree = ast.parse(text)
    except SyntaxError:
        tree = ast.Module(body=[], type_ignores=[])
    package = module if path.name == "__init__.py" else module.rpartition(".")[0]
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            base = node.module or ""
            if node.level:
                anchor = package.split(".")[: len(package.split(".")) - node.level + 1]
                base = ".".join(anchor + ([base] if base else []))
            found.add(base)
            # `from clawbio.common import checksums` imports a module, not a name
            found.update(f"{base}.{alias.name}" for alias in node.names)
    if "clawbio.py" in text:
        found.add("clawbio.cli")
    return {name for name in found if name in known}


def _path_segment(name: str) -> re.Pattern[str]:
    n = re.escape(name)
    return re.compile(rf"/{n}/|[\"']{n}[\"']\s*/|/\s*[\"']{n}[\"']")


def _python_files(folder: Path) -> list[Path]:
    return [p for p in folder.rglob("*.py") if "__pycache__" not in p.parts]


def select_skills(root: Path, changed: list[str]) -> list[str]:
    root = Path(root)
    skills_dir = root / "skills"
    skills = sorted(p.name for p in skills_dir.iterdir() if p.is_dir()) if skills_dir.is_dir() else []

    if any(path in ENVIRONMENT_FILES for path in changed):
        return skills

    modules = {_module_name(root, p): p for p in _python_files(root / "clawbio")}
    known = set(modules)

    # Changed clawbio modules. A non-Python file counts as a change to every
    # module in its package, its likely readers; clawbio.py is the CLI's entry point.
    affected: set[str] = set()
    for rel in changed:
        parts = Path(rel).parts
        if rel == "clawbio.py":
            affected.add("clawbio.cli")
        elif parts and parts[0] == "clawbio":
            if rel.endswith(".py"):
                affected.add(_module_name(root, root / rel))
                continue
            package = ".".join(parts[:-1])
            while package and package not in known:
                package = package.rpartition(".")[0]
            affected.update(m for m in known if m == package or m.startswith(f"{package}."))

    # Close over clawbio's own imports: a module is affected if anything it imports is.
    deps = {name: _imports(path, name, known) for name, path in modules.items()}
    grew = True
    while grew:
        grew = False
        for name, imported in deps.items():
            if name not in affected and imported & affected:
                affected.add(name)
                grew = True

    selected = {
        parts[1] for rel in changed
        if len(parts := Path(rel).parts) > 2 and parts[0] == "skills" and parts[1] in skills
    }
    texts: dict[str, str] = {}
    for skill in skills:
        files = _python_files(skills_dir / skill)
        # A path built on tmp_path is a throwaway copy, not a reference to the skill.
        texts[skill] = "\n".join(
            line for p in files
            for line in p.read_text(encoding="utf-8", errors="replace").splitlines()
            if "tmp_path" not in line
        )
        imported: set[str] = set()
        for path in files:
            imported |= _imports(path, "", known)
        if imported & affected:
            selected.add(skill)

    # Skills that reach a selected skill's folder by path: skills/<name>/... or
    # SKILLS_DIR / "<name>". A bare name in a routing table does not count.
    grew = True
    while grew:
        grew = False
        for skill in skills:
            if skill in selected:
                continue
            for other in selected:
                if _path_segment(other).search(texts[skill]):
                    selected.add(skill)
                    grew = True
                    break

    return sorted(selected)


def main() -> None:
    changed = [line.strip() for line in sys.stdin if line.strip()]
    root = Path(__file__).resolve().parents[1]
    print("\n".join(select_skills(root, changed)))


if __name__ == "__main__":
    main()
