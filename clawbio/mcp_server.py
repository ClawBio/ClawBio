"""Model Context Protocol server for ClawBio.

Exposes the ClawBio skill library to any MCP-capable client (Cursor, Zed,
VS Code agent mode, Claude Desktop, and others) through three tools:

    list_skills(query)      discover skills
    describe_skill(name)    read a skill's SKILL.md contract
    run_skill(...)          execute a skill and return its structured result

Claude Code users should prefer the plugin marketplace, which ships the skills
themselves rather than a remote-call shim:

    /plugin marketplace add ClawBio/ClawBio

Safety
------
`run_skill` runs demo data freely. Pointing a skill at a local file is refused
unless CLAWBIO_MCP_ALLOW_LOCAL_FILES=1 is set, so that connecting the server
cannot, on its own, hand an agent a patient genome.

The functions below are deliberately free of any `mcp` SDK import so they stay
testable without the optional dependency. `serve()` holds the only binding.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

from clawbio.cli import run_skill as _cli_run_skill

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILLS_DIR = REPO_ROOT / "skills"
CATALOG_PATH = SKILLS_DIR / "catalog.json"

ALLOW_LOCAL_FILES_ENV = "CLAWBIO_MCP_ALLOW_LOCAL_FILES"
MAX_OUTPUT_CHARS = 20_000

_LIST_FIELDS = (
    "name",
    "cli_alias",
    "description",
    "status",
    "maturity_tier",
    "has_demo",
    "demo_command",
    "tags",
    "runnable",
)


class SkillNotFoundError(KeyError):
    """Raised when a skill name or CLI alias is not in the catalog."""


class SkillNotRunnableError(RuntimeError):
    """Raised for catalog skills that have no CLI entry point (spec-only skills)."""


class LocalFileAccessDenied(PermissionError):
    """Raised when a local path is requested without the opt-in env var."""


def _registry_key(entry: dict) -> str | None:
    """Return the key `clawbio run` dispatches on, or None if not executable.

    The CLI registry is keyed by `cli_alias`, not by the catalog `name`.
    """
    from clawbio.cli import SKILLS

    for candidate in (entry.get("cli_alias"), entry.get("name")):
        if candidate and candidate in SKILLS:
            return candidate
    return None


def catalog_skills() -> list[dict]:
    """Return every catalog entry. The catalog is the authoritative skill list."""
    return json.loads(CATALOG_PATH.read_text())["skills"]


def _haystack(entry: dict) -> str:
    parts = [
        entry.get("name", ""),
        entry.get("cli_alias") or "",
        entry.get("description", ""),
        " ".join(entry.get("tags") or []),
        " ".join(entry.get("trigger_keywords") or []),
    ]
    return " ".join(parts).lower()


def list_skills(query: str = "") -> list[dict]:
    """Return compact catalog entries, optionally filtered by a free-text query."""
    entries = catalog_skills()
    if query:
        needle = query.lower().strip()
        entries = [e for e in entries if needle in _haystack(e)]
    out = []
    for e in entries:
        item = {k: e[k] for k in _LIST_FIELDS if k in e}
        item["runnable"] = _registry_key(e) is not None
        out.append(item)
    return out


def _resolve(name: str) -> dict:
    needle = name.strip().lower()
    for entry in catalog_skills():
        if needle in (entry.get("name", "").lower(), (entry.get("cli_alias") or "").lower()):
            return entry
    raise SkillNotFoundError(
        f"Unknown skill {name!r}. Call list_skills to see the {len(catalog_skills())} available skills."
    )


def describe_skill(name: str) -> dict:
    """Return a skill's catalog metadata plus the full text of its SKILL.md contract."""
    entry = _resolve(name)
    spec_path = SKILLS_DIR / entry["name"] / "SKILL.md"
    spec = spec_path.read_text() if spec_path.exists() else ""
    return {**entry, "spec": spec}


def _truncate(text: str) -> str:
    if text is None or len(text) <= MAX_OUTPUT_CHARS:
        return text or ""
    dropped = len(text) - MAX_OUTPUT_CHARS
    return text[:MAX_OUTPUT_CHARS] + f"\n... [truncated, {dropped} more characters; see output_dir for full results]"


def run_skill(
    skill: str,
    demo: bool = False,
    input_path: str | None = None,
    output_dir: str | None = None,
    extra_args: list[str] | None = None,
    timeout: int = 300,
) -> dict:
    """Run a ClawBio skill and return its structured result.

    Raises SkillNotFoundError for unknown skills and LocalFileAccessDenied when a
    local path is supplied without CLAWBIO_MCP_ALLOW_LOCAL_FILES=1.
    """
    entry = _resolve(skill)

    registry_key = _registry_key(entry)
    if registry_key is None:
        raise SkillNotRunnableError(
            f"Skill {entry['name']!r} is agent-readable only: it ships a SKILL.md contract but no "
            "CLI entry point, so it cannot be executed. Use describe_skill to read its spec."
        )

    if (input_path or output_dir) and os.environ.get(ALLOW_LOCAL_FILES_ENV) not in {"1", "true", "yes"}:
        raise LocalFileAccessDenied(
            "This ClawBio MCP server is restricted to demo data. To let it read and write "
            f"local files, restart it with {ALLOW_LOCAL_FILES_ENV}=1."
        )

    result = _cli_run_skill(
        registry_key,
        input_path=input_path,
        output_dir=output_dir,
        demo=demo,
        extra_args=list(extra_args or []),
        timeout=timeout,
    )
    result = dict(result)
    result["stdout"] = _truncate(result.get("stdout", ""))
    result["stderr"] = _truncate(result.get("stderr", ""))
    return result


def serve() -> None:  # pragma: no cover - thin SDK binding
    """Run the MCP server over stdio. Requires the optional `mcp` dependency."""
    try:
        import mcp  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise SystemExit(
            "The MCP server needs the optional 'mcp' dependency:\n"
            "    pip install 'clawbio[mcp]'"
        ) from exc

    try:
        from mcp.server.fastmcp import FastMCP
    except ImportError as exc:  # pragma: no cover
        from importlib.metadata import version

        raise SystemExit(
            f"Installed mcp {version('mcp')} does not provide mcp.server.fastmcp "
            "(removed in mcp 2.0). ClawBio needs mcp>=1.9,<2:\n"
            "    pip install 'clawbio[mcp]'"
        ) from exc

    # Note: FastMCP reports the mcp SDK version in serverInfo, not ClawBio's.
    # The ClawBio version is available to clients via describe_skill metadata.
    app = FastMCP("clawbio")

    @app.tool()
    def clawbio_list_skills(query: str = "") -> list[dict]:
        """Search the ClawBio bioinformatics skill library. Empty query lists everything."""
        return list_skills(query)

    @app.tool()
    def clawbio_describe_skill(name: str) -> dict:
        """Read a ClawBio skill's SKILL.md contract: inputs, outputs, safety rules, demo command."""
        return describe_skill(name)

    @app.tool()
    def clawbio_run_skill(
        skill: str,
        demo: bool = False,
        input_path: str | None = None,
        output_dir: str | None = None,
        extra_args: list[str] | None = None,
    ) -> dict:
        """Run a ClawBio skill. Demo data always works; local files need CLAWBIO_MCP_ALLOW_LOCAL_FILES=1."""
        return run_skill(
            skill,
            demo=demo,
            input_path=input_path,
            output_dir=output_dir,
            extra_args=extra_args,
        )

    app.run()
