"""Embed one maintained standard-library runtime in portable remap scripts."""

from __future__ import annotations

import ast
import shutil
from pathlib import Path

from clawbio.common.textio import write_text_lf

_RUNTIME_START = "# >>> SHARED NFCORE RUNTIME >>>"
_RUNTIME_END = "# <<< SHARED NFCORE RUNTIME <<<"
_RUNTIME_PATH = Path(__file__).with_name("bundle_runtime.py")
# The runtime functions adapters call, bound under private aliases in both modes.
_EXPORTS = (
    "write_text_lf",
    "prefix_matches",
    "remap_csv",
    "run_actions",
    "repair_bundle",
    "format_output_value",
    "escape_double_quoted",
)


def _bindings(statements: list[ast.stmt]) -> dict[str, object]:
    """Top-level names bound by ``statements``; an import maps to what it imports."""
    names: dict[str, object] = {}
    for node in statements:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names[node.name] = node
        elif isinstance(node, ast.Import):
            for alias in node.names:
                bound = alias.asname or alias.name.split(".")[0]
                names[bound] = ("import", alias.name if alias.asname else bound)
        elif isinstance(node, ast.ImportFrom):
            for alias in node.names:
                names[alias.asname or alias.name] = ("import", f"{node.module}.{alias.name}")
        elif isinstance(node, (ast.Assign, ast.AnnAssign, ast.AugAssign, ast.For, ast.With)):
            targets = (
                node.targets
                if isinstance(node, ast.Assign)
                else [item.optional_vars for item in node.items if item.optional_vars]
                if isinstance(node, ast.With)
                else [node.target]
            )
            for target in targets:
                for sub in ast.walk(target):
                    if isinstance(sub, ast.Name):
                        names[sub.id] = node
        for field in ("body", "orelse", "finalbody"):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                names.update(_bindings(getattr(node, field, None) or []))
        for handler in getattr(node, "handlers", None) or []:
            names.update(_bindings(handler.body))
    return names


def _same_import(left: object, right: object) -> bool:
    return isinstance(left, tuple) and left == right


def _check_no_shadowing(adapter_source: str, runtime_source: str) -> None:
    """Refuse adapters whose names would change the pasted runtime's behaviour.

    In-repo, the adapter imports the runtime as a module, so the runtime resolves
    its own globals. Pasted into the adapter, both share one namespace: a name the
    adapter binds before the block is overwritten by the runtime, and one bound
    after it replaces whatever the runtime's functions call by that name.
    """
    runtime = ast.parse(runtime_source)
    runtime_names = _bindings(runtime.body)
    dependencies = {
        node.id
        for function in runtime.body
        if isinstance(function, ast.FunctionDef)
        for node in ast.walk(function)
        if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load)
    } & set(runtime_names)
    aliases = {f"_nfcore_{name}" for name in _EXPORTS}

    adapter = ast.parse(adapter_source)
    start_line = adapter_source[: adapter_source.index(_RUNTIME_START)].count("\n") + 1
    end_line = adapter_source[: adapter_source.index(_RUNTIME_END)].count("\n") + 1
    before = _bindings([n for n in adapter.body if n.lineno < start_line])
    after = _bindings([n for n in adapter.body if n.lineno > end_line])

    conflicts = {
        name
        for name, binding in before.items()
        if name in runtime_names and not _same_import(binding, runtime_names[name])
    } | {
        name
        for name, binding in after.items()
        if name in dependencies and not _same_import(binding, runtime_names[name])
    } | ((set(before) | set(after)) & aliases)
    if conflicts:
        raise ValueError(
            "Adapter names would shadow the pasted nf-core runtime: "
            + ", ".join(sorted(conflicts))
        )


def render_standalone_remapper(adapter_source: str) -> str:
    """Replace the adapter's runtime import block with the runtime source itself."""
    if (
        adapter_source.count(_RUNTIME_START) != 1
        or adapter_source.count(_RUNTIME_END) != 1
    ):
        raise ValueError("Expected exactly one shared nf-core runtime block.")
    runtime_source = _RUNTIME_PATH.read_text(encoding="utf-8")
    _check_no_shadowing(adapter_source, runtime_source)
    first = adapter_source.index(_RUNTIME_START)
    last = adapter_source.index(_RUNTIME_END, first) + len(_RUNTIME_END)
    aliases = "\n".join(f"_nfcore_{name} = {name}" for name in _EXPORTS)
    generated = (
        adapter_source[:first] + runtime_source + "\n" + aliases + adapter_source[last:]
    )
    compile(generated, "remap_paths.py", "exec")
    return generated


def write_standalone_remapper(destination: Path, *, adapter_path: Path) -> Path:
    """Write the single-file remapper for a bundle, keeping the adapter's mode bits."""
    generated = render_standalone_remapper(adapter_path.read_text(encoding="utf-8"))
    destination.parent.mkdir(parents=True, exist_ok=True)
    write_text_lf(destination, generated)
    try:
        shutil.copymode(adapter_path, destination)
    except OSError:
        pass  # Mode bits are a convenience; replay runs `python remap_paths.py`.
    return destination
