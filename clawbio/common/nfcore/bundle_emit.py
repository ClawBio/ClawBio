"""Embed one maintained standard-library runtime in portable remap scripts."""

from __future__ import annotations

import shutil
from collections.abc import Mapping
from pathlib import Path

from clawbio.common.textio import write_text_lf

_RUNTIME_START = "# >>> SHARED NFCORE RUNTIME >>>"
_RUNTIME_END = "# <<< SHARED NFCORE RUNTIME <<<"


def render_standalone_remapper(
    *, adapter_source: str, runtime_source: str, substitutions: Mapping[str, str]
) -> str:
    if (
        adapter_source.count(_RUNTIME_START) != 1
        or adapter_source.count(_RUNTIME_END) != 1
    ):
        raise ValueError("Expected exactly one shared nf-core runtime block.")
    first = adapter_source.index(_RUNTIME_START)
    last = adapter_source.index(_RUNTIME_END, first) + len(_RUNTIME_END)
    aliases = "\n".join(
        f"_nfcore_{name} = {name}"
        for name in (
            "write_text_lf",
            "prefix_matches",
            "remap_csv",
            "run_actions",
            "repair_bundle",
            "format_output_value",
            "escape_double_quoted",
        )
    )
    generated = (
        adapter_source[:first] + runtime_source + "\n" + aliases + adapter_source[last:]
    )
    for marker, replacement in substitutions.items():
        if generated.count(marker) != 1:
            raise ValueError(f"Expected exactly one substitution marker: {marker}")
        generated = generated.replace(marker, replacement, 1)
    compile(generated, "remap_paths.py", "exec")
    return generated


def write_standalone_remapper(
    destination: Path, *, adapter_path: Path, substitutions: Mapping[str, str]
) -> Path:
    runtime_source = (
        Path(__file__).with_name("bundle_runtime.py").read_text(encoding="utf-8")
    )
    generated = render_standalone_remapper(
        adapter_source=adapter_path.read_text(encoding="utf-8"),
        runtime_source=runtime_source,
        substitutions=substitutions,
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    write_text_lf(destination, generated)
    shutil.copystat(adapter_path, destination)
    return destination
