#!/usr/bin/env python3
"""
snakemake_bio_scaffold.py — generate a Snakemake bioinformatics project in the
config-by-concern layout.

Layout produced (one stage = one rule file + one standalone script):

    <project>/
    ├── Snakefile                 # load config once, wire the DAG, rule all
    ├── config/
    │   ├── data.yaml             # WHAT to process (items, paths, per-item overrides)
    │   ├── analysis.yaml         # HOW to process (per-stage settings, output_dir, targets)
    │   ├── software.yaml         # WHERE tools live (binaries, threads)
    │   └── config_loader.py      # the only reader of the 3 yaml files
    ├── rules/<stage>.smk         # DAG wiring only, no logic
    ├── scripts/<stage>.py        # standalone CLI per stage (never imports snakemake)
    ├── scripts/lib/pipeline_io.py
    ├── resources/  input/  utilities/
    └── results/                  # <stage>/, done/ sentinels, logs/ (gitignored)

Usage:
    python snakemake_bio_scaffold.py --input spec.yaml --output out/
    python snakemake_bio_scaffold.py --name my_pipe --stages align,count,report --output out/
    python snakemake_bio_scaffold.py --demo --output /tmp/scaffold_demo [--check]
"""

from __future__ import annotations

import argparse
import json
import keyword
import random
import re
import shlex
import subprocess
import sys
import unicodedata
from datetime import date
from pathlib import Path, PurePosixPath, PureWindowsPath

import yaml

_SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SKILL_DIR.parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)

SKILL_NAME = "snakemake-bio-scaffold"
SKILL_VERSION = "0.1.0"
DEMO_SPEC = _SKILL_DIR / "examples" / "demo_spec.yaml"

_IDENT_RE = re.compile(r"^[a-z][a-z0-9_]*$")
_PROJECT_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_-]*$")
_ITEM_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")
_RESERVED_STAGES = {
    "all", "lib", "input", "output", "params", "log", "rule", "shell",
    "config", "workflow", "rules", "checkpoints", "scatter", "gather",
}
_RESERVED_PARAMS = {"input", "out", "summary_json"}


class SpecError(ValueError):
    """Raised when a scaffold spec is invalid."""


# ── Spec handling ─────────────────────────────────────────────────────────────


def _check_single_line(value: str, what: str) -> str:
    """Values that become command-line arguments must not hold control characters:
    on Windows cmd a newline ends the command, whatever the quoting."""
    if any(unicodedata.category(ch) == "Cc" for ch in value):
        raise SpecError(f"{what} contains a control character (newline, tab, NUL, ...); use plain text")
    return value


def _check_ident(value, what: str) -> str:
    if not isinstance(value, str) or not _IDENT_RE.match(value) or keyword.iskeyword(value):
        raise SpecError(
            f"{what} '{value}' must be lower_snake_case (letters, digits, underscore; "
            "start with a letter) and not a Python keyword"
        )
    return value


def validate_spec(raw: dict) -> dict:
    """Validate a scaffold spec and return a normalised copy with defaults filled."""
    if not isinstance(raw, dict):
        raise SpecError("spec must be a YAML mapping")

    project = raw.get("project")
    if not isinstance(project, str) or not _PROJECT_RE.match(project):
        raise SpecError(
            f"project name '{project}' must start with a letter and contain only "
            "letters, digits, '_' or '-'"
        )

    items_raw = raw.get("items") or {}
    items_key = _check_ident(items_raw.get("key", "datasets"), "items.key")
    wildcard = _check_ident(items_raw.get("wildcard", "dataset"), "items.wildcard")
    entries_raw = items_raw.get("entries") or {
        "example": {"path": "input/example.tsv", "label": "Example item (replace me)"}
    }
    entries: dict[str, dict] = {}
    for name, entry in entries_raw.items():
        if not isinstance(name, str) or not _ITEM_RE.match(name):
            raise SpecError(f"item name '{name}' may contain only letters, digits, '_' or '-'")
        entry = dict(entry or {})
        if not entry.get("path"):
            raise SpecError(f"item '{name}' is missing required key 'path'")
        entries[name] = {
            "path": _check_single_line(str(entry["path"]), f"item '{name}' path").replace("\\", "/"),
            "label": entry.get("label", name),
            "description": entry.get("description", ""),
            "columns": entry.get("columns") or {},
            "overrides": entry.get("overrides"),
        }

    stages_raw = raw.get("stages") or []
    if not stages_raw:
        raise SpecError("spec must define at least one stage")
    stages = []
    seen: set[str] = set()
    for s in stages_raw:
        s = {"name": s} if isinstance(s, str) else dict(s or {})
        name = s.get("name")
        if name in _RESERVED_STAGES:
            raise SpecError(f"stage name '{name}' is reserved (clashes with a Snakemake keyword)")
        _check_ident(name, "stage name")
        if name in seen:
            raise SpecError(f"duplicate stage name '{name}'")
        seen.add(name)
        params = dict(s.get("params") or {})
        for key, value in params.items():
            _check_ident(key, f"stage '{name}' param")
            if key in _RESERVED_PARAMS:
                raise SpecError(f"stage '{name}' param '{key}' clashes with a built-in script flag")
            if value is None:
                raise SpecError(
                    f"stage '{name}' param '{key}' is null; give it a default value"
                )
            if not isinstance(value, (bool, int, float, str)):
                raise SpecError(f"stage '{name}' param '{key}' must be a scalar")
            if isinstance(value, str):
                _check_single_line(value, f"stage '{name}' param '{key}'")
        ext = str(s.get("ext", "tsv")).lstrip(".")
        if not re.match(r"^[A-Za-z0-9][A-Za-z0-9.]*$", ext):
            raise SpecError(f"stage '{name}' ext '{ext}' is not a valid file extension")
        stages.append({
            "name": name,
            "description": s.get("description") or f"{name} stage (describe me).",
            "ext": ext,
            "params": params,
        })

    stage_names = {s["name"] for s in stages}
    for item, entry in entries.items():
        for stage, values in (entry["overrides"] or {}).items():
            if stage not in stage_names:
                raise SpecError(f"item '{item}' overrides unknown stage '{stage}'")
            params = next(s["params"] for s in stages if s["name"] == stage)
            for key in values or {}:
                if key not in params:
                    raise SpecError(f"item '{item}' overrides unknown param '{stage}.{key}'")

    software = dict(raw.get("software") or {})
    for key in software:
        _check_ident(key, "software key")

    return {
        "project": project,
        "description": str(raw.get("description") or "").strip(),
        "items": {"key": items_key, "wildcard": wildcard, "entries": entries},
        "stages": stages,
        "software": software,
        "output_dir": _check_output_dir(raw.get("output_dir", "results")),
    }


def _check_output_dir(value) -> str:
    """output_dir must be a relative path that stays inside the project.

    Checked with both path flavours, so 'C:/x', '\\\\server\\share' and '/x'
    are all rejected whichever OS generates or later runs the project.
    """
    text = str(value if value is not None else "").strip()
    win = PureWindowsPath(text)
    if not text or PurePosixPath(text).is_absolute() or win.drive or win.root or ".." in win.parts:
        raise SpecError(
            f"output_dir '{text}' must be a relative path inside the project (no absolute path, drive or '..')"
        )
    return win.as_posix()


def spec_from_flags(name: str, stages: str, items: str | None) -> dict:
    entries = None
    if items:
        entries = {}
        for pair in items.split(","):
            if "=" not in pair:
                raise SpecError(f"--items entry '{pair}' must look like name=path")
            k, v = pair.split("=", 1)
            entries[k.strip()] = {"path": v.strip()}
    raw = {
        "project": name,
        "stages": [s.strip() for s in stages.split(",") if s.strip()],
        "items": {"entries": entries} if entries else {},
    }
    return raw


# ── Templates ─────────────────────────────────────────────────────────────────


_PLACEHOLDER_RE = re.compile(r"@@([A-Z_]+)@@")


def _fill(template: str, **values: str) -> str:
    """Substitute @@NAME@@ placeholders in ONE pass.

    Single pass matters: spec text that happens to contain "@@ARGPARSE@@" must
    stay text, not pull generated code into a string literal.
    """
    missing = set(_PLACEHOLDER_RE.findall(template)) - set(values)
    if missing:
        raise RuntimeError(f"unfilled template placeholders: {sorted(missing)}")
    return _PLACEHOLDER_RE.sub(lambda m: values[m.group(1)], template)


# Spec text reaches generated Python only in two forms, never pasted into a
# string literal or docstring: as comment lines (_comment), or as a literal
# built by repr() (_py_literal).


def _comment(text: str, indent: str = "") -> str:
    """Free text as '#' lines. str.splitlines() splits on every line break
    Python recognises (and more), so no line can end the comment early."""
    lines = str(text).replace("\x00", "").splitlines() or [""]
    return "\n".join(f"{indent}# {line}".rstrip() for line in lines)


def _py_literal(text: str) -> str:
    return repr(str(text))


SNAKEFILE_T = '''"""
Snakefile — @@PROJECT@@ pipeline entry point (description in the comments below).

Pipeline stages (one rule each, in rules/*.smk):
  @@CHAIN@@

Run:
  snakemake -n --cores 1          # dry run: show what would execute
  snakemake --cores 4             # build every target in analysis.yaml
"""
@@DESCRIPTION@@
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), "config"))
from config_loader import load_config

CFG = load_config()

PYTHON = sys.executable  # same interpreter that runs Snakemake, never a bare `python`
OUT = CFG["analysis"]["output_dir"]
TARGETS = CFG["analysis"]["targets"]
ITEMS = CFG["data"]["@@ITEMS_KEY@@"]
SOFTWARE = CFG["software"]


def _item_path(name):
    return ITEMS[name]["path"]


def _item_analysis(name, stage):
    """Per-item stage settings: analysis.yaml values with data.yaml `overrides` applied."""
    return ITEMS[name]["analysis"][stage]


wildcard_constraints:
    @@WILDCARD@@="|".join(re.escape(n) for n in ITEMS),


@@INCLUDES@@


rule all:
    input:
        [f"{OUT}/done/@@LAST_STAGE@@_{n}.done" for n in TARGETS],
'''

RULE_T = '''@@DESC@@
rule @@STAGE@@:
    """Stage @@N@@: @@STAGE@@ (description in the comments above)."""
    input:
        data=@@INPUT_EXPR@@,
@@DONE_INPUT@@    output:
        result=f"{OUT}/@@STAGE@@/{{@@WILDCARD@@}}.@@EXT@@",
        summary=f"{OUT}/@@STAGE@@/{{@@WILDCARD@@}}.summary.json",
        done=touch(f"{OUT}/done/@@STAGE@@_{{@@WILDCARD@@}}.done"),
@@PARAMS_BLOCK@@    log:
        f"{OUT}/logs/@@STAGE@@/{{@@WILDCARD@@}}.log",
    shell:
        # :q quotes each value for the shell, so $(...), backticks and quotes
        # in paths or settings stay literal text.
        '{PYTHON:q} scripts/@@STAGE@@.py '
        '--input {input.data:q} --out {output.result:q} --summary-json {output.summary:q} '
@@PARAM_FLAGS@@        '> {log:q} 2>&1'
'''

SCRIPT_T = '''#!/usr/bin/env python
"""
scripts/@@STAGE@@.py — stage @@N@@ of the pipeline (description in the comments below).

SCAFFOLD STUB: process() currently passes every row through unchanged so the
pipeline runs end to end. Replace its body with the real logic; keep the CLI
(flags mirror config/analysis.yaml) so the rule in rules/@@STAGE@@.smk keeps
working. Run with --help for the stage description and an example command.
"""
@@DESC@@
import argparse
import json
import logging
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "lib"))
from pipeline_io import read_table, str2bool, write_table  # noqa: E402,F401

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

PARAMS = [@@PARAM_NAMES@@]
USAGE = @@USAGE@@


def parse_args():
    p = argparse.ArgumentParser(description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="Input table (tab-separated, optional .gz)")
    p.add_argument("--out", required=True, help="Output table path")
    p.add_argument("--summary-json", required=True, help="Where to write the per-run summary JSON")
@@ARGPARSE@@    return p.parse_args()


def process(header, rows, args):
    """TODO: implement @@STAGE@@. Return (header, rows) for the output table."""
    return header, rows


def main():
    args = parse_args()
    header, rows = read_table(args.input)
    logger.info("Read %d rows from %s", len(rows), args.input)

    out_header, out_rows = process(header, rows, args)
    logger.info("%d in, %d dropped, %d kept", len(rows), len(rows) - len(out_rows), len(out_rows))

    write_table(args.out, out_header, out_rows)
    summary = {
        "stage": "@@STAGE@@",
        "input": args.input,
        "rows_in": len(rows),
        "rows_out": len(out_rows),
        "params": {name: getattr(args, name) for name in PARAMS},
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.summary_json)), exist_ok=True)
    with open(args.summary_json, "w") as fh:
        json.dump(summary, fh, indent=2)
    logger.info("Wrote %s and %s", args.out, args.summary_json)


if __name__ == "__main__":
    main()
'''

PIPELINE_IO = '''"""
scripts/lib/pipeline_io.py

Plain-python table helpers shared by every stage script. Standard library only,
so the scaffold runs in any Python 3 environment; switch a stage to
pandas/polars when it needs them.
"""
import csv
import gzip
import os


def _open(path, mode):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode + "t", newline="")
    return open(path, mode, newline="")


def read_table(path, sep="\\t"):
    """Return (header, rows) for a delimited text file."""
    with _open(path, "r") as fh:
        reader = csv.reader(fh, delimiter=sep)
        header = next(reader, [])
        rows = list(reader)
    return header, rows


def write_table(path, header, rows, sep="\\t"):
    """Write a delimited text file with LF line endings, creating parent dirs."""
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)
    with _open(path, "w") as fh:
        writer = csv.writer(fh, delimiter=sep, lineterminator="\\n")
        writer.writerow(header)
        writer.writerows(rows)


def str2bool(value):
    """argparse type for true/false flags rendered from YAML booleans."""
    if isinstance(value, bool):
        return value
    v = str(value).strip().lower()
    if v in ("true", "1", "yes", "y"):
        return True
    if v in ("false", "0", "no", "n"):
        return False
    raise ValueError(f"not a boolean: {value!r}")
'''

CONFIG_LOADER_T = '''"""
config/config_loader.py

The single place that reads data.yaml, analysis.yaml and software.yaml,
validates them, resolves every relative path against the project root, merges
per-item overrides, and returns one plain nested dict (CFG). No other file in
this project should call yaml.safe_load on these three files.

Example:
    from config_loader import load_config
    CFG = load_config()
    CFG["data"]["@@ITEMS_KEY@@"]["<name>"]["path"]                 # -> absolute path
    CFG["data"]["@@ITEMS_KEY@@"]["<name>"]["analysis"]["<stage>"]  # -> merged settings
    CFG["analysis"]["output_dir"]                                  # -> absolute path
"""
import os

import yaml


class ConfigError(Exception):
    """Raised when data.yaml / analysis.yaml / software.yaml fail validation."""


PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ITEMS_KEY = "@@ITEMS_KEY@@"
STAGES = (@@STAGE_TUPLE@@)


def _resolve(path):
    if path is None:
        return None
    path = str(path)
    resolved = path if os.path.isabs(path) else os.path.normpath(os.path.join(PROJECT_ROOT, path))
    # Forward slashes everywhere: on Windows, backslash-separated absolute paths
    # mixed with the rules' forward-slash suffixes (f"{OUT}/<stage>/...") can stop
    # Snakemake from matching an input to the rule that produces it, giving a
    # spurious MissingInputException on a first build. Forward slashes work for
    # Python, R, PLINK and the shell on every OS.
    return resolved.replace("\\\\", "/")


def _load_yaml(path):
    if not os.path.exists(path):
        raise ConfigError(f"Config file not found: {path}")
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def _require(d, key, where):
    if key not in d or d[key] is None:
        raise ConfigError(f"Missing required key '{key}' in {where}")
    return d[key]


def _merge(overrides, defaults):
    merged = dict(defaults)
    if overrides:
        merged.update(overrides)
    return merged


def _validate_item(name, entry, analysis):
    """Validate one data.yaml item and attach its per-stage merged settings."""
    entry = dict(entry or {})
    if not entry.get("path"):
        raise ConfigError(f"data.yaml {ITEMS_KEY} '{name}' is missing required key 'path'")
    entry["path"] = _resolve(entry["path"])
    if not os.path.exists(entry["path"]):
        raise ConfigError(f"data.yaml {ITEMS_KEY} '{name}': input file not found: {entry['path']}")
    entry.setdefault("label", name)
    entry.setdefault("description", "")
    entry["columns"] = entry.get("columns") or {}
    overrides = entry.get("overrides") or {}
    for stage, values in overrides.items():
        if stage not in STAGES:
            raise ConfigError(f"data.yaml {ITEMS_KEY} '{name}' overrides unknown stage '{stage}'")
        unknown = set(values or {}) - set(analysis[stage])
        if unknown:
            raise ConfigError(
                f"data.yaml {ITEMS_KEY} '{name}' overrides unknown {stage} setting(s): {sorted(unknown)}"
            )
    entry["overrides"] = overrides or None
    entry["analysis"] = {stage: _merge(overrides.get(stage), analysis[stage]) for stage in STAGES}
    return entry


def _build_analysis(analysis_cfg):
    analysis = {
        "output_dir": _resolve(_require(analysis_cfg, "output_dir", "analysis.yaml")),
        "targets": list(_require(analysis_cfg, "targets", "analysis.yaml")),
        "resources": {k: _resolve(v) for k, v in (analysis_cfg.get("resources") or {}).items()},
    }
    for stage in STAGES:
        analysis[stage] = dict(analysis_cfg.get(stage) or {})
    return analysis


def _build_software(software_cfg):
    software = {}
    for key, value in software_cfg.items():
        # Bare command names (found on PATH) stay as-is; anything path-like is resolved.
        if isinstance(value, str) and ("/" in value or "\\\\" in value):
            value = _resolve(value)
        software[key] = value
    software.setdefault("threads", 1)
    return software


def load_config(config_dir=None):
    config_dir = config_dir or os.path.join(PROJECT_ROOT, "config")

    data_cfg = _load_yaml(os.path.join(config_dir, "data.yaml"))
    analysis_cfg = _load_yaml(os.path.join(config_dir, "analysis.yaml"))
    software_cfg = _load_yaml(os.path.join(config_dir, "software.yaml"))

    analysis = _build_analysis(analysis_cfg)

    raw_items = _require(data_cfg, ITEMS_KEY, "data.yaml")
    items = {name: _validate_item(name, entry, analysis) for name, entry in raw_items.items()}

    missing = [t for t in analysis["targets"] if t not in items]
    if missing:
        raise ConfigError(f"analysis.yaml targets not defined in data.yaml {ITEMS_KEY}: {missing}")

    return {
        "project_root": _resolve(PROJECT_ROOT),
        "data": {ITEMS_KEY: items},
        "analysis": analysis,
        "software": _build_software(software_cfg),
    }
'''

README_T = '''# @@PROJECT@@

@@DESCRIPTION@@

Generated by ClawBio `snakemake-bio-scaffold` (config-by-concern layout).

## Pipeline

```
@@CHAIN@@
```

## Layout

```
Snakefile                 entry point: load config once, include rules, define rule all
config/data.yaml          WHAT to process: @@ITEMS_KEY@@, file paths, column maps, per-item overrides
config/analysis.yaml      HOW to process: per-stage settings, output_dir, targets, resources
config/software.yaml      WHERE tools live: binaries, threads
config/config_loader.py   the only reader of the three yaml files (validates, resolves paths)
rules/<stage>.smk         DAG wiring only: paths, params, one shell line calling a script
scripts/<stage>.py        one standalone CLI per stage; never imports snakemake
scripts/lib/              helpers shared by several scripts
resources/                reference data the pipeline reads but never writes
input/                    raw input data
utilities/                one-off scripts kept OUTSIDE the DAG
results/                  all output (gitignored): <stage>/, done/, logs/
```

## Run

```bash
snakemake -n --cores 1                       # dry run
snakemake --cores 4                          # build every target in analysis.yaml
snakemake --cores 4 -R <stage>               # rerun a stage after editing its script
python scripts/<stage>.py --help             # run or debug one stage by hand
```

## Conventions

- Every rule writes `results/done/<stage>_<@@WILDCARD@@>.done`; downstream rules depend on it.
- Every rule launches scripts with `"{PYTHON}"` (= `sys.executable`), never a bare `python`.
- `{@@WILDCARD@@}` is constrained to the closed set of names in `data.yaml`.
- Thresholds live in `analysis.yaml`; per-item exceptions go under that item's `overrides:`.
- Scripts log "N in, M dropped, K kept" for every filtering step.

## Gotchas

- Snakemake tracks the rule's text, not the script it calls. After editing
  `scripts/<stage>.py`, force a rerun with `-R <stage>`.
- Deleting a stage's output but not its `results/done/` sentinel leaves the DAG
  thinking that stage is finished. Delete both.
- Editing a file many rules share as input (a sample list, a region file)
  invalidates every completed job that reads it. When adding one new item, build
  just its targets: `snakemake --cores 4 results/done/@@LAST_STAGE@@_<name>.done`.
- On Windows, pass absolute target paths on the command line if Snakemake
  reports `MissingRuleException` for a relative target that clearly exists.
- On Windows, `.snakemake/metadata` file names are the base64 of each output's
  full path, so a deeply nested project can exceed MAX_PATH ("Error recording
  metadata for finished job"). Keep the project near the drive root, enable
  long paths, or run with `--drop-metadata`.
- Close Excel outputs before rerunning the rule that writes them.
'''

GITIGNORE = """results/
.snakemake/
__pycache__/
*.pyc
logs/
"""

UTILITIES_README = """# utilities/

Standalone, run-by-hand scripts (one-off audits, ad hoc QC, exploratory
checks). They follow the same CLI and config conventions as `scripts/` but are
NOT wired into any rule and NOT part of `rule all`.
"""

RESOURCES_README = """# resources/

Static reference data the pipeline reads but never writes (reference panels,
gene annotations, LD scores, ...). Point to these from `analysis.yaml`
`resources:` so every path is resolved once by `config_loader.py`.
"""

INPUT_README = """# input/

Raw input data. Register each file in `config/data.yaml`.
"""


# ── Rendering ─────────────────────────────────────────────────────────────────


def _render_rule(spec: dict, idx: int) -> str:
    stage = spec["stages"][idx]
    wc = spec["items"]["wildcard"]
    if idx == 0:
        input_expr = f"lambda wc: _item_path(wc.{wc})"
        done_input = ""
    else:
        prev = spec["stages"][idx - 1]
        input_expr = f'f"{{OUT}}/{prev["name"]}/{{{{{wc}}}}}.{prev["ext"]}"'
        done_input = f'        done=f"{{OUT}}/done/{prev["name"]}_{{{{{wc}}}}}.done",\n'
    params = stage["params"]
    if params:
        lines = ["    params:"]
        for key in params:
            lines.append(f'        {key}=lambda wc: _item_analysis(wc.{wc}, "{stage["name"]}")["{key}"],')
        params_block = "\n".join(lines) + "\n"
        flags = " ".join(f'--{k.replace("_", "-")} {{params.{k}:q}}' for k in params)
        param_flags = f"        '{flags} '\n"
    else:
        params_block = ""
        param_flags = ""
    return _fill(
        RULE_T,
        STAGE=stage["name"],
        N=str(idx + 1),
        DESC=_comment(f"Stage {idx + 1}: {stage['description']}"),
        INPUT_EXPR=input_expr,
        DONE_INPUT=done_input,
        WILDCARD=wc,
        EXT=stage["ext"],
        PARAMS_BLOCK=params_block,
        PARAM_FLAGS=param_flags,
    )


def _argparse_line(key: str, value) -> str:
    flag = "--" + key.replace("_", "-")
    if isinstance(value, bool):
        typ = "str2bool"
    elif isinstance(value, int):
        typ = "int"
    elif isinstance(value, float):
        typ = "float"
    else:
        typ = "str"
    # argparse %-formats help strings, so a literal % must be doubled.
    help_text = _py_literal("analysis.yaml default: " + str(value).replace("%", "%%"))
    return f'    p.add_argument("{flag}", type={typ}, required=True, help={help_text})\n'


def _usage(spec: dict, idx: int) -> str:
    """The stage's --help text: description plus an example command line."""
    stage = spec["stages"][idx]
    entries = spec["items"]["entries"]
    ex_item = next(iter(entries))
    if idx == 0:
        ex_input = entries[ex_item]["path"]
    else:
        prev = spec["stages"][idx - 1]
        ex_input = f"results/{prev['name']}/{ex_item}.{prev['ext']}"
    name, ext = stage["name"], stage["ext"]
    out_path = f"results/{name}/{ex_item}.{ext}"
    summary_path = f"results/{name}/{ex_item}.summary.json"
    lines = [
        f"python scripts/{name}.py",
        "--input " + shlex.quote(ex_input),
        "--out " + shlex.quote(out_path),
        "--summary-json " + shlex.quote(summary_path),
    ] + ["--" + k.replace("_", "-") + " " + shlex.quote(str(v)) for k, v in stage["params"].items()]
    text = f"Stage {idx + 1}: {stage['description']}\n\nExample:\n    " + " \\\n        ".join(lines)
    # RawDescriptionHelpFormatter %-formats the description only when it holds "%(prog)".
    return text.replace("%", "%%") if "%(prog)" in text else text


def _render_script(spec: dict, idx: int) -> str:
    stage = spec["stages"][idx]
    return _fill(
        SCRIPT_T,
        STAGE=stage["name"],
        N=str(idx + 1),
        DESC=_comment(f"Stage {idx + 1}: {stage['description']}"),
        USAGE=_py_literal(_usage(spec, idx)),
        PARAM_NAMES=", ".join(f'"{k}"' for k in stage["params"]),
        ARGPARSE="".join(_argparse_line(k, v) for k, v in stage["params"].items()),
    )


def _chain(spec: dict) -> str:
    return " -> ".join(s["name"] for s in spec["stages"])


def _dump_yaml(header: str, payload: dict) -> str:
    return header + yaml.safe_dump(payload, sort_keys=False, default_flow_style=False)


def render_project(spec: dict) -> dict[str, str]:
    """Return {relative_path: file_content} for the whole generated project."""
    items_key = spec["items"]["key"]
    wc = spec["items"]["wildcard"]
    stages = spec["stages"]
    description = spec["description"] or f"{spec['project']} Snakemake pipeline."
    files: dict[str, str] = {}

    files["Snakefile"] = _fill(
        SNAKEFILE_T,
        PROJECT=spec["project"],
        DESCRIPTION=_comment(description),
        CHAIN=_chain(spec),
        ITEMS_KEY=items_key,
        WILDCARD=wc,
        INCLUDES="\n".join(f'include: "rules/{s["name"]}.smk"' for s in stages),
        LAST_STAGE=stages[-1]["name"],
    )

    data_payload = {items_key: {}}
    for name, entry in spec["items"]["entries"].items():
        data_payload[items_key][name] = {
            "path": entry["path"],
            "label": entry["label"],
            "description": entry["description"],
            "columns": entry["columns"],
            "overrides": entry["overrides"],
        }
    files["config/data.yaml"] = _dump_yaml(
        "# data.yaml — WHAT to process.\n"
        f"# One entry per item under `{items_key}`: its input `path` (relative to the\n"
        "# project root or absolute), a display `label`, an optional `columns` map\n"
        "# (canonical name -> this file's real column name) and optional `overrides`\n"
        "# of analysis.yaml stage settings for this item only, e.g.\n"
        "#   overrides:\n"
        "#     <stage>: {<setting>: <value>}\n\n",
        data_payload,
    )

    analysis_payload = {
        "output_dir": spec["output_dir"],
        "targets": list(spec["items"]["entries"]),
        "resources": {},
    }
    for s in stages:
        analysis_payload[s["name"]] = dict(s["params"])
    files["config/analysis.yaml"] = _dump_yaml(
        "# analysis.yaml — HOW to process.\n"
        "# output_dir: where results/ goes. targets: which data.yaml items `rule all` builds.\n"
        "# resources: shared reference files (resolved to absolute paths by config_loader.py).\n"
        "# One block per stage holds that stage's settings; each becomes a --flag on its script.\n\n",
        analysis_payload,
    )

    software_payload = dict(spec["software"])
    software_payload.setdefault("threads", 1)
    files["config/software.yaml"] = _dump_yaml(
        "# software.yaml — WHERE tools live.\n"
        "# Bare names are looked up on PATH; paths are resolved against the project root.\n"
        "# Python is always the interpreter running Snakemake (sys.executable), not set here.\n\n",
        software_payload,
    )

    files["config/config_loader.py"] = _fill(
        CONFIG_LOADER_T,
        ITEMS_KEY=items_key,
        STAGE_TUPLE="".join(f'"{s["name"]}", ' for s in stages).rstrip(" "),
    )

    for i, s in enumerate(stages):
        files[f"rules/{s['name']}.smk"] = _render_rule(spec, i)
        files[f"scripts/{s['name']}.py"] = _render_script(spec, i)
    files["scripts/lib/pipeline_io.py"] = PIPELINE_IO

    files["README.md"] = _fill(
        README_T,
        PROJECT=spec["project"],
        DESCRIPTION=description,
        CHAIN=_chain(spec),
        ITEMS_KEY=items_key,
        WILDCARD=wc,
        LAST_STAGE=stages[-1]["name"],
    )
    files[".gitignore"] = GITIGNORE
    files["utilities/README.md"] = UTILITIES_README
    files["resources/README.md"] = RESOURCES_README
    files["input/README.md"] = INPUT_README
    return files


def synthetic_sumstats(seed: int, n: int = 200) -> str:
    """Deterministic synthetic GWAS summary statistics (no real data)."""
    rng = random.Random(seed)
    lines = ["\t".join(["SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N"])]
    alleles = ["A", "C", "G", "T"]
    for i in range(n):
        chrom = 1 + i * 22 // n
        bp = 1_000_000 + i * 50_000 + rng.randint(0, 40_000)
        a1, a2 = rng.sample(alleles, 2)
        frq = round(rng.uniform(0.001, 0.5), 4)
        se = round(rng.uniform(0.01, 0.08), 4)
        beta = round(rng.gauss(0, se * (8 if i % 40 == 0 else 1)), 4)
        p = 10 ** -rng.uniform(9, 12) if i % 40 == 0 else rng.uniform(1e-6, 1.0)
        lines.append("\t".join([
            f"synth_rs{i + 1}", str(chrom), str(bp), a1, a2,
            str(frq), str(beta), str(se), f"{p:.4g}", "5000",
        ]))
    return "\n".join(lines) + "\n"


def write_project(files: dict[str, str], project_dir: Path, force: bool) -> list[str]:
    if project_dir.exists() and any(project_dir.iterdir()) and not force:
        raise SpecError(
            f"{project_dir} already exists and is not empty; refusing to overwrite. "
            "Pass --force to overwrite the scaffold files (other files are left alone)."
        )
    for rel, content in files.items():
        path = project_dir / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(content)
    return sorted(files)


def dry_run_check(project_dir: Path) -> dict:
    try:
        import importlib.util

        if importlib.util.find_spec("snakemake") is None:
            return {"status": "skipped", "detail": "snakemake is not installed in this Python"}
    except Exception as exc:  # pragma: no cover
        return {"status": "skipped", "detail": str(exc)}
    proc = subprocess.run(
        [sys.executable, "-m", "snakemake", "-n", "--cores", "1"],
        cwd=project_dir,
        capture_output=True,
        text=True,
    )
    tail = (proc.stdout + proc.stderr).strip().splitlines()[-15:]
    return {
        "status": "passed" if proc.returncode == 0 else "failed",
        "detail": "\n".join(tail),
    }


# ── Report ────────────────────────────────────────────────────────────────────


_ROLES = {
    "Snakefile": "entry point: load config, include rules, rule all",
    "config/data.yaml": "WHAT: items, paths, per-item overrides",
    "config/analysis.yaml": "HOW: per-stage settings, output_dir, targets",
    "config/software.yaml": "WHERE: tool paths, threads",
    "config/config_loader.py": "only reader of the yaml files",
    "scripts/lib/pipeline_io.py": "shared table I/O helpers",
    "README.md": "project conventions and run commands",
}


def write_report(output_dir: Path, spec: dict, files: list[str], dry_run: dict) -> Path:
    project = spec["project"]
    wc = spec["items"]["wildcard"]
    stages = spec["stages"]
    lines = [
        "# Snakemake Bio Scaffold Report",
        "",
        f"**Project**: `{project}`  ",
        f"**Generated**: {date.today().isoformat()}  ",
        f"**Stages**: {len(stages)} · **{spec['items']['key']}**: {len(spec['items']['entries'])} · "
        f"**Wildcard**: `{{{wc}}}`",
        "",
        "## Pipeline",
        "",
        "```",
        _chain(spec),
        "```",
        "",
        "| # | Stage | Output | Settings (analysis.yaml) |",
        "|---|-------|--------|--------------------------|",
    ]
    for i, s in enumerate(stages, 1):
        settings = ", ".join(f"`{k}={v}`" for k, v in s["params"].items()) or "none"
        lines.append(f"| {i} | `{s['name']}` | `results/{s['name']}/{{{wc}}}.{s['ext']}` | {settings} |")
    lines += [
        "",
        "## Files generated",
        "",
        "| Path | Role |",
        "|------|------|",
    ]
    for rel in files:
        role = _ROLES.get(rel)
        if role is None and rel.startswith("rules/"):
            role = "DAG wiring for one stage"
        elif role is None and rel.startswith("scripts/"):
            role = "standalone CLI stub for one stage (edit `process()`)"
        elif role is None and rel.startswith("input/") and not rel.endswith(".md"):
            role = "synthetic demo input"
        elif role is None:
            role = "placeholder / guidance"
        lines.append(f"| `{project}/{rel}` | {role} |")
    lines += [
        "",
        "## Next steps",
        "",
        "```bash",
        f"cd {project}",
        "snakemake -n --cores 1        # dry run",
        "snakemake --cores 4           # build every target",
        "```",
        "",
        "1. Point each item in `config/data.yaml` at its real input file.",
        "2. Replace the pass-through `process()` in each `scripts/<stage>.py` with real logic.",
        "3. After editing a script, rerun that stage with `snakemake --cores 4 -R <stage>`.",
        "",
        "## Dry-run check",
        "",
        f"**Status**: {dry_run['status']}",
    ]
    if dry_run.get("detail"):
        lines += ["", "```", dry_run["detail"], "```"]
    lines += ["", "---", "", f"*{DISCLAIMER}*", ""]
    path = output_dir / "report.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


# ── CLI ───────────────────────────────────────────────────────────────────────


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Generate a Snakemake bioinformatics project (config-by-concern layout)."
    )
    p.add_argument("--input", help="Scaffold spec YAML (see examples/demo_spec.yaml)")
    p.add_argument("--output", required=True, help="Output directory; project goes in <output>/<project>")
    p.add_argument("--demo", action="store_true", help="Use the bundled demo spec and synthetic inputs")
    p.add_argument("--name", help="Quick mode: project name (instead of --input)")
    p.add_argument("--stages", help="Quick mode: comma-separated stage names, in order")
    p.add_argument("--items", help="Quick mode: comma-separated name=path items")
    p.add_argument("--force", action="store_true", help="Overwrite scaffold files in a non-empty project dir and an earlier report in --output")
    p.add_argument("--check", action="store_true", help="Run `snakemake -n` on the generated project")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        if args.demo:
            raw = yaml.safe_load(DEMO_SPEC.read_text(encoding="utf-8"))
        elif args.input:
            raw = yaml.safe_load(Path(args.input).read_text(encoding="utf-8"))
        elif args.name and args.stages:
            raw = spec_from_flags(args.name, args.stages, args.items)
        else:
            raise SpecError("give --input <spec.yaml>, --name + --stages, or --demo")
        spec = validate_spec(raw)

        output_dir = Path(args.output).resolve()
        clobbered = [n for n in ("report.md", "result.json", "reproducibility") if (output_dir / n).exists()]
        if clobbered and not args.force:
            raise SpecError(
                f"{output_dir} already holds {', '.join(clobbered)} from an earlier run; refusing to "
                "overwrite. Pass --force to overwrite, or choose another --output."
            )
        output_dir.mkdir(parents=True, exist_ok=True)
        project_dir = output_dir / spec["project"]
        files = render_project(spec)
        if args.demo:
            for i, name in enumerate(spec["items"]["entries"]):
                files[spec["items"]["entries"][name]["path"]] = synthetic_sumstats(seed=42 + i)
        written = write_project(files, project_dir, args.force)
    except (SpecError, OSError, yaml.YAMLError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    dry_run = dry_run_check(project_dir) if args.check else {
        "status": "skipped", "detail": "pass --check to run `snakemake -n` on the project",
    }

    report = write_report(output_dir, spec, written, dry_run)
    result = {
        "skill": SKILL_NAME,
        "version": SKILL_VERSION,
        "project": spec["project"],
        "project_dir": project_dir.as_posix(),
        "items_key": spec["items"]["key"],
        "wildcard": spec["items"]["wildcard"],
        "targets": list(spec["items"]["entries"]),
        "stages": spec["stages"],
        "rule_count": len(spec["stages"]),
        "files": written,
        "dry_run": dry_run,
    }
    result_path = output_dir / "result.json"
    result_path.write_text(json.dumps(result, indent=2), encoding="utf-8")

    cmd_args = shlex.join(sys.argv[1:] if argv is None else list(argv))
    write_commands_sh(output_dir, f"python skills/{SKILL_NAME}/snakemake_bio_scaffold.py {cmd_args}")
    write_environment_yml(
        output_dir, env_name="clawbio-snakemake-bio-scaffold",
        pip_deps=["pyyaml>=6.0", "snakemake>=8"],
    )
    write_checksums(
        [report, result_path] + [project_dir / rel for rel in written],
        output_dir,
        anchor=output_dir,
    )
    print(f"Scaffolded {spec['project']} ({len(spec['stages'])} stages) -> {project_dir}")
    print(f"Report: {report}")
    if dry_run["status"] == "failed":
        print("warning: snakemake dry run failed; see report.md", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
