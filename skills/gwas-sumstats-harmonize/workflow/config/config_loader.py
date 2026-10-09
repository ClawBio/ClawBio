"""
config/config_loader.py

The single place that reads data.yaml, analysis.yaml and software.yaml for the
harmonization workflow, validates them, resolves relative paths against the
config directory's parent, and returns one plain nested dict (CFG).

The workflow lives inside the skill; the three yaml files live with each run
(the ClawBio wrapper writes them to <output>/config/). Point Snakemake at them
with `--config config_dir=<dir>`.

Example:
    from config_loader import load_config
    CFG = load_config("/runs/my_run/config")
    CFG["data"]["datasets"]["cohort_a"]["path"]        # -> absolute path
    CFG["analysis"]["qc_filter"]["palindromic"]        # -> "ambiguous"
    CFG["analysis"]["resources"]["reference"]          # -> absolute path or None
"""
import os

import yaml


class ConfigError(Exception):
    """Raised when data.yaml / analysis.yaml / software.yaml fail validation."""


STAGES = ("map_columns", "derive_effects", "qc_filter", "align_reference")
CANONICAL = ("SNP", "CHR", "BP", "EA", "NEA", "EAF", "BETA", "SE", "P", "N", "OR", "LOG10P", "Z")

_DEFAULTS = {
    "map_columns": {},
    "derive_effects": {},
    "qc_filter": {"palindromic": "ambiguous", "min_maf": 0.0, "keep_indels": True},
    # With a reference, unmatched rows are dropped by default so every output
    # row really has EA = reference ALT. Without one the stage passes through.
    "align_reference": {"drop_unmatched": True},
}


def _resolve(path, root):
    if path is None:
        return None
    path = str(path)
    resolved = path if os.path.isabs(path) else os.path.normpath(os.path.join(root, path))
    # Forward slashes: backslash paths mixed with the rules' forward-slash
    # suffixes can stop Snakemake matching an input to its producing rule.
    return resolved.replace("\\", "/")


def _load_yaml(path):
    if not os.path.exists(path):
        raise ConfigError(f"Config file not found: {path}")
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def _require(d, key, where):
    if key not in d or d[key] is None:
        raise ConfigError(f"Missing required key '{key}' in {where}")
    return d[key]


def _validate_dataset(name, entry, root):
    entry = dict(entry or {})
    entry["path"] = _resolve(_require(entry, "path", f"data.yaml datasets.{name}"), root)
    if not os.path.exists(entry["path"]):
        raise ConfigError(f"data.yaml datasets '{name}': input file not found: {entry['path']}")
    entry.setdefault("label", name)
    columns = entry.get("columns") or {}
    unknown = sorted(set(k.upper() for k in columns) - set(CANONICAL))
    if unknown:
        raise ConfigError(f"data.yaml datasets '{name}': unknown canonical column(s) {unknown}")
    entry["columns"] = {k.upper(): v for k, v in columns.items()}
    return entry


def _build_analysis(cfg, root):
    analysis = {
        "output_dir": _resolve(_require(cfg, "output_dir", "analysis.yaml"), root),
        "targets": list(_require(cfg, "targets", "analysis.yaml")),
    }
    resources = cfg.get("resources") or {}
    reference = _resolve(resources.get("reference"), root)
    if reference and not os.path.exists(reference):
        raise ConfigError(f"analysis.yaml resources.reference not found: {reference}")
    analysis["resources"] = {"reference": reference}
    for stage in STAGES:
        merged = dict(_DEFAULTS[stage])
        merged.update(cfg.get(stage) or {})
        analysis[stage] = merged
    if analysis["qc_filter"]["palindromic"] not in ("ambiguous", "all", "none"):
        raise ConfigError("analysis.yaml qc_filter.palindromic must be ambiguous, all or none")
    return analysis


def load_config(config_dir):
    if not config_dir:
        raise ConfigError("no config_dir given (snakemake --config config_dir=<dir>)")
    config_dir = os.path.abspath(config_dir)
    root = os.path.dirname(config_dir)

    data_cfg = _load_yaml(os.path.join(config_dir, "data.yaml"))
    analysis_cfg = _load_yaml(os.path.join(config_dir, "analysis.yaml"))
    software_cfg = _load_yaml(os.path.join(config_dir, "software.yaml"))

    raw = _require(data_cfg, "datasets", "data.yaml")
    datasets = {name: _validate_dataset(name, entry, root) for name, entry in raw.items()}
    analysis = _build_analysis(analysis_cfg, root)
    missing = [t for t in analysis["targets"] if t not in datasets]
    if missing:
        raise ConfigError(f"analysis.yaml targets not defined in data.yaml datasets: {missing}")

    software = dict(software_cfg)
    software.setdefault("threads", 1)
    return {"data": {"datasets": datasets}, "analysis": analysis, "software": software}
