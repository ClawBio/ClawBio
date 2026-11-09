"""LD reference execution skill (single-client on-demand mode).

Computes r² between a lead variant and partner variants in a chromosomal
window using the 1000 Genomes Phase 3 GRCh38 reference panel, ancestry-
stratified by super-population.

The skill ships ONE client, `OnDemand1000GLDClient` (in `ondemand_client.py`),
which tabix-fetches the region VCF from EBI 1000G FTP per-call (~5-50 MB),
super-pop-filters, and shells out to plink 1.9 for the r² compute. No
multi-GB pre-baked panel download is required; a fresh install renders an
LD-coloured regional plot in seconds.

License posture:

- 1000G Phase 3 GRCh38 data: open access with attribution (Auton 2015
  Nature, Clarke 2017 NAR). Treated as Green-with-attribution.
- plink 1.9 binary: GPL-3 standalone. Subprocess invocation is FSF mere
  aggregation (https://www.gnu.org/licenses/gpl-faq.html#MereAggregation),
  not linkage. Codebase remains MIT. Do NOT bundle the binary in our
  wheels. Users install via brew / apt / conda.

This module exposes:
- Public schema dataclasses: `LDPair`, `LDResult`, `SuperPop` enum.
- `LDComputeError` for upstream catch-blocks.
- A CLI entry point (`main`) that wires the on-demand client and writes
  `ld_pairs.tsv`, `manifest.yaml`, and `report.md` outputs.

Caching is the orchestrator's job. The wrapper itself only caches region
VCFs under `~/.clawbio/locuscompare_cache/1000g/`.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path


class SuperPop(str, Enum):
    """1000G Phase 3 super-populations.

    Per-study ancestry tag drives the choice. EUR is the explicit fallback
    for unspecified-ancestry studies.
    """
    EUR = "EUR"
    AFR = "AFR"
    AMR = "AMR"
    EAS = "EAS"
    SAS = "SAS"


@dataclass
class LDPair:
    """r² between the lead and one partner variant."""

    partner_variant_id: str  # chr_pos_ref_alt (matches OT join key)
    r2: float
    dprime: float | None = None  # plink also reports D'; optional


@dataclass
class LDResult:
    panel_id: str  # e.g. "1000g_phase3_v5b_grch38_basic"
    panel_version: str  # e.g. "5b_remote_2019_03_12"
    super_pop: SuperPop
    plink_version: str  # binary version string returned by `plink --version`
    chromosome: str
    lead_variant_id: str
    window_bp: int
    n_partners_requested: int
    n_partners_returned: int
    pairs: list[LDPair]
    fetched_at_utc: str
    notes: list[str] = field(default_factory=list)


class LDComputeError(Exception):
    """Raised when r² computation cannot proceed (plink missing, fetch fails, etc.)."""


__all__ = [
    "LDComputeError",
    "LDPair",
    "LDResult",
    "SuperPop",
]


# ---------------------------------------------------------------------------
# CLI entry point: --input <config> --output <dir> --demo.
# ---------------------------------------------------------------------------

import argparse  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402

# When run directly as a script, ensure the script's own directory is on
# sys.path so the sibling `ondemand_client` module resolves.
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

# Local cache for LD-compute outputs (region VCFs are cached separately by
# OnDemand1000GLDClient under ~/.clawbio/locuscompare_cache/1000g/).
DEFAULT_LD_RESULT_CACHE_DIR = Path(
    os.environ.get(
        "LD_1000G_RESULT_CACHE_DIR",
        Path.home() / ".clawbio" / "ld_1000g_region_compute_cache",
    )
).expanduser()


def main(argv: list[str] | None = None) -> int:
    """Standard skill CLI: --input <config> --output <dir> --demo.

    Config schema (JSON or YAML):
        lead: 1_109274968_G_T
        partners:
          - 1_109270398_G_A
          - 1_109274570_A_G
        chromosome: "1"
        window_bp: 1000000
        super_pop: EUR
        # optional plink binary override:
        plink_bin: /usr/local/bin/plink

    Writes <output>/{ld_pairs.tsv, manifest.yaml, report.md}.
    """
    parser = argparse.ArgumentParser(
        prog="ld-1000g-region-compute",
        description="Compute pairwise r² between a lead and partner variants using 1000G Phase 3 GRCh38, ancestry-stratified.",
    )
    parser.add_argument("--input", type=Path, help="JSON or YAML config (see docstring).")
    parser.add_argument("--output", type=Path,
                        help="Output directory; created if missing. Required unless --list-demos.")
    parser.add_argument("--demo", nargs="?", const="__default__", default=None,
                        metavar="NAME",
                        help="Run a bundled demo. Bare --demo runs the default; "
                             "pass a name to choose a specific one. See --list-demos.")
    parser.add_argument("--list-demos", action="store_true",
                        help="List bundled demo configs in this skill's examples/ directory.")
    parser.add_argument("--no-cache", action="store_true",
                        help="Bypass the LD-result cache; re-compute every time.")
    args = parser.parse_args(argv)

    if args.list_demos:
        _print_available_demos()
        return 0
    if args.demo is None and args.input is None:
        parser.error("either --input <config> or --demo [NAME] or --list-demos is required")
    if args.output is None:
        parser.error("--output is required")
    args.output.mkdir(parents=True, exist_ok=True)

    if args.demo is not None:
        cfg_path = _resolve_demo_path(args.demo)
        cfg = _load_config(cfg_path)
        print(f"info: using bundled demo {cfg_path.name}", file=sys.stderr)
    else:
        cfg = _load_config(args.input)

    from ondemand_client import (
        DEFAULT_PLINK_BIN,
        OnDemand1000GLDClient,
    )

    super_pop_str = cfg.get("super_pop", "EUR")
    plink_bin = cfg.get("plink_bin") or DEFAULT_PLINK_BIN
    lead = cfg["lead"]
    partners = [p for p in cfg["partners"] if p != lead]
    chromosome = str(cfg["chromosome"])
    window_bp = int(cfg["window_bp"])

    cache_dir = None if args.no_cache else DEFAULT_LD_RESULT_CACHE_DIR
    cache_path = None
    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        import hashlib
        partners_hash = hashlib.sha1(("\n".join(sorted(partners))).encode()).hexdigest()[:12]
        cache_path = cache_dir / f"{lead}__{super_pop_str}__win{window_bp}__{partners_hash}.json"
        if cache_path.is_file():
            cached = json.loads(cache_path.read_text())
            _write_ld_outputs(cached, args.output)
            print(f"ld-1000g-region-compute: {cached['n_partners_returned']} pairs (cache hit) -> {args.output / 'ld_pairs.tsv'}")
            return 0

    client = OnDemand1000GLDClient(
        super_pop=super_pop_str,
        plink_bin=plink_bin,
    )

    result = client.r2_with_lead(
        lead=lead,
        partners=partners,
        chromosome=chromosome,
        window_bp=window_bp,
    )
    payload = {
        "lead": lead,
        "chromosome": chromosome,
        "window_bp": window_bp,
        "super_pop": getattr(result.super_pop, "value", str(result.super_pop)),
        "panel_id": result.panel_id,
        "panel_version": result.panel_version,
        "plink_version": result.plink_version,
        "n_partners_requested": result.n_partners_requested,
        "n_partners_returned": result.n_partners_returned,
        "fetched_at_utc": result.fetched_at_utc,
        "notes": list(result.notes),
        "pairs": [{"partner_variant_id": p.partner_variant_id,
                   "r2": p.r2,
                   "dprime": getattr(p, "dprime", None)} for p in result.pairs],
    }
    if cache_path is not None:
        cache_path.write_text(json.dumps(payload, default=str))
    _write_ld_outputs(payload, args.output)
    print(f"ld-1000g-region-compute: {result.n_partners_returned} pairs -> {args.output / 'ld_pairs.tsv'}")
    print(f"  panel: {result.panel_id} ({payload['super_pop']}) | plink {result.plink_version}")
    return 0


def _load_config(path: Path) -> dict:
    text = path.read_text()
    if path.suffix.lower() in (".yaml", ".yml"):
        import yaml as _yaml
        return _yaml.safe_load(text) or {}
    if path.suffix.lower() == ".json":
        return json.loads(text)
    raise ValueError(f"unsupported config extension: {path.suffix}")


def _examples_dir() -> Path:
    return Path(__file__).resolve().parent / "examples"


def _list_demos() -> list[Path]:
    out: list[Path] = []
    for ext in ("*.json", "*.yaml", "*.yml"):
        out.extend(sorted(_examples_dir().glob(ext)))
    return [p for p in out if p.name not in {"expected_output.md", "README.md"}]


def _resolve_demo_path(name: str) -> Path:
    examples = _examples_dir()
    if name == "__default__":
        for cand in ("default.json", "default.yaml", "default.yml", "input.json"):
            p = examples / cand
            if p.is_file():
                return p
        files = _list_demos()
        if not files:
            raise FileNotFoundError(f"no bundled demo configs found in {examples}")
        return files[0]
    for ext in (".json", ".yaml", ".yml", ""):
        p = examples / (name if ext == "" else f"{name}{ext}")
        if p.is_file():
            return p
    available = ", ".join(p.stem for p in _list_demos())
    raise FileNotFoundError(
        f"no bundled demo named {name!r} in {examples}. Available: {available}"
    )


def _print_available_demos() -> None:
    paths = _list_demos()
    if not paths:
        print(f"no bundled demos in {_examples_dir()}")
        return
    try:
        default_path = _resolve_demo_path("__default__")
    except FileNotFoundError:
        default_path = None
    print(f"Bundled demos in {_examples_dir()}:")
    for p in paths:
        marker = " (default)" if default_path is not None and p == default_path else ""
        print(f"  {p.stem}{marker}    [{p.name}]")


def _write_ld_outputs(payload: dict, output: Path) -> None:
    tsv_path = output / "ld_pairs.tsv"
    cols = ["lead_variant_id", "partner_variant_id", "r2", "dprime", "panel_id", "super_pop"]
    with tsv_path.open("w") as f:
        f.write("# ld-1000g-region-compute v0.1.1\n")
        f.write(f"# panel: {payload['panel_id']}\n")
        f.write(f"# super_pop: {payload['super_pop']}\n")
        f.write(f"# plink: {payload['plink_version']}\n")
        f.write("\t".join(cols) + "\n")
        for p in payload["pairs"]:
            row = [
                payload["lead"],
                p["partner_variant_id"],
                f"{p['r2']:.6f}",
                "" if p.get("dprime") is None else f"{p['dprime']:.6f}",
                payload["panel_id"],
                str(payload["super_pop"]),
            ]
            f.write("\t".join(row) + "\n")

    manifest = {
        "skill": "ld-1000g-region-compute",
        "version": "0.1.1",
        "lead": payload["lead"],
        "chromosome": payload["chromosome"],
        "window_bp": payload["window_bp"],
        "super_pop": payload["super_pop"],
        "panel_id": payload["panel_id"],
        "panel_version": payload["panel_version"],
        "plink_version": payload["plink_version"],
        "n_partners_requested": payload["n_partners_requested"],
        "n_partners_returned": payload["n_partners_returned"],
        "fetched_at_utc": payload["fetched_at_utc"],
        "outputs": {"ld_pairs_tsv": "ld_pairs.tsv"},
        "notes": payload.get("notes") or [],
    }
    try:
        import yaml as _yaml
        (output / "manifest.yaml").write_text(_yaml.safe_dump(manifest, sort_keys=False))
    except ImportError:
        (output / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str))

    report = [
        "# ld-1000g-region-compute report",
        "",
        f"- **Lead:** `{payload['lead']}`",
        f"- **Region:** chr{payload['chromosome']} ±{payload['window_bp']//2//1000} kb",
        f"- **Reference panel:** {payload['panel_id']} ({payload['super_pop']})",
        f"- **plink:** {payload['plink_version']}",
        f"- **Partners requested / returned:** {payload['n_partners_requested']} / {payload['n_partners_returned']}",
        f"- **Output TSV:** ld_pairs.tsv",
    ]
    (output / "report.md").write_text("\n".join(report) + "\n")


if __name__ == "__main__":
    sys.exit(main())
