"""gwas-lookup API modules — also provides the importable run() entry point.

Usage:
    import importlib, sys, pathlib
    _skill_dir = pathlib.Path("<project_root>/skills/gwas-lookup")
    if str(_skill_dir) not in sys.path:
        sys.path.insert(0, str(_skill_dir))
    from gwas_lookup_api import run

    result = run(
        genotypes=None,
        options={"rsid": "rs3798220", "output_dir": "/tmp/gwas_demo"},
    )
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

# Ensure project root is importable
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

# Ensure skill directory is importable
_SKILL_DIR = Path(__file__).resolve().parent.parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

from clawbio.common.report import write_result_json

import gwas_lookup as _engine  # noqa: E402  (sibling module in skill dir)


def run(genotypes: dict[str, str] | None = None, options: dict | None = None) -> dict:
    """Run GWAS Lookup for an rsID.

    Args:
        genotypes: Not used directly (GWAS Lookup queries by rsID, not
                   genotype data). Reserved for future use.
        options: Dict with keys:
            - rsid (str): rsID to look up (e.g. "rs3798220"). Required
              unless demo is True.
            - demo (bool): If True, use pre-fetched demo data (default False).
            - skip_apis (list[str]): API names to skip.
            - max_hits (int): Max GWAS associations (default 100).
            - make_figures (bool): Generate matplotlib figures (default True).
            - cache_dir (str): Cache directory path.
            - use_cache (bool): Whether to use cache (default True).
            - output_dir (str): If provided, write report files and result.json.

    Returns:
        Dict with keys: merged (full merged results), summary (key stats).
    """
    options = options or {}

    rsid = options.get("rsid", _engine.DEMO_RSID)
    demo = options.get("demo", False)
    skip_apis = options.get("skip_apis", [])
    max_hits = options.get("max_hits", 100)
    make_figures = options.get("make_figures", True)
    cache_dir = Path(options.get("cache_dir", _engine.DEFAULT_CACHE_DIR))
    use_cache = options.get("use_cache", True)
    output_dir = options.get("output_dir")

    # Load demo data if requested
    demo_data = None
    if demo:
        if _engine.DEMO_DATA_PATH.exists():
            demo_data = json.loads(_engine.DEMO_DATA_PATH.read_text())

    # Determine output directory
    if output_dir:
        output_path = Path(output_dir)
    else:
        import tempfile
        output_path = Path(tempfile.mkdtemp(prefix="gwas_lookup_"))

    merged = _engine.run_lookup(
        rsid=rsid,
        output_dir=output_path,
        cache_dir=cache_dir,
        use_cache=use_cache,
        skip_apis=skip_apis,
        max_hits=max_hits,
        make_figures=make_figures,
        demo_data=demo_data,
    )

    summary = merged.get("summary", {})

    result = {
        "merged": merged,
        "summary": {
            "rsid": rsid,
            "total_gwas": summary.get("total_gwas", 0),
            "total_gwas_significant": summary.get("total_gwas_significant", 0),
            "total_phewas_ukb": summary.get("total_phewas_ukb", 0),
            "total_phewas_finngen": summary.get("total_phewas_finngen", 0),
            "total_phewas_bbj": summary.get("total_phewas_bbj", 0),
            "total_eqtls": summary.get("total_eqtls", 0),
            "total_credible_sets": summary.get("total_credible_sets", 0),
        },
        "output_dir": str(output_path),
    }

    # write_result_json is already called by _engine.run_lookup,
    # so we don't need to call it again here.

    return result
