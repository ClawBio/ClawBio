from __future__ import annotations

import sys
from pathlib import Path

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

from errors import ErrorCode, SkillError
from schemas import PREFERRED_H5AD_ORDER


def _first_match(root: Path, pattern: str) -> str:
    matches = sorted(root.glob(pattern))
    return str(matches[0]) if matches else ""


def parse_outputs(output_dir: Path) -> dict[str, object]:
    upstream_dir = _require_upstream_results_dir(output_dir)
    h5ad_candidates = find_h5ad_candidates(upstream_dir)
    rds_candidates = find_rds_outputs(upstream_dir)
    preferred_h5ad = select_preferred_h5ad(h5ad_candidates)

    return {
        "upstream_dir": str(upstream_dir),
        "multiqc_report": find_multiqc_report(upstream_dir),
        "pipeline_info_dir": find_pipeline_info_dir(upstream_dir),
        "h5ad_candidates": h5ad_candidates,
        "rds_candidates": rds_candidates,
        "preferred_h5ad": preferred_h5ad,
        "preferred_h5ad_selection_log": build_selection_log(h5ad_candidates, preferred_h5ad),
        "handoff_available": bool(preferred_h5ad),
        "samples_detected": detect_sample_names(h5ad_candidates),
        "cellbender_used": detect_cellbender_outputs(h5ad_candidates),
    }


def build_selection_log(h5ad_candidates: list[str], preferred_h5ad: str) -> str:
    if not h5ad_candidates:
        return "No h5ad output files found."
    if not preferred_h5ad:
        return (
            f"No canonical h5ad selected — {len(h5ad_candidates)} candidate(s) found "
            "but selection is ambiguous. Inspect h5ad_candidates in result.json."
        )
    name = Path(preferred_h5ad).name
    all_names = [Path(c).name for c in h5ad_candidates]
    return f"Selected '{name}' from {len(h5ad_candidates)} candidate(s). All candidates: {all_names}."


def _require_upstream_results_dir(output_dir: Path) -> Path:
    upstream_dir = output_dir / "upstream" / "results"
    if upstream_dir.exists():
        return upstream_dir
    raise SkillError(
        stage="parsing",
        error_code=ErrorCode.EXPECTED_OUTPUTS_NOT_FOUND,
        message="Pipeline output directory was not created.",
        fix="Re-run the wrapper after checking the Nextflow logs.",
        details={"expected_dir": str(upstream_dir)},
    )


def find_h5ad_candidates(upstream_dir: Path) -> list[str]:
    return sorted(str(path) for path in upstream_dir.rglob("*.h5ad"))


def find_rds_outputs(upstream_dir: Path) -> list[str]:
    return sorted(str(path) for path in upstream_dir.rglob("*.rds"))


def find_multiqc_report(upstream_dir: Path) -> str:
    return _first_match(upstream_dir, "multiqc/**/multiqc_report.html") or _first_match(
        upstream_dir, "**/multiqc_report.html"
    )


def find_pipeline_info_dir(upstream_dir: Path) -> str:
    pipeline_info_dir = upstream_dir / "pipeline_info"
    return str(pipeline_info_dir) if pipeline_info_dir.exists() else ""


def select_preferred_h5ad(h5ad_candidates: list[str]) -> str:
    for preferred_name in PREFERRED_H5AD_ORDER:
        preferred_matches = _find_candidates_named(h5ad_candidates, preferred_name)
        if len(preferred_matches) == 1:
            return preferred_matches[0]
        if len(preferred_matches) > 1:
            return ""
    # A single sample output is safe to hand off; multiple sample outputs are ambiguous.
    return h5ad_candidates[0] if len(h5ad_candidates) == 1 else ""


def _find_candidates_named(h5ad_candidates: list[str], filename: str) -> list[str]:
    return [candidate for candidate in h5ad_candidates if Path(candidate).name == filename]


def detect_sample_names(h5ad_candidates: list[str]) -> list[str]:
    return sorted(
        {
            _sample_name_from_h5ad(candidate)
            for candidate in h5ad_candidates
            if _sample_name_from_h5ad(candidate)
        }
    )


def _sample_name_from_h5ad(candidate: str) -> str:
    sample_name = (
        Path(candidate).name
        .replace("_cellbender_filter_matrix.h5ad", "")
        .replace("_filtered_matrix.h5ad", "")
        .replace("_raw_matrix.h5ad", "")
    )
    if sample_name == "combined" or sample_name.endswith(".h5ad"):
        return ""
    return sample_name


def detect_cellbender_outputs(h5ad_candidates: list[str]) -> bool:
    return any("cellbender" in path.lower() for path in h5ad_candidates)
