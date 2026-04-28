#!/usr/bin/env python3
"""ClawBio skill API — data-extractor.

Importable run() interface following ClawBio conventions.
"""

from __future__ import annotations

import asyncio
import sys
from pathlib import Path

# Ensure skill root is on sys.path for core imports
_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))


def run(genotypes=None, options=None) -> dict:
    """Run the data extractor.

    Parameters
    ----------
    genotypes : ignored (this skill does not consume genotype data)
    options : dict with keys:
        - image_path (str): Path to figure image
        - image_bytes (bytes): Raw image bytes (alternative to path)
        - output_dir (str): Where to write CSV/JSON results (default: ./output)
        - web (bool): If True, launch web UI (default False)
        - port (int): Web UI port (default 8765)
        - plot_type (str): Force plot type (optional, auto-detected)

    Returns
    -------
    dict with keys:
        - success (bool)
        - results (list[dict]): extracted data per panel
        - summary (dict): {n_panels, n_series, total_points, confidence}
        - output_files (list[str]): paths to generated files
        - output_dir (str)
    """
    options = options or {}

    # Web UI mode
    if options.get("web"):
        from web.server import launch
        port = options.get("port", 8765)
        launch(port=port)
        return {"success": True, "mode": "web", "port": port}

    # Extraction mode
    image_path = options.get("image_path")
    image_bytes = options.get("image_bytes")
    output_dir = options.get("output_dir", str(_SKILL_DIR / "output"))
    plot_type = options.get("plot_type")

    if not image_path and not image_bytes:
        return {
            "success": False,
            "error": "No image provided. Set image_path or image_bytes in options.",
            "results": [],
            "summary": {},
            "output_files": [],
            "output_dir": output_dir,
        }

    try:
        from data_extractor_core.digitizer import extract_from_image, export_csv, export_json
        results, _panel_figs = asyncio.run(
            extract_from_image(
                image_path=image_path,
                image_bytes=image_bytes,
                plot_type=plot_type,
            )
        )
    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "results": [],
            "summary": {},
            "output_files": [],
            "output_dir": output_dir,
        }

    # Export results
    output_files = []
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    csv_path = export_csv(results, str(out / "extracted_data.csv"))
    output_files.append(csv_path)

    json_path = export_json(results, str(out / "extracted_data.json"))
    output_files.append(json_path)

    # Summary
    total_series = sum(len(r.series) for r in results)
    total_points = sum(len(s.y_values) for r in results for s in r.series)
    confidences = [r.confidence.value for r in results]

    return {
        "success": True,
        "results": [r.model_dump() for r in results],
        "summary": {
            "n_panels": len(results),
            "n_series": total_series,
            "total_points": total_points,
            "confidences": confidences,
        },
        "output_files": output_files,
        "output_dir": output_dir,
    }
