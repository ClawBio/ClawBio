#!/usr/bin/env python3
"""ClawBio data-extractor CLI.

Usage:
    python data_extractor.py --image figure.png --output results/
    python data_extractor.py --web --port 8765
    python data_extractor.py --demo
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Ensure skill root is on sys.path
_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))


def main():
    parser = argparse.ArgumentParser(
        description="ClawBio Data Extractor — digitize scientific figures",
    )
    parser.add_argument("--image", dest="image_path", help="Path to figure image (PNG/JPG)")
    parser.add_argument("--output", "-o", dest="output_dir", help="Output directory for CSV/JSON")
    parser.add_argument("--web", action="store_true", help="Launch interactive web UI")
    parser.add_argument("--port", type=int, default=8765, help="Web UI port (default: 8765)")
    parser.add_argument("--plot-type", dest="plot_type", help="Force plot type (skip auto-detection)")
    parser.add_argument("--demo", action="store_true", help="Run on bundled demo figure")
    parser.add_argument("--json", action="store_true", help="Output results as JSON to stdout")

    args = parser.parse_args()

    from api import run

    if args.demo:
        demo_fig = _SKILL_DIR / "data" / "demo_figure.png"
        if not demo_fig.exists():
            print("Demo figure not found. Place a figure at data/demo_figure.png")
            sys.exit(1)
        result = run(options={
            "image_path": str(demo_fig),
            "output_dir": args.output_dir or str(_SKILL_DIR / "output" / "demo"),
        })
    elif args.web:
        result = run(options={"web": True, "port": args.port})
    elif args.image_path:
        result = run(options={
            "image_path": args.image_path,
            "output_dir": args.output_dir or str(_SKILL_DIR / "output"),
            "plot_type": args.plot_type,
        })
    else:
        parser.print_help()
        sys.exit(1)

    if args.json:
        print(json.dumps(result, indent=2, default=str))
    elif result.get("success"):
        summary = result.get("summary", {})
        print(f"\n  Extracted {summary.get('n_panels', 0)} panel(s), "
              f"{summary.get('n_series', 0)} series, "
              f"{summary.get('total_points', 0)} data points")
        for f in result.get("output_files", []):
            print(f"  -> {f}")
        if result.get("mode") == "web":
            print(f"  Web UI running on port {result.get('port', 8765)}")
    else:
        print(f"\n  Error: {result.get('error', 'Unknown error')}")
        sys.exit(1)


if __name__ == "__main__":
    main()
