#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    ReproCommand,
    ReproPath,
    write_checksums,
    write_environment_yml,
    write_portable_commands_sh,
)

SKILL_DIR = Path(__file__).resolve().parent
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)
REQUIRED_COLUMNS = {"spot_id", "x", "y", "total_counts", "EPCAM", "PTPRC", "COL1A1", "MKI67"}
MARKERS = ["EPCAM", "PTPRC", "COL1A1", "MKI67"]


def load_spots(path: Path) -> list[dict[str, str | float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Input is missing required columns: {sorted(missing)}")
        spots = []
        for raw in reader:
            row: dict[str, str | float] = {"spot_id": raw["spot_id"]}
            for column in sorted(REQUIRED_COLUMNS - {"spot_id"}):
                try:
                    row[column] = float(raw[column])
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"Column {column} must be numeric") from exc
            spots.append(row)
    if not spots:
        raise ValueError("Input contains no spots")
    return spots


def _region_for(row: dict[str, str | float]) -> str:
    values = {marker: float(row[marker]) for marker in MARKERS}
    marker = max(values, key=values.get)
    return {"PTPRC": "immune_edge", "EPCAM": "tumor_core", "COL1A1": "stromal_zone", "MKI67": "proliferative_core"}[marker]


def map_spots(spots: list[dict[str, str | float]]) -> dict:
    mapped = []
    region_stats: dict[str, dict[str, float | int | str]] = {}
    for spot in spots:
        region = _region_for(spot)
        hotspot = region in {"tumor_core", "proliferative_core"}
        mapped_spot = {**spot, "region": region, "hotspot": hotspot}
        mapped.append(mapped_spot)
        stats = region_stats.setdefault(region, {"region": region, "spot_count": 0, "mean_total_counts": 0.0})
        stats["spot_count"] = int(stats["spot_count"]) + 1
        stats["mean_total_counts"] = float(stats["mean_total_counts"]) + float(spot["total_counts"])
    for stats in region_stats.values():
        stats["mean_total_counts"] = round(float(stats["mean_total_counts"]) / int(stats["spot_count"]), 2)
    priority = {"tumor_core": 0, "immune_edge": 1, "stromal_zone": 2, "proliferative_core": 3}
    regions = sorted(region_stats.values(), key=lambda item: (-int(item["spot_count"]), priority.get(str(item["region"]), 99)))
    return {
        "skill": "marker-dominance-mapper",
        "summary": {"spot_count": len(mapped), "region_count": len(regions), "hotspot_count": sum(1 for spot in mapped if spot["hotspot"])},
        "spots": mapped,
        "regions": regions,
        "disclaimer": DISCLAIMER,
    }


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["spot_id", "x", "y", "total_counts", "EPCAM", "PTPRC", "COL1A1", "MKI67", "region", "hotspot"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def _write_region_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["region", "spot_count", "mean_total_counts"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def _write_svg(path: Path, spots: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    colors = {
        "immune_edge": "#2b6cb0",
        "tumor_core": "#c53030",
        "stromal_zone": "#2f855a",
        "proliferative_core": "#b7791f",
    }
    max_x = max(float(spot["x"]) for spot in spots) or 1.0
    max_y = max(float(spot["y"]) for spot in spots) or 1.0
    elements = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="420" height="320" viewBox="0 0 420 320">',
        '<rect width="420" height="320" fill="#ffffff"/>',
        '<text x="16" y="24" font-family="Arial" font-size="16" font-weight="bold">Marker dominance map</text>',
    ]
    for spot in spots:
        x = 40 + (float(spot["x"]) / max_x) * 300
        y = 60 + (float(spot["y"]) / max_y) * 210
        region = str(spot["region"])
        radius = 13 if spot["hotspot"] else 10
        elements.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{radius}" fill="{colors.get(region, "#718096")}" opacity="0.85"/>')
        elements.append(f'<text x="{x + 12:.1f}" y="{y + 4:.1f}" font-family="Arial" font-size="10">{spot["spot_id"]}</text>')
    legend_y = 70
    for region, color in colors.items():
        elements.append(f'<circle cx="360" cy="{legend_y}" r="6" fill="{color}"/>')
        elements.append(f'<text x="372" y="{legend_y + 4}" font-family="Arial" font-size="10">{region}</text>')
        legend_y += 18
    elements.append("</svg>")
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def repro_command(input_path: Path, output_dir: Path, demo: bool) -> ReproCommand:
    args: list[str | ReproPath] = ["--demo"] if demo else [
        "--input",
        ReproPath(input_path, "repo_root" if input_path.is_relative_to(_PROJECT_ROOT) else "auto"),
    ]
    args += ["--output", ReproPath(output_dir, "output_dir")]
    return ReproCommand(
        script_path=Path("skills/marker-dominance-mapper/marker_dominance_mapper.py"),
        args=args,
        comment="Reproduce this marker-dominance-mapper run",
    )


def write_repro_bundle(input_path: Path, output_dir: Path, demo: bool) -> None:
    write_environment_yml(
        output_dir,
        env_name="clawbio-marker-dominance-mapper",
        pip_deps=[],
        python_version="3.10",
    )
    write_portable_commands_sh(
        output_dir,
        repro_command(input_path, output_dir, demo),
        repo_root=_PROJECT_ROOT,
    )
    write_checksums(
        [
            output_dir / "report.md",
            output_dir / "result.json",
            output_dir / "tables" / "mapped_spots.csv",
            output_dir / "tables" / "region_summary.csv",
            output_dir / "figures" / "marker_map.svg",
        ],
        output_dir,
        anchor=output_dir,
    )


def write_outputs(result: dict, input_path: Path, output_dir: Path, demo: bool) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        print(f"WARNING: output directory already exists and files may be overwritten: {output_dir}", file=sys.stderr)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)
    (output_dir / "reproducibility").mkdir(exist_ok=True)
    _write_csv(output_dir / "tables" / "mapped_spots.csv", result["spots"])
    _write_region_csv(output_dir / "tables" / "region_summary.csv", result["regions"])
    _write_svg(output_dir / "figures" / "marker_map.svg", result["spots"])
    (output_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    rows = [
        "# Marker Dominance Mapper Report",
        "",
        f"**Input**: `{input_path}`",
        f"**Mode**: {'Synthetic demo data' if demo else 'User-provided local data'}",
        f"**Spots mapped**: {result['summary']['spot_count']}",
        f"**Hotspots**: {result['summary']['hotspot_count']}",
        "",
        "| Spot | X | Y | Region | Hotspot | EPCAM | PTPRC | COL1A1 | MKI67 |",
        "|---|---:|---:|---|---|---:|---:|---:|---:|",
    ]
    for spot in result["spots"]:
        rows.append(f"| {spot['spot_id']} | {spot['x']} | {spot['y']} | {spot['region']} | {spot['hotspot']} | {spot['EPCAM']} | {spot['PTPRC']} | {spot['COL1A1']} | {spot['MKI67']} |")
    rows.extend([
        "",
        "## Interpretation",
        "",
        "Regions are assigned by dominant marker expression. Coordinates are used for SVG placement only, not for region assignment.",
        "Hotspots are tumor-core or MKI67-dominant proliferative-core spots.",
        "",
        DISCLAIMER,
        "",
    ])
    (output_dir / "report.md").write_text("\n".join(rows), encoding="utf-8")
    write_repro_bundle(input_path, output_dir, demo)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Marker Dominance Mapper")
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path, default=Path("marker_dominance_mapper_out"))
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args(argv)
    input_path = SKILL_DIR / "demo_marker_counts.csv" if args.demo else args.input
    if input_path is None:
        parser.error("--input is required unless --demo is used")
    try:
        result = map_spots(load_spots(input_path))
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    write_outputs(result, input_path, args.output, args.demo)
    print(f"Marker Dominance Mapper wrote {args.output / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
