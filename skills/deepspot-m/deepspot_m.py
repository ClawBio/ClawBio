"""
deepspot_m.py, the DeepSpot-M virtual spatial transcriptomics skill for ClawBio.

Scores a 224x224 H&E tile with DeepSpot-M and writes per-gene log1p-CPM values.

Usage:
    python skills/deepspot-m/deepspot_m.py --input tile.png --output /tmp/deepspot_out
    python skills/deepspot-m/deepspot_m.py --input tile.png --genes BRAF,CD37,COL1A1 \
        --source evo2 --output /tmp/deepspot_out
    python skills/deepspot-m/deepspot_m.py --demo --output /tmp/deepspot_demo
"""
from __future__ import annotations

import argparse
import csv
import datetime
import json
import sys
from pathlib import Path

# Add project root so clawbio.common is importable when running as a script.
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

SKILL_DIR = Path(__file__).resolve().parent
EXAMPLES_DIR = SKILL_DIR / "examples"
DEMO_TILE = EXAMPLES_DIR / "demo_tile.png"
DEMO_EXPRESSION = EXAMPLES_DIR / "demo_expression.json"

DISCLAIMER = (
    "ClawBio is a research and educational tool. "
    "It is not a medical device and does not provide clinical diagnoses. "
    "Consult a healthcare professional before making any medical decisions."
)

MODEL_REPO = "ratschlab/DeepSpotM"
EXPRESSION_UNIT = "log1p-CPM"

# Frozen biological embedding spaces the gene router hypernetwork can draw from.
EMBEDDING_SOURCES = ("evo2", "orthrus", "prott5", "scgpt", "apertus")
DEFAULT_SOURCE = "scgpt"

# DeepSpot-M reads one tile of exactly this edge length, cut at roughly 20x
# magnification (about 0.5 microns per pixel).
TILE_SIZE = 224
TARGET_MPP = 0.5

DEFAULT_GENES = (
    "EPCAM", "KRT19", "COL1A1", "VIM", "ACTA2",
    "PTPRC", "CD68", "CD3D", "CD8A", "MKI67",
)

PIP_DEPS = ["deepspotm>=1.0", "Pillow>=9.0"]


# ---------------------------------------------------------------------------
# Input handling
# ---------------------------------------------------------------------------


def parse_genes(value: str | None) -> list[str]:
    """Split a comma separated gene list into unique uppercase HGNC symbols."""
    if value is None:
        return list(DEFAULT_GENES)

    genes: list[str] = []
    for token in value.replace(";", ",").split(","):
        symbol = token.strip().upper()
        if symbol and symbol not in genes:
            genes.append(symbol)
    if not genes:
        raise ValueError("--genes was given but contained no gene symbols.")
    return genes


def validate_tile_size(width: int, height: int) -> None:
    """Accept only a square tile of exactly TILE_SIZE pixels per side."""
    if width != TILE_SIZE or height != TILE_SIZE:
        raise ValueError(
            f"DeepSpot-M reads {TILE_SIZE}x{TILE_SIZE} tiles cut at roughly 20x "
            f"({TARGET_MPP} microns per pixel). Got {width}x{height}. "
            f"Re-tile the slide on a {TILE_SIZE}-pixel grid at native 20x resolution."
        )


def load_tile(path: Path | str):
    """Load an H&E tile as an RGB PIL image after checking its dimensions."""
    from PIL import Image

    tile = Image.open(str(path)).convert("RGB")
    validate_tile_size(tile.width, tile.height)
    return tile


# ---------------------------------------------------------------------------
# Prediction
# ---------------------------------------------------------------------------


def predict_expression(tile, genes: list[str], source: str) -> dict[str, float]:
    """Score one tile with DeepSpot-M and return per-gene log1p-CPM values.

    `deepspotm` and `torch` are imported here rather than at module scope so the
    skill loads, self-documents and runs its demo without the model stack.
    """
    try:
        import torch
        from deepspotm import DeepSpotM
    except ImportError as exc:
        raise RuntimeError(
            f"{exc}. Install the model stack with "
            "'pip install -r skills/deepspot-m/requirements.txt', then request "
            f"access to https://huggingface.co/{MODEL_REPO} and run "
            "'huggingface-cli login'. Use --demo to inspect the output format "
            "without any of that."
        ) from exc

    model, image_processor = DeepSpotM.from_pretrained(MODEL_REPO, source=source)
    model.eval()

    batch = image_processor(tile).unsqueeze(0)
    with torch.no_grad():
        values = model.predict_genes(batch, genes)

    scores = [float(v) for v in values.squeeze(0).tolist()]
    return dict(zip(genes, scores))


def load_demo_expression(genes: list[str]) -> dict[str, float]:
    """Read the bundled offline fixture for the demo gene panel."""
    fixture = json.loads(DEMO_EXPRESSION.read_text(encoding="utf-8"))
    table = fixture["expression"]

    missing = [gene for gene in genes if gene not in table]
    if missing:
        raise ValueError(
            f"Demo mode covers the bundled panel {sorted(table)}. "
            f"No fixture values for: {missing}. "
            "Run with --input and an installed deepspotm package to score other genes."
        )
    return {gene: float(table[gene]) for gene in genes}


def rank_genes(expression: dict[str, float]) -> list[dict[str, object]]:
    """Order genes by descending expression and attach a 1-based rank."""
    ordered = sorted(expression.items(), key=lambda item: (-item[1], item[0]))
    return [
        {"rank": index, "gene": gene, "expression": round(value, 4), "unit": EXPRESSION_UNIT}
        for index, (gene, value) in enumerate(ordered, start=1)
    ]


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------


def write_gene_table(rows: list[dict[str, object]], output_dir: Path) -> Path:
    """Write tables/gene_expression.csv."""
    tables_dir = output_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    path = tables_dir / "gene_expression.csv"
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["rank", "gene", "expression", "unit"])
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_result_json(
    output_dir: Path,
    rows: list[dict[str, object]],
    meta: dict[str, object],
) -> Path:
    """Write result.json."""
    payload = {
        "skill": "deepspot-m",
        "model": MODEL_REPO,
        "unit": EXPRESSION_UNIT,
        **meta,
        "genes": [row["gene"] for row in rows],
        "expression": {row["gene"]: row["expression"] for row in rows},
        "ranked": rows,
    }
    path = output_dir / "result.json"
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return path


def write_report(
    output_dir: Path,
    rows: list[dict[str, object]],
    meta: dict[str, object],
) -> Path:
    """Write report.md."""
    timestamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    demo_tag = " (demo)" if meta["demo"] else ""

    lines = [
        f"# DeepSpot-M Virtual Spatial Transcriptomics Report{demo_tag}",
        "",
        f"**Date**: {timestamp}",
        f"**Tile**: {meta['tile']}",
        f"**Tile size**: {TILE_SIZE}x{TILE_SIZE} px at roughly 20x ({TARGET_MPP} microns per pixel)",
        f"**Model**: {MODEL_REPO}",
        f"**Gene embedding source**: {meta['source']}",
        f"**Unit**: {EXPRESSION_UNIT}",
        f"**Genes scored**: {len(rows)}",
        "",
    ]

    if meta["demo"]:
        lines += [
            "> Demo mode. The values below come from the bundled offline fixture "
            "`examples/demo_expression.json`, not from a model run. They exist so "
            "the report format, the CSV schema and the reproducibility bundle can "
            "be inspected without the model weights.",
            "",
        ]

    heading = "Fixture Expression" if meta["demo"] else "Predicted Expression"
    lines += [
        f"## {heading}",
        "",
        f"| Rank | Gene | Expression ({EXPRESSION_UNIT}) |",
        "|------|------|------------------------|",
    ]
    for row in rows:
        lines.append(f"| {row['rank']} | {row['gene']} | {row['expression']:.2f} |")

    top = rows[0]
    subject = "this fixture" if meta["demo"] else "this tile"
    lines += [
        "",
        "## Summary",
        "",
        f"{top['gene']} holds the highest value in {subject} "
        f"({top['expression']:.2f} {EXPRESSION_UNIT}). "
        f"Values sit on the {EXPRESSION_UNIT} scale, the same scale as the "
        f"TCGA virtual spatial transcriptomics atlas of 28,664 slides "
        f"across 32 cancer types.",
        "",
        "## Output Files",
        "",
        "| File | Description |",
        "|------|-------------|",
        "| `result.json` | Machine-readable per-gene values and run parameters |",
        "| `tables/gene_expression.csv` | Ranked gene table, one row per gene |",
        "| `reproducibility/commands.sh` | Exact command that produced this run |",
        "| `reproducibility/environment.yml` | Conda and pip environment snapshot |",
        "| `reproducibility/checksums.sha256` | SHA-256 digests of the outputs |",
        "",
        "---",
        "",
        f"*{DISCLAIMER}*",
        "",
    ]

    path = output_dir / "report.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_reproducibility(
    output_dir: Path,
    gene_table: Path,
    meta: dict[str, object],
) -> None:
    """Write the reproducibility bundle."""
    script = "skills/deepspot-m/deepspot_m.py"
    gene_arg = ",".join(str(gene) for gene in meta["requested_genes"])
    if meta["demo"]:
        command = f"python {script} --demo --output {output_dir}"
    else:
        command = (
            f"python {script} --input {meta['tile']} --genes {gene_arg} "
            f"--source {meta['source']} --output {output_dir}"
        )
    write_commands_sh(output_dir, command)

    write_environment_yml(
        output_dir,
        env_name="clawbio-deepspot-m",
        pip_deps=PIP_DEPS,
        python_version="3.11",
        channels=["conda-forge", "nodefaults"],
    )

    write_checksums(
        [output_dir / "report.md", output_dir / "result.json", gene_table],
        output_dir,
        anchor=output_dir,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "DeepSpot-M: transcriptome-wide virtual spatial transcriptomics "
            "from a 224x224 H&E tile"
        )
    )
    parser.add_argument("--input", help=f"H&E tile, {TILE_SIZE}x{TILE_SIZE} px (PNG, JPEG or TIFF)")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument(
        "--genes",
        default=None,
        help="Comma separated HGNC symbols to score (default: a 10 gene marker panel)",
    )
    parser.add_argument(
        "--source",
        choices=EMBEDDING_SOURCES,
        default=DEFAULT_SOURCE,
        help="Frozen gene embedding space used by the gene router hypernetwork",
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Score the bundled demo tile from an offline fixture (no weights needed)",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.demo and not args.input:
        parser.error("--input is required unless --demo is used")

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        genes = parse_genes(args.genes)

        if args.demo:
            tile_label = str(DEMO_TILE.relative_to(_PROJECT_ROOT))
            # Check the bundled tile really is 224x224 when Pillow is available,
            # so demo runs exercise the same guard as real runs.
            try:
                load_tile(DEMO_TILE)
            except ImportError:
                pass
            expression = load_demo_expression(genes)
            source = DEFAULT_SOURCE
            if args.source != DEFAULT_SOURCE:
                print(
                    f"[deepspot-m] Note: the bundled fixture was recorded under "
                    f"'{DEFAULT_SOURCE}', so demo mode reports that source rather "
                    f"than '{args.source}'.",
                    file=sys.stderr,
                )
        else:
            tile_label = args.input
            tile = load_tile(args.input)
            source = args.source
            print(f"[deepspot-m] Scoring {len(genes)} genes with source '{source}'...")
            expression = predict_expression(tile, genes, source)
    except (ValueError, FileNotFoundError, OSError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    rows = rank_genes(expression)
    meta: dict[str, object] = {
        "demo": bool(args.demo),
        "tile": tile_label,
        "source": source,
        "tile_size_px": TILE_SIZE,
        "microns_per_pixel": TARGET_MPP,
        "requested_genes": genes,
    }

    gene_table = write_gene_table(rows, output_dir)
    write_result_json(output_dir, rows, meta)
    write_report(output_dir, rows, meta)
    write_reproducibility(output_dir, gene_table, meta)

    print(f"[deepspot-m] Done. {len(rows)} genes written.")
    print(f"[deepspot-m] Report: {output_dir / 'report.md'}")


if __name__ == "__main__":
    main()
