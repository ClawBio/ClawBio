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
import os
import sys
from pathlib import Path

# Add project root so clawbio.common is importable when running as a script.
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (
    ReproCommand,
    ReproPath,
    write_checksums,
    write_environment_yml,
    write_portable_commands_sh,
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
# Pinned so a run today and a run next year read the same weights. The gated HF
# repo moves independently of the PyPI package, so a floating fetch would change
# results silently. Bump deliberately; see ## Maintenance in SKILL.md.
MODEL_REVISION = "86113ee431248c892d25cf55e1f8017cccec2926"

# The three files upstream's loader reads, in the order it reads them. Named here
# so the skill can resolve them itself and keep the download gate honest; see
# _resolve_checkpoint_dir.
MODEL_FILES = ("config.json", "model.safetensors", "tokens.csv")

EXPRESSION_UNIT = "log1p-CPM"

# Upstream's WEIGHTS_LICENSE.md puts the NonCommercial term on "the weights or
# their outputs", so the restriction rides along into the numbers this skill
# writes. Real runs stamp it on the report; demo runs do not, because fixture
# values never touched the weights.
WEIGHTS_LICENSE = "CC-BY-NC-SA-4.0"
OUTPUT_LICENSE_NOTE = (
    "These values are outputs of CC-BY-NC-SA-4.0 model weights. Upstream applies "
    "the NonCommercial term to the weights and their outputs, so this file is for "
    "non-commercial use, and attribution is required: cite the DeepSpot-M paper "
    "(doi:10.64898/2026.06.19.26356060)."
)

# Quoted verbatim from the "Limitations and biases" section of the model card at
# https://huggingface.co/ratschlab/DeepSpotM. They travel into every report
# because they qualify every number in it.
UPSTREAM_LIMITATIONS = (
    "Trained on a finite set of cancer indications.",
    "Performance on unseen tissue types, stains, scanners or resolutions may degrade.",
    "Predicts relative expression rather than absolute counts.",
    "Under-sequenced genes are predicted less reliably.",
    "Trained on oncology cohorts, so it is not representative of healthy tissue "
    "or non-oncology contexts.",
    "Not for clinical or diagnostic use.",
)

# Frozen biological embedding spaces the gene router hypernetwork can draw from.
EMBEDDING_SOURCES = ("evo2", "orthrus", "prott5", "scgpt", "apertus")
DEFAULT_SOURCE = "scgpt"

# DeepSpot-M reads one tile of exactly this edge length. The upstream README cuts
# tiles on a 224-pixel grid at native (~20x) resolution; no microns-per-pixel
# figure is published anywhere upstream, so this skill never assumes one.
TILE_SIZE = 224

# Mean-pixel threshold above which a tile counts as background. Upstream's own
# default, from examples/predict.py --white-mean.
WHITE_MEAN_DEFAULT = 220.0
# Mean HSV saturation on a 0-1 scale below which a tile is essentially greyscale.
# ClawBio-side heuristic, not an upstream check: H&E is a two-dye stain, so a
# colourless tile is not H&E at all.
MIN_SATURATION_DEFAULT = 0.05

# Sanity band for resolution metadata, in microns per pixel. This filters
# decoding mistakes, not biology: a 72-dpi print resolution decodes to 353
# microns per pixel, which is no microscope. It makes no claim about what
# magnification the model wants.
PLAUSIBLE_MPP_RANGE = (0.05, 5.0)

# What "~20x" comes out as on a slide scanner, in microns per pixel. This is a
# scanner convention, not a figure from the model card: upstream publishes no
# pixel size at all, only the "~20x" above. It exists to warn, never to resize,
# refuse or rescale, because 224x224 is a pixel count and not a field of view. A
# tile cut at 40x has identical pixel dimensions and covers a quarter of the
# tissue, and validate_tile_size cannot tell the two apart.
TWENTY_X_MPP_RANGE = (0.4, 0.6)

DEFAULT_GENES = (
    "EPCAM", "KRT19", "COL1A1", "VIM", "ACTA2",
    "PTPRC", "CD68", "CD3D", "CD8A", "MKI67",
)

# What environment.yml pins. Kept in step with requirements.txt, which lists
# huggingface_hub and torch because this file imports both directly rather than
# leaning on them arriving as deepspotm dependencies.
PIP_DEPS = ["deepspotm>=1.0,<2", "Pillow>=9.0", "huggingface_hub>=0.30", "torch>=2.0"]

# Everything main() turns into a one-line message rather than a traceback.
# KeyError is in here because upstream's genes_to_indices raises it for a symbol
# outside the released panel.
INPUT_ERRORS = (ValueError, KeyError, FileNotFoundError, OSError, RuntimeError)


# ---------------------------------------------------------------------------
# Input handling
# ---------------------------------------------------------------------------


def parse_genes(value: str | None) -> list[str]:
    """Split a comma separated gene list, keeping each symbol as it was typed.

    Deduplication is case-insensitive but the spelling survives, because HGNC
    keeps 'orf' lower case in roughly 200 symbols (C9orf72 among them) and
    upper-casing them makes a correct request unresolvable. resolve_genes turns
    whatever was typed into the panel's own spelling later on.
    """
    if value is None:
        return list(DEFAULT_GENES)

    genes: list[str] = []
    seen: set[str] = set()
    for token in value.replace(";", ",").split(","):
        symbol = token.strip()
        if not symbol or symbol.casefold() in seen:
            continue
        seen.add(symbol.casefold())
        genes.append(symbol)
    if not genes:
        raise ValueError("--genes was given but contained no gene symbols.")
    return genes


def resolve_genes(requested: list[str], panel) -> list[str]:
    """Match requested symbols against the model's own panel, case-insensitively.

    Returns the panel's spelling, so a request for 'c9orf72' or 'C9ORF72' both
    come back as 'C9orf72' and the report names the gene the way HGNC does.
    """
    lookup: dict[str, str] = {}
    for symbol in panel:
        lookup.setdefault(str(symbol).casefold(), str(symbol))

    resolved: list[str] = []
    missing: list[str] = []
    for symbol in requested:
        canonical = lookup.get(symbol.casefold())
        if canonical is None:
            missing.append(symbol)
        else:
            resolved.append(canonical)

    if missing:
        raise ValueError(
            f"{len(missing)} symbol(s) not in the released panel: {missing}. "
            "The checkpoint scores the roughly 19,000 genes in tokens.csv. "
            "Use a current HGNC symbol; a retired alias will not resolve, and "
            "genes outside the panel need regenerated source gene embeddings "
            "that are not part of the release."
        )
    return resolved


def panel_order(canonical: list[str], panel) -> list[str]:
    """Sort resolved symbols into the model panel's own index order.

    `predict_genes` hands back a bare vector with no symbols attached, so mapping
    it onto names rests on a convention upstream never states: values in the order
    the genes were requested, or in panel-index order. Its documented entry point
    is `genes_to_indices`, and resolving to indices then gathering is the natural
    vectorised implementation, which yields the second. `requirements.txt` admits
    any `deepspotm` 1.x, so a minor release could change which one holds, and the
    failure is silent: every gene keeps a plausible magnitude and a clean rank
    while carrying another gene's value.

    Requesting the genes already in panel-index order collapses the two
    conventions onto the same list, so the mapping is right under either. The
    caller's order is restored afterwards.
    """
    index: dict[str, int] = {}
    for position, symbol in enumerate(panel):
        index.setdefault(str(symbol), position)
    return sorted(canonical, key=lambda gene: index[gene])


def validate_tile_size(width: int, height: int) -> None:
    """Accept only a square tile of exactly TILE_SIZE pixels per side."""
    if width != TILE_SIZE or height != TILE_SIZE:
        raise ValueError(
            f"DeepSpot-M reads {TILE_SIZE}x{TILE_SIZE} tiles cut at native (~20x) "
            f"resolution. Got {width}x{height}. "
            f"Re-tile the slide on a {TILE_SIZE}-pixel grid at native 20x resolution."
        )


def load_tile(path: Path | str):
    """Load an H&E tile as a PIL image after checking its dimensions.

    The image is returned unconverted so resolution metadata survives; callers
    that need pixels use assess_tile or hand it to the image processor.
    """
    from PIL import Image

    tile = Image.open(str(path))
    validate_tile_size(tile.width, tile.height)
    return tile


# ---------------------------------------------------------------------------
# Physical scale: measured or declared, never assumed
# ---------------------------------------------------------------------------


def validate_mpp(value: float) -> float:
    """Reject a declared pixel size that cannot describe an image."""
    if value is None:
        return value
    if not (value > 0) or value != value or value in (float("inf"), float("-inf")):
        raise ValueError(f"--mpp must be a positive number of microns per pixel, got {value}.")
    return value


def field_of_view_warning(mpp: float | None) -> str | None:
    """Say so when a declared pixel size implies the wrong magnification.

    The 224x224 check is strict about pixel dimensions and blind to how much
    tissue they cover. The pixel size is the only datum in the run that can see
    it, and until now it was recorded and never read.
    """
    if mpp is None:
        return None

    low, high = TWENTY_X_MPP_RANGE
    if low <= mpp <= high:
        return None

    covered = TILE_SIZE * mpp
    expected = TILE_SIZE * (low + high) / 2
    magnification = "higher" if mpp < low else "lower"
    return (
        f"{mpp:.4g} microns per pixel is outside the {low}-{high} that a ~20x scan "
        f"typically produces, so this tile covers {covered:.0f} x {covered:.0f} "
        f"microns rather than roughly {expected:.0f} x {expected:.0f}. "
        f"That reads as a {magnification} magnification than the model was trained "
        "on. The pixel dimensions are identical either way, so nothing else in "
        "this run can catch it. Scored anyway; the band is a scanner convention "
        "and upstream publishes no pixel size."
    )


def read_pixel_size_um(image) -> float | None:
    """Read microns per pixel from an image's own resolution metadata.

    Returns None when the file carries no resolution, or carries one that is not
    a microscopy scale. Nothing here is inferred from magnification.
    """
    candidate: float | None = None

    # TIFF resolution tags come with an explicit unit, so prefer them.
    tags = getattr(image, "tag_v2", None)
    if tags:
        try:
            x_res = tags.get(282)  # XResolution
            unit = tags.get(296, 2)  # ResolutionUnit: 2 = inch, 3 = centimetre
            if x_res:
                x_res = float(x_res)
                if x_res > 0 and unit in (2, 3):
                    candidate = (25400.0 if unit == 2 else 10000.0) / x_res
        except (TypeError, ValueError, ZeroDivisionError):
            candidate = None

    if candidate is None:
        dpi = (getattr(image, "info", None) or {}).get("dpi")
        try:
            if dpi and float(dpi[0]) > 0:
                candidate = 25400.0 / float(dpi[0])
        except (TypeError, ValueError, IndexError, ZeroDivisionError):
            candidate = None

    if candidate is None:
        return None
    low, high = PLAUSIBLE_MPP_RANGE
    return candidate if low <= candidate <= high else None


def resolve_pixel_size(image, declared: float | None) -> tuple[float | None, str | None]:
    """Settle on a pixel size and say where it came from.

    An explicit --mpp always wins, but when the file disagrees the report says so
    rather than quietly picking one.
    """
    measured = read_pixel_size_um(image)

    if declared is not None:
        validate_mpp(declared)
        if measured is not None and abs(measured - declared) / declared > 0.05:
            return declared, (
                f"declared with --mpp; the image resolution tags report "
                f"{measured:.3f}, which disagrees"
            )
        return declared, "declared with --mpp"

    if measured is not None:
        return measured, "read from the image resolution tags"

    return None, None


# ---------------------------------------------------------------------------
# Does this tile look like H&E at all?
# ---------------------------------------------------------------------------


def assess_tile(
    image,
    white_mean: float = WHITE_MEAN_DEFAULT,
    min_saturation: float = MIN_SATURATION_DEFAULT,
) -> dict[str, object]:
    """Flag tiles that plainly are not stained tissue.

    Two cheap checks. The near-white test is upstream's own background filter
    from examples/predict.py, threshold and all. The saturation test is
    ClawBio-side. Neither decides whether the tissue is in distribution; they
    only catch tiles the model was never meant to read, which it would otherwise
    score and return with no hint that anything was wrong.
    """
    from PIL import ImageStat

    rgb = image.convert("RGB")
    mean_pixel = sum(ImageStat.Stat(rgb).mean) / 3.0
    mean_saturation = ImageStat.Stat(rgb.convert("HSV")).mean[1] / 255.0

    warnings: list[str] = []
    if mean_pixel > white_mean:
        warnings.append(
            f"The tile looks like background: mean pixel {mean_pixel:.1f} is above the "
            f"--white-mean threshold of {white_mean:.0f}. Blank glass still returns "
            "numbers, but not meaningful ones."
        )
    if mean_saturation < min_saturation:
        warnings.append(
            f"The tile does not look like H&E: mean saturation {mean_saturation:.3f} is "
            f"below {min_saturation:.2f}, so it is close to greyscale. Upstream notes "
            "performance may degrade on unseen stains."
        )

    return {
        "mean_pixel": round(mean_pixel, 2),
        "mean_saturation": round(mean_saturation, 4),
        "warnings": warnings,
    }


# ---------------------------------------------------------------------------
# Prediction
# ---------------------------------------------------------------------------


def _resolve_checkpoint_dir(allow_download: bool) -> str:
    """Put the pinned checkpoint on disk and return the directory holding it.

    Upstream's loader calls `hf_hub_download(repo_id, fn, revision=revision)` with
    no `local_files_only`, so handing it a repo id fetches whatever is missing.
    `HF_HUB_OFFLINE` cannot be used to stop that from here: huggingface_hub reads
    it once, at import, into `constants.HF_HUB_OFFLINE`, and `deepspotm` imports
    the library before this function ever runs. Resolving the three files here,
    with the flag passed explicitly, is what makes --allow-download mean anything.

    `DeepSpotM.from_pretrained` takes a local directory as well as a repo id, and
    builds the vision backbone from a config bundled in the package
    (`backbone_pretrained=False`), so once these three files are cached the load
    itself touches no network.
    """
    from huggingface_hub import hf_hub_download

    resolved = [
        hf_hub_download(
            MODEL_REPO,
            filename,
            revision=MODEL_REVISION,
            local_files_only=not allow_download,
        )
        for filename in MODEL_FILES
    ]

    # One repo at one pinned revision lands in one snapshot directory. Checked
    # rather than assumed, because the directory is what gets loaded.
    directories = {os.path.dirname(path) for path in resolved}
    if len(directories) != 1:
        raise RuntimeError(
            f"Expected {', '.join(MODEL_FILES)} to resolve into a single directory, "
            f"got {sorted(directories)}."
        )
    return directories.pop()


def predict_expression(
    tile,
    genes: list[str],
    source: str,
    allow_download: bool = False,
) -> dict[str, float]:
    """Score one tile with DeepSpot-M and return per-gene log1p-CPM values.

    `deepspotm` and `torch` are imported here rather than at module scope so the
    skill loads, self-documents and runs its demo without the model stack.
    Keys come back as the model's panel spells them, not as the caller typed them.
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

    try:
        checkpoint_dir = _resolve_checkpoint_dir(allow_download)
    except Exception as exc:
        if allow_download:
            raise
        raise RuntimeError(
            f"Could not read {MODEL_REPO} at revision {MODEL_REVISION[:12]} from the "
            f"local Hugging Face cache ({exc}). The weights are gated: request access "
            f"on https://huggingface.co/{MODEL_REPO}, run 'huggingface-cli login', then "
            "re-run with --allow-download to fetch them once. Use --demo to inspect "
            "the output format without any of that."
        ) from exc

    # The revision is pinned above, on the download; a directory has no revision
    # to pass on.
    model, image_processor = DeepSpotM.from_pretrained(checkpoint_dir, source=source)

    model.eval()

    canonical = resolve_genes(genes, model.gene_names)
    ordered = panel_order(canonical, model.gene_names)

    batch = image_processor(tile).unsqueeze(0)
    with torch.no_grad():
        values = model.predict_genes(batch, ordered)

    scores = [float(v) for v in values.squeeze(0).tolist()]
    if len(scores) != len(ordered):
        raise RuntimeError(
            f"deepspotm returned {len(scores)} value(s) for {len(ordered)} requested "
            f"gene(s) ({', '.join(ordered)}). Refusing to zip them: the shorter list "
            "would silently truncate the longer one and every gene past the cut would "
            "either vanish or take another gene's value."
        )

    scored = dict(zip(ordered, scores))
    # Back into the order the caller asked for; `ordered` was the panel's order.
    return {gene: scored[gene] for gene in canonical}


def load_demo_expression(genes: list[str]) -> dict[str, float]:
    """Read the bundled offline fixture for the demo gene panel."""
    fixture = json.loads(DEMO_EXPRESSION.read_text(encoding="utf-8"))
    table = fixture["expression"]

    lookup = {symbol.casefold(): symbol for symbol in table}
    missing = [gene for gene in genes if gene.casefold() not in lookup]
    if missing:
        raise ValueError(
            f"Demo mode covers the bundled panel {sorted(table)}. "
            f"No fixture values for: {missing}. "
            "Run with --input and an installed deepspotm package to score other genes."
        )
    return {lookup[g.casefold()]: float(table[lookup[g.casefold()]]) for g in genes}


def build_rows(expression: dict[str, float]) -> list[dict[str, object]]:
    """One row per gene, in the order it was requested.

    `rank` rides along as a secondary column for convenience, but the report
    never leads with it: DeepSpot-M predicts relative expression, so ordering
    different genes within one tile mostly recovers each gene's training-set
    mean rather than anything about this tile.
    """
    order = sorted(expression, key=lambda gene: (-expression[gene], gene))
    ranks = {gene: index for index, gene in enumerate(order, start=1)}
    return [
        {
            "gene": gene,
            "expression": round(value, 4),
            "unit": EXPRESSION_UNIT,
            "rank": ranks[gene],
        }
        for gene, value in expression.items()
    ]


# ---------------------------------------------------------------------------
# Run metadata
# ---------------------------------------------------------------------------


def tile_label(path: str) -> str:
    """The tile's file name, without the directory it came from.

    report.md and result.json get shared; the directory a slide sat in can name
    a patient, a cohort or a case number, and none of that belongs in a file
    someone forwards. reproducibility/commands.sh keeps the full path, because
    replaying the run is the one thing that needs it.
    """
    return Path(path).name


def build_meta(
    *,
    demo: bool,
    tile_label: str,
    source: str,
    genes: list[str],
    mpp: float | None,
    mpp_source: str | None,
    assessment: dict[str, object],
) -> dict[str, object]:
    """Collect everything the writers report about a run."""
    return {
        "demo": bool(demo),
        "tile": tile_label,
        "source": source,
        "tile_size_px": TILE_SIZE,
        "microns_per_pixel": mpp,
        "microns_per_pixel_source": mpp_source,
        "microns_per_pixel_warning": field_of_view_warning(mpp),
        "requested_genes": list(genes),
        "tile_mean_pixel": assessment["mean_pixel"],
        "tile_mean_saturation": assessment["mean_saturation"],
        "tile_warnings": list(assessment["warnings"]),
    }


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------


GENE_TABLE_FIELDS = (
    "gene",
    "expression_log1p_cpm",
    "unit",
    "rank",
    "provenance",
    "model",
    "model_revision",
)


def write_gene_table(
    rows: list[dict[str, object]],
    meta: dict[str, object],
    output_dir: Path,
) -> Path:
    """Write tables/gene_expression.csv, provenance on every row.

    This is the output built to be detached: SKILL.md names diff-visualizer as a
    consumer, and a CSV that says only `expression` cannot tell a heatmap whether
    its numbers came from the checkpoint or from the ten hand-written values in
    examples/demo_expression.json. report.md and result.json both carry that
    distinction; the file most likely to travel without them has to carry it too.
    Repeating the stamp on every row rather than in a header comment keeps it
    through a read_csv that would drop a comment line.
    """
    tables_dir = output_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    path = tables_dir / "gene_expression.csv"
    provenance = "demo_fixture" if meta["demo"] else "model_prediction"

    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(GENE_TABLE_FIELDS))
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "gene": row["gene"],
                    "expression_log1p_cpm": row["expression"],
                    "unit": row["unit"],
                    "rank": row["rank"],
                    "provenance": provenance,
                    "model": MODEL_REPO,
                    "model_revision": MODEL_REVISION,
                }
            )
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
        "model_revision": MODEL_REVISION,
        "unit": EXPRESSION_UNIT,
        **meta,
        "genes": [row["gene"] for row in rows],
        "expression": {row["gene"]: row["expression"] for row in rows},
        "ranked": sorted(rows, key=lambda row: row["rank"]),
        "interpretation": {
            "cross_gene_comparison": (
                "Values are comparable for one gene across tiles. Comparing "
                "different genes within one tile largely recovers each gene's "
                "average abundance in the training data."
            ),
            # Stated rather than omitted. A consumer that finds no uncertainty
            # key reads the point estimates as confident ones, and a tile from
            # an organ the model never trained on passes both tile checks and
            # comes back looking like any other.
            "per_gene_uncertainty": None,
            "uncertainty_note": (
                "The checkpoint returns one point estimate per gene and no "
                "confidence, interval, variance or out-of-distribution score. "
                "null here means none was computed, not that the estimate is "
                "certain."
            ),
            "upstream_limitations": list(UPSTREAM_LIMITATIONS),
        },
    }
    if not meta["demo"]:
        payload["weights_license"] = WEIGHTS_LICENSE
        payload["output_license_note"] = OUTPUT_LICENSE_NOTE
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

    if meta["microns_per_pixel"] is None:
        mpp_line = (
            "**Microns per pixel**: not declared "
            "(pass --mpp, or use a tile whose resolution tags carry it)"
        )
    else:
        mpp_line = (
            f"**Microns per pixel**: {meta['microns_per_pixel']:.4g} "
            f"({meta['microns_per_pixel_source']})"
        )

    lines = [
        f"# DeepSpot-M Virtual Spatial Transcriptomics Report{demo_tag}",
        "",
        f"**Date**: {timestamp}",
        f"**Tile**: {meta['tile']}",
        f"**Tile size**: {TILE_SIZE}x{TILE_SIZE} px, cut at native (~20x) resolution",
        mpp_line,
        f"**Model**: {MODEL_REPO} @ {MODEL_REVISION[:12]}",
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

    tile_checks = list(meta["tile_warnings"])
    # The field-of-view warning is kept out of tile_warnings on purpose:
    # --skip-background refuses to score on those, and a pixel size that reads
    # as 40x is a caveat on the result, not a reason to withhold it.
    if meta.get("microns_per_pixel_warning"):
        tile_checks.append(meta["microns_per_pixel_warning"])

    if tile_checks:
        lines += ["## Tile Checks", ""]
        lines += [f"- {warning}" for warning in tile_checks]
        lines += [""]

    heading = "Fixture Expression" if meta["demo"] else "Predicted Expression"
    lines += [
        f"## {heading}",
        "",
        "Genes appear in the order they were requested.",
        "",
        f"| Gene | Expression ({EXPRESSION_UNIT}) |",
        "|------|------------------------|",
    ]
    for row in rows:
        lines.append(f"| {row['gene']} | {row['expression']:.2f} |")

    lines += [
        "",
        "## How to Read These Values",
        "",
        "DeepSpot-M predicts relative expression, so a value means something next "
        "to the same gene in another tile, not next to a different gene in this "
        "one. Ordering the genes in this table by value would largely recover each "
        "gene's average abundance in the training data rather than anything "
        f"specific to this {'fixture' if meta['demo'] else 'tile'}. "
        "`tables/gene_expression.csv` carries a `rank` column for convenience; it "
        "inherits that caveat.",
        "",
        "Upstream states the following limitations, quoted from the "
        "\"Limitations and biases\" section of the model card:",
        "",
    ]
    lines += [f"- {limitation}" for limitation in UPSTREAM_LIMITATIONS]

    lines += [
        "",
        "## Output Files",
        "",
        "| File | Description |",
        "|------|-------------|",
        "| `result.json` | Machine-readable per-gene values and run parameters |",
        "| `tables/gene_expression.csv` | Gene table, one row per gene |",
        "| `reproducibility/commands.sh` | Exact command that produced this run |",
        "| `reproducibility/environment.yml` | Conda and pip environment snapshot |",
        "| `reproducibility/checksums.sha256` | SHA-256 digests of the outputs |",
        "",
        "---",
        "",
    ]

    if not meta["demo"]:
        lines += [f"*Licence: {OUTPUT_LICENSE_NOTE}*", ""]

    lines += [
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
    tile_path: str | None,
    args: argparse.Namespace,
) -> None:
    """Write the reproducibility bundle."""
    script = Path("skills/deepspot-m/deepspot_m.py")
    cmd_args: list[str | ReproPath] = []

    if meta["demo"]:
        cmd_args.append("--demo")
    else:
        cmd_args += ["--input", ReproPath(Path(tile_path).resolve(), anchor="auto")]
        cmd_args += ["--genes", ",".join(str(g) for g in meta["requested_genes"])]
        cmd_args += ["--source", str(meta["source"])]
        if meta["microns_per_pixel"] is not None and args.mpp is not None:
            cmd_args += ["--mpp", str(args.mpp)]
        # Both thresholds feed the tile assessment, so a bundle that drops a
        # non-default one can refuse to score a tile the original run scored.
        if args.white_mean != WHITE_MEAN_DEFAULT:
            cmd_args += ["--white-mean", str(args.white_mean)]
        if args.min_saturation != MIN_SATURATION_DEFAULT:
            cmd_args += ["--min-saturation", str(args.min_saturation)]
        if args.skip_background:
            cmd_args.append("--skip-background")
        if args.allow_download:
            cmd_args.append("--allow-download")

    cmd_args += ["--output", ReproPath(output_dir, anchor="output_dir")]

    write_portable_commands_sh(
        output_dir,
        ReproCommand(
            script_path=script,
            args=cmd_args,
            comment=(
                "Replays this deepspot-m run. The --input path is the only place "
                "the bundle records where the tile came from."
            ),
        ),
        repo_root=_PROJECT_ROOT,
    )

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
        "--mpp",
        type=float,
        default=None,
        help="Microns per pixel of the tile. Recorded as declared; nothing is "
             "assumed when it is absent and the file carries no resolution tags.",
    )
    parser.add_argument(
        "--skip-background",
        action="store_true",
        help="Refuse to score a tile that looks like background instead of warning",
    )
    parser.add_argument(
        "--white-mean",
        type=float,
        default=WHITE_MEAN_DEFAULT,
        help="Mean pixel value above which a tile counts as background "
             f"(default: {WHITE_MEAN_DEFAULT:.0f}, upstream's own threshold)",
    )
    parser.add_argument(
        "--min-saturation",
        type=float,
        default=MIN_SATURATION_DEFAULT,
        help="Mean HSV saturation (0-1) below which a tile is flagged as not H&E "
             f"(default: {MIN_SATURATION_DEFAULT})",
    )
    parser.add_argument(
        "--allow-download",
        action="store_true",
        help="Permit the one-time gated weight download from Hugging Face. "
             "Without it the model is loaded from the local cache only.",
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
        if args.mpp is not None:
            validate_mpp(args.mpp)

        tile_source = str(DEMO_TILE) if args.demo else args.input
        tile = load_tile(tile_source)
        assessment = assess_tile(tile, args.white_mean, args.min_saturation)
        mpp, mpp_source = resolve_pixel_size(tile, args.mpp)

        for warning in assessment["warnings"]:
            print(f"[deepspot-m] Warning: {warning}", file=sys.stderr)

        # Not folded into assessment["warnings"]: those decide --skip-background,
        # and a pixel size reading as 40x is a caveat, not grounds for refusing.
        scale_warning = field_of_view_warning(mpp)
        if scale_warning:
            print(f"[deepspot-m] Warning: {scale_warning}", file=sys.stderr)
        if args.skip_background and assessment["warnings"]:
            raise ValueError(
                "Refusing to score this tile: "
                + " ".join(assessment["warnings"])
                + " Drop --skip-background to score it anyway, or raise --white-mean."
            )

        if args.demo:
            label = tile_label(str(DEMO_TILE))
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
            label = tile_label(args.input)
            source = args.source
            print(f"[deepspot-m] Scoring {len(genes)} genes with source '{source}'...")
            expression = predict_expression(
                tile.convert("RGB"), genes, source, allow_download=args.allow_download
            )
    except INPUT_ERRORS as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    rows = build_rows(expression)
    meta = build_meta(
        demo=args.demo,
        tile_label=label,
        source=source,
        genes=list(expression),
        mpp=mpp,
        mpp_source=mpp_source,
        assessment=assessment,
    )

    gene_table = write_gene_table(rows, meta, output_dir)
    write_result_json(output_dir, rows, meta)
    write_report(output_dir, rows, meta)
    write_reproducibility(output_dir, gene_table, meta, args.input, args)

    print(f"[deepspot-m] Done. {len(rows)} genes written.")
    print(f"[deepspot-m] Report: {output_dir / 'report.md'}")


if __name__ == "__main__":
    main()
