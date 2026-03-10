"""Get Me The Data — lightweight FastAPI server for the data-extractor skill."""

from __future__ import annotations

import base64 as b64mod
import hashlib
import io
import logging
import sys
from pathlib import Path

# Add skill root to sys.path so `from core.models import ...` works
_WEB_DIR = Path(__file__).resolve().parent
_SKILL_ROOT = _WEB_DIR.parent
if str(_SKILL_ROOT) not in sys.path:
    sys.path.insert(0, str(_SKILL_ROOT))

from fastapi import FastAPI, HTTPException, UploadFile
from fastapi.responses import FileResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from PIL import Image as PILImage
from pydantic import BaseModel

from core.models import (
    Confidence,
    ExtractedData,
    Figure,
    PlotType,
)
from core.digitizer import digitize_figure

logger = logging.getLogger(__name__)

# In-memory stores
images_store: dict[str, bytes] = {}  # image_id -> PNG bytes
images_meta: dict[str, dict] = {}  # image_id -> {width, height}
figures_store: dict[str, list[Figure]] = {}  # image_id -> figures
extracted_store: dict[str, list[ExtractedData]] = {}  # image_id -> results

app = FastAPI(title="Get Me The Data", version="0.2.0")

# Serve static files from the web/ directory (CSS, JS)
app.mount("/static", StaticFiles(directory=str(_WEB_DIR)), name="static")

# Also serve CSS/JS at root level (index.html uses relative paths)
@app.get("/styles.css")
async def serve_css():
    return FileResponse(str(_WEB_DIR / "styles.css"), media_type="text/css")

@app.get("/app.js")
async def serve_js():
    return FileResponse(str(_WEB_DIR / "app.js"), media_type="application/javascript")


@app.get("/")
async def index():
    return FileResponse(str(_WEB_DIR / "index.html"))


# --- Image Upload & Serving ---


@app.post("/api/upload-image")
async def upload_image(file: UploadFile):
    """Upload an image (PNG/JPG/etc) and return an image_id."""
    if not file.content_type or not file.content_type.startswith("image/"):
        raise HTTPException(status_code=400, detail="File must be an image")

    content = await file.read()
    image_id = hashlib.sha256(content).hexdigest()[:12]

    # Convert to PNG and store
    img = PILImage.open(io.BytesIO(content))
    if img.mode == "RGBA":
        bg = PILImage.new("RGB", img.size, (255, 255, 255))
        bg.paste(img, mask=img.split()[3])
        img = bg
    elif img.mode != "RGB":
        img = img.convert("RGB")

    buf = io.BytesIO()
    img.save(buf, format="PNG", optimize=True)
    png_bytes = buf.getvalue()

    images_store[image_id] = png_bytes
    images_meta[image_id] = {"width": img.width, "height": img.height}

    return {"image_id": image_id, "width": img.width, "height": img.height}


@app.get("/api/image/{image_id}")
async def get_image(image_id: str):
    """Serve an uploaded image as PNG."""
    if image_id not in images_store:
        raise HTTPException(status_code=404, detail="Image not found")
    return StreamingResponse(
        io.BytesIO(images_store[image_id]),
        media_type="image/png",
        headers={"Cache-Control": "public, max-age=3600"},
    )


# --- Plot Detection ---


@app.get("/api/detect-plots/{image_id}")
async def detect_plots(image_id: str):
    """Use Claude vision to detect plot regions on an uploaded image."""
    import anthropic

    if image_id not in images_store:
        raise HTTPException(status_code=404, detail="Image not found")

    b64 = b64mod.b64encode(images_store[image_id]).decode()

    client = anthropic.AsyncAnthropic()
    response = await client.messages.create(
        model="claude-sonnet-4-20250514",
        max_tokens=2048,
        messages=[{
            "role": "user",
            "content": [
                {
                    "type": "image",
                    "source": {"type": "base64", "media_type": "image/png", "data": b64},
                },
                {
                    "type": "text",
                    "text": (
                        "This is an image (likely a screenshot from a scientific paper). "
                        "Identify ALL distinct quantitative plots/charts "
                        "(bar charts, scatter plots, line graphs, box plots, heatmaps, "
                        "forest plots, Kaplan-Meier curves, histograms, violin plots, etc.).\n\n"
                        "Do NOT include: text-only areas, photographs, Western blots, gel images, "
                        "microscopy images, schematics, flowcharts, or diagrams without axes.\n"
                        "DO include: tables with numerical data.\n\n"
                        "For each plot found, return its bounding box as percentages (0-100) "
                        "of the image. Include the axis labels and tick marks in the box.\n\n"
                        "Return JSON:\n"
                        '{"plots": [\n'
                        '  {"label": "a", "type": "bar", "title": "short description", '
                        '"x_pct": 5, "y_pct": 10, "w_pct": 45, "h_pct": 40},\n'
                        "  ...\n"
                        "]}\n\n"
                        "If no quantitative plots are found, return: {\"plots\": []}\n"
                        "Plot type must be one of: scatter, bar, line, box, violin, histogram, "
                        "heatmap, forest, kaplan_meier, dot_strip, stacked_bar, funnel, roc, "
                        "volcano, waterfall, bland_altman, paired, bubble, area, dose_response, "
                        "manhattan, correlation_matrix, error_bar, table, other.\n"
                        "Return ONLY the JSON."
                    ),
                },
            ],
        }],
    )

    import json
    raw = response.content[0].text.strip()
    if raw.startswith("```"):
        raw = raw.split("\n", 1)[1]
    if raw.endswith("```"):
        raw = raw[:-3]
    raw = raw.strip()

    try:
        data = json.loads(raw)
        return data
    except json.JSONDecodeError:
        return {"plots": []}


# --- Extraction ---


class ImageRegionRequest(BaseModel):
    image_id: str
    crop_x_pct: float | None = None
    crop_y_pct: float | None = None
    crop_w_pct: float | None = None
    crop_h_pct: float | None = None


@app.post("/api/extract-image-region")
async def extract_image_region(req: ImageRegionRequest):
    """Extract data from a region of an uploaded image."""
    if req.image_id not in images_store:
        raise HTTPException(status_code=404, detail="Image not found")

    img = PILImage.open(io.BytesIO(images_store[req.image_id]))

    # Crop if region specified
    if req.crop_x_pct is not None:
        x = int(req.crop_x_pct / 100 * img.width)
        y = int(req.crop_y_pct / 100 * img.height)
        w = int(req.crop_w_pct / 100 * img.width)
        h = int(req.crop_h_pct / 100 * img.height)
        img = img.crop((x, y, x + w, y + h))

    # Resize if too large
    if img.width > 1500 or img.height > 1500:
        img.thumbnail((1500, 1500), PILImage.LANCZOS)

    buf = io.BytesIO()
    img.save(buf, format="PNG", optimize=True)
    b64 = b64mod.b64encode(buf.getvalue()).decode()

    figure_id = req.image_id + "_region"

    fig = Figure(
        figure_id=figure_id,
        paper_id=req.image_id,
        page_number=1,
        image_index=0,
        width=img.width,
        height=img.height,
        image_base64=b64,
        plot_type=PlotType.OTHER,
        plot_type_confidence=Confidence.MEDIUM,
    )

    # Store so image endpoint can serve panels
    if req.image_id not in figures_store:
        figures_store[req.image_id] = []
    figures_store[req.image_id].append(fig)

    try:
        results, panel_figures = await digitize_figure(fig)
        if req.image_id not in extracted_store:
            extracted_store[req.image_id] = []
        extracted_store[req.image_id].extend(results)
        if panel_figures:
            figures_store[req.image_id].extend(panel_figures)
        return [r.model_dump() for r in results]
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Extraction failed: {str(e)}")


# --- Figure image serving (for panel crops) ---


@app.get("/api/figure-image/{figure_id}")
async def get_figure_image(figure_id: str):
    """Serve a figure/panel image as PNG."""
    for figs in figures_store.values():
        for fig in figs:
            if fig.figure_id == figure_id:
                img_bytes = b64mod.b64decode(fig.image_base64)
                return StreamingResponse(
                    io.BytesIO(img_bytes),
                    media_type="image/png",
                    headers={"Cache-Control": "public, max-age=3600"},
                )
    raise HTTPException(status_code=404, detail="Figure not found")


# --- Edit extracted data ---


class EditCellRequest(BaseModel):
    image_id: str
    result_index: int
    series_index: int
    field: str  # "x_values" or "y_values"
    point_index: int
    value: float | str


@app.patch("/api/edit-cell")
async def edit_cell(req: EditCellRequest):
    """Edit a single extracted data point (user correction)."""
    if req.image_id not in extracted_store:
        raise HTTPException(status_code=404, detail="No results for this image")
    results = extracted_store[req.image_id]
    if req.result_index >= len(results):
        raise HTTPException(status_code=404, detail="Result index out of range")

    result = results[req.result_index]
    if req.series_index >= len(result.series):
        raise HTTPException(status_code=404, detail="Series index out of range")

    series = result.series[req.series_index]
    arr = getattr(series, req.field, None)
    if arr is None or req.field not in ("x_values", "y_values", "error_bars_lower", "error_bars_upper"):
        raise HTTPException(status_code=400, detail=f"Invalid field: {req.field}")
    if req.point_index >= len(arr):
        raise HTTPException(status_code=404, detail="Point index out of range")

    # Apply the edit
    if req.field == "y_values":
        arr[req.point_index] = float(req.value)
    elif req.field in ("error_bars_lower", "error_bars_upper"):
        arr[req.point_index] = float(req.value) if req.value is not None else None
    else:
        arr[req.point_index] = req.value

    return {"ok": True}


def launch(port: int = 8765):
    """Launch the server on the given port."""
    import uvicorn
    uvicorn.run(
        app,
        host="0.0.0.0",
        port=port,
    )


if __name__ == "__main__":
    launch()
