"""cell-detection skill — cell segmentation using cpsam (Cellpose 4.0) and future backends."""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

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

DISCLAIMER = (
    "*ClawBio is a research and educational tool. "
    "It is not a medical device and does not provide clinical diagnoses. "
    "Consult a healthcare professional before making any medical decisions.*"
)


def make_demo_image(seed: int = 42) -> np.ndarray:
    """Generate a synthetic fluorescence nuclei image (no network required).

    Produces a 512×512 greyscale uint8 image with ~80 randomly placed interior
    nuclei plus 8 cells that are deliberately cut off at each image border —
    ideal for testing --exclude_on_edges.  Similar in character to a
    DAPI-stained mitosis field.
    """
    from scipy.ndimage import gaussian_filter

    rng = np.random.default_rng(seed)
    h, w = 512, 512
    img = np.zeros((h, w), dtype=np.float32)

    def _draw_nucleus(cy: int, cx: int, radius: float, intensity: float) -> None:
        a = radius * rng.uniform(0.7, 1.3)
        b = radius * rng.uniform(0.7, 1.3)
        angle = rng.uniform(0, np.pi)
        ys, xs = np.ogrid[:h, :w]
        dy = ys - cy
        dx = xs - cx
        cos_a, sin_a = np.cos(angle), np.sin(angle)
        ell = ((dx * cos_a + dy * sin_a) / a) ** 2 + ((-dx * sin_a + dy * cos_a) / b) ** 2
        mask = ell <= 1.0
        img[mask] = np.maximum(img[mask], intensity * rng.uniform(0.85, 1.0))

    # Interior cells
    n_cells = 80
    for _ in range(n_cells):
        cy = rng.integers(20, h - 20)
        cx = rng.integers(20, w - 20)
        _draw_nucleus(cy, cx, rng.uniform(8, 18), rng.uniform(180, 255))

    # Edge cells — centres placed just outside the border so the nucleus is
    # partially clipped, giving --exclude_on_edges something to act on.
    edge_positions = [
        (0, w // 5),           # top
        (0, 3 * w // 5),       # top
        (h - 1, w // 4),       # bottom
        (h - 1, 3 * w // 4),   # bottom
        (h // 5, 0),           # left
        (3 * h // 5, 0),       # left
        (h // 4, w - 1),       # right
        (3 * h // 4, w - 1),   # right
    ]
    for cy, cx in edge_positions:
        _draw_nucleus(cy, cx, rng.uniform(10, 16), rng.uniform(180, 240))

    # Soft glow around nuclei
    img = gaussian_filter(img, sigma=1.5)
    # Background noise
    img += rng.normal(8, 4, img.shape).astype(np.float32)
    img = np.clip(img, 0, 255).astype(np.uint8)
    return img


def _normalize_axes(axes: str) -> str:
    """Normalize vendor axis labels while preserving positional metadata.

    Notes:
    - Keep alphanumeric axis markers so trailing placeholders like "0" in
      CZI axes strings (e.g. "VBTCZYX0") remain aligned with array ndim.
    - Map TIFF sample axis "S" to channel axis "C" for consistent handling.
    """
    normalized: list[str] = []
    for ch in axes.upper():
        if not ch.isalnum():
            continue
        normalized.append("C" if ch == "S" else ch)
    return "".join(normalized)


def _drop_singleton_non_czyx(arr: np.ndarray, axes: str) -> tuple[np.ndarray, str]:
    for idx in reversed(range(len(axes))):
        if axes[idx] in {"C", "Z", "Y", "X"}:
            continue
        if arr.shape[idx] != 1:
            raise ValueError(f"Unsupported non-singleton axis '{axes[idx]}' with shape {arr.shape}")
        arr = np.take(arr, 0, axis=idx)
        axes = axes[:idx] + axes[idx + 1:]
    return arr, axes


def _transpose_to_axes(arr: np.ndarray, axes: str, target_axes: str) -> np.ndarray:
    order = [axes.index(axis) for axis in target_axes]
    return np.transpose(arr, order)


def _finalize_loaded_array(arr: np.ndarray, axes: str, path: str, z_projection: str) -> tuple[np.ndarray, int] | None:
    axes = _normalize_axes(axes)
    if not axes or arr.ndim != len(axes):
        return None

    arr, axes = _drop_singleton_non_czyx(arr, axes)

    if set(axes) == {"C", "Z", "Y", "X"} and len(axes) == 4:
        arr = _transpose_to_axes(arr, axes, "CZYX")
        if z_projection == "max":
            cyx = arr.max(axis=1)
            return cyx.transpose(1, 2, 0), cyx.shape[0]
        if z_projection == "none":
            zcyx = arr.transpose(1, 0, 2, 3)
            return zcyx, zcyx.shape[1]
        raise ValueError(f"Unsupported z_projection '{z_projection}'")

    if set(axes) == {"C", "Y", "X"} and len(axes) == 3:
        cyx = _transpose_to_axes(arr, axes, "CYX")
        return cyx.transpose(1, 2, 0), cyx.shape[0]

    if set(axes) == {"Z", "Y", "X"} and len(axes) == 3:
        zyx = _transpose_to_axes(arr, axes, "ZYX")
        return zyx, 1

    if set(axes) == {"Y", "X"} and len(axes) == 2:
        yx = _transpose_to_axes(arr, axes, "YX")
        return yx, 1

    return None


def load_image(path: str, z_projection: str = "max") -> tuple[np.ndarray, int]:
    """Load image and return (array, n_channels).

    Supported inputs: TIFF/TIF, CZI, ND2, PNG, JPG/JPEG.

    Output layouts:
      - 2D greyscale: H×W, n_channels=1
      - 2D multi-channel: H×W×C, n_channels=C
      - 3D volume from explicit Z input: Z×H×W, n_channels=1
      - 4D volume when z_projection='none': Z×C×H×W, n_channels=C

    z_projection:
      - "max": default, max-project 4D C/Z stacks over Z while preserving channels.
      - "none": keep Z; 4D stacks return Z×C×H×W (including C=1).

    Axis metadata is used when available (including mapping TIFF sample axis
    "S" to channel axis "C"), with shape-based fallback for common
    microscopy stack layouts.
    """
    p = Path(path)
    extension = p.suffix.lower()
    arr: np.ndarray
    axes = ""
    if extension in (".tif", ".tiff"):
        import tifffile
        with tifffile.TiffFile(str(p)) as tf:
            arr = tf.asarray()
            # Try to read axes from OME/series metadata
            if tf.series:
                axes = tf.series[0].axes.upper()
    elif extension == ".czi":
        import czifile
        with czifile.CziFile(str(p)) as czi:
            arr = np.asarray(czi.asarray())
            axes = getattr(czi, "axes", "") or ""
    elif extension == ".nd2":
        import nd2
        arr = np.asarray(nd2.imread(str(p)))
        with nd2.ND2File(str(p)) as nd2f:
            sizes = getattr(nd2f, "sizes", {}) or {}
            if hasattr(sizes, "keys"):
                axes = "".join(str(k).upper() for k in sizes.keys())
    else:
        from PIL import Image
        arr = np.array(Image.open(str(p)))

    finalized = _finalize_loaded_array(arr, axes, path, z_projection)
    if finalized is not None:
        return finalized

    if arr.ndim == 4:
        # Microscopy stacks can come as C×Z×H×W (CZI) or Z×C×H×W (ND2).
        # For 2D segmentation, reduce Z by max projection while preserving channels.
        if arr.shape[0] <= 20 and arr.shape[1] <= 20:
            if arr.shape[0] <= arr.shape[1]:
                raise ValueError(
                    f"Cannot infer channel vs Z axis for 4D image shape {arr.shape} at {path}: "
                    "both leading dimensions are <= 20 and the first is not larger than the second. "
                    "Provide CZI/ND2/TIFF axis metadata, "
                    "or use a stack where Z and C differ enough for the size heuristic "
                    "(Z > C for Z×C×Y×X, or C > Z for C×Z×Y×X)."
                )
            channel_axis = 1
        elif arr.shape[0] <= 20:
            channel_axis = 0
        elif arr.shape[1] <= 20:
            channel_axis = 1
        else:
            raise ValueError(f"Unsupported 4D image shape {arr.shape} for {path}")
        z_axis = 1 if channel_axis == 0 else 0
        if z_projection == "max":
            arr = arr.max(axis=z_axis)  # C×H×W
        elif z_projection == "none":
            if channel_axis == 0:
                arr = np.moveaxis(arr, 0, 1)  # CZYX -> ZCYX
            # If channel_axis == 1, already ZCYX
        else:
            raise ValueError(f"Unsupported z_projection '{z_projection}'")
    if arr.ndim == 2:
        return arr, 1
    if arr.ndim == 4:
        if z_projection == "none":
            # 4D no-projection path uses volumetric layout Z×C×Y×X.
            return arr, arr.shape[1]
        raise ValueError(f"Unexpected image ndim=4 for {path}")
    if arr.ndim == 3:
        c, h, w = arr.shape
        # CYX: either explicit from metadata or heuristic (C small, H/W much larger)
        is_cyx = ("C" in axes and axes.index("C") == 0) or (
            c <= 20 and h >= 4 * c and w >= 4 * c
        )
        if is_cyx:
            arr = arr.transpose(1, 2, 0)  # → H×W×C
        return arr, arr.shape[2]
    raise ValueError(f"Unexpected image ndim={arr.ndim} for {path}")


def prepare_image(img: np.ndarray, n_channels: int) -> np.ndarray:
    """Prepare image for cpsam segmentation.

    cpsam is channel-order invariant and uses up to 3 channels.
    - 1 channel (greyscale): pass as-is (H×W)
    - 2–3 channels: pass as-is (H×W×C) — cpsam handles any ordering
    - >3 channels: truncate to first 3 with a warning
    """
    if n_channels <= 3:
        return img
    print(
        f"[cell-detection] Warning: image has {n_channels} channels; "
        "cpsam uses the first 3. Truncating."
    )
    return img[:, :, :3]


def compute_metrics(masks: np.ndarray) -> list[dict]:
    """Compute per-cell morphology metrics from a label image.

    Args:
        masks: H×W uint16 array. 0=background, 1..N=cell labels.

    Returns:
        List of dicts with keys: id, area, diameter, centroid_x, centroid_y, eccentricity.
    """
    from skimage.measure import regionprops

    rows = []
    is_3d = masks.ndim == 3
    for prop in regionprops(masks.astype(int)):
        # skimage does not implement eccentricity for 3D regions.
        eccentricity = float("nan") if is_3d else prop.eccentricity
        centroid_x = prop.centroid[-1]
        centroid_y = prop.centroid[-2]
        rows.append(
            {
                "id": prop.label,
                "area": prop.area,
                "diameter": prop.equivalent_diameter_area,
                "centroid_x": centroid_x,
                "centroid_y": centroid_y,
                "eccentricity": eccentricity,
            }
        )
    return rows


def write_report(metrics: list[dict], meta: dict, output_dir: Path | str, outlines_filename: str = "image_cp_outlines.png", masks_filename: str = "image_cp_masks.tif", seg_filename: str = "image_seg.npy", csv_filename: str = "image_measurements.csv", histogram_filename: str = "image_histogram.png") -> None:
    """Write report.md to output_dir."""
    import statistics

    output_dir = Path(output_dir)
    n = len(metrics)
    areas = [m["area"] for m in metrics]
    diameters = [m["diameter"] for m in metrics]

    def _stats(vals: list) -> str:
        if not vals:
            return "N/A"
        if len(vals) == 1:
            return f"median={vals[0]:.1f}, mean={vals[0]:.1f}, SD=N/A"
        return (
            f"median={statistics.median(vals):.1f}, "
            f"mean={statistics.mean(vals):.1f}, "
            f"SD={statistics.stdev(vals):.1f}"
        )

    diameter_used = meta.get("diameter") or "auto-estimated"
    device = "GPU" if meta.get("use_gpu") else "CPU"
    edge_excluded = "yes" if meta.get("exclude_on_edges") else "no"
    flow_threshold = meta.get("flow_threshold", 0.4)
    cellprob_threshold = meta.get("cellprob_threshold", 0.0)
    z_projection = meta.get("z_projection", "max")

    lines = [
        "# Cell Segmentation Report",
        "",
        f"**Image:** {meta.get('image_path', 'demo')}",
        "**Backend:** cpsam (Cellpose 4.0)",
        f"**Device:** {device}",
        f"**Diameter used:** {diameter_used} px",
        f"**Exclude edge cells:** {edge_excluded}",
        f"**Flow threshold:** {flow_threshold}",
        f"**Cellprob threshold:** {cellprob_threshold}",
        f"**Z projection policy:** {z_projection}",
        "",
        "## Results",
        "",
        f"- **Cells detected:** {n}",
        f"- **Area (px²):** {_stats(areas)}",
        f"- **Diameter (px):** {_stats(diameters)}",
        "",
        "## Output Files",
        "",
        "| File | Description |",
        "|------|-------------|",
        f"| `{csv_filename}` | Per-cell metrics (area, diameter, centroid, eccentricity) |",
        f"| `{masks_filename}` | Label image (uint16, each cell a unique integer) |",
        f"| `{seg_filename}` | Cellpose seg dict (masks + flows, reload with `np.load(..., allow_pickle=True)`) |",
        f"| `figures/{outlines_filename}` | Original image with cell outlines (cellpose --save_outlines) |",
        f"| `figures/{histogram_filename}` | Histogram of cell equivalent diameters |",
        "",
        "---",
        "",
        f"> {DISCLAIMER}",
    ]
    (output_dir / "report.md").write_text("\n".join(lines))


def save_outlines(img: np.ndarray, masks: np.ndarray, flows: list, output_dir: Path | str, stem: str = "image") -> str:
    """Save outlines via cellpose's built-in --save_outlines mechanism.

    Produces <stem>_cp_outlines.png in output_dir/figures/.
    Returns the filename produced so callers can report it.
    """
    from cellpose import io as cp_io

    fig_dir = Path(output_dir) / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    if masks.ndim != 2:
        # Cellpose does not support saving 3D outlines as PNG. Emit a note so
        # report generation remains consistent for volumetric runs.
        note_name = f"{stem}_cp_outlines_unavailable.txt"
        (fig_dir / note_name).write_text(
            "Outlines PNG not generated for volumetric segmentation; "
            "inspect label volume in napari instead."
        )
        return note_name

    cp_io.save_masks(
        images=[img],
        masks=[masks],
        flows=[flows],
        file_names=[stem],
        save_outlines=True,
        savedir=str(fig_dir),
        png=False,
        tif=False,
    )
    # Cellpose saves {stem}_outlines_cp_masks.png into fig_dir; rename to {stem}_cp_outlines.png
    src = fig_dir / f"{stem}_outlines_cp_masks.png"
    final_name = f"{stem}_cp_outlines.png"
    if src.exists():
        target = fig_dir / final_name
        if target.exists():
            target.unlink()
        src.rename(target)
    return final_name


def save_histogram(metrics: list[dict], output_dir: Path | str, stem: str = "image") -> str:
    """Save {stem}_histogram.png to output_dir/figures/. Returns the filename."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig_dir = Path(output_dir) / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    diameters = [m["diameter"] for m in metrics]
    fig, ax = plt.subplots(figsize=(6, 4))
    if diameters:
        ax.hist(diameters, bins=min(20, len(diameters)), color="steelblue", edgecolor="white")
    ax.set_xlabel("Equivalent diameter (px)")
    ax.set_ylabel("Cell count")
    ax.set_title("Cell size distribution")
    fig.tight_layout()
    filename = f"{stem}_histogram.png"
    fig.savefig(fig_dir / filename, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return filename


def save_masks(masks: np.ndarray, output_dir: Path | str, stem: str = "image") -> str:
    """Write masks as a uint16 TIFF label image. Returns the filename produced."""
    import tifffile
    filename = f"{stem}_cp_masks.tif"
    tifffile.imwrite(str(Path(output_dir) / filename), masks.astype(np.uint16))
    return filename


def save_seg_npy(masks: np.ndarray, flows: list, output_dir: Path | str, stem: str = "image") -> str:
    """Save cellpose seg.npy (masks + flows) for downstream reuse. Returns the filename."""
    filename = f"{stem}_seg.npy"
    np.save(str(Path(output_dir) / filename), {"masks": masks, "flows": flows}, allow_pickle=True)
    return filename


def run_segmentation(
    img: np.ndarray,
    diameter: float | None,
    use_gpu: bool,
    do_3D: bool = False,
    exclude_on_edges: bool = False,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
) -> tuple[np.ndarray, list]:
    """Run cpsam on an image (greyscale H×W, multi-channel H×W×C, z-stack Z×H×W, or Z×C×H×W).

    For 2D mode, cpsam is channel-order invariant and no explicit channel
    mapping is passed.
    Pass do_3D=True for volumetric segmentation of (Z, H, W) or (Z, C, H, W) stacks.
    For 3D multi-channel inputs (Z, C, H, W), channel_axis=1 is passed.
    Pass exclude_on_edges=True to remove any cell whose mask touches a border
    pixel (mirrors cellpose CLI --exclude_on_edges).
    flow_threshold: max allowed error of the flows for each mask (default 0.4).
      Increase to accept more masks; decrease to be more stringent.
    cellprob_threshold: threshold on the cell probability output (default 0.0).
      Decrease (e.g. -6) to detect more/dimmer cells; increase to be stricter.
    Returns (uint16 label mask, flows) — flows are needed for save_outlines.
    """
    from cellpose.models import CellposeModel
    from cellpose import utils as cp_utils

    model = CellposeModel(gpu=use_gpu)
    eval_kwargs = {
        "diameter": diameter,
        "do_3D": do_3D,
        "z_axis": 0 if do_3D else None,
        "flow_threshold": flow_threshold,
        "cellprob_threshold": cellprob_threshold,
    }
    if do_3D and img.ndim == 4:
        # For volumetric multi-channel inputs, keep image as Z×C×Y×X.
        eval_kwargs["channel_axis"] = 1

    masks, flows, _ = model.eval(img, **eval_kwargs)
    masks = masks.astype(np.uint16)
    if exclude_on_edges:
        masks = cp_utils.remove_edge_masks(masks).astype(np.uint16)
    return masks, flows


def _detect_gpu(requested: bool) -> bool:
    """Return True if GPU should be used (requested AND available)."""
    if not requested:
        return False
    try:
        import torch

        return torch.cuda.is_available() or torch.backends.mps.is_available()
    except ImportError:
        return False


def _repro_path(value: Path, *, output_dir: Path) -> ReproPath:
    """Classify a path for portable reproducibility command rendering."""
    try:
        value.relative_to(_PROJECT_ROOT)
        return ReproPath(value, anchor="repo_root")
    except ValueError:
        pass

    try:
        value.relative_to(output_dir)
        return ReproPath(value, anchor="output_dir")
    except ValueError:
        return ReproPath(value, anchor="auto")


def repro_command_for_bundle(args: argparse.Namespace, output_dir: Path) -> ReproCommand:
    """Build a structured reproducibility command for this skill."""
    cmd_args: list[str | ReproPath] = []
    if args.demo:
        cmd_args.append("--demo")
    else:
        cmd_args.extend(["--input", _repro_path(Path(args.input), output_dir=output_dir)])

    if args.diameter is not None:
        cmd_args.extend(["--diameter", str(args.diameter)])
    if args.gpu is True:
        cmd_args.append("--use_gpu")
    if args.use_cpu:
        cmd_args.append("--use_cpu")
    if args.do_3D:
        cmd_args.append("--do_3D")
    if args.exclude_on_edges:
        cmd_args.append("--exclude_on_edges")
    if args.flow_threshold != 0.4:
        cmd_args.extend(["--flow_threshold", str(args.flow_threshold)])
    if args.cellprob_threshold != 0.0:
        cmd_args.extend(["--cellprob_threshold", str(args.cellprob_threshold)])
    if args.z_projection != "max":
        cmd_args.extend(["--z_projection", args.z_projection])

    cmd_args.extend(["--output", ReproPath(output_dir, anchor="output_dir")])

    return ReproCommand(
        script_path=Path("skills/cell-detection/cell_detection.py"),
        args=cmd_args,
        comment="Replay this ClawBio cell-detection run",
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="cell-detection — cell segmentation using cpsam (Cellpose 4.0)"
    )
    parser.add_argument("--input", help="Input image (TIFF, CZI, ND2, PNG, JPG)")
    parser.add_argument("--diameter", type=float, default=None, help="Cell diameter in pixels (default: auto)")
    parser.add_argument(
        "--use_gpu",
        dest="gpu",
        action="store_true",
        default=None,
        help="Prefer GPU acceleration (default behavior: auto-use GPU when available).",
    )
    parser.add_argument(
        "--use_cpu",
        action="store_true",
        default=False,
        help="Force CPU execution even when a GPU backend is available.",
    )
    parser.add_argument("--do_3D", action="store_true", default=False, help="Volumetric 3D segmentation for Z×H×W or Z×C×H×W input")
    parser.add_argument("--exclude_on_edges", action="store_true", default=False, help="Remove cells touching the image border (mirrors cellpose CLI --exclude_on_edges)")
    parser.add_argument("--flow_threshold", type=float, default=0.4, help="Max flow error per mask (default: 0.4). Increase to accept more masks, decrease for stricter filtering.")
    parser.add_argument("--cellprob_threshold", type=float, default=0.0, help="Cell probability threshold (default: 0.0). Decrease (e.g. -6) to detect dimmer/more cells, increase for stricter detection.")
    parser.add_argument("--z_projection", choices=("max", "none"), default="max", help="How to reduce 4D stacks (C/Z/Y/X). 'max' performs Z max projection (default). 'none' preserves Z; multi-channel stacks remain Z×C×H×W.")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run on a synthetic demo image")
    args = parser.parse_args()

    # Guardrails for mode consistency:
    # - 2D mode must use z_projection=max
    # - 3D mode must use z_projection=none
    # Exception: if --do_3D is requested for a truly 2D single-channel input,
    # fall back to 2D mode so the pipeline still runs.
    if not args.do_3D and args.z_projection != "max":
        parser.error("--z_projection none is only valid with --do_3D")

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    effective_do_3d = args.do_3D

    # Build input image
    if args.demo:
        if args.do_3D:
            print(
                "[cell-detection] Warning: --do_3D requested for 2D demo image; "
                "falling back to 2D mode."
            )
            effective_do_3d = False
            args.do_3D = False
        img_prep = make_demo_image()
        image_path = "demo (synthetic fluorescence nuclei — offline)"
        stem = "demo"
    else:
        if not args.input:
            parser.error("--input is required unless --demo is used")

        if args.do_3D and args.z_projection != "none":
            img_none, n_channels_none = load_image(args.input, z_projection="none")
            if img_none.ndim == 2 and n_channels_none == 1:
                print(
                    "[cell-detection] Warning: --do_3D requested for 2D single-channel "
                    "input with z_projection=max; falling back to 2D mode."
                )
                effective_do_3d = False
                args.do_3D = False
            else:
                parser.error("--do_3D requires --z_projection none for volumetric inputs")

        img, n_channels = load_image(args.input, z_projection=args.z_projection)
        if args.do_3D and img.ndim == 2 and n_channels == 1:
            print(
                "[cell-detection] Warning: --do_3D requested for 2D single-channel "
                "input; falling back to 2D mode."
            )
            effective_do_3d = False
            args.do_3D = False
        if effective_do_3d and img.ndim not in (3, 4):
            parser.error(
                "--do_3D requires volumetric input with shape Z×Y×X or Z×C×Y×X"
            )
        if effective_do_3d and img.ndim == 4 and img.shape[0] < 2:
            parser.error(
                "--do_3D requires a true Z stack; got 4D input with fewer than 2 Z slices"
            )
        image_path = args.input
        # Strip all compound suffixes (e.g. .ome.tif → base stem)
        stem = Path(args.input).name.split(".")[0]
        if effective_do_3d:
            img_prep = img  # pass Z×H×W directly to cellpose
        else:
            img_prep = prepare_image(img, n_channels)

    # Segmentation defaults to GPU auto-detection unless explicitly forced to CPU.
    gpu_requested = True if args.demo else (True if args.gpu is None else args.gpu)
    if args.use_cpu:
        gpu_requested = False
    use_gpu = _detect_gpu(gpu_requested)
    device_label = "GPU" if use_gpu else "CPU"
    print(f"[cell-detection] Segmenting with cpsam on {device_label}...")
    masks, flows = run_segmentation(img_prep, args.diameter, use_gpu, do_3D=effective_do_3d, exclude_on_edges=args.exclude_on_edges, flow_threshold=args.flow_threshold, cellprob_threshold=args.cellprob_threshold)

    # Save masks + seg.npy
    masks_filename = save_masks(masks, output_dir, stem=stem)
    seg_filename = save_seg_npy(masks, flows, output_dir, stem=stem)

    # Metrics
    metrics = compute_metrics(masks)

    # Save CSV
    csv_path = output_dir / f"{stem}_measurements.csv"
    with open(csv_path, "w", newline="") as f:
        fieldnames = ["id", "area", "diameter", "centroid_x", "centroid_y", "eccentricity"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(metrics)

    # Outlines image (cellpose --save_outlines)
    outlines_filename = save_outlines(img_prep, masks, flows, output_dir, stem=stem)

    # Size distribution
    save_histogram(metrics, output_dir, stem=stem)

    # Report
    effective_z_projection = "none" if effective_do_3d else "max"
    meta = {"image_path": image_path, "use_gpu": use_gpu, "diameter": args.diameter, "exclude_on_edges": args.exclude_on_edges, "flow_threshold": args.flow_threshold, "cellprob_threshold": args.cellprob_threshold, "z_projection": effective_z_projection}
    write_report(metrics, meta, output_dir, outlines_filename=outlines_filename, masks_filename=masks_filename, seg_filename=seg_filename, csv_filename=f"{stem}_measurements.csv", histogram_filename=f"{stem}_histogram.png")

    # Reproducibility — environment.yml must be written first so REPLAY.md
    # can include the env name from the file.
    write_environment_yml(
        output_dir,
        env_name="clawbio-cell-detection",
        pip_deps=["cellpose>=4.0", "tifffile", "czifile>=2019.7.2.2", "nd2>=0.11.1", "Pillow", "numpy", "matplotlib", "scikit-image", "scipy"],
        python_version="3.10",
    )
    write_portable_commands_sh(
        output_dir,
        repro_command_for_bundle(args, output_dir),
        repo_root=_PROJECT_ROOT,
    )
    write_checksums(
        [
            output_dir / "report.md",
            output_dir / masks_filename,
            output_dir / seg_filename,
            csv_path,
            output_dir / "figures" / outlines_filename,
            output_dir / "figures" / f"{stem}_histogram.png",
        ],
        output_dir,
        anchor=output_dir,
    )

    print(f"[cell-detection] Done — {len(metrics)} cells detected.")
    print(f"[cell-detection] Report: {output_dir / 'report.md'}")


if __name__ == "__main__":
    main()
