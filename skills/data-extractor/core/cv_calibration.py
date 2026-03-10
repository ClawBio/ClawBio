"""Computer vision calibration for plot data extraction.

Uses OpenCV to detect:
1. Plot area boundaries (axis lines)
2. Colored marker positions (dots, squares, diamonds)
3. Bar top positions

Returns pixel-level calibration data that gets injected into the
Claude extraction prompt, so Claude maps pixel positions to axis
values instead of eyeballing.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import cv2
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class PlotRegion:
    """Detected plot area boundaries in pixel coordinates."""
    x_left: int  # left edge of plot area (y-axis line)
    x_right: int  # right edge of plot area
    y_top: int  # top edge of plot area
    y_bottom: int  # bottom edge (x-axis line)

    @property
    def width(self) -> int:
        return self.x_right - self.x_left

    @property
    def height(self) -> int:
        return self.y_bottom - self.y_top


@dataclass
class DetectedMarker:
    """A detected data marker with pixel position and color."""
    px: int  # pixel x
    py: int  # pixel y
    color_name: str  # e.g., "red", "blue", "green"
    color_bgr: tuple[int, int, int]  # average BGR color
    area: float  # contour area (for size filtering)


@dataclass
class DetectedBar:
    """A detected bar with position and extent."""
    px_center: int  # bar center x pixel
    py_top: int  # bar top y pixel
    px_left: int
    px_right: int
    color_name: str
    color_bgr: tuple[int, int, int]


@dataclass
class CalibrationResult:
    """Full calibration data for a plot image."""
    plot_region: PlotRegion | None = None
    markers: list[DetectedMarker] = field(default_factory=list)
    bars: list[DetectedBar] = field(default_factory=list)
    image_width: int = 0
    image_height: int = 0


# --- Color definitions for common plot markers ---
# HSV ranges: (H_low, S_low, V_low), (H_high, S_high, V_high)
# OpenCV HSV: H=0-179, S=0-255, V=0-255
COLOR_RANGES = {
    "red": [
        ((0, 70, 70), (10, 255, 255)),      # red low hue
        ((170, 70, 70), (179, 255, 255)),    # red high hue (wraps around)
    ],
    "blue": [
        ((100, 70, 50), (130, 255, 255)),
    ],
    "green": [
        ((35, 70, 50), (85, 255, 255)),
    ],
    "orange": [
        ((10, 100, 100), (25, 255, 255)),
    ],
    "purple": [
        ((130, 50, 50), (160, 255, 255)),
    ],
    "cyan": [
        ((80, 70, 50), (100, 255, 255)),
    ],
    "yellow": [
        ((25, 100, 100), (35, 255, 255)),
    ],
    "pink": [
        ((160, 50, 70), (170, 255, 255)),
    ],
}


def _load_image(image_bytes: bytes) -> np.ndarray:
    """Load image from bytes into OpenCV BGR array."""
    arr = np.frombuffer(image_bytes, dtype=np.uint8)
    img = cv2.imdecode(arr, cv2.IMREAD_COLOR)
    if img is None:
        raise ValueError("Failed to decode image")
    return img


def detect_plot_region(img: np.ndarray) -> PlotRegion | None:
    """Detect the plot area by finding axis lines.

    Looks for the L-shaped axis frame (vertical left line + horizontal bottom line).
    Uses Hough line detection on edges.
    """
    h, w = img.shape[:2]
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    edges = cv2.Canny(gray, 50, 150, apertureSize=3)

    # Detect lines
    lines = cv2.HoughLinesP(edges, 1, np.pi / 180, threshold=80,
                            minLineLength=min(w, h) * 0.2, maxLineGap=10)
    if lines is None:
        return None

    # Separate horizontal and vertical lines
    h_lines = []  # (y, x1, x2)
    v_lines = []  # (x, y1, y2)

    for line in lines:
        x1, y1, x2, y2 = line[0]
        dx = abs(x2 - x1)
        dy = abs(y2 - y1)

        if dy < 5 and dx > w * 0.15:  # horizontal
            h_lines.append((min(y1, y2), min(x1, x2), max(x1, x2)))
        elif dx < 5 and dy > h * 0.15:  # vertical
            v_lines.append((min(x1, x2), min(y1, y2), max(y1, y2)))

    if not h_lines or not v_lines:
        # Fallback: estimate plot region as inner 80% of image
        margin_x = int(w * 0.12)
        margin_y = int(h * 0.08)
        return PlotRegion(margin_x, w - int(w * 0.05), margin_y, h - int(h * 0.15))

    # Bottom x-axis: horizontal line with largest y (closest to bottom)
    h_lines.sort(key=lambda l: l[0], reverse=True)
    x_axis_y = h_lines[0][0]

    # Left y-axis: vertical line with smallest x (closest to left)
    v_lines.sort(key=lambda l: l[0])
    y_axis_x = v_lines[0][0]

    # Right boundary: rightmost extent of horizontal lines near bottom
    right_x = max(l[2] for l in h_lines if abs(l[0] - x_axis_y) < 20)

    # Top boundary: topmost extent of vertical lines near left axis
    top_y = min(l[1] for l in v_lines if abs(l[0] - y_axis_x) < 20)

    return PlotRegion(
        x_left=y_axis_x,
        x_right=min(right_x, w - 5),
        y_top=max(top_y, 5),
        y_bottom=x_axis_y,
    )


def detect_markers(img: np.ndarray, plot_region: PlotRegion | None = None,
                   min_area: float = 15, max_area: float = 2000) -> list[DetectedMarker]:
    """Detect colored data point markers using HSV color segmentation.

    Returns centroids of detected markers within the plot region.
    """
    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
    markers: list[DetectedMarker] = []

    # Define region of interest (with padding to catch edge markers)
    if plot_region:
        pad = 10
        roi_mask = np.zeros(img.shape[:2], dtype=np.uint8)
        y1 = max(0, plot_region.y_top - pad)
        y2 = min(img.shape[0], plot_region.y_bottom + pad)
        x1 = max(0, plot_region.x_left - pad)
        x2 = min(img.shape[1], plot_region.x_right + pad)
        roi_mask[y1:y2, x1:x2] = 255
    else:
        roi_mask = np.ones(img.shape[:2], dtype=np.uint8) * 255

    for color_name, ranges in COLOR_RANGES.items():
        # Combine all HSV ranges for this color
        combined_mask = np.zeros(img.shape[:2], dtype=np.uint8)
        for (low, high) in ranges:
            mask = cv2.inRange(hsv, np.array(low), np.array(high))
            combined_mask = cv2.bitwise_or(combined_mask, mask)

        # Apply ROI and clean up
        combined_mask = cv2.bitwise_and(combined_mask, roi_mask)
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
        combined_mask = cv2.morphologyEx(combined_mask, cv2.MORPH_CLOSE, kernel)
        combined_mask = cv2.morphologyEx(combined_mask, cv2.MORPH_OPEN, kernel)

        # Find contours
        contours, _ = cv2.findContours(combined_mask, cv2.RETR_EXTERNAL,
                                       cv2.CHAIN_APPROX_SIMPLE)

        for contour in contours:
            area = cv2.contourArea(contour)
            if area < min_area or area > max_area:
                continue

            M = cv2.moments(contour)
            if M["m00"] == 0:
                continue

            cx = int(M["m10"] / M["m00"])
            cy = int(M["m01"] / M["m00"])

            # Get average color at this point
            mask_point = np.zeros(img.shape[:2], dtype=np.uint8)
            cv2.drawContours(mask_point, [contour], -1, 255, -1)
            mean_color = cv2.mean(img, mask=mask_point)[:3]

            markers.append(DetectedMarker(
                px=cx, py=cy,
                color_name=color_name,
                color_bgr=(int(mean_color[0]), int(mean_color[1]), int(mean_color[2])),
                area=area,
            ))

    # Sort by x position
    markers.sort(key=lambda m: (m.px, m.py))
    return markers


def detect_bars(img: np.ndarray, plot_region: PlotRegion | None = None) -> list[DetectedBar]:
    """Detect bars in a bar chart by finding colored rectangular regions.

    Returns bar center x, top y, and color for each detected bar.
    """
    if plot_region is None:
        return []

    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
    bars: list[DetectedBar] = []

    for color_name, ranges in COLOR_RANGES.items():
        combined_mask = np.zeros(img.shape[:2], dtype=np.uint8)
        for (low, high) in ranges:
            mask = cv2.inRange(hsv, np.array(low), np.array(high))
            combined_mask = cv2.bitwise_or(combined_mask, mask)

        # Only look in plot region
        roi_mask = np.zeros(img.shape[:2], dtype=np.uint8)
        roi_mask[plot_region.y_top:plot_region.y_bottom,
                 plot_region.x_left:plot_region.x_right] = 255
        combined_mask = cv2.bitwise_and(combined_mask, roi_mask)

        # Clean up — use larger kernel for bars
        kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (5, 5))
        combined_mask = cv2.morphologyEx(combined_mask, cv2.MORPH_CLOSE, kernel)

        contours, _ = cv2.findContours(combined_mask, cv2.RETR_EXTERNAL,
                                       cv2.CHAIN_APPROX_SIMPLE)

        for contour in contours:
            area = cv2.contourArea(contour)
            if area < 200:  # bars should be sizable
                continue

            x, y, bw, bh = cv2.boundingRect(contour)
            aspect = bh / max(bw, 1)

            # Bars are tall rectangles (aspect > 0.5) or wide if horizontal
            if aspect < 0.3 and bw / max(bh, 1) < 0.3:
                continue  # not bar-shaped

            # Get average color
            mask_bar = np.zeros(img.shape[:2], dtype=np.uint8)
            cv2.drawContours(mask_bar, [contour], -1, 255, -1)
            mean_color = cv2.mean(img, mask=mask_bar)[:3]

            bars.append(DetectedBar(
                px_center=x + bw // 2,
                py_top=y,
                px_left=x,
                px_right=x + bw,
                color_name=color_name,
                color_bgr=(int(mean_color[0]), int(mean_color[1]), int(mean_color[2])),
            ))

    # Sort by x position
    bars.sort(key=lambda b: b.px_center)
    return bars


def calibrate_image(image_bytes: bytes) -> CalibrationResult:
    """Run full calibration pipeline on a plot image.

    Returns detected plot region, markers, and bars with pixel coordinates.
    """
    try:
        img = _load_image(image_bytes)
    except Exception as e:
        logger.warning(f"Failed to load image for calibration: {e}")
        return CalibrationResult()

    h, w = img.shape[:2]
    result = CalibrationResult(image_width=w, image_height=h)

    # Step 1: Detect plot region
    result.plot_region = detect_plot_region(img)
    if result.plot_region:
        logger.info(
            f"Plot region: ({result.plot_region.x_left},{result.plot_region.y_top}) "
            f"to ({result.plot_region.x_right},{result.plot_region.y_bottom})"
        )

    # Step 2: Detect markers
    result.markers = detect_markers(img, result.plot_region)
    logger.info(f"Detected {len(result.markers)} markers")

    # Step 3: Detect bars
    result.bars = detect_bars(img, result.plot_region)
    logger.info(f"Detected {len(result.bars)} bars")

    return result


def format_calibration_prompt(cal: CalibrationResult) -> str:
    """Format calibration data as text to inject into the extraction prompt.

    Tells Claude the exact pixel positions of detected features so it can
    map them to axis values instead of eyeballing.
    """
    if not cal.plot_region and not cal.markers and not cal.bars:
        return ""

    parts: list[str] = []
    parts.append("COMPUTER VISION CALIBRATION DATA (use this to improve accuracy):")

    if cal.plot_region:
        pr = cal.plot_region
        parts.append(
            f"Plot area spans pixels: x=[{pr.x_left}..{pr.x_right}] ({pr.width}px wide), "
            f"y=[{pr.y_top}..{pr.y_bottom}] ({pr.height}px tall). "
            f"Image size: {cal.image_width}x{cal.image_height}px."
        )
        parts.append(
            "AXIS CALIBRATION INSTRUCTIONS:\n"
            "Step 1: Read ALL tick labels on the x-axis and y-axis.\n"
            "Step 2: The leftmost x-axis tick is near pixel x=" + str(pr.x_left) +
            ", the rightmost near x=" + str(pr.x_right) + ".\n"
            "Step 3: The bottom y-axis tick is near pixel y=" + str(pr.y_bottom) +
            ", the top near y=" + str(pr.y_top) + " (y increases downward in pixels).\n"
            "Step 4: For each data point below, interpolate its value using "
            "the tick positions as reference."
        )

    if cal.markers:
        # Group markers by color
        by_color: dict[str, list[DetectedMarker]] = {}
        for m in cal.markers:
            by_color.setdefault(m.color_name, []).append(m)

        parts.append(f"\nDETECTED DATA POINTS ({len(cal.markers)} total):")
        for color, mks in sorted(by_color.items()):
            coords = ", ".join(f"({m.px},{m.py})" for m in mks)
            parts.append(f"  {color} markers ({len(mks)}): {coords}")

        parts.append(
            "Each (px,py) is a pixel coordinate. Map each to (x_value, y_value) "
            "by interpolating between the nearest axis ticks. "
            "Group markers by color into series using the legend."
        )

    if cal.bars:
        parts.append(f"\nDETECTED BARS ({len(cal.bars)} total):")
        by_color: dict[str, list[DetectedBar]] = {}
        for b in cal.bars:
            by_color.setdefault(b.color_name, []).append(b)

        for color, bs in sorted(by_color.items()):
            coords = ", ".join(
                f"(center_x={b.px_center}, top_y={b.py_top})" for b in bs
            )
            parts.append(f"  {color} bars ({len(bs)}): {coords}")

        parts.append(
            "Bar top y-pixel tells you the bar height — map to y-axis value. "
            "Bar center x-pixel tells you which category — map to x-axis label."
        )

    return "\n".join(parts)
