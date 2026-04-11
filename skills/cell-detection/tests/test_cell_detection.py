"""Tests for the CellposeSAM cell segmentation skill."""
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

import cell_detection


# ---------------------------------------------------------------------------
# TestDemoImage
# ---------------------------------------------------------------------------


class TestDemoImage:
    """Tests for the synthetic demo image generator."""

    def test_shape(self):
        img = cell_detection.make_demo_image()
        assert img.shape == (512, 512)

    def test_dtype(self):
        img = cell_detection.make_demo_image()
        assert img.dtype == np.uint8

    def test_value_range(self):
        img = cell_detection.make_demo_image()
        assert img.min() >= 0
        assert img.max() <= 255

    def test_reproducible(self):
        img1 = cell_detection.make_demo_image()
        img2 = cell_detection.make_demo_image()
        np.testing.assert_array_equal(img1, img2)

    def test_not_blank(self):
        """Demo image has non-zero pixels (blobs are present)."""
        img = cell_detection.make_demo_image()
        assert img.max() > 0

    def test_different_seeds_differ(self):
        img1 = cell_detection.make_demo_image(seed=1)
        img2 = cell_detection.make_demo_image(seed=2)
        assert not np.array_equal(img1, img2)


# ---------------------------------------------------------------------------
# TestPrepareImage
# ---------------------------------------------------------------------------


class TestPrepareImage:
    """Tests for image preparation / channel handling."""

    def test_greyscale_passthrough(self):
        img = np.zeros((64, 64), dtype=np.uint8)
        out = cell_detection.prepare_image(img, n_channels=1)
        assert out is img

    def test_two_channel_passthrough(self):
        img = np.zeros((64, 64, 2), dtype=np.uint8)
        out = cell_detection.prepare_image(img, n_channels=2)
        assert out is img

    def test_three_channel_passthrough(self):
        img = np.zeros((64, 64, 3), dtype=np.uint8)
        out = cell_detection.prepare_image(img, n_channels=3)
        assert out is img

    def test_four_channel_truncated(self, capsys):
        img = np.zeros((64, 64, 4), dtype=np.uint8)
        out = cell_detection.prepare_image(img, n_channels=4)
        assert out.shape == (64, 64, 3)
        captured = capsys.readouterr()
        assert "Truncating" in captured.out

    def test_many_channels_truncated(self):
        img = np.zeros((32, 32, 10), dtype=np.uint8)
        out = cell_detection.prepare_image(img, n_channels=10)
        assert out.shape == (32, 32, 3)

    def test_truncation_keeps_first_three_channels(self):
        img = np.zeros((64, 64, 5), dtype=np.uint8)
        img[:, :, 0] = 10
        img[:, :, 1] = 20
        img[:, :, 2] = 30
        out = cell_detection.prepare_image(img, n_channels=5)
        assert out[:, :, 0].mean() == 10
        assert out[:, :, 1].mean() == 20
        assert out[:, :, 2].mean() == 30


# ---------------------------------------------------------------------------
# TestComputeMetrics
# ---------------------------------------------------------------------------


class TestComputeMetrics:
    """Tests for per-cell morphology metric extraction."""

    def _single_cell_mask(self):
        """10x10 mask with one 4x4 square cell at label=1."""
        masks = np.zeros((10, 10), dtype=np.uint16)
        masks[3:7, 3:7] = 1
        return masks

    def test_empty_mask_returns_empty(self):
        masks = np.zeros((64, 64), dtype=np.uint16)
        rows = cell_detection.compute_metrics(masks)
        assert rows == []

    def test_single_cell_count(self):
        masks = self._single_cell_mask()
        rows = cell_detection.compute_metrics(masks)
        assert len(rows) == 1

    def test_single_cell_keys(self):
        masks = self._single_cell_mask()
        row = cell_detection.compute_metrics(masks)[0]
        for key in ("id", "area", "diameter", "centroid_x", "centroid_y", "eccentricity"):
            assert key in row, f"Missing key: {key}"

    def test_single_cell_area(self):
        masks = self._single_cell_mask()
        row = cell_detection.compute_metrics(masks)[0]
        assert row["area"] == 16  # 4x4 = 16 px

    def test_single_cell_id(self):
        masks = self._single_cell_mask()
        row = cell_detection.compute_metrics(masks)[0]
        assert row["id"] == 1

    def test_two_cells(self):
        masks = np.zeros((20, 20), dtype=np.uint16)
        masks[2:5, 2:5] = 1
        masks[12:16, 12:16] = 2
        rows = cell_detection.compute_metrics(masks)
        assert len(rows) == 2
        ids = {r["id"] for r in rows}
        assert ids == {1, 2}

    def test_eccentricity_range(self):
        masks = self._single_cell_mask()
        row = cell_detection.compute_metrics(masks)[0]
        assert 0.0 <= row["eccentricity"] <= 1.0

    def test_diameter_positive(self):
        masks = self._single_cell_mask()
        row = cell_detection.compute_metrics(masks)[0]
        assert row["diameter"] > 0


# ---------------------------------------------------------------------------
# TestWriteReport
# ---------------------------------------------------------------------------


class TestWriteReport:
    """Tests for the Markdown report writer."""

    def _sample_metrics(self):
        return [
            {"id": 1, "area": 200, "diameter": 16.0, "centroid_x": 10.0, "centroid_y": 10.0, "eccentricity": 0.2},
            {"id": 2, "area": 150, "diameter": 13.8, "centroid_x": 30.0, "centroid_y": 30.0, "eccentricity": 0.5},
        ]

    def test_report_created(self, tmp_path):
        cell_detection.write_report(self._sample_metrics(), {}, tmp_path)
        assert (tmp_path / "report.md").exists()

    def test_report_heading(self, tmp_path):
        cell_detection.write_report(self._sample_metrics(), {}, tmp_path)
        text = (tmp_path / "report.md").read_text()
        assert "# Cell Segmentation Report" in text

    def test_report_cell_count(self, tmp_path):
        cell_detection.write_report(self._sample_metrics(), {}, tmp_path)
        text = (tmp_path / "report.md").read_text()
        assert "Cells detected:** 2" in text

    def test_report_disclaimer(self, tmp_path):
        cell_detection.write_report(self._sample_metrics(), {}, tmp_path)
        text = (tmp_path / "report.md").read_text()
        assert "ClawBio is a research and educational tool" in text

    def test_report_gpu_device(self, tmp_path):
        cell_detection.write_report(self._sample_metrics(), {"use_gpu": True}, tmp_path)
        text = (tmp_path / "report.md").read_text()
        assert "GPU" in text

    def test_report_cpu_device(self, tmp_path):
        cell_detection.write_report(self._sample_metrics(), {"use_gpu": False}, tmp_path)
        text = (tmp_path / "report.md").read_text()
        assert "CPU" in text

    def test_report_outlines_filename(self, tmp_path):
        cell_detection.write_report(self._sample_metrics(), {}, tmp_path, outlines_filename="demo_cp_outlines.png")
        text = (tmp_path / "report.md").read_text()
        assert "demo_cp_outlines.png" in text

    def test_report_empty_metrics(self, tmp_path):
        """Empty metrics list produces N/A stats and 0 cells."""
        cell_detection.write_report([], {}, tmp_path)
        text = (tmp_path / "report.md").read_text()
        assert "Cells detected:** 0" in text
        assert "N/A" in text

    def test_report_single_cell_no_sd(self, tmp_path):
        """Single cell produces SD=N/A (stdev undefined for n=1)."""
        metrics = [{"id": 1, "area": 100, "diameter": 11.3, "centroid_x": 5.0, "centroid_y": 5.0, "eccentricity": 0.0}]
        cell_detection.write_report(metrics, {}, tmp_path)
        text = (tmp_path / "report.md").read_text()
        assert "SD=N/A" in text

    def test_report_diameter_auto(self, tmp_path):
        """None diameter renders as auto-estimated."""
        cell_detection.write_report(self._sample_metrics(), {"diameter": None}, tmp_path)
        text = (tmp_path / "report.md").read_text()
        assert "auto-estimated" in text

    def test_report_diameter_explicit(self, tmp_path):
        """Explicit diameter value is shown in report."""
        cell_detection.write_report(self._sample_metrics(), {"diameter": 30.0}, tmp_path)
        text = (tmp_path / "report.md").read_text()
        assert "30.0" in text


# ---------------------------------------------------------------------------
# TestSaveMasks
# ---------------------------------------------------------------------------


class TestSaveMasks:
    """Tests for TIFF mask output."""

    def test_masks_tif_created(self, tmp_path):
        import tifffile
        masks = np.zeros((64, 64), dtype=np.uint16)
        masks[10:20, 10:20] = 1
        cell_detection.save_masks(masks, tmp_path)
        assert (tmp_path / "image_cp_masks.tif").exists()

    def test_masks_tif_roundtrip(self, tmp_path):
        import tifffile
        masks = np.zeros((32, 32), dtype=np.uint16)
        masks[5:10, 5:10] = 3
        cell_detection.save_masks(masks, tmp_path)
        loaded = tifffile.imread(str(tmp_path / "image_cp_masks.tif"))
        np.testing.assert_array_equal(loaded, masks)

    def test_masks_tif_dtype(self, tmp_path):
        import tifffile
        masks = np.ones((16, 16), dtype=np.int32)  # non-uint16 input
        cell_detection.save_masks(masks, tmp_path)
        loaded = tifffile.imread(str(tmp_path / "image_cp_masks.tif"))
        assert loaded.dtype == np.uint16


# ---------------------------------------------------------------------------
# TestSaveSizeDistribution
# ---------------------------------------------------------------------------


class TestSaveSizeDistribution:
    """Tests for the cell size histogram figure."""

    def test_png_created(self, tmp_path):
        metrics = [
            {"diameter": 12.0}, {"diameter": 15.5}, {"diameter": 9.3},
        ]
        cell_detection.save_histogram(metrics, tmp_path)
        assert (tmp_path / "figures" / "image_histogram.png").exists()

    def test_empty_metrics_no_error(self, tmp_path):
        """Empty metrics produces a blank histogram without raising."""
        cell_detection.save_histogram([], tmp_path)
        assert (tmp_path / "figures" / "image_histogram.png").exists()


# ---------------------------------------------------------------------------
# TestDetectGpu
# ---------------------------------------------------------------------------


class TestDetectGpu:
    """Tests for GPU auto-detection."""

    def test_false_when_not_requested(self):
        assert cell_detection._detect_gpu(requested=False) is False

    def test_false_when_torch_missing(self):
        with patch.dict("sys.modules", {"torch": None}):
            result = cell_detection._detect_gpu(requested=True)
        assert result is False

    def test_false_when_no_cuda_no_mps(self):
        import types
        fake_torch = types.ModuleType("torch")
        fake_torch.cuda = types.SimpleNamespace(is_available=lambda: False)
        fake_torch.backends = types.SimpleNamespace(
            mps=types.SimpleNamespace(is_available=lambda: False)
        )
        with patch.dict("sys.modules", {"torch": fake_torch}):
            result = cell_detection._detect_gpu(requested=True)
        assert result is False


# ---------------------------------------------------------------------------
# TestLoadImage
# ---------------------------------------------------------------------------


class TestLoadImage:
    """Tests for image loading from disk."""

    def test_load_greyscale_tiff(self, tmp_path):
        import tifffile
        img = np.random.randint(0, 255, (64, 64), dtype=np.uint8)
        path = tmp_path / "grey.tif"
        tifffile.imwrite(str(path), img)
        arr, n = cell_detection.load_image(str(path))
        assert arr.shape == (64, 64)
        assert n == 1

    def test_load_rgb_tiff(self, tmp_path):
        import tifffile
        img = np.random.randint(0, 255, (32, 32, 3), dtype=np.uint8)
        path = tmp_path / "rgb.tif"
        tifffile.imwrite(str(path), img)
        arr, n = cell_detection.load_image(str(path))
        assert arr.shape == (32, 32, 3)
        assert n == 3

    def test_load_png(self, tmp_path):
        from PIL import Image as PILImage
        img = PILImage.fromarray(np.zeros((16, 16), dtype=np.uint8))
        path = tmp_path / "img.png"
        img.save(str(path))
        arr, n = cell_detection.load_image(str(path))
        assert n == 1


# ---------------------------------------------------------------------------
# TestDemoImageEdgeCells
# ---------------------------------------------------------------------------


class TestDemoImageEdgeCells:
    """Verify that the demo image contains cells cut off at the border."""

    def test_border_pixels_lit(self):
        """At least one border pixel must be non-zero (edge nuclei touch the frame)."""
        img = cell_detection.make_demo_image()
        h, w = img.shape
        border = np.concatenate([
            img[0, :],     # top row
            img[h - 1, :], # bottom row
            img[:, 0],     # left col
            img[:, w - 1], # right col
        ])
        assert border.max() > 0, "Expected edge-touching nuclei but no lit border pixels found"

    def test_interior_cells_still_present(self):
        """Interior of the image should also have bright pixels (interior cells remain)."""
        img = cell_detection.make_demo_image()
        interior = img[20:-20, 20:-20]
        assert interior.max() > 100


# ---------------------------------------------------------------------------
# TestExcludeOnEdges
# ---------------------------------------------------------------------------


class TestExcludeOnEdges:
    """Tests for the --exclude_on_edges post-processing step."""

    def _masks_with_edge_cell(self) -> np.ndarray:
        """64×64 mask: cell 1 interior, cell 2 touches top border."""
        masks = np.zeros((64, 64), dtype=np.uint16)
        masks[20:30, 20:30] = 1   # fully interior
        masks[0:5, 30:40] = 2     # touches top edge
        return masks

    def test_remove_edge_masks_removes_border_cell(self):
        from cellpose import utils as cp_utils
        masks = self._masks_with_edge_cell()
        cleaned = cp_utils.remove_edge_masks(masks)
        assert 2 not in cleaned, "Edge cell (label 2) should have been removed"

    def test_remove_edge_masks_keeps_interior_cell(self):
        from cellpose import utils as cp_utils
        masks = self._masks_with_edge_cell()
        cleaned = cp_utils.remove_edge_masks(masks)
        assert 1 in cleaned, "Interior cell (label 1) should be retained"

    def test_report_shows_exclude_on_edges_yes(self, tmp_path):
        metrics = [{"id": 1, "area": 200, "diameter": 16.0, "centroid_x": 10.0, "centroid_y": 10.0, "eccentricity": 0.2}]
        cell_detection.write_report(metrics, {"exclude_on_edges": True}, tmp_path)
        text = (tmp_path / "report.md").read_text()
        assert "Exclude edge cells:** yes" in text

    def test_report_shows_exclude_on_edges_no(self, tmp_path):
        metrics = [{"id": 1, "area": 200, "diameter": 16.0, "centroid_x": 10.0, "centroid_y": 10.0, "eccentricity": 0.2}]
        cell_detection.write_report(metrics, {"exclude_on_edges": False}, tmp_path)
        text = (tmp_path / "report.md").read_text()
        assert "Exclude edge cells:** no" in text


# ---------------------------------------------------------------------------
# TestFlowCellprobThresholds
# ---------------------------------------------------------------------------


class TestFlowCellprobThresholds:
    """Tests for --flow_threshold and --cellprob_threshold parameters."""

    def _metrics(self):
        return [{"id": 1, "area": 200, "diameter": 16.0, "centroid_x": 10.0, "centroid_y": 10.0, "eccentricity": 0.2}]

    # Report rendering

    def test_report_shows_flow_threshold(self, tmp_path):
        cell_detection.write_report(self._metrics(), {"flow_threshold": 0.4}, tmp_path)
        text = (tmp_path / "report.md").read_text()
        assert "0.4" in text

    def test_report_shows_cellprob_threshold(self, tmp_path):
        cell_detection.write_report(self._metrics(), {"cellprob_threshold": 0.0}, tmp_path)
        text = (tmp_path / "report.md").read_text()
        assert "cellprob" in text.lower()

    def test_report_shows_custom_flow_threshold(self, tmp_path):
        cell_detection.write_report(self._metrics(), {"flow_threshold": 0.8}, tmp_path)
        text = (tmp_path / "report.md").read_text()
        assert "0.8" in text

    def test_report_shows_custom_cellprob_threshold(self, tmp_path):
        cell_detection.write_report(self._metrics(), {"cellprob_threshold": -2.0}, tmp_path)
        text = (tmp_path / "report.md").read_text()
        assert "-2.0" in text

    # run_segmentation passes params to model.eval

    def test_run_segmentation_passes_flow_threshold(self):
        import types
        calls = {}

        def fake_eval(img, diameter=None, do_3D=False, z_axis=None,
                      flow_threshold=0.4, cellprob_threshold=0.0):
            calls["flow_threshold"] = flow_threshold
            calls["cellprob_threshold"] = cellprob_threshold
            h, w = img.shape[:2]
            return np.zeros((h, w), dtype=np.uint16), [[]], None

        fake_model = types.SimpleNamespace(eval=fake_eval)

        with patch("cellpose.models.CellposeModel", return_value=fake_model):
            img = np.zeros((32, 32), dtype=np.uint8)
            cell_detection.run_segmentation(img, diameter=None, use_gpu=False,
                                            flow_threshold=0.9, cellprob_threshold=0.0)
        assert calls["flow_threshold"] == 0.9

    def test_run_segmentation_passes_cellprob_threshold(self):
        import types
        calls = {}

        def fake_eval(img, diameter=None, do_3D=False, z_axis=None,
                      flow_threshold=0.4, cellprob_threshold=0.0):
            calls["flow_threshold"] = flow_threshold
            calls["cellprob_threshold"] = cellprob_threshold
            h, w = img.shape[:2]
            return np.zeros((h, w), dtype=np.uint16), [[]], None

        fake_model = types.SimpleNamespace(eval=fake_eval)

        with patch("cellpose.models.CellposeModel", return_value=fake_model):
            img = np.zeros((32, 32), dtype=np.uint8)
            cell_detection.run_segmentation(img, diameter=None, use_gpu=False,
                                            flow_threshold=0.4, cellprob_threshold=-3.0)
        assert calls["cellprob_threshold"] == -3.0

    # Argparse defaults

    def test_argparse_flow_threshold_default(self):
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("--flow_threshold", type=float, default=0.4)
        parser.add_argument("--cellprob_threshold", type=float, default=0.0)
        args = parser.parse_args([])
        assert args.flow_threshold == 0.4
        assert args.cellprob_threshold == 0.0
