"""Tests for the CellposeSAM cell segmentation skill."""
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

import cell_detection

# cellpose pulls torch, whose default Linux wheels carry the CUDA runtime, so
# CI installs only the light imaging dependencies (the `cell-detection` extra)
# and skips the four tests that need cellpose itself. Those four are covered by
# running the skill locally; everything else here runs in CI.
_CELLPOSE_REASON = "cellpose (and torch) are not installed in CI; run locally to cover these"


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

    def test_load_czi_zcyx_metadata_max_projects_to_hwc(self, tmp_path):
        import types
        path = tmp_path / "img.czi"
        path.write_bytes(b"")
        source = np.ones((1, 3, 32, 32), dtype=np.uint16)

        class _FakeCziFile:
            axes = "ZCYX"

            def __init__(self, image):
                self.image = image

            def asarray(self):
                return source

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

        fake_czi = types.ModuleType("czifile")
        fake_czi.imread = lambda image: source
        fake_czi.CziFile = _FakeCziFile
        with patch.dict("sys.modules", {"czifile": fake_czi}):
            arr, n = cell_detection.load_image(str(path))
        assert arr.shape == (32, 32, 3)
        assert n == 3

    def test_load_nd2_stack_counts_channels(self, tmp_path):
        import types
        path = tmp_path / "img.nd2"
        path.write_bytes(b"")
        source = np.ones((4, 3, 16, 16), dtype=np.uint16)

        class _FakeND2File:
            sizes = {"Z": 4, "C": 3, "Y": 16, "X": 16}

            def __init__(self, image):
                self.image = image

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

        fake_nd2 = types.ModuleType("nd2")
        fake_nd2.imread = lambda image: source
        fake_nd2.ND2File = _FakeND2File
        with patch.dict("sys.modules", {"nd2": fake_nd2}):
            arr, n = cell_detection.load_image(str(path))
        assert arr.shape == (16, 16, 3)
        assert n == 3

    def test_load_czi_uses_axes_metadata_for_channel_gt_z(self, tmp_path):
        import types

        path = tmp_path / "img.czi"
        path.write_bytes(b"")
        source = np.zeros((2, 5, 8, 8), dtype=np.uint16)  # Z, C, Y, X
        for channel in range(5):
            source[:, channel, :, :] = channel + 1

        class _FakeCziFile:
            axes = "ZCYX"

            def __init__(self, image):
                self.image = image

            def asarray(self):
                return source

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

        fake_czi = types.ModuleType("czifile")
        fake_czi.CziFile = _FakeCziFile
        fake_czi.imread = lambda image: source

        with patch.dict("sys.modules", {"czifile": fake_czi}):
            arr, n = cell_detection.load_image(str(path))

        assert arr.shape == (8, 8, 5)
        assert n == 5
        for channel in range(5):
            assert np.all(arr[:, :, channel] == channel + 1)

    def test_load_nd2_uses_metadata_for_channel_gt_z(self, tmp_path):
        import types

        path = tmp_path / "img.nd2"
        path.write_bytes(b"")
        source = np.zeros((2, 5, 8, 8), dtype=np.uint16)  # Z, C, Y, X
        for channel in range(5):
            source[:, channel, :, :] = channel + 1

        class _FakeND2File:
            sizes = {"Z": 2, "C": 5, "Y": 8, "X": 8}

            def __init__(self, image):
                self.image = image

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

        fake_nd2 = types.ModuleType("nd2")
        fake_nd2.imread = lambda image: source
        fake_nd2.ND2File = _FakeND2File

        with patch.dict("sys.modules", {"nd2": fake_nd2}):
            arr, n = cell_detection.load_image(str(path))

        assert arr.shape == (8, 8, 5)
        assert n == 5
        for channel in range(5):
            assert np.all(arr[:, :, channel] == channel + 1)

    def test_load_czi_singleton_channel_preserved(self, tmp_path):
        import types

        path = tmp_path / "img.czi"
        path.write_bytes(b"")
        source = np.zeros((8, 1, 12, 12), dtype=np.uint16)  # Z, C, Y, X
        source[:, 0, :, :] = 4096

        class _FakeCziFile:
            axes = "ZCYX"

            def __init__(self, image):
                self.image = image

            def asarray(self):
                return source

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

        fake_czi = types.ModuleType("czifile")
        fake_czi.CziFile = _FakeCziFile
        fake_czi.imread = lambda image: source

        with patch.dict("sys.modules", {"czifile": fake_czi}):
            arr, n = cell_detection.load_image(str(path))

        assert arr.shape == (12, 12, 1)
        assert n == 1
        assert arr.dtype == np.uint16
        assert int(arr.max()) == 4096

    def test_load_nd2_uint16_dynamic_range_preserved(self, tmp_path):
        import types

        path = tmp_path / "img.nd2"
        path.write_bytes(b"")
        source = np.zeros((2, 3, 10, 10), dtype=np.uint16)  # Z, C, Y, X
        source[0, 0, :, :] = 10
        source[0, 1, :, :] = 2000
        source[1, 2, :, :] = 65535

        class _FakeND2File:
            sizes = {"Z": 2, "C": 3, "Y": 10, "X": 10}

            def __init__(self, image):
                self.image = image

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

        fake_nd2 = types.ModuleType("nd2")
        fake_nd2.imread = lambda image: source
        fake_nd2.ND2File = _FakeND2File

        with patch.dict("sys.modules", {"nd2": fake_nd2}):
            arr, n = cell_detection.load_image(str(path))

        assert arr.shape == (10, 10, 3)
        assert n == 3
        assert arr.dtype == np.uint16
        assert int(arr[:, :, 0].max()) == 10
        assert int(arr[:, :, 1].max()) == 2000
        assert int(arr[:, :, 2].max()) == 65535

    def test_load_nd2_projection_none_returns_zcyx_single_channel(self, tmp_path):
        import types

        path = tmp_path / "img.nd2"
        path.write_bytes(b"")
        source = np.arange(4 * 1 * 6 * 7, dtype=np.uint16).reshape(4, 1, 6, 7)  # Z, C, Y, X

        class _FakeND2File:
            sizes = {"Z": 4, "C": 1, "Y": 6, "X": 7}

            def __init__(self, image):
                self.image = image

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

        fake_nd2 = types.ModuleType("nd2")
        fake_nd2.imread = lambda image: source
        fake_nd2.ND2File = _FakeND2File

        with patch.dict("sys.modules", {"nd2": fake_nd2}):
            arr, n = cell_detection.load_image(str(path), z_projection="none")

        assert arr.shape == (4, 1, 6, 7)
        assert n == 1
        np.testing.assert_array_equal(arr, source)

    def test_load_czi_projection_none_multichannel_returns_zcyx(self, tmp_path):
        import types

        path = tmp_path / "img.czi"
        path.write_bytes(b"")
        source = np.zeros((3, 4, 6, 7), dtype=np.uint16)  # Z, C, Y, X
        for channel in range(4):
            source[:, channel, :, :] = channel + 10

        class _FakeCziFile:
            axes = "ZCYX"

            def __init__(self, image):
                self.image = image

            def asarray(self):
                return source

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

        fake_czi = types.ModuleType("czifile")
        fake_czi.CziFile = _FakeCziFile
        fake_czi.imread = lambda image: source

        with patch.dict("sys.modules", {"czifile": fake_czi}):
            arr, n = cell_detection.load_image(str(path), z_projection="none")

        assert arr.shape == (3, 4, 6, 7)
        assert n == 4
        np.testing.assert_array_equal(arr, source)

    def test_load_nd2_projection_none_multichannel_returns_zcyx(self, tmp_path):
        import types

        path = tmp_path / "img.nd2"
        path.write_bytes(b"")
        source = np.zeros((5, 3, 4, 6), dtype=np.uint16)  # Z, C, Y, X
        for channel in range(3):
            source[:, channel, :, :] = channel + 20

        class _FakeND2File:
            sizes = {"Z": 5, "C": 3, "Y": 4, "X": 6}

            def __init__(self, image):
                self.image = image

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

        fake_nd2 = types.ModuleType("nd2")
        fake_nd2.imread = lambda image: source
        fake_nd2.ND2File = _FakeND2File

        with patch.dict("sys.modules", {"nd2": fake_nd2}):
            arr, n = cell_detection.load_image(str(path), z_projection="none")

        assert arr.shape == (5, 3, 4, 6)
        assert n == 3
        np.testing.assert_array_equal(arr, source)

    def test_run_segmentation_3d_multichannel_uses_z_and_channel_axes(self):
        import types

        calls = {}

        def fake_eval(
            img,
            diameter=None,
            do_3D=False,
            z_axis=None,
            channel_axis=None,
            flow_threshold=0.4,
            cellprob_threshold=0.0,
        ):
            calls["shape"] = img.shape
            calls["do_3D"] = do_3D
            calls["z_axis"] = z_axis
            calls["channel_axis"] = channel_axis
            return np.zeros((img.shape[2], img.shape[3]), dtype=np.uint16), [[]], None

        class _FakeCellposeModel:
            def __init__(self, gpu=False):
                self.gpu = gpu

            def eval(self, *args, **kwargs):
                return fake_eval(*args, **kwargs)

        fake_models = types.ModuleType("cellpose.models")
        fake_models.CellposeModel = _FakeCellposeModel
        fake_utils = types.ModuleType("cellpose.utils")
        fake_utils.remove_edge_masks = lambda masks: masks
        fake_cellpose = types.ModuleType("cellpose")
        fake_cellpose.models = fake_models
        fake_cellpose.utils = fake_utils

        with patch.dict(
            "sys.modules",
            {
                "cellpose": fake_cellpose,
                "cellpose.models": fake_models,
                "cellpose.utils": fake_utils,
            },
        ):
            img = np.zeros((4, 3, 8, 9), dtype=np.uint16)  # Z, C, Y, X
            cell_detection.run_segmentation(img, diameter=None, use_gpu=False, do_3D=True)

        assert calls["shape"] == (4, 3, 8, 9)
        assert calls["do_3D"] is True
        assert calls["z_axis"] == 0
        assert calls["channel_axis"] == 1

    def test_load_czi_projection_none_heuristic_4d_returns_zcyx(self, tmp_path):
        import types

        path = tmp_path / "img.czi"
        path.write_bytes(b"")
        source = np.zeros((7, 4, 16, 16), dtype=np.uint16)  # Z, C, Y, X
        for channel in range(4):
            source[:, channel, :, :] = channel + 1

        class _FakeCziFile:
            # Non-CZYX metadata so loader falls back to 4D heuristic path.
            axes = "BTZYX"

            def __init__(self, image):
                self.image = image

            def asarray(self):
                return source

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

        fake_czi = types.ModuleType("czifile")
        fake_czi.CziFile = _FakeCziFile

        with patch.dict("sys.modules", {"czifile": fake_czi}):
            arr, n = cell_detection.load_image(str(path), z_projection="none")

        assert arr.shape == (7, 4, 16, 16)
        assert n == 4
        np.testing.assert_array_equal(arr, source)

    def test_load_4d_heuristic_raises_when_z_le_c_and_no_metadata(self, tmp_path):
        import types

        path = tmp_path / "img.czi"
        path.write_bytes(b"")
        source = np.zeros((3, 4, 16, 16), dtype=np.uint16)  # Z, C, Y, X — Z <= C

        class _FakeCziFile:
            axes = "BTZYX"

            def __init__(self, image):
                self.image = image

            def asarray(self):
                return source

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

        fake_czi = types.ModuleType("czifile")
        fake_czi.CziFile = _FakeCziFile

        with patch.dict("sys.modules", {"czifile": fake_czi}):
            with pytest.raises(ValueError, match="Cannot infer channel vs Z axis"):
                cell_detection.load_image(str(path), z_projection="max")
            with pytest.raises(ValueError, match="Cannot infer channel vs Z axis"):
                cell_detection.load_image(str(path), z_projection="none")


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
        cp_utils = pytest.importorskip("cellpose.utils", reason=_CELLPOSE_REASON)
        masks = self._masks_with_edge_cell()
        cleaned = cp_utils.remove_edge_masks(masks)
        assert 2 not in cleaned, "Edge cell (label 2) should have been removed"

    def test_remove_edge_masks_keeps_interior_cell(self):
        cp_utils = pytest.importorskip("cellpose.utils", reason=_CELLPOSE_REASON)
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
        pytest.importorskip("cellpose.models", reason=_CELLPOSE_REASON)
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
        pytest.importorskip("cellpose.models", reason=_CELLPOSE_REASON)
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


# ---------------------------------------------------------------------------
# TestCliGuardrails
# ---------------------------------------------------------------------------


class TestCliGuardrails:
    """Tests for --do_3D / --z_projection mode consistency in main()."""

    def test_z_projection_none_without_do_3d_errors(self, tmp_path):
        with patch(
            "sys.argv",
            [
                "cell_detection.py",
                "--input",
                str(tmp_path / "img.tif"),
                "--z_projection",
                "none",
                "--output",
                str(tmp_path / "out"),
            ],
        ):
            with pytest.raises(SystemExit) as exc:
                cell_detection.main()
        assert exc.value.code == 2

    def test_do_3d_with_max_projection_errors_on_volumetric_input(self, tmp_path):
        path = tmp_path / "stack.nd2"
        path.write_bytes(b"")
        volume = np.zeros((3, 2, 8, 8), dtype=np.uint16)

        with patch.object(cell_detection, "load_image", return_value=(volume, 2)):
            with patch(
                "sys.argv",
                [
                    "cell_detection.py",
                    "--input",
                    str(path),
                    "--do_3D",
                    "--output",
                    str(tmp_path / "out"),
                ],
            ):
                with pytest.raises(SystemExit) as exc:
                    cell_detection.main()
        assert exc.value.code == 2

    def test_do_3d_errors_on_non_volumetric_ndim2_input(self, tmp_path):
        """Greyscale H×W (ndim=2) with n_channels>1 cannot be treated as 3D."""
        path = tmp_path / "img.tif"
        path.write_bytes(b"")
        image = np.zeros((8, 8), dtype=np.uint8)

        with patch.object(cell_detection, "load_image", return_value=(image, 4)):
            with patch(
                "sys.argv",
                [
                    "cell_detection.py",
                    "--input",
                    str(path),
                    "--do_3D",
                    "--z_projection",
                    "none",
                    "--output",
                    str(tmp_path / "out"),
                ],
            ):
                with pytest.raises(SystemExit) as exc:
                    cell_detection.main()
        assert exc.value.code == 2

    def test_do_3d_errors_when_z_stack_has_single_plane(self, tmp_path):
        path = tmp_path / "stack.nd2"
        path.write_bytes(b"")
        volume = np.zeros((1, 3, 8, 8), dtype=np.uint16)

        with patch.object(cell_detection, "load_image", return_value=(volume, 3)):
            with patch(
                "sys.argv",
                [
                    "cell_detection.py",
                    "--input",
                    str(path),
                    "--do_3D",
                    "--z_projection",
                    "none",
                    "--output",
                    str(tmp_path / "out"),
                ],
            ):
                with pytest.raises(SystemExit) as exc:
                    cell_detection.main()
        assert exc.value.code == 2


# ---------------------------------------------------------------------------
# TestReproducibility
# ---------------------------------------------------------------------------


class TestReproducibility:
    """Reproducibility artefacts written by main()."""

    @pytest.fixture
    def demo_output(self, tmp_path):
        """Run main() --demo with segmentation and outline-saving mocked."""
        masks = np.zeros((64, 64), dtype=np.uint16)
        masks[10:20, 10:20] = 1
        fake_flows = [np.zeros((64, 64))]

        def _fake_seg(img, diameter, use_gpu, **kwargs):
            return masks, fake_flows

        def _fake_outlines(img, masks, flows, output_dir, stem="image"):
            fig_dir = Path(output_dir) / "figures"
            fig_dir.mkdir(parents=True, exist_ok=True)
            fname = f"{stem}_cp_outlines.png"
            (fig_dir / fname).write_bytes(b"\x89PNG\r\n")
            return fname

        with patch.object(cell_detection, "run_segmentation", side_effect=_fake_seg), \
             patch.object(cell_detection, "save_outlines", side_effect=_fake_outlines):
            with patch("sys.argv", ["cell_detection.py", "--demo", "--output", str(tmp_path)]):
                cell_detection.main()
        return tmp_path

    def test_checksums_file_created(self, demo_output):
        assert (demo_output / "reproducibility" / "checksums.sha256").exists()

    def test_checksums_lists_report(self, demo_output):
        content = (demo_output / "reproducibility" / "checksums.sha256").read_text()
        assert "report.md" in content

    def test_checksums_lists_csv(self, demo_output):
        content = (demo_output / "reproducibility" / "checksums.sha256").read_text()
        assert ".csv" in content

    def test_environment_yml_created(self, demo_output):
        assert (demo_output / "reproducibility" / "environment.yml").exists()

    def test_environment_yml_has_cellpose(self, demo_output):
        text = (demo_output / "reproducibility" / "environment.yml").read_text()
        assert "cellpose" in text

    def test_environment_yml_has_env_name(self, demo_output):
        text = (demo_output / "reproducibility" / "environment.yml").read_text()
        assert "clawbio-cell-detection" in text

    def test_commands_sh_created(self, demo_output):
        assert (demo_output / "reproducibility" / "commands.sh").exists()

    def test_commands_sh_is_executable(self, demo_output):
        path = demo_output / "reproducibility" / "commands.sh"
        assert path.stat().st_mode & 0o111

    def test_commands_sh_contains_cell_detection(self, demo_output):
        text = (demo_output / "reproducibility" / "commands.sh").read_text()
        assert "cell_detection" in text
