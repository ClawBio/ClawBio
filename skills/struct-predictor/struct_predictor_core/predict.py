"""
predict.py — subprocess wrappers for `run_openfold predict` and `boltz predict`.

Calls the backend CLI, streams output, and locates the resulting CIF and
confidence JSON in the output directory.
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path


def run_boltz(
    input_path: Path,
    boltz_output_dir: Path,
) -> dict:
    """Run `boltz predict` and return paths to outputs.

    Runs fully offline — no MSA server is used so no data leaves the machine.
    The input YAML must have ``msa: empty`` injected for each protein chain
    (handled by ``validate_and_prepare``).

    Args:
        input_path: Path to the Boltz input YAML file.
        boltz_output_dir: Directory where Boltz writes its results.

    Returns:
        {
            "cif_path": Path,
            "confidence_json_path": Path | None,
            "boltz_output_dir": Path,
        }

    Raises:
        RuntimeError: If boltz exits with a non-zero return code.
        FileNotFoundError: If no CIF is found in the output after a successful run.
    """
    boltz_output_dir = Path(boltz_output_dir)
    boltz_output_dir.mkdir(parents=True, exist_ok=True)

    cmd = _build_boltz_cmd(input_path, boltz_output_dir)
    print(f"  Running: {' '.join(str(c) for c in cmd)}")

    proc = subprocess.run(
        cmd,
        capture_output=False,   # let stdout/stderr stream to terminal
        text=True,
    )

    if proc.returncode != 0:
        stderr = getattr(proc, "stderr", "") or ""
        raise RuntimeError(
            f"Boltz exited with code {proc.returncode}.\n"
            f"stderr: {stderr.strip()}"
        )

    return _find_cif(boltz_output_dir)


def _build_boltz_cmd(
    input_path: Path,
    boltz_output_dir: Path,
) -> list[str]:
    """Build the boltz predict command list (fully offline, no MSA server)."""
    return [
        "boltz", "predict",
        str(input_path),
        "--out_dir", str(boltz_output_dir),
    ]


def _find_cif(boltz_output_dir: Path) -> dict:
    """Locate outputs written by Boltz-2.

    Boltz-2 writes to:
        <boltz_output_dir>/predictions/<name>/<name>_model_0.cif
        <boltz_output_dir>/predictions/<name>/confidence_<name>_model_0.json
    """
    boltz_output_dir = Path(boltz_output_dir)
    cifs = list(boltz_output_dir.rglob("*_model_0.cif"))
    if not cifs:
        raise FileNotFoundError(
            f"No CIF file found under {boltz_output_dir}. "
            "Boltz may not have produced output — check the logs above."
        )

    cif_path = cifs[0]
    pred_dir = cif_path.parent

    conf_candidates = list(pred_dir.glob("confidence_*_model_0.json"))

    return {
        "cif_path": cif_path,
        "confidence_json_path": conf_candidates[0] if conf_candidates else None,
        "boltz_output_dir": boltz_output_dir,
    }


def run_openfold3(query_json_path: Path, output_dir: Path) -> dict:
    """Run `run_openfold predict` and return paths to the best-ranked sample.

    Fully offline: no MSA server, no templates. Needs a CUDA GPU.

    Returns:
        {"cif_path": Path, "confidence_json_path": Path | None, "output_dir": Path}

    Raises:
        RuntimeError: If openfold3 is not installed or exits non-zero.
        FileNotFoundError: If no CIF is found after a successful run.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = _build_openfold3_cmd(query_json_path, output_dir)
    print(f"  Running: {' '.join(str(c) for c in cmd)}")

    try:
        proc = subprocess.run(cmd, capture_output=False, text=True)
    except FileNotFoundError:
        raise RuntimeError(
            "run_openfold not found. Install with: pip install openfold3 "
            "&& setup_openfold  (or omit --backend to use Boltz-2)"
        )

    if proc.returncode != 0:
        raise RuntimeError(f"OpenFold3 exited with code {proc.returncode}.")

    return _find_openfold3_output(output_dir)


def _build_openfold3_cmd(query_json_path: Path, output_dir: Path) -> list[str]:
    """Build the run_openfold predict command (fully offline)."""
    return [
        "run_openfold", "predict",
        f"--query-json={query_json_path}",
        f"--output-dir={output_dir}",
        "--use-msa-server=false",
        "--use-templates=false",
    ]


def _find_openfold3_output(output_dir: Path) -> dict:
    """Locate the best-ranked sample written by OpenFold3.

    OpenFold3 writes, per sample:
        <out>/<name>/seed_<n>/<name>_seed_<n>_sample_<k>_model.cif
        <out>/<name>/seed_<n>/<name>_seed_<n>_sample_<k>_confidences.json
        <out>/<name>/seed_<n>/<name>_seed_<n>_sample_<k>_confidences_aggregated.json
    The sample with the highest ``sample_ranking_score`` is returned.
    """
    output_dir = Path(output_dir)
    cifs = list(output_dir.rglob("*_model.cif"))
    if not cifs:
        raise FileNotFoundError(
            f"No CIF file found under {output_dir}. "
            "OpenFold3 may not have produced output — check the logs above."
        )

    def score(cif: Path) -> float:
        agg = cif.with_name(cif.name.replace("_model.cif", "_confidences_aggregated.json"))
        try:
            return float(json.loads(agg.read_text())["sample_ranking_score"])
        except (OSError, KeyError, ValueError):
            return float("-inf")

    best = max(cifs, key=score)
    conf = best.with_name(best.name.replace("_model.cif", "_confidences.json"))
    return {
        "cif_path": best,
        "confidence_json_path": conf if conf.exists() else None,
        "output_dir": output_dir,
    }
