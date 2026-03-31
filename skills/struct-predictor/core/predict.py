"""
predict.py — subprocess wrapper for `boltz predict`.

Calls Boltz-2 CLI, streams output, and locates the resulting CIF and
confidence JSON in the output directory.
"""
from __future__ import annotations

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
