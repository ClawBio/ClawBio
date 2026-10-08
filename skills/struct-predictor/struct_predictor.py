#!/usr/bin/env python3
"""
Struct Predictor — Protein structure prediction with Boltz-2 (default) or OpenFold3.

Usage:
    # Single protein or multi-chain complex (YAML)
    python skills/struct-predictor/struct_predictor.py \
        --input complex.yaml --output /tmp/struct_out

    # Same input, OpenFold3 instead of Boltz-2 (needs a CUDA GPU)
    python skills/struct-predictor/struct_predictor.py \
        --input complex.yaml --output /tmp/struct_out --backend openfold3

    # Demo (Trp-cage miniprotein, no input needed)
    python skills/struct-predictor/struct_predictor.py \
        --demo --output /tmp/struct_demo
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

_SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SKILL_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

from struct_predictor_core.io import validate_and_prepare, write_openfold3_query
from struct_predictor_core.predict import run_boltz, run_openfold3
from struct_predictor_core.confidence import extract_confidence
from struct_predictor_core.report import generate_report

# ---------------------------------------------------------------------------
# Demo data — Trp-cage miniprotein (PDB: 1L2Y)
# ---------------------------------------------------------------------------

_DEMO_YAML = _SKILL_DIR / "demo_data" / "trpcage.yaml"
DEMO_NAME  = "Trpcage"


# ---------------------------------------------------------------------------
# Top-level pipeline
# ---------------------------------------------------------------------------


def run_struct_prediction(
    input_path: Path | None,
    output_dir: Path,
    demo: bool = False,
    backend: str = "boltz",
) -> dict:
    """Run the full struct-predictor pipeline (fully offline, no MSA server).

    Args:
        input_path: Path to YAML input. Required unless demo=True.
        output_dir: Where to write the final report and artefacts.
        demo: Run with Trp-cage miniprotein (PDB 1L2Y, 20 residues).
        backend: "boltz" (default) or "openfold3".

    Returns:
        result dict (same content as result.json).

    Raises:
        ValueError: If neither input_path nor demo is supplied.
    """
    if not demo and input_path is None:
        raise ValueError("Provide --input or --demo.")

    if backend not in ("openfold3", "boltz"):
        raise ValueError(f"Unknown backend {backend!r}. Use 'openfold3' or 'boltz'.")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Validate input and prepare it in a temp dir. The Boltz YAML
    # (msa: empty, offline) is written either way; OpenFold3 reads the
    # validated sequences from it and gets its own query JSON.
    with tempfile.TemporaryDirectory(prefix="struct_predictor_input_") as _tmpdir:
        work_dir = Path(_tmpdir)

        if demo:
            print(f"  Demo mode: {DEMO_NAME} (20 residues, PDB 1L2Y)")
            prepared = validate_and_prepare(_DEMO_YAML, work_dir / "boltz_input")
            input_label = "Trp-cage miniprotein (demo)"
        else:
            print(f"  Input: {input_path}")
            prepared = validate_and_prepare(Path(input_path), work_dir / "boltz_input")
            input_label = str(input_path)

        print(f"  Chains: {len(prepared['sequences'])}, type: {prepared['input_type']}")

        # Step 2: the backend writes directly to output_dir, preserving its native layout:
        #   boltz:     output_dir/predictions/<name>/<name>_model_0.cif
        #   openfold3: output_dir/<name>/seed_<n>/<name>_seed_<n>_sample_<k>_model.cif
        if backend == "openfold3":
            name = DEMO_NAME if demo else Path(input_path).stem
            query = write_openfold3_query(prepared["sequences"], name, work_dir / "openfold3_input")
            print("  Running OpenFold3 prediction...")
            predict_result = run_openfold3(query, output_dir, name)
        else:
            print("  Running Boltz-2 prediction...")
            predict_result = run_boltz(
                input_path=prepared["boltz_input_path"],
                boltz_output_dir=output_dir,
            )

    cif_path = predict_result["cif_path"]
    conf_json_path = predict_result["confidence_json_path"]

    cmd = _build_cmd(input_label=input_label, output_dir=output_dir, demo=demo, backend=backend)

    # Step 3: Extract confidence
    print("  Extracting pLDDT and PAE...")
    conf = extract_confidence(cif_path, conf_json_path)
    plddt = conf["plddt"]
    pae   = conf["pae"]
    chain_boundaries = conf["chain_boundaries"]

    print(f"  Mean pLDDT: {plddt.mean():.1f}")
    print(f"  Residues: {len(plddt)}")

    # Step 4: Write ClawBio artifacts alongside Boltz output at output_dir root
    print(f"  Writing report to {output_dir}/")
    generate_report(
        output_dir=output_dir,
        sequences_info=prepared["sequences"],
        plddt=plddt,
        pae=pae,
        chain_boundaries=chain_boundaries,
        cif_path=cif_path,
        cmd=cmd,
        input_label=input_label,
        demo=demo,
        backend=backend,
    )

    result = json.loads((output_dir / "result.json").read_text())
    print(f"\n  Report: {output_dir / 'report.md'}")
    print(f"  Full output: {output_dir}/")
    print(
        "\n  ClawBio is a research and educational tool. It is not a medical device "
        "and does not provide clinical diagnoses. Consult a healthcare professional "
        "before making any medical decisions."
    )
    return result


def _build_cmd(input_label: str, output_dir: Path, demo: bool, backend: str) -> str:
    parts = ["python skills/struct-predictor/struct_predictor.py"]
    if demo:
        parts.append("--demo")
    else:
        parts.append(f"--input {input_label}")
    parts.append(f"--output {output_dir}")
    parts.append(f"--backend {backend}")
    return " \\\n  ".join(parts)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Struct Predictor — protein structure prediction with Boltz-2 or OpenFold3"
    )
    parser.add_argument("--input", "-i", help="Input YAML file")
    parser.add_argument("--output", "-o", required=True, help="Output directory")
    parser.add_argument("--backend", choices=["boltz", "openfold3"], default="boltz",
                        help="Prediction engine (default: boltz)")
    parser.add_argument("--demo", action="store_true",
                        help="Run demo with Trp-cage miniprotein, PDB 1L2Y (no input needed)")
    return parser


def main():
    parser = _build_parser()
    args = parser.parse_args()

    if not args.input and not args.demo:
        parser.print_help()
        print("\nError: provide --input or --demo")
        sys.exit(1)

    print(f"Struct Predictor — {'OpenFold3' if args.backend == 'openfold3' else 'Boltz-2'}")
    print("=" * 60)
    print()

    run_struct_prediction(
        input_path=Path(args.input) if args.input else None,
        output_dir=Path(args.output),
        demo=args.demo,
        backend=args.backend,
    )


if __name__ == "__main__":
    main()
