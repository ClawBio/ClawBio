#!/usr/bin/env python3
"""Generate the bundled toy repurposing screen dataset.

Synthetic but PRISM-shaped: 10 samples × 20 compounds, two contexts
(context_A target, context_B reference), three context-selective hits.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
DEMO_DIR = SCRIPT_DIR / "demo"

RNG = np.random.default_rng(42)
SAMPLES = [f"SAMPLE_{i:02d}" for i in range(1, 11)]
CONTEXT_A = SAMPLES[:5]
CONTEXT_B = SAMPLES[5:]
COMPOUNDS = [f"BRD-{i:04d}" for i in range(1, 21)]
SELECTIVE = {"BRD-0003", "BRD-0007", "BRD-0015"}
PLATES = ["PLATE_01", "PLATE_02"]


def _wells(n: int = 96) -> list[str]:
    out = []
    for row in "ABCDEFGH":
        for col in range(1, 13):
            out.append(f"{row}{col:02d}")
            if len(out) >= n:
                return out
    return out


def build_treatment_info() -> pd.DataFrame:
    rows = []
    wells = _wells(120)
    wi = 0
    for plate in PLATES:
        for _ in range(6):
            rows.append({
                "well_id": wells[wi], "detection_plate": plate,
                "perturbation_type": "vehicle_control", "compound_id": "DMSO",
                "compound_name": "DMSO", "dose_um": 0.0,
                "clinical_phase": None, "moa": "vehicle",
            })
            wi += 1
        for _ in range(6):
            rows.append({
                "well_id": wells[wi], "detection_plate": plate,
                "perturbation_type": "positive_control", "compound_id": "POS",
                "compound_name": "Bortezomib", "dose_um": 10.0,
                "clinical_phase": "Launched", "moa": "proteasome inhibitor",
            })
            wi += 1
        for compound_id in COMPOUNDS:
            name = f"Drug_{compound_id.split('-')[1]}"
            phase = "Launched" if int(compound_id.split("-")[1]) % 3 == 0 else "Phase 2"
            rows.append({
                "well_id": wells[wi], "detection_plate": plate,
                "perturbation_type": "experimental_treatment",
                "compound_id": compound_id, "compound_name": name,
                "dose_um": 2.5, "clinical_phase": phase,
                "moa": "kinase inhibitor" if compound_id in SELECTIVE else "unknown",
            })
            wi += 1
    return pd.DataFrame(rows)


def build_sample_info() -> pd.DataFrame:
    rows = []
    for sid in SAMPLES:
        context = "context_A" if sid in CONTEXT_A else "context_B"
        disease = "IBD" if sid in CONTEXT_A else "fibrosis"
        row = {
            "sample_id": sid,
            "context": context,
            "disease": disease,
            "lineage": "organoid" if sid in CONTEXT_A else "fibroblast",
        }
        rows.append(row)
    df = pd.DataFrame(rows)
    for compound_id in sorted(SELECTIVE):
        sens = []
        for sid in SAMPLES:
            base = RNG.normal(0.8, 0.05)
            if sid in CONTEXT_A:
                base = RNG.normal(0.15, 0.05)
            elif compound_id == "BRD-0007":
                base = RNG.normal(0.55, 0.08)
            sens.append(max(0, min(1.2, base)))
        df[f"sensitivity_{compound_id}"] = sens
    return df


def build_readout(treatment: pd.DataFrame, sample_info: pd.DataFrame) -> pd.DataFrame:
    trt = treatment.set_index("well_id")
    data = {}
    for sample_id in SAMPLES:
        vals = []
        context = "context_A" if sample_id in CONTEXT_A else "context_B"
        for well_id in trt.index:
            row = trt.loc[well_id]
            if row["perturbation_type"] == "vehicle_control":
                mfi = RNG.normal(15000, 800)
            elif row["perturbation_type"] == "positive_control":
                mfi = RNG.normal(3500, 400)
            else:
                compound_id = row["compound_id"]
                if compound_id in SELECTIVE and context == "context_A":
                    mfi = RNG.normal(5000, 600)
                elif compound_id in SELECTIVE:
                    mfi = RNG.normal(14000, 700)
                else:
                    mfi = RNG.normal(14500, 900)
            vals.append(max(500, mfi))
        data[sample_id] = vals
    wide = pd.DataFrame(data, index=trt.index).T
    wide.index.name = "sample_id"
    return wide


def build_features(sample_info: pd.DataFrame, demo_dir: Path) -> None:
    feat_dir = demo_dir / "features"
    feat_dir.mkdir(parents=True, exist_ok=True)
    expr = pd.DataFrame(
        RNG.normal(5, 1, size=(len(SAMPLES), 30)),
        index=SAMPLES,
        columns=[f"GENE_{i}" for i in range(1, 31)],
    )
    expr["MT1A"] = sample_info.set_index("sample_id").loc[SAMPLES, "sensitivity_BRD-0003"].values * -2 + RNG.normal(8, 0.3, len(SAMPLES))
    meth = pd.DataFrame(
        RNG.normal(0.5, 0.1, size=(len(SAMPLES), 20)),
        index=SAMPLES,
        columns=[f"cg_{i:05d}" for i in range(1, 21)],
    )
    meth["cg_context_A"] = (sample_info.set_index("sample_id").loc[SAMPLES, "context"] == "context_A").astype(float) + RNG.normal(0, 0.05, len(SAMPLES))
    expr.to_csv(feat_dir / "expression.csv")
    meth.to_csv(feat_dir / "methylation.csv")


def build_configs(demo_dir: Path) -> None:
    schema = {
        "name": "demo_toy_screen",
        "screen_topology": "primary_only",
        "paths": {
            "primary_readout": "readouts/primary.csv",
            "treatment_info": "metadata/treatment_info.csv",
            "sample_info": "metadata/sample_info.csv",
            "features_dir": "features",
        },
        "columns": {
            "sample_id": "sample_id",
            "well_id": "well_id",
            "perturbation_type": "perturbation_type",
            "detection_plate": "detection_plate",
            "compound_id": "compound_id",
            "compound_name": "compound_name",
        },
        "controls": {"negative": "vehicle_control", "positive": "positive_control"},
        "qc": {"ssmd_cutoff": 1.5},
        "hit_calling": {"viability_cutoff": 0.5, "robust_z_cutoff": -2.0, "min_samples": 2},
    }
    objective = {
        "name": "Approved compounds selective in IBD organoid context",
        "compound_filters": {"clinical_phase": ["Launched", "Phase 2", "Phase 3"]},
        "target_context": {"sample_info_query": "context == 'context_A'"},
        "off_target_context": {"sample_info_query": "context == 'context_B'"},
        "selectivity": {"prefer": "context_dependent"},
        "priority_weights": {
            "selectivity": 0.30,
            "biomarker_strength": 0.25,
            "clinical_phase": 0.20,
            "mechanism_novelty": 0.15,
            "phenocopy_support": 0.10,
        },
    }
    (demo_dir / "schema.yaml").write_text(
        __import__("yaml").dump(schema, sort_keys=False)
    )
    (demo_dir / "objective.yaml").write_text(
        __import__("yaml").dump(objective, sort_keys=False)
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--output-dir", type=Path, default=DEMO_DIR,
        help="Where to write the bundle (default: the skill's checked-in demo/).",
    )
    demo_dir = parser.parse_args(argv).output_dir
    demo_dir.mkdir(parents=True, exist_ok=True)
    (demo_dir / "readouts").mkdir(exist_ok=True)
    (demo_dir / "metadata").mkdir(exist_ok=True)

    treatment = build_treatment_info()
    sample_info = build_sample_info()
    primary = build_readout(treatment, sample_info)

    treatment.to_csv(demo_dir / "metadata/treatment_info.csv", index=False)
    sample_info.to_csv(demo_dir / "metadata/sample_info.csv", index=False)
    primary.to_csv(demo_dir / "readouts/primary.csv")
    build_features(sample_info, demo_dir)
    build_configs(demo_dir)

    manifest = {
        "description": "Synthetic toy repurposing screen (10 samples, 20 compounds, 2 plates)",
        "selective_compounds": sorted(SELECTIVE),
        "target_context": "context_A (IBD organoid)",
        "reference_context": "context_B (fibroblast)",
    }
    (demo_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"Demo bundle written to {demo_dir}")


if __name__ == "__main__":
    main()
