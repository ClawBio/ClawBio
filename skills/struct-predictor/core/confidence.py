"""
confidence.py — Extract pLDDT and PAE from Boltz-2 output files.

pLDDT is stored as B-factors in the output CIF (CA atoms only).
PAE is stored in the companion confidence JSON under key "pae".
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def extract_confidence(
    cif_path: Path,
    confidence_json_path: Path | None,
) -> dict:
    """Extract pLDDT, PAE, and chain boundaries from Boltz output.

    Returns:
        {
            "plddt": np.ndarray,           # shape [n_residues], float32
            "pae": np.ndarray,             # shape [n_residues, n_residues], float32
            "chain_boundaries": [
                {"chain_id": str, "start": int, "end": int}
            ]
        }
    """
    plddt, chain_boundaries = _parse_cif_atoms(cif_path)
    pae = _parse_pae_from_json(confidence_json_path, len(plddt))
    return {
        "plddt": plddt,
        "pae": pae,
        "chain_boundaries": chain_boundaries,
    }


def _read_atom_site_columns(cif_path: Path) -> dict[str, list[str]]:
    """Minimal mmCIF loop reader — extracts _atom_site columns without biopython.

    Finds the loop_ block that declares _atom_site fields, records column
    indices, then collects the data rows that follow.  Only the four columns
    needed for pLDDT/chain extraction are returned; other columns are ignored.
    """
    # Use lowercase keys throughout to avoid case-sensitivity bugs
    WANTED = {
        "_atom_site.label_atom_id",
        "_atom_site.label_seq_id",
        "_atom_site.label_asym_id",
        "_atom_site.b_iso_or_equiv",
    }
    # Map lowercase key → original key for the returned dict
    WANTED_ORIG = {
        "_atom_site.label_atom_id": "_atom_site.label_atom_id",
        "_atom_site.label_seq_id": "_atom_site.label_seq_id",
        "_atom_site.label_asym_id": "_atom_site.label_asym_id",
        "_atom_site.b_iso_or_equiv": "_atom_site.B_iso_or_equiv",
    }

    lines = Path(cif_path).read_text(errors="replace").splitlines()

    # --- locate the _atom_site loop ---
    # Only peek at consecutive field-declaration lines (starting with "_"),
    # so we don't accidentally see fields from a later loop block.
    header_start = None
    for i, line in enumerate(lines):
        if line.strip() != "loop_":
            continue
        for j in range(i + 1, len(lines)):
            s = lines[j].strip()
            if not s or s.startswith("#"):
                continue  # skip blank/comment lines between loop_ and fields
            if s.startswith("_"):
                if s.lower().startswith("_atom_site."):
                    header_start = i
                break  # first non-blank non-comment non-field line → wrong loop
        if header_start is not None:
            break

    if header_start is None:
        return {v: [] for v in WANTED_ORIG.values()}

    # --- collect column names (lowercased) ---
    col_names: list[str] = []
    data_start = header_start + 1
    for i in range(header_start + 1, len(lines)):
        s = lines[i].strip()
        if s.startswith("_"):
            col_names.append(s.split()[0].lower())
            data_start = i + 1
        else:
            data_start = i
            break

    wanted_idx = {name: idx for idx, name in enumerate(col_names) if name in WANTED}

    # Build result with original-case keys
    result: dict[str, list[str]] = {WANTED_ORIG[k]: [] for k in WANTED}

    # --- collect data rows until next loop_ or data_ block ---
    for line in lines[data_start:]:
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("loop_") or s.startswith("data_"):
            if s.startswith("loop_") or s.startswith("data_"):
                break
            continue
        tokens = s.split()
        for lc_name, idx in wanted_idx.items():
            if idx < len(tokens):
                result[WANTED_ORIG[lc_name]].append(tokens[idx])

    return result


def _parse_cif_atoms(cif_path: Path) -> tuple[np.ndarray, list[dict]]:
    """Single-pass CIF parse returning (plddt, chain_boundaries).

    Reads _atom_site records once, extracting CA B-factors (pLDDT) and
    chain boundary information in the same iteration.
    """
    mmcif = _read_atom_site_columns(cif_path)

    atom_ids  = mmcif.get("_atom_site.label_atom_id", [])
    seq_ids   = mmcif.get("_atom_site.label_seq_id", [])
    chain_ids = mmcif.get("_atom_site.label_asym_id", [])
    bfacs     = mmcif.get("_atom_site.B_iso_or_equiv", [])

    if not atom_ids:
        raise ValueError(f"No _atom_site records found in {cif_path}")

    seen: dict[tuple[str, str], float] = {}
    order: list[tuple[str, str]] = []

    for atom, seq, chain, bfac in zip(atom_ids, seq_ids, chain_ids, bfacs):
        if atom.strip() != "CA":
            continue
        key = (chain.strip(), seq.strip())
        if key not in seen:
            seen[key] = float(bfac)
            order.append(key)

    if not seen:
        raise ValueError(
            f"No CA atoms found in {cif_path}. "
            "Check that the CIF file is a valid Boltz output."
        )

    plddt = np.array([seen[k] for k in order], dtype=np.float32)

    boundaries: list[dict] = []
    current_chain = order[0][0]
    start_idx = 0
    for i, (chain, _) in enumerate(order):
        if chain != current_chain:
            boundaries.append({"chain_id": current_chain, "start": start_idx, "end": i - 1})
            current_chain = chain
            start_idx = i
    boundaries.append({"chain_id": current_chain, "start": start_idx, "end": len(order) - 1})

    return plddt, boundaries


def _parse_plddt_from_cif(cif_path: Path) -> np.ndarray:
    """Extract per-residue pLDDT from CIF B-factors (CA atoms only)."""
    plddt, _ = _parse_cif_atoms(cif_path)
    return plddt


def _parse_pae_from_json(
    confidence_json_path: Path | None,
    n_residues: int,
) -> np.ndarray:
    """Load PAE matrix from Boltz confidence JSON.

    If the JSON is absent or lacks a 'pae' key, returns a zero matrix.
    """
    if confidence_json_path is None or not Path(confidence_json_path).exists():
        return np.zeros((n_residues, n_residues), dtype=np.float32)

    data = json.loads(Path(confidence_json_path).read_text())

    if "pae" not in data:
        return np.zeros((n_residues, n_residues), dtype=np.float32)

    return np.array(data["pae"], dtype=np.float32)
