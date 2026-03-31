"""
io.py — Input validation and Boltz-ready file preparation.

Accepts .yaml/.yml (single protein or multi-chain complex).
Validates sequences against the standard amino acid alphabet.
Writes a normalised input file that Boltz can consume.
"""
from __future__ import annotations

from pathlib import Path

# Standard 20 amino acids + selenocysteine (U) + pyrrolysine (O)
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYUO")


def validate_and_prepare(
    input_path: Path,
    work_dir: Path,
) -> dict:
    """Validate and prepare YAML input for Boltz.

    Injects ``msa: empty`` into protein entries so that Boltz runs fully
    offline (no MSA server, no data egress).

    Returns:
        {
            "boltz_input_path": Path,
            "input_type": "yaml",
            "sequences": [{"name": str, "sequence": str, "chain_id": str}]
        }
    """
    try:
        import yaml
    except ImportError:
        raise ImportError("PyYAML is required. Install with: uv pip install pyyaml")

    input_path = Path(input_path)
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    if input_path.suffix.lower() not in (".yaml", ".yml"):
        raise ValueError(
            f"Expected a .yaml/.yml input file, got: {input_path.suffix!r}. "
            "FASTA input is no longer supported."
        )

    data = yaml.safe_load(input_path.read_text())
    sequences = _parse_sequences_from_data(data, input_path)
    _validate_sequences(sequences)

    if len(sequences) > 26:
        raise ValueError(
            f"Too many chains ({len(sequences)}). Boltz-2 supports at most 26 chains (A–Z)."
        )
    for i, seq in enumerate(sequences):
        if "chain_id" not in seq:
            seq["chain_id"] = chr(ord("A") + i)

    out_path = (work_dir / input_path.stem).with_suffix(input_path.suffix)
    _inject_msa_empty(data)
    out_path.write_text(yaml.dump(data, default_flow_style=False, allow_unicode=True))

    return {
        "boltz_input_path": out_path,
        "input_type": "yaml",
        "sequences": sequences,
    }


def _parse_sequences_from_data(data: dict, input_path: Path | None = None) -> list[dict]:
    """Extract a normalised sequence list from an already-parsed YAML dict."""
    if not isinstance(data, dict) or "sequences" not in data:
        raise ValueError(
            f"YAML input must have a top-level 'sequences' key. "
            f"Got: {list(data.keys()) if isinstance(data, dict) else type(data)}"
        )

    sequences = []
    for i, entry in enumerate(data["sequences"]):
        # Native Boltz format: {protein: {id: A, sequence: ...}} or {ligand: {id: B, smiles: ...}}
        if "sequence" not in entry and "entity_type" not in entry:
            entity_type = None
            inner = None
            for key in ("protein", "rna", "dna", "ligand"):
                if key in entry:
                    entity_type = key
                    inner = entry[key]
                    break
            if inner is None:
                raise ValueError(
                    f"YAML entry {i} is not in a recognised format. "
                    "Expected a 'protein', 'rna', 'dna', or 'ligand' key, "
                    "or a flat entry with 'sequence' and 'entity_type'."
                )
            if entity_type == "ligand":
                raw_ids = inner.get("id", f"lig_{i}")
                ids = raw_ids if isinstance(raw_ids, list) else [raw_ids]
                for chain_id in ids:
                    seq_entry: dict = {
                        "name": str(chain_id),
                        "sequence": inner.get("smiles") or inner.get("ccd") or "",
                        "entity_type": "ligand",
                    }
                    if inner.get("smiles"):
                        seq_entry["smiles"] = inner["smiles"]
                    if inner.get("ccd"):
                        seq_entry["ccd"] = inner["ccd"]
                    sequences.append(seq_entry)
                continue
            raw_ids = inner.get("id", f"chain_{i}")
            ids = raw_ids if isinstance(raw_ids, list) else [raw_ids]
            seq_str = str(inner.get("sequence", "")).upper()
            for chain_id in ids:
                sequences.append({
                    "name": str(chain_id),
                    "sequence": seq_str,
                    "entity_type": entity_type,
                })
            continue

        # Flat legacy format: {id: A, sequence: ..., entity_type: protein}
        if "sequence" not in entry:
            raise ValueError(
                f"YAML entry {i} (id={entry.get('id', '?')}) is missing 'sequence' field."
            )
        sequences.append({
            "name": entry.get("id", f"chain_{i}"),
            "sequence": str(entry["sequence"]).upper(),
            "entity_type": entry.get("entity_type", "protein"),
        })

    return sequences


# Keep _parse_yaml as an alias for callers outside this module
def _parse_yaml(path: Path) -> list[dict]:
    try:
        import yaml
    except ImportError:
        raise ImportError("PyYAML is required. Install with: uv pip install pyyaml")
    data = yaml.safe_load(path.read_text())
    return _parse_sequences_from_data(data, path)


def _inject_msa_empty(data: dict) -> None:
    """Mutate parsed YAML data to add ``msa: empty`` to protein/rna/dna entries lacking it.

    Boltz requires either an MSA path or ``msa: empty`` for every polymer chain
    when not using the MSA server.
    """
    for entry in data.get("sequences", []):
        for polymer_key in ("protein", "rna", "dna"):
            if polymer_key in entry:
                inner = entry[polymer_key]
                if isinstance(inner, dict) and "msa" not in inner:
                    inner["msa"] = "empty"
        # flat legacy format
        if "entity_type" in entry and entry["entity_type"] == "protein" and "msa" not in entry:
            entry["msa"] = "empty"


def _validate_sequences(sequences: list[dict]) -> None:
    """Validate all sequences against the amino acid alphabet."""
    for seq_info in sequences:
        seq = seq_info["sequence"]
        entity_type = seq_info.get("entity_type", "protein")
        if entity_type != "protein":
            continue  # skip ligands / nucleic acids
        if len(seq) == 0:
            raise ValueError(f"Sequence '{seq_info['name']}' is empty.")
        invalid = sorted(set(seq) - VALID_AA)
        if invalid:
            raise ValueError(
                f"Invalid amino acid character(s) {invalid} in sequence '{seq_info['name']}'. "
                f"Valid characters: {sorted(VALID_AA)}"
            )
