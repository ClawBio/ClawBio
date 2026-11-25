"""PatientProfile — upload once, query many times.

Stores parsed genotypes and accumulated skill results in a JSON file
so multiple skills can reuse the same parsed data.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from clawbio.common.checksums import sha256_file
from clawbio.common.parsers import (
    GenotypeRecord,
    parse_genetic_file,
)


class PatientProfile:
    """In-memory patient profile with JSON persistence."""

    def __init__(
        self,
        patient_id: str = "",
        input_file: str = "",
        checksum: str = "",
        upload_date: str = "",
        genotypes: dict[str, dict] | None = None,
        ancestry: dict | None = None,
        skill_results: dict[str, Any] | None = None,
    ):
        self.metadata = {
            "patient_id": patient_id,
            "input_file": input_file,
            "checksum": checksum,
            "upload_date": upload_date or datetime.now(timezone.utc).isoformat(),
        }
        self._genotypes: dict[str, dict] = genotypes or {}
        self.ancestry: dict | None = ancestry
        self.skill_results: dict[str, Any] = skill_results or {}

    # --- Construction from file ---

    @classmethod
    def from_genetic_file(
        cls,
        filepath: str | Path,
        patient_id: str = "",
        fmt: str = "auto",
    ) -> "PatientProfile":
        """Parse a genetic file and create a new profile."""
        filepath = Path(filepath)
        records = parse_genetic_file(filepath, fmt=fmt)

        # Store as serializable dicts
        genotypes = {rsid: rec.to_dict() for rsid, rec in records.items()}
        checksum = sha256_file(filepath)

        if not patient_id:
            patient_id = filepath.stem.replace(" ", "_")[:32]

        return cls(
            patient_id=patient_id,
            input_file=str(filepath.resolve()),
            checksum=checksum,
            genotypes=genotypes,
        )

    # --- Genotype access ---

    @property
    def genotype_count(self) -> int:
        return len(self._genotypes)

    def get_genotypes(self, rsids: list[str] | None = None) -> dict[str, str]:
        """Get simple {rsid: genotype_str} dict, optionally filtered to a subset."""
        if rsids is None:
            return {rsid: rec.get("genotype", "") for rsid, rec in self._genotypes.items()}
        return {
            rsid: self._genotypes[rsid].get("genotype", "")
            for rsid in rsids
            if rsid in self._genotypes
        }

    def get_records(self, rsids: list[str] | None = None) -> dict[str, GenotypeRecord]:
        """Get full GenotypeRecord objects, optionally filtered."""
        if rsids is None:
            targets = self._genotypes.items()
        else:
            targets = ((r, self._genotypes[r]) for r in rsids if r in self._genotypes)
        return {
            rsid: GenotypeRecord(**rec)
            for rsid, rec in targets
        }

    # --- Skill results ---

    def add_skill_result(self, skill_name: str, result_dict: dict) -> None:
        """Store the result of a skill run."""
        self.skill_results[skill_name] = {
            "run_at": datetime.now(timezone.utc).isoformat(),
            "data": result_dict,
        }

    def get_skill_result(self, skill_name: str) -> dict | None:
        """Retrieve a previous skill result."""
        entry = self.skill_results.get(skill_name)
        if entry:
            return entry.get("data")
        return None

    # --- Persistence ---

    def save(self, path: str | Path) -> Path:
        """Save profile to a JSON file."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        data = {
            "metadata": self.metadata,
            "genotypes": self._genotypes,
            "ancestry": self.ancestry,
            "skill_results": self.skill_results,
        }
        path.write_text(json.dumps(data, indent=2, default=str))
        return path

    @classmethod
    def load(cls, path: str | Path) -> "PatientProfile":
        """Load profile from a JSON file."""
        path = Path(path)
        data = json.loads(path.read_text())
        meta = data.get("metadata", {})
        return cls(
            patient_id=meta.get("patient_id", ""),
            input_file=meta.get("input_file", ""),
            checksum=meta.get("checksum", ""),
            upload_date=meta.get("upload_date", ""),
            genotypes=data.get("genotypes"),
            ancestry=data.get("ancestry"),
            skill_results=data.get("skill_results"),
        )

    def __repr__(self) -> str:
        pid = self.metadata.get("patient_id", "unknown")
        n = self.genotype_count
        skills = list(self.skill_results.keys())
        return f"PatientProfile(id={pid!r}, genotypes={n}, skills={skills})"
