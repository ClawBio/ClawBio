"""Attach the cross-biobank outcome harmonization to lookup rows (thin; the crossings live
in `snp_phewas_crosswalk` and its committed table).

Each row gets `harmonized_outcome_id`, `harmonized_label`, `harmonization_basis`; an
outcome with no counterpart keeps its row with basis `unmapped` (or `icd10_broad_only`)
and stays out of the cross-dataset views. When several outcomes of ONE dataset share a
harmonized id (several FinnGen endpoints FinnGen matched to one Pan-UKBB outcome), one is
the representative the replication view compares (`harmonized_representative`): the
strongest basis, then the most cases, then the id.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any

from snp_phewas_crosswalk import (
    BASIS_PRECEDENCE,
    BASIS_UNMAPPED,
    HARMONIZATION_BASIS_LABELS,
    load_crosswalk,
)

from snp_phewas_catalog import get_dataset


def attach_harmonization(result: Any) -> None:
    cw = load_crosswalk()
    groups: dict[tuple[str, str], list] = defaultdict(list)
    family = {d: get_dataset(d).family for d in result.dataset_ids}
    for r, in_catalog in [*((r, True) for r in result.rows),
                          *((r, False) for r in result.uncatalogued)]:
        fam = family[r.dataset_id]
        x = cw.get((fam, r.outcome_id))
        if x is None:
            r.harmonization_basis = BASIS_UNMAPPED
            r.harmonized_representative = None
            continue
        r.harmonized_outcome_id = x.harmonized_outcome_id
        r.harmonized_label = x.harmonized_label
        r.harmonization_basis = x.harmonization_basis
        if x.harmonized_outcome_id and in_catalog:
            groups[(r.dataset_id, x.harmonized_outcome_id)].append(r)
    for members in groups.values():
        members.sort(key=lambda r: (BASIS_PRECEDENCE.get(r.harmonization_basis, 9),
                                    -(r.n_cases or 0), r.outcome_id))
        for i, r in enumerate(members):
            r.harmonized_representative = i == 0


def basis_label(basis: str | None) -> str:
    return HARMONIZATION_BASIS_LABELS.get(basis or "", basis or "")
