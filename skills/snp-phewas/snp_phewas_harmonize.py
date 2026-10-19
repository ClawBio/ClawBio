"""Attach the cross-biobank outcome harmonization to lookup rows (thin; the crossings live
in `snp_phewas_crosswalk` and its committed table).

Each row gets `harmonized_outcome_id`, `harmonized_label`, `harmonization_basis`; an
outcome with no counterpart keeps its row with basis `unmapped` (or `icd10_broad_only`)
and stays out of the cross-dataset views. When several outcomes of ONE dataset share a
harmonized id (several FinnGen endpoints FinnGen matched to one Pan-UKBB outcome), one is
the representative the replication view compares (`harmonized_representative`): one the
source tested, then the strongest basis, then the most cases, then the id. A row whose
basis is never compared (the maximum / minimum summary of an MVP measurement) keeps its id
and is never the representative.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any

from snp_phewas_crosswalk import (
    BASIS_PRECEDENCE,
    BASIS_UNMAPPED,
    HARMONIZATION_BASIS_LABELS,
    NON_REPRESENTATIVE_BASES,
    load_crosswalk,
)

import snp_phewas_schema as S
from snp_phewas_catalog import get_dataset

_COMPARABLE = S.TESTED_STATUSES | {S.ST_API_CENSORED}


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
        if x.harmonization_basis in NON_REPRESENTATIVE_BASES:
            # the maximum / minimum of an MVP measurement whose mean is the one compared
            r.harmonized_representative = False
            continue
        if x.harmonized_outcome_id and in_catalog:
            groups[(r.dataset_id, x.harmonized_outcome_id)].append(r)
    for members in groups.values():
        # A member the source tested comes first, so an outage or absence in one endpoint
        # does not hide a sibling's comparable value. The rule looks at whether a value
        # exists, never at how strong it is. Next, a study known unreadable from the catalog
        # never outranks a readable one: the threshold does not count it, so letting an
        # outage hand it the place would cut n_shared and loosen 0.05 / n.
        members.sort(key=lambda r: (r.status not in _COMPARABLE,
                                    r.status == S.ST_NOT_READABLE,
                                    BASIS_PRECEDENCE.get(r.harmonization_basis, 9),
                                    -(r.n_cases or 0), r.outcome_id))
        for i, r in enumerate(members):
            r.harmonized_representative = i == 0


def basis_label(basis: str | None) -> str:
    return HARMONIZATION_BASIS_LABELS.get(basis or "", basis or "")
