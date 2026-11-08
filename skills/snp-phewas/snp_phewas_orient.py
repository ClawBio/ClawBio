"""Re-express a biobank record on the Open Targets effect allele (the Open Targets `alt`).

Every source here is reference-anchored: FinnGen reports GRCh38 reference / alternative
alleles, Pan-UKBB GRCh37 reference / alternative, Open Targets GRCh38 reference /
alternative, and the GRCh37 placement used to read Pan-UKBB carries the mapping strand
(alleles are reverse-complemented on a minus-strand mapping before comparison). So the
strand is known and allele ORDER decides orientation:

  * same order                    -> keep as is
  * swapped order, not A/T or C/G -> flip the sign of beta and take 1 - EAF
  * swapped order, A/T or C/G     -> withhold the direction. Under the
    reference-anchored palindromic rule a strand-ambiguous pair is kept only when both
    sides list it in the same (reference, alternative) order, so a swap cannot be
    resolved by allele identity. The p-value and SE stay (they do not depend on
    orientation); beta and EAF are withheld with `direction_withheld_reason` set.
  * any other alleles             -> `allele_mismatch` (a different variant at the position)

No allele-frequency rescue is attempted for a swapped palindromic pair: that rescue is
for two studies on unknown strands, which is not the situation here, and Open Targets
carries no allele frequency for the target side to compare against.
"""

from __future__ import annotations

from dataclasses import dataclass

from snp_phewas_schema import WITHHELD_PALINDROMIC_ORDER

# The four strand-ambiguous (palindromic) SNV allele pairs.
PALINDROMIC_PAIRS = frozenset({("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")})


def reference_anchored_in_order(src_ref: str, src_alt: str, tgt_ref: str, tgt_alt: str) -> bool:
    """True when both sides list the alleles in the same (reference, alternative) order.

    The reference-anchored palindromic rule: when both studies report alleles against
    the forward-strand reference, a palindromic pair is trustworthy only in matching
    order; in swapped order the strand cannot be told from the alleles alone.
    """
    return (src_ref.upper(), src_alt.upper()) == (tgt_ref.upper(), tgt_alt.upper())


_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def complement(allele: str) -> str:
    return allele.translate(_COMPLEMENT).upper()


@dataclass(frozen=True)
class Orientation:
    matches: bool
    flip: bool = False
    palindromic: bool = False
    withheld_reason: str | None = None


def orient(src_ref: str, src_alt: str, tgt_ref: str, tgt_alt: str) -> Orientation:
    """How the source record (src_ref, src_alt) relates to the target (tgt_ref, tgt_alt)."""
    s = (src_ref.upper(), src_alt.upper())
    t = (tgt_ref.upper(), tgt_alt.upper())
    palin = t in PALINDROMIC_PAIRS
    if s == t:
        return Orientation(matches=True, palindromic=palin)
    if s == (t[1], t[0]):
        if palin and not reference_anchored_in_order(s[0], s[1], t[0], t[1]):
            return Orientation(matches=True, palindromic=True,
                               withheld_reason=WITHHELD_PALINDROMIC_ORDER)
        return Orientation(matches=True, flip=True, palindromic=palin)
    return Orientation(matches=False, palindromic=palin)


def apply(o: Orientation, beta: float | None, eaf: float | None) -> tuple[float | None, float | None]:
    """(beta, eaf) on the target effect allele."""
    if o.withheld_reason:
        return None, None
    if o.flip:
        return (None if beta is None else -beta), (None if eaf is None else 1.0 - eaf)
    return beta, eaf
