"""Every ORCID iD in the citation and maintainer metadata must be a valid iD.

From v0.5.0 to v0.7.1 CITATION.cff, .zenodo.json and MAINTAINERS.md carried
0000-0002-5765-9827 for the lead maintainer. That iD fails the ORCID checksum
and resolves to nothing, so every Zenodo release minted from these files was
attributed to a non-existent ORCID record. A link-presence check cannot see
this; the ISO 7064 MOD 11-2 check digit can, offline.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METADATA_FILES = ["CITATION.cff", ".zenodo.json", "MAINTAINERS.md"]
ORCID_RE = re.compile(r"\b(\d{4}-\d{4}-\d{4}-\d{3}[\dX])\b")


def orcid_checksum_ok(orcid: str) -> bool:
    digits = orcid.replace("-", "")
    total = 0
    for ch in digits[:-1]:
        total = (total + int(ch)) * 2
    remainder = (12 - total % 11) % 11
    return digits[-1] == ("X" if remainder == 10 else str(remainder))


def _ids():
    found = []
    for name in METADATA_FILES:
        path = ROOT / name
        if path.exists():
            found += [(name, i) for i in ORCID_RE.findall(path.read_text())]
    return found


def test_checksum_accepts_orcids_published_example():
    # The worked example in ORCID's own documentation.
    assert orcid_checksum_ok("0000-0002-1825-0097")
    assert not orcid_checksum_ok("0000-0002-1825-0098")


def test_metadata_declares_at_least_one_orcid():
    assert _ids(), "no ORCID iD found; the checksum test below would pass vacuously"


def test_every_metadata_orcid_has_a_valid_check_digit():
    bad = [f"{name}: {orcid}" for name, orcid in _ids() if not orcid_checksum_ok(orcid)]
    assert bad == [], f"invalid ORCID iDs: {bad}"
