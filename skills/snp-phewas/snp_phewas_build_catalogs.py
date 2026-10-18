"""Rebuild the in-repo outcome catalogs from each source's published manifest.

    python skills/snp-phewas/snp_phewas_build_catalogs.py [--pan-ukbb-manifest PATH]
                                                          [--finngen-manifest PATH]

Maintenance only (network; not run by the skill or its tests). It rewrites the two
outcome catalogs under `catalogs/`.

With no paths, the manifests are downloaded from the URLs pinned in `datasets.yaml`.
Either way the file's md5 must equal the pin, so a source that changed under us fails
loudly instead of silently producing a different universe (a release bump edits the pin
and re-derives everything keyed on it).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import os
import sys

import requests

from snp_phewas_catalog import PAN_UKBB_POPS, _SKILL_DIR, _load_yaml

PAN_UKBB_COLUMNS: tuple[str, ...] = (
    "outcome_id", "filename", "filename_tabix", "trait_type", "phenocode", "pheno_sex",
    "coding", "modifier", "description", "coding_description", "category", "pops",
    *(f"n_cases_{p}" for p in PAN_UKBB_POPS), *(f"n_controls_{p}" for p in PAN_UKBB_POPS),
    "md5_hex", "size_in_bytes",
)
FINNGEN_COLUMNS: tuple[str, ...] = (
    "outcome_id", "phenotype", "category", "num_cases", "num_controls", "filename",
)


class PinMismatchError(RuntimeError):
    """The manifest read does not match the md5 pinned in datasets.yaml."""


def pan_ukbb_outcome_id(row: dict[str, str]) -> str:
    """Compose the Pan-UKBB outcome id: trait_type-phenocode-pheno_sex[-coding][-modifier]."""
    parts = [row["trait_type"], row["phenocode"], row["pheno_sex"]]
    parts += [x for x in (row.get("coding") or "", row.get("modifier") or "") if x]
    return "-".join(parts)


def _read_bytes(path_or_url: str) -> bytes:
    if path_or_url.startswith(("http://", "https://")):
        r = requests.get(path_or_url, timeout=120)
        r.raise_for_status()
        return r.content
    with open(path_or_url, "rb") as fh:
        return fh.read()


def _check_pin(data: bytes, expected_md5: str, what: str) -> None:
    got = hashlib.md5(data).hexdigest()
    if got != expected_md5:
        raise PinMismatchError(
            f"{what} md5 {got} != pinned {expected_md5}; the source changed. Re-pin in "
            f"datasets.yaml and rebuild everything keyed on the release before regenerating."
        )


def build_pan_ukbb(manifest: bytes) -> list[dict[str, str]]:
    text = gzip.decompress(manifest).decode("utf-8")
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    out = []
    for r in rows:
        rec = {c: r.get(c, "") for c in PAN_UKBB_COLUMNS if c != "outcome_id"}
        rec["outcome_id"] = pan_ukbb_outcome_id(r)
        out.append(rec)
    if len({r["outcome_id"] for r in out}) != len(out):
        raise ValueError("Pan-UKBB manifest composes duplicate outcome ids")
    return out


def build_finngen(manifest: bytes) -> list[dict[str, str]]:
    rows = list(csv.DictReader(io.StringIO(manifest.decode("utf-8")), delimiter="\t"))
    out = []
    for r in rows:
        out.append({
            "outcome_id": r["phenocode"],
            "phenotype": r["phenotype"],
            "category": r["category"],
            "num_cases": r["num_cases"],
            "num_controls": r["num_controls"],
            "filename": r["path_https"].rsplit("/", 1)[-1],
        })
    if len({r["outcome_id"] for r in out}) != len(out):
        raise ValueError("FinnGen manifest has duplicate phenocodes")
    return out


def _write(rows: list[dict[str, str]], columns: tuple[str, ...], path: str) -> None:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow(r)
    # mtime=0 so a rebuild from the same manifest is byte-identical.
    with open(path, "wb") as fh:
        with gzip.GzipFile(fileobj=fh, mode="wb", mtime=0) as gz:
            gz.write(buf.getvalue().encode("utf-8"))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pan-ukbb-manifest")
    ap.add_argument("--finngen-manifest")
    args = ap.parse_args(argv)
    fams = _load_yaml()["families"]

    pin = fams["pan_ukbb"]["release_pin"]
    data = _read_bytes(args.pan_ukbb_manifest or pin["manifest_url"])
    _check_pin(data, pin["manifest_md5"], "Pan-UKBB phenotype manifest")
    pan = build_pan_ukbb(data)
    _write(pan, PAN_UKBB_COLUMNS, os.path.join(_SKILL_DIR, fams["pan_ukbb"]["outcome_catalog"]))

    pin = fams["finngen"]["release_pin"]
    data = _read_bytes(args.finngen_manifest or pin["manifest_url"])
    _check_pin(data, pin["manifest_md5"], "FinnGen R12 manifest")
    fg = build_finngen(data)
    _write(fg, FINNGEN_COLUMNS, os.path.join(_SKILL_DIR, fams["finngen"]["outcome_catalog"]))
    print(f"wrote {len(pan)} Pan-UKBB outcomes and {len(fg)} FinnGen endpoints")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
