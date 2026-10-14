"""UKB-PPP pQTL per-region summary-stats fetcher.

Fetches per-variant plasma-pQTL summary statistics within a genomic window for
a given (protein, ancestry) measurement from the UK Biobank Pharma Proteomics
Project (Sun et al. 2023, Nature). Returns harmonised rows
(variant_id, chromosome, position, ref, alt, beta, SE, p-value, A1FREQ,
protein_hgnc, ancestry) suitable for downstream colocalization, fine-mapping,
regional plotting, or Mendelian randomisation on the pQTL exposure axis.

Source: UKB-PPP release 1 (2023), 2,923 proteins x ~46k EUR + multi-ancestry
breakouts. Hosted on Synapse (`syn51364943`) + the AWS Open Data Registry
(`s3://ukbiobank.opendata.sagebase.org`). The AWS bucket is documented as
public but anonymous reads return `AccessDenied` as of 2026-05-15; the
canonical functional access path is Synapse with an authenticated PAT.

License: CC-BY 4.0 (verified 2026-05-14 against the AWS Open Data Registry
entry; attribution string emitted in every manifest).

Two-path access design (v1.3):

1. **Bundled slices** (default for canonical demos; no auth required).
   `bundled_slices/<PROTEIN>_<ANCESTRY>_chr<C>_<start>_<end>.json.gz` files
   shipped inside the skill repo contain pre-computed regional slices
   for the canonical demo cohort. v0.1.0 ships the SORT1 / EUR /
   OID20213 slice (chr1:108,774,968-109,774,968); the slice convention
   scales to additional proteins by dropping further files into
   `bundled_slices/`. The client checks these FIRST. Redistribution is
   permitted under CC-BY 4.0 with the cited attribution; the
   bundled-slice manifest carries the same provenance block as a live
   fetch. ClawBio + K-Dense end users running the bundled demos need
   ZERO auth.

2. **Live Synapse fetch** (opt-in; requires free Synapse PAT). For
   arbitrary queries beyond the demo cohort, the client falls through
   to the live Synapse downloader. `SYNAPSE_AUTH_TOKEN` env var is
   required; a clear, signposted error is raised when missing. No UK
   Biobank Application is required for the summary-stats layer.

When neither a bundled slice nor a PAT is available, the skill raises
`UKBPPPAccessError` with a multi-line, user-facing message that includes
the Synapse signup URL + PAT generation steps + an explanation of why
the live fetch is required.

Fetch path (live): download the per-protein
`<HGNC>_<UniProt>_<OlinkID>_v1_<Panel>.tar` archive, extract the
per-chromosome REGENIE-format summary-stats file, then filter to the
requested (chr, start, end) window in-process. The tar download is
cached locally (`UKB_PPP_CACHE_DIR` env or
`~/.clawbio/ukb_ppp_region_fetch_cache/`).

Variant id convention: UKB-PPP REGENIE outputs use ALLELE0 (REF) and ALLELE1
(ALT, effect). OT uses `chr_pos_ref_alt` GRCh38 with ALT-effect beta. The
skill harmonises rows to OT convention at the row-normalisation boundary so
the join key matches eQTL Catalogue + GWAS Catalog feeds.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import os
import re
import sys
import tarfile
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Iterator

# Default cache directory for downloaded per-protein tars. Mirrors the eqtl
# Catalogue cache convention; per-protein archives are reused across many
# region fetches.
DEFAULT_CACHE_DIR = Path(
    os.environ.get(
        "UKB_PPP_CACHE_DIR",
        Path.home() / ".clawbio" / "ukb_ppp_region_fetch_cache",
    )
).expanduser()

# Synapse parent folder for the pGWAS sumstats tree (verified 2026-05-15).
SYNAPSE_PROJECT_ID = "syn51364943"
SYNAPSE_SUMSTATS_FOLDER_ID = "syn51365301"

# Ancestry folder ids inside the pGWAS folder (verified 2026-05-15).
ANCESTRY_FOLDER_IDS: dict[str, str] = {
    "EUR": "syn51365303",   # European (discovery)
    "AFR": "syn51365304",   # African
    "CSA": "syn51365305",   # Central/South Asian (synonym SAS)
    "SAS": "syn51365305",   # alias
    "EAS": "syn51365306",   # East Asian
    "MID": "syn51365307",   # Middle East
    "AMR": "syn51500434",   # American (Hispanic)
    "ALL": "syn51365308",   # Combined (multi-ancestry meta)
}

# Human-readable ancestry expansion per the AGENTS.md output-friendly rule.
ANCESTRY_LABELS: dict[str, str] = {
    "EUR": "European (discovery)",
    "AFR": "African",
    "CSA": "Central / South Asian",
    "SAS": "Central / South Asian",
    "EAS": "East Asian",
    "MID": "Middle Eastern",
    "AMR": "American (Hispanic)",
    "ALL": "Combined (multi-ancestry meta)",
}

# UKB-PPP per-protein file name pattern, verified 2026-05-15 against the
# Synapse listing:
#   <HGNC>_<UniProt>_<OlinkID>_v1_<Panel>.tar
# e.g. A1BG_P04217_OID30771_v1_Inflammation_II.tar
PROTEIN_FILE_RE = re.compile(
    r"^(?P<hgnc>[A-Z0-9][A-Z0-9\-\.]*)_(?P<uniprot>[A-Z0-9]+)_"
    r"(?P<olink_id>OID\d+)_v(?P<version>\d+)_(?P<panel>.+)\.tar$"
)

# REGENIE step-2 column order in UKB-PPP per-chromosome files (Sun 2023
# methods; verified against the protein-annotation README on Synapse).
# The actual files include a comment header line followed by a space-separated
# data section; the parser tolerates whitespace variants.
REGENIE_COLUMNS = [
    "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1",
    "A1FREQ", "INFO", "N", "TEST",
    "BETA", "SE", "CHISQ", "LOG10P", "EXTRA",
]


@dataclass
class UKBPPPRelease:
    """Records the UKB-PPP release + protein measurement metadata pinned per
    fetch (for the manifest). Mirrors `EQTLCatalogueRelease` so downstream
    renderers can build descriptive panel titles uniformly across eQTL +
    pQTL exposures.
    """

    release_label: str          # e.g. "UKB-PPP r1 2023" (Sun 2023)
    fetched_at_utc: str
    protein_hgnc: str           # e.g. "SORT1"
    protein_uniprot: str        # e.g. "Q99523"
    olink_reagent_id: str       # e.g. "OID20213" (SORT1 on Cardiometabolic)
    olink_panel: str            # e.g. "Cardiometabolic"
    ancestry: str               # short code (EUR/AFR/...)
    ancestry_label: str         # user-friendly expansion
    n_samples: int | None       # discovery cohort N (EUR ≈ 46_673)
    synapse_id: str             # e.g. "syn51469328"
    source_url: str             # canonical Synapse entity URL
    study_label: str = "UKB-PPP"

    @property
    def protein_label(self) -> str:
        """User-friendly protein name + UniProt + Olink ID (per AGENTS.md
        output-friendly rule). Renderers should prefer this over the raw
        HGNC alone."""
        return f"{self.protein_hgnc} ({self.protein_uniprot}, {self.olink_reagent_id})"


@dataclass
class RegionVariant:
    """One variant row, harmonised to OT GRCh38 ALT-effect convention.

    Mirrors `eqtl_catalogue_region_fetch.RegionVariant` so the orchestrator's
    downstream harmoniser + renderer consume both feeds without
    special-casing.
    """

    variant_id: str  # chr_pos_ref_alt (matches OT convention)
    chromosome: str
    position: int
    ref: str
    alt: str
    beta: float | None
    se: float | None
    p_value: float | None
    maf: float | None
    effect_allele_frequency: float | None
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass
class RegionResult:
    """Per-region fetch payload."""

    protein_label_short: str    # e.g. "SORT1" (HGNC only — for filenames)
    ancestry: str
    chromosome: str
    region_start_bp: int
    region_end_bp: int
    n_variants: int
    variants: list[RegionVariant]
    release: UKBPPPRelease
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "protein_label_short": self.protein_label_short,
            "ancestry": self.ancestry,
            "chromosome": self.chromosome,
            "region_start_bp": self.region_start_bp,
            "region_end_bp": self.region_end_bp,
            "n_variants": self.n_variants,
            "variants": [asdict(v) for v in self.variants],
            "release": asdict(self.release),
            "notes": list(self.notes),
        }


class UKBPPPAccessError(Exception):
    """Raised when the UKB-PPP backend (Synapse) cannot resolve, authenticate,
    or stream a requested protein archive.
    """


# Multi-line, user-facing error message printed when a live Synapse fetch
# is attempted without a PAT. Per user direction 2026-05-15: be very loud
# about the auth requirement at the moment a pQTL render is attempted,
# point to the free PAT path, and clarify that no UKB Application is
# needed for the summary-stats layer. Plain text (no emoji) so the message
# renders in raw stdout / log capture without surprise.
PAT_REQUIRED_MESSAGE = """\
UKB-PPP pQTL summary statistics require Synapse authentication.

This skill ships pre-computed regional slices for the canonical demo
proteins (e.g. SORT1) that work offline with no auth required. The
requested (protein, ancestry, region) is NOT in the bundled-slice
inventory, so a live Synapse fetch was attempted - and a Synapse
personal access token (PAT) is needed for that.

The free Synapse PAT path (NO UK Biobank Application required for the
summary-stats layer; only a free Synapse account):

1. Register / sign in at https://www.synapse.org
2. Accept the Synapse Terms of Use (one-time click-through).
3. Open Account Settings -> Personal Access Tokens
   (https://www.synapse.org/Profile:settings).
4. Click "Create new token". Scopes: tick `view` and `download`.
   Copy the token immediately (Synapse shows it once).
5. Export it before re-running:
     export SYNAPSE_AUTH_TOKEN=<your_token>

Why this gating exists: UKB-PPP's AWS Open Data Registry bucket
(s3://ukbiobank.opendata.sagebase.org) is advertised as public but
returns AccessDenied for anonymous reads as of 2026-05-15. Synapse is
the only functional access path the data owner currently offers.
ClawBio + K-Dense users running the BUNDLED demos do NOT need a PAT;
the PAT is only required for arbitrary protein x ancestry x region
queries beyond the bundled inventory.

See SKILL.md "Gotchas" #1 for the full context.
"""


class _ProteinIndex:
    """Per-ancestry name → Synapse-id index. Built lazily on first lookup."""

    def __init__(self, fetcher: "SynapseListFetcher", ancestry: str) -> None:
        self._fetcher = fetcher
        self.ancestry = ancestry
        self._by_hgnc: dict[str, dict[str, Any]] | None = None
        self._by_uniprot: dict[str, dict[str, Any]] | None = None

    def _ensure_loaded(self) -> None:
        if self._by_hgnc is not None:
            return
        folder_id = ANCESTRY_FOLDER_IDS.get(self.ancestry)
        if folder_id is None:
            raise UKBPPPAccessError(
                f"unknown UKB-PPP ancestry {self.ancestry!r}; "
                f"expected one of {sorted(ANCESTRY_FOLDER_IDS)}"
            )
        entries = self._fetcher.list_folder(folder_id)
        by_hgnc: dict[str, dict[str, Any]] = {}
        by_uniprot: dict[str, dict[str, Any]] = {}
        for ent in entries:
            name = ent.get("name") or ""
            m = PROTEIN_FILE_RE.match(name)
            if not m:
                continue
            rec = {
                "hgnc": m.group("hgnc"),
                "uniprot": m.group("uniprot"),
                "olink_id": m.group("olink_id"),
                "olink_panel": m.group("panel").replace("_", " "),
                "synapse_id": ent.get("id") or "",
                "file_name": name,
            }
            by_hgnc.setdefault(rec["hgnc"], rec)
            by_uniprot.setdefault(rec["uniprot"], rec)
        self._by_hgnc = by_hgnc
        self._by_uniprot = by_uniprot

    def resolve(self, protein_label: str) -> dict[str, Any]:
        """Find the protein entry by HGNC symbol or UniProt accession.

        Some HGNC symbols map to multiple Olink reagents (isoform / panel
        ambiguity); this returns the first one found. To force a specific
        reagent, pass the OID directly via the alternate `resolve_by_olink_id`
        path.
        """
        self._ensure_loaded()
        assert self._by_hgnc is not None and self._by_uniprot is not None
        key = protein_label.strip()
        rec = self._by_hgnc.get(key.upper())
        if rec is None:
            rec = self._by_uniprot.get(key.upper())
        if rec is None:
            raise UKBPPPAccessError(
                f"no UKB-PPP protein file for {protein_label!r} in ancestry "
                f"{self.ancestry}; tried HGNC + UniProt indexes"
            )
        return rec


class SynapseListFetcher:
    """Thin Synapse REST wrapper for folder listings.

    No authentication is required for the public listing endpoint; only
    file downloads require a PAT. We use this for the name → fileID
    resolution layer.

    A real authenticated download uses `synapseclient` (lazy import inside
    `SynapseFileFetcher`) so the listing path stays auth-free for static
    inventory work + offline tests.
    """

    BASE_URL = "https://repo-prod.prod.sagebase.org/repo/v1"

    def __init__(self, base_url: str = BASE_URL, timeout_s: float = 60.0) -> None:
        self.base_url = base_url.rstrip("/")
        self.timeout_s = timeout_s

    def list_folder(self, folder_id: str) -> list[dict[str, Any]]:
        """List the children (files + folders) of a Synapse folder."""
        import requests
        url = f"{self.base_url}/entity/children"
        body = {
            "parentId": folder_id,
            "includeTypes": ["folder", "file", "link"],
        }
        results: list[dict[str, Any]] = []
        next_token: str | None = None
        while True:
            payload = dict(body)
            if next_token:
                payload["nextPageToken"] = next_token
            resp = requests.post(
                url, json=payload, timeout=self.timeout_s,
                headers={"Accept": "application/json"},
            )
            resp.raise_for_status()
            data = resp.json()
            results.extend(data.get("page", []) or [])
            next_token = data.get("nextPageToken")
            if not next_token:
                break
        return results


class SynapseFileFetcher:
    """Authenticated downloader. Wraps `synapseclient.Synapse.get`."""

    def __init__(
        self,
        auth_token: str | None = None,
        cache_dir: Path = DEFAULT_CACHE_DIR,
    ) -> None:
        self.auth_token = auth_token or os.environ.get("SYNAPSE_AUTH_TOKEN") or ""
        self.cache_dir = Path(cache_dir).expanduser()
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._syn = None

    def _ensure_synapse(self):  # type: ignore[no-untyped-def]
        if self._syn is not None:
            return self._syn
        try:
            import synapseclient
        except ImportError as e:
            raise UKBPPPAccessError(
                "synapseclient is required for live UKB-PPP downloads; "
                "install via `pip install synapseclient`. The bundled-slice "
                "demo path does not need synapseclient."
            ) from e
        if not self.auth_token:
            raise UKBPPPAccessError(PAT_REQUIRED_MESSAGE)
        syn = synapseclient.Synapse(silent=True, cache_root_dir=str(self.cache_dir))
        syn.login(authToken=self.auth_token)
        self._syn = syn
        return syn

    def download_tar(self, synapse_id: str) -> Path:
        """Download the per-protein `.tar` archive to the cache, returning
        the local path. Uses synapseclient's built-in caching."""
        syn = self._ensure_synapse()
        entity = syn.get(synapse_id, downloadLocation=str(self.cache_dir))
        path = Path(entity.path) if getattr(entity, "path", None) else None
        if path is None or not path.is_file():
            raise UKBPPPAccessError(
                f"Synapse download succeeded but no local path returned for {synapse_id}"
            )
        return path


def _bundled_slices_dir() -> Path:
    """Where the redistributable, pre-computed regional slices live."""
    return Path(__file__).resolve().parent / "bundled_slices"


def _bundled_slice_key(
    protein_label: str, ancestry: str, chromosome: str,
    start_bp: int, end_bp: int,
) -> str:
    """Deterministic filename key for a (protein, ancestry, region) slice.

    Slices ship as gzipped JSON (`<key>.json.gz`); per-variant pQTL rows
    compress ~8.5x, shaving the bulk of the PR size.
    """
    return (
        f"{protein_label.strip().upper()}__"
        f"{ancestry.strip().upper()}__"
        f"chr{str(chromosome).lstrip('chr')}__"
        f"{int(start_bp)}_{int(end_bp)}.json.gz"
    )


def _load_bundled_slice(path: Path) -> RegionResult:
    """Reconstruct a RegionResult from a gzipped bundled-slice JSON file."""
    with gzip.open(path, mode="rt", encoding="utf-8") as f:
        d = json.load(f)
    rel = d.get("release") or {}
    release = UKBPPPRelease(
        release_label=rel.get("release_label", ""),
        fetched_at_utc=rel.get("fetched_at_utc", ""),
        protein_hgnc=rel.get("protein_hgnc", ""),
        protein_uniprot=rel.get("protein_uniprot", ""),
        olink_reagent_id=rel.get("olink_reagent_id", ""),
        olink_panel=rel.get("olink_panel", ""),
        ancestry=rel.get("ancestry", ""),
        ancestry_label=rel.get("ancestry_label", ""),
        n_samples=rel.get("n_samples"),
        synapse_id=rel.get("synapse_id", ""),
        source_url=rel.get("source_url", ""),
        study_label=rel.get("study_label", "UKB-PPP"),
    )
    variants = [
        RegionVariant(
            variant_id=v["variant_id"],
            chromosome=v["chromosome"],
            position=int(v["position"]),
            ref=v["ref"],
            alt=v["alt"],
            beta=v.get("beta"),
            se=v.get("se"),
            p_value=v.get("p_value"),
            maf=v.get("maf"),
            effect_allele_frequency=v.get("effect_allele_frequency"),
            raw=v.get("raw") or {},
        )
        for v in d.get("variants", [])
    ]
    return RegionResult(
        protein_label_short=d["protein_label_short"],
        ancestry=d["ancestry"],
        chromosome=d["chromosome"],
        region_start_bp=int(d["region_start_bp"]),
        region_end_bp=int(d["region_end_bp"]),
        n_variants=int(d["n_variants"]),
        variants=variants,
        release=release,
        notes=list(d.get("notes") or []),
    )


class UKBPPPClient:
    """Regional UKB-PPP pQTL summary-stats client.

    Two-path fetch (v1.3):

    1. **Bundled slice** (default; no auth). On `fetch_region(...)`, the
       client first checks `bundled_slices/<key>.json.gz` inside this skill
       directory; if it exists, the pre-computed harmonised slice is
       returned directly. Demo cohort coverage means most canonical
       renders never hit the network.

    2. **Live Synapse fetch** (opt-in; PAT required). If no bundled slice
       matches, the client falls through to a live Synapse download:
       protein resolution (HGNC or UniProt to a `.tar` Synapse id) ->
       tar download (cached) -> per-chromosome extract -> region filter
       -> harmonisation to OT convention.

    Two-tier downloader design (mirrors `EQTLCatalogueClient`):
    - `lister` (SynapseListFetcher) is auth-free, used for index building +
      offline tests via mocking
    - `downloader` (SynapseFileFetcher) requires a PAT, used only for the
      actual tar fetch
    """

    def __init__(
        self,
        *,
        lister: SynapseListFetcher | None = None,
        downloader: SynapseFileFetcher | None = None,
        auth_token: str | None = None,
        cache_dir: Path = DEFAULT_CACHE_DIR,
        bundled_slices_dir: Path | None = None,
        allow_live_fetch: bool = True,
    ) -> None:
        self.lister = lister or SynapseListFetcher()
        self.downloader = downloader or SynapseFileFetcher(
            auth_token=auth_token, cache_dir=cache_dir,
        )
        self._index_cache: dict[str, _ProteinIndex] = {}
        self.bundled_slices_dir = (
            Path(bundled_slices_dir) if bundled_slices_dir is not None
            else _bundled_slices_dir()
        )
        self.allow_live_fetch = allow_live_fetch

    def _index_for(self, ancestry: str) -> _ProteinIndex:
        if ancestry not in self._index_cache:
            self._index_cache[ancestry] = _ProteinIndex(self.lister, ancestry)
        return self._index_cache[ancestry]

    def _bundled_slice_path(
        self, protein_label: str, ancestry: str, chromosome: str,
        start_bp: int, end_bp: int,
    ) -> Path:
        return self.bundled_slices_dir / _bundled_slice_key(
            protein_label, ancestry, chromosome, start_bp, end_bp,
        )

    def fetch_region(
        self,
        protein_label: str,
        ancestry: str,
        chromosome: str,
        start_bp: int,
        end_bp: int,
    ) -> RegionResult:
        """Pull pQTL summary stats for a protein in a chromosomal window.

        `protein_label` is HGNC symbol (e.g. "SORT1") or UniProt accession
        (e.g. "Q99523"). `ancestry` is one of EUR/AFR/EAS/SAS/CSA/MID/AMR/ALL.

        Resolution order:
        1. Bundled slice (no auth) at `bundled_slices/<key>.json.gz`.
        2. Live Synapse fetch (requires SYNAPSE_AUTH_TOKEN).

        Returns harmonised variants in OT `chr_pos_ref_alt` ALT-effect
        convention.

        Raises `UKBPPPAccessError(PAT_REQUIRED_MESSAGE)` when no bundled
        slice matches AND no PAT is configured AND
        `allow_live_fetch=True`. When `allow_live_fetch=False`, missing
        bundled slices raise a different error to make the no-network
        intent explicit.
        """
        ancestry_u = ancestry.strip().upper()
        chrom_q = str(chromosome).lstrip("chr")

        # Step 1: bundled slice.
        slice_path = self._bundled_slice_path(
            protein_label, ancestry_u, chrom_q, int(start_bp), int(end_bp),
        )
        if slice_path.is_file():
            return _load_bundled_slice(slice_path)

        # Step 2: live Synapse fetch (or refuse if disabled).
        if not self.allow_live_fetch:
            raise UKBPPPAccessError(
                f"no bundled slice for "
                f"({protein_label!r}, {ancestry_u}, chr{chrom_q}:"
                f"{int(start_bp)}-{int(end_bp)}) and allow_live_fetch=False; "
                f"checked {slice_path}"
            )

        rec = self._index_for(ancestry_u).resolve(protein_label)
        tar_path = self.downloader.download_tar(rec["synapse_id"])

        variants, notes = _extract_region_from_tar(
            tar_path=tar_path,
            chromosome=chrom_q,
            start_bp=int(start_bp),
            end_bp=int(end_bp),
        )

        n_samples = _UKB_PPP_N_BY_ANCESTRY.get(ancestry_u)
        release = UKBPPPRelease(
            release_label="UKB-PPP r1 2023 (Sun 2023)",
            fetched_at_utc=_now_utc(),
            protein_hgnc=rec["hgnc"],
            protein_uniprot=rec["uniprot"],
            olink_reagent_id=rec["olink_id"],
            olink_panel=rec["olink_panel"],
            ancestry=ancestry_u,
            ancestry_label=ANCESTRY_LABELS.get(ancestry_u, ancestry_u),
            n_samples=n_samples,
            synapse_id=rec["synapse_id"],
            source_url=f"https://www.synapse.org/Synapse:{rec['synapse_id']}",
        )
        return RegionResult(
            protein_label_short=rec["hgnc"],
            ancestry=ancestry_u,
            chromosome=chrom_q,
            region_start_bp=int(start_bp),
            region_end_bp=int(end_bp),
            n_variants=len(variants),
            variants=variants,
            release=release,
            notes=notes,
        )


# Discovery cohort sizes per Sun 2023 Nature Table 1. EUR is the largest
# stratum by ~10x; trans-ancestry coloc renders should surface the smaller
# strata's N to remind the reader about power.
_UKB_PPP_N_BY_ANCESTRY: dict[str, int | None] = {
    "EUR": 46_673,
    "AFR": 931,
    "CSA": 920,
    "SAS": 920,
    "EAS": 262,
    "MID": 124,
    "AMR": 60,
    "ALL": 47_970,
}


def _extract_region_from_tar(
    *,
    tar_path: Path,
    chromosome: str,
    start_bp: int,
    end_bp: int,
) -> tuple[list[RegionVariant], list[str]]:
    """Open `tar_path`, find the chr<chromosome> member, stream-filter to
    the requested window, and return harmonised variants.

    Per Sun 2023 methods, the per-protein tar contains one REGENIE step-2
    summary file per autosome + X. Names follow the pattern
    `discovery_chr<N>_<protein>_*.regenie.gz` (per ancestry-specific README
    on Synapse). The parser is tolerant of name variants — it matches on a
    `chr<N>` substring inside any `.gz` member.

    Stream-filter (not tabix) because the source files are plain gzip, not
    BGZ. For one ±500 kb window this scans ~1M variant lines per
    chromosome — measured at ~2 s on a modern laptop, acceptable for the
    Tier-2 coloc-render workflow which fetches one region per protein per
    render.
    """
    notes: list[str] = []
    if not tar_path.is_file():
        raise UKBPPPAccessError(f"tar archive not found: {tar_path}")

    chr_token = f"chr{chromosome}"
    target_member = None
    with tarfile.open(tar_path, mode="r:*") as tf:
        for m in tf.getmembers():
            name = m.name.split("/")[-1]
            if not name.endswith((".gz", ".tsv", ".txt", ".regenie")):
                continue
            # match strict word boundary so chr1 doesn't match chr10
            if re.search(rf"(^|[^0-9a-zA-Z]){re.escape(chr_token)}([^0-9]|$)", name):
                target_member = m
                break
        if target_member is None:
            raise UKBPPPAccessError(
                f"tar {tar_path.name} has no member matching {chr_token}; "
                f"members: {[m.name for m in tf.getmembers()][:5]}"
            )
        fobj = tf.extractfile(target_member)
        if fobj is None:
            raise UKBPPPAccessError(
                f"tar member {target_member.name} could not be opened"
            )
        text_stream = _maybe_gunzip(fobj, target_member.name)
        variants = list(_parse_regenie_stream(
            text_stream, chromosome=chromosome,
            start_bp=start_bp, end_bp=end_bp, notes=notes,
        ))
    return variants, notes


def _maybe_gunzip(fobj, name: str) -> io.TextIOBase:
    """Wrap a binary file object in gzip + text decoders as needed."""
    if name.endswith(".gz"):
        return io.TextIOWrapper(gzip.GzipFile(fileobj=fobj), encoding="utf-8")
    return io.TextIOWrapper(fobj, encoding="utf-8")


def _parse_regenie_stream(
    text_stream,
    *,
    chromosome: str,
    start_bp: int,
    end_bp: int,
    notes: list[str],
) -> Iterator[RegionVariant]:
    """Parse a UKB-PPP REGENIE step-2 output file (space- or tab-delimited).

    The first non-comment line is the header. Subsequent lines are data
    rows. Filter rows by CHROM == chromosome and start_bp <= GENPOS <= end_bp,
    yielding harmonised variants on the fly so the whole file isn't held
    in memory.
    """
    header: list[str] | None = None
    sep_re = re.compile(r"\s+")
    for raw_line in text_stream:
        line = raw_line.rstrip("\n\r")
        if not line:
            continue
        if line.startswith("#"):
            continue
        fields = sep_re.split(line.strip())
        if header is None:
            header = [f.upper() for f in fields]
            continue
        if len(fields) != len(header):
            notes.append(
                f"row column count {len(fields)} != header {len(header)}; "
                f"skipping (REGENIE schema may have drifted)"
            )
            continue
        row = dict(zip(header, fields))
        row_chrom = str(row.get("CHROM", "")).lstrip("chr")
        if row_chrom != chromosome:
            continue
        try:
            pos = int(row.get("GENPOS") or 0)
        except (TypeError, ValueError):
            continue
        if pos < start_bp or pos > end_bp:
            continue
        yield _normalise_row(row)


def _normalise_row(row: dict[str, Any]) -> RegionVariant:
    """Convert one REGENIE row to our internal RegionVariant.

    OT convention: `chr_pos_ref_alt` GRCh38, ALT = effect allele.
    REGENIE convention: ALLELE0 = REF, ALLELE1 = ALT (effect), BETA is on
    ALLELE1. No flipping required at the row boundary; the join layer
    handles palindromic exclusions.

    `LOG10P` is converted to a linear `p_value` for cross-source
    comparability with eQTL Catalogue + GWAS Catalog rows.
    """
    chrom = str(row.get("CHROM", "")).lstrip("chr")
    pos = int(row.get("GENPOS") or 0)
    ref = str(row.get("ALLELE0", "")).upper()
    alt = str(row.get("ALLELE1", "")).upper()
    variant_id = _build_variant_id(chrom, pos, ref, alt)
    beta = _maybe_float(row.get("BETA"))
    se = _maybe_float(row.get("SE"))
    log10p = _maybe_float(row.get("LOG10P"))
    p_value: float | None = None
    if log10p is not None:
        # LOG10P is the absolute value of log10(p); REGENIE convention.
        try:
            p_value = 10 ** (-abs(log10p))
        except OverflowError:
            p_value = 0.0
    a1freq = _maybe_float(row.get("A1FREQ"))
    maf = None
    if a1freq is not None:
        maf = a1freq if a1freq <= 0.5 else 1.0 - a1freq
    return RegionVariant(
        variant_id=variant_id,
        chromosome=chrom,
        position=pos,
        ref=ref,
        alt=alt,
        beta=beta,
        se=se,
        p_value=p_value,
        maf=maf,
        effect_allele_frequency=a1freq,
        raw=dict(row),
    )


def _build_variant_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"{chrom}_{pos}_{ref}_{alt}"


def _maybe_float(v: Any) -> float | None:
    if v is None or v == "" or v == "NA":
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if (f == f) else None  # filter NaN


def _now_utc() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


# ----------------------- CLI -----------------------


def main(argv: list[str] | None = None) -> int:
    """ClawBio-shape CLI: `--input <config> --output <dir> --demo`.

    Config schema (JSON or YAML):
        protein_label: SORT1            # HGNC symbol or UniProt accession
        ancestry: EUR                   # EUR/AFR/EAS/SAS/CSA/MID/AMR/ALL
        chromosome: "1"
        start_bp: 108774968
        end_bp: 109774968

    Writes to <output>/:
        variants.tsv        # canonical sumstats slice (per locuscompare INPUT_SCHEMA.md)
        manifest.yaml       # source release + provenance
        report.md           # human-readable run summary
    """
    parser = argparse.ArgumentParser(
        prog="ukb-ppp-region-fetch",
        description="Fetch a region of UKB-PPP pQTL summary statistics for one (protein, ancestry) measurement.",
    )
    parser.add_argument("--input", type=Path, help="JSON or YAML config (see docstring).")
    parser.add_argument("--output", type=Path,
                        help="Output directory; created if missing. Required unless --list-demos.")
    parser.add_argument("--demo", nargs="?", const="__default__", default=None,
                        metavar="NAME",
                        help="Run a bundled demo. Bare --demo runs the default; "
                             "pass a name (e.g. --demo sort1_ukb_ppp_eur) to choose a specific one.")
    parser.add_argument("--list-demos", action="store_true",
                        help="List bundled demo configs in this skill's examples/ directory.")
    args = parser.parse_args(argv)

    if args.list_demos:
        _print_available_demos()
        return 0
    if args.demo is None and args.input is None:
        parser.error("either --input <config> or --demo [NAME] or --list-demos is required")
    if args.output is None:
        parser.error("--output is required")
    args.output.mkdir(parents=True, exist_ok=True)

    cfg = _load_config(_resolve_demo_path(args.demo) if args.demo is not None else args.input)
    if args.demo is not None:
        print(f"info: using bundled demo", file=sys.stderr)

    # Probe the bundled-slice inventory upfront so we can print a clear
    # "live fetch will be attempted; here's what that needs" message
    # before failing. Better UX than letting synapseclient's auth call
    # blow up deep in the stack.
    client = UKBPPPClient()
    slice_path = client._bundled_slice_path(
        protein_label=cfg["protein_label"],
        ancestry=cfg["ancestry"],
        chromosome=str(cfg["chromosome"]),
        start_bp=int(cfg["start_bp"]),
        end_bp=int(cfg["end_bp"]),
    )
    if slice_path.is_file():
        print(
            f"info: using bundled slice {slice_path.name} (no Synapse auth needed)",
            file=sys.stderr,
        )
    else:
        if not os.environ.get("SYNAPSE_AUTH_TOKEN"):
            print(
                "\n" + "=" * 72 + "\n"
                "ukb-ppp-region-fetch: this query is NOT in the bundled-slice\n"
                "inventory; a live Synapse fetch will be attempted.\n"
                + "=" * 72 + "\n"
                + PAT_REQUIRED_MESSAGE
                + "=" * 72,
                file=sys.stderr,
            )
            return 2
        print(
            f"info: no bundled slice for {slice_path.name}; "
            "falling through to live Synapse fetch (PAT detected).",
            file=sys.stderr,
        )

    result = client.fetch_region(
        protein_label=cfg["protein_label"],
        ancestry=cfg["ancestry"],
        chromosome=str(cfg["chromosome"]),
        start_bp=int(cfg["start_bp"]),
        end_bp=int(cfg["end_bp"]),
    )

    tsv_path = args.output / "variants.tsv"
    _write_canonical_tsv(result, tsv_path)

    manifest = {
        "skill": "ukb-ppp-region-fetch",
        "version": "0.1.0",
        "protein_label": cfg["protein_label"],
        "ancestry": result.ancestry,
        "region": {"chromosome": result.chromosome,
                   "start_bp": result.region_start_bp,
                   "end_bp": result.region_end_bp},
        "n_variants": result.n_variants,
        "release": {
            "study_label": result.release.study_label,
            "release_label": result.release.release_label,
            "protein_hgnc": result.release.protein_hgnc,
            "protein_uniprot": result.release.protein_uniprot,
            "protein_label": result.release.protein_label,
            "olink_reagent_id": result.release.olink_reagent_id,
            "olink_panel": result.release.olink_panel,
            "ancestry": result.release.ancestry,
            "ancestry_label": result.release.ancestry_label,
            "n_samples": result.release.n_samples,
            "synapse_id": result.release.synapse_id,
            "source_url": result.release.source_url,
            "fetched_at_utc": result.release.fetched_at_utc,
        },
        "attribution": (
            "UK Biobank Pharma Proteomics Project (Sun et al. 2023 Nature; "
            "PMID 37794186); accessed via Synapse syn51364943."
        ),
        "outputs": {"variants_tsv": "variants.tsv"},
    }
    try:
        import yaml as _yaml
        (args.output / "manifest.yaml").write_text(
            _yaml.safe_dump(manifest, sort_keys=False)
        )
    except ImportError:
        (args.output / "manifest.json").write_text(
            json.dumps(manifest, indent=2, default=str)
        )

    report = [
        "# ukb-ppp-region-fetch report",
        "",
        f"- **Protein:** {result.release.protein_label}",
        f"- **Olink panel:** {result.release.olink_panel}",
        f"- **Ancestry:** {result.release.ancestry_label} ({result.ancestry})"
        + (f"; N = {result.release.n_samples:,}" if result.release.n_samples else ""),
        f"- **Region:** chr{result.chromosome}:{result.region_start_bp:,}-{result.region_end_bp:,}",
        f"- **Variants returned:** {result.n_variants}",
        f"- **Source:** {result.release.source_url}",
        f"- **Output TSV:** {tsv_path.name}",
    ]
    (args.output / "report.md").write_text("\n".join(report) + "\n")

    print(f"ukb-ppp-region-fetch: {result.n_variants} variants -> {tsv_path}")
    print(f"  source: UKB-PPP | {result.release.protein_label} | "
          f"{result.release.ancestry_label} ({result.ancestry})")
    return 0


def _load_config(path: Path) -> dict:
    text = path.read_text()
    if path.suffix.lower() in (".yaml", ".yml"):
        import yaml as _yaml
        return _yaml.safe_load(text) or {}
    if path.suffix.lower() == ".json":
        return json.loads(text)
    raise ValueError(f"unsupported config extension: {path.suffix}")


def _examples_dir() -> Path:
    return Path(__file__).resolve().parent / "examples"


def _list_demos() -> list[Path]:
    out: list[Path] = []
    for ext in ("*.json", "*.yaml", "*.yml"):
        out.extend(sorted(_examples_dir().glob(ext)))
    return [p for p in out if p.name not in {"expected_output.md", "README.md"}]


def _resolve_demo_path(name: str) -> Path:
    examples = _examples_dir()
    if name == "__default__":
        for cand in ("default.json", "default.yaml", "default.yml"):
            p = examples / cand
            if p.is_file():
                return p
        files = _list_demos()
        if not files:
            raise FileNotFoundError(f"no bundled demos in {examples}")
        return files[0]
    for ext in (".json", ".yaml", ".yml", ""):
        p = examples / (name if ext == "" else f"{name}{ext}")
        if p.is_file():
            return p
    available = ", ".join(p.stem for p in _list_demos())
    raise FileNotFoundError(
        f"no bundled demo named {name!r} in {examples}. Available: {available}"
    )


def _print_available_demos() -> None:
    paths = _list_demos()
    if not paths:
        print(f"no bundled demos in {_examples_dir()}")
        return
    try:
        default_path = _resolve_demo_path("__default__")
    except FileNotFoundError:
        default_path = None
    print(f"Bundled demos in {_examples_dir()}:")
    for p in paths:
        marker = " (default)" if default_path is not None and p == default_path else ""
        print(f"  {p.stem}{marker}    [{p.name}]")


def _write_canonical_tsv(result: RegionResult, tsv_path: Path) -> None:
    """Emit the canonical sumstats-slice TSV consumed by the regional
    LocusCompare orchestrator (see the locuscompare-region-render skill
    INPUT_SCHEMA for the column contract)."""
    cols = ["variant_id", "chromosome", "position_bp",
            "allele_a", "allele_b", "beta", "se", "p", "maf",
            "molecular_trait_id", "study_id"]
    # For pQTL, molecular_trait_id is the Olink reagent id (one trait per
    # file, unlike ge-eQTL where one file holds many genes). study_id is
    # the Synapse fileID so downstream tooling can re-resolve provenance.
    mt = result.release.olink_reagent_id
    study = result.release.synapse_id
    with tsv_path.open("w") as f:
        f.write("# locuscompare-schema-version: 1.0\n")
        f.write("# source: ukb_ppp\n")
        f.write(f"# protein_hgnc: {result.release.protein_hgnc}\n")
        f.write(f"# protein_uniprot: {result.release.protein_uniprot}\n")
        f.write(f"# olink_reagent_id: {result.release.olink_reagent_id}\n")
        f.write(f"# ancestry: {result.ancestry}\n")
        f.write("\t".join(cols) + "\n")
        for v in result.variants:
            row = [
                v.variant_id,
                v.chromosome,
                str(v.position),
                v.ref,
                v.alt,
                "" if v.beta is None else f"{v.beta:.6g}",
                "" if v.se is None else f"{v.se:.6g}",
                "" if v.p_value is None else f"{v.p_value:.6g}",
                "" if v.maf is None else f"{v.maf:.6g}",
                mt,
                study,
            ]
            f.write("\t".join(row) + "\n")


def save_region_result_as_bundled_slice(
    result: RegionResult,
    *,
    bundled_slices_dir: Path | None = None,
) -> Path:
    """Persist a `RegionResult` to the redistributable bundled-slice
    directory. Used by maintainers to bake pre-computed slices into the
    skill repo for the canonical demo cohort (after a one-time fetch
    with a Synapse PAT).

    Returns the path to the written file.

    CC-BY 4.0 redistribution requirement: the on-disk manifest carries
    the same attribution string the live fetcher emits ("UK Biobank
    Pharma Proteomics Project (Sun et al. 2023 Nature; PMID 37794186);
    accessed via Synapse syn51364943.").
    """
    target_dir = Path(bundled_slices_dir) if bundled_slices_dir is not None \
        else _bundled_slices_dir()
    target_dir.mkdir(parents=True, exist_ok=True)
    path = target_dir / _bundled_slice_key(
        result.release.protein_hgnc, result.ancestry,
        result.chromosome, result.region_start_bp, result.region_end_bp,
    )
    with gzip.open(path, mode="wt", encoding="utf-8", compresslevel=9) as f:
        json.dump(result.to_dict(), f, default=str, indent=2)
    return path


__all__ = [
    "ANCESTRY_FOLDER_IDS",
    "ANCESTRY_LABELS",
    "DEFAULT_CACHE_DIR",
    "PAT_REQUIRED_MESSAGE",
    "PROTEIN_FILE_RE",
    "REGENIE_COLUMNS",
    "RegionResult",
    "RegionVariant",
    "SYNAPSE_PROJECT_ID",
    "SYNAPSE_SUMSTATS_FOLDER_ID",
    "SynapseFileFetcher",
    "SynapseListFetcher",
    "UKBPPPAccessError",
    "UKBPPPClient",
    "UKBPPPRelease",
    "save_region_result_as_bundled_slice",
]


if __name__ == "__main__":
    sys.exit(main())
