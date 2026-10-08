"""GWAS Catalog harmonised per-region summary-stats execution skill.

Fetches per-variant outcome-side summary statistics within a genomic window
for a given GWAS Catalog accession. Substrate for regional LocusCompare
panels and downstream coloc / fine-mapping / Mendelian randomisation.

Source layout:
- Per-study harmonised TSV.bgz at
  https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/<GCST_PREFIX>/<GCST>/harmonised/
  with a sibling .tbi tabix index.
- Forward-strand-aligned to GRCh38 (the harmoniser's `hm_*` columns).

We range-fetch via pysam's tabix support, which streams only the bytes
covering the requested chromosome+window. pysam is MIT-licensed and is the
bioinformatics-standard tabix client.

License: open access per Sollis 2023. Per-study attribution = the original GWAS
publication. See allowlist key `gwas_catalog_summary_stats`.

Variant id mapping: harmonised columns provide `hm_chrom`, `hm_pos`, `hm_other_allele`,
`hm_effect_allele` (the effect allele after forward-strand alignment), `hm_beta`,
`hm_odds_ratio`, `standard_error`, `p_value`, `hm_effect_allele_frequency`. We
emit `variant_id` in OT's `chr_pos_ref_alt` form using `hm_other_allele` as ref
and `hm_effect_allele` as alt; downstream renderer joins on this id.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import re
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any

import requests

DEFAULT_FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
DEFAULT_TIMEOUT_S = 120.0


@dataclass
class GWASCatalogRelease:
    """Records the GWAS Catalog mirror release pinned per fetch (for the manifest)."""

    accession: str
    harmonised_url: str
    fetched_at_utc: str
    harmoniser_version: str | None = None  # if discoverable from sibling metadata


@dataclass
class RegionVariant:
    """One harmonised variant row, GRCh38 ALT-effect convention (matches OT)."""

    variant_id: str  # chr_pos_ref_alt (matches OT join key)
    chromosome: str
    position: int
    ref: str  # = hm_other_allele
    alt: str  # = hm_effect_allele
    beta: float | None
    se: float | None
    p_value: float | None
    odds_ratio: float | None
    effect_allele_frequency: float | None
    # Odds-ratio 95% CI bounds, present for binary-trait studies that publish them.
    # Passed through as published so a consumer can derive beta = ln(OR) and a
    # CI-based SE; the fetcher does not derive them itself.
    ci_lower: float | None = None
    ci_upper: float | None = None
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass
class RegionResult:
    accession: str
    chromosome: str
    region_start_bp: int
    region_end_bp: int
    n_variants: int
    variants: list[RegionVariant]
    release: GWASCatalogRelease
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "accession": self.accession,
            "chromosome": self.chromosome,
            "region_start_bp": self.region_start_bp,
            "region_end_bp": self.region_end_bp,
            "n_variants": self.n_variants,
            "variants": [asdict(v) for v in self.variants],
            "release": asdict(self.release),
            "notes": list(self.notes),
        }


class GWASCatalogFetchError(Exception):
    """Raised when the harmonised file cannot be fetched or parsed."""


_GCST_RE = re.compile(r"^GCST\d+$")


def gcst_url_base(accession: str, ftp_base: str = DEFAULT_FTP_BASE) -> str:
    """Compute the harmonised dir URL.

    GWAS Catalog organises by GCST prefix bucket of 1000:
    GCST90475990 → GCST90475001-GCST90476000/GCST90475990/harmonised/
    """
    if not _GCST_RE.fullmatch(accession):
        raise ValueError(f"invalid GCST accession (expected ^GCST\\d+$): {accession}")
    digits = accession[4:]
    if not digits:
        raise ValueError(f"unparseable GCST accession: {accession}")
    n = int(digits)
    bucket_lo = ((n - 1) // 1000) * 1000 + 1
    bucket_hi = bucket_lo + 999
    bucket = f"GCST{bucket_lo:>06d}-GCST{bucket_hi:>06d}".replace(" ", "0")
    return f"{ftp_base.rstrip('/')}/{bucket}/{accession}/harmonised"


def harmonised_file_url(
    accession: str,
    ftp_base: str = DEFAULT_FTP_BASE,
    session: requests.Session | None = None,
    timeout_s: float = 30.0,
) -> str:
    """Resolve the harmonised TSV URL for a GCST accession.

    GWAS Catalog uses two file-naming conventions for harmonised files:
      1. Simple: `<GCST>.h.tsv.gz`
      2. Prefixed: `<PMID>-<GCST>-<EFO_id>.h.tsv.gz`
    Try the simple form first via HEAD; fall back to listing the directory
    and picking the first `*.h.tsv.gz` file.
    """
    base = gcst_url_base(accession, ftp_base=ftp_base)
    sess = session or requests.Session()
    simple = f"{base}/{accession}.h.tsv.gz"
    try:
        r = sess.head(simple, timeout=timeout_s, allow_redirects=True)
        if r.status_code == 200:
            return simple
    except requests.RequestException:
        pass
    # Fallback: list the harmonised/ directory and find any .h.tsv.gz file.
    try:
        r = sess.get(f"{base}/", timeout=timeout_s)
        if r.status_code == 200:
            import re
            matches = re.findall(r'href="([^"]+\.h\.tsv\.gz)"', r.text)
            if matches:
                return f"{base}/{matches[0]}"
    except requests.RequestException as e:
        raise GWASCatalogFetchError(
            f"could not list harmonised dir for {accession}: {e!s}"
        ) from e
    raise GWASCatalogFetchError(
        f"no harmonised .h.tsv.gz found in {base}/ for accession {accession}"
    )


class GWASCatalogClient:
    """Range-fetcher for GWAS Catalog harmonised summary stats.

    Uses pysam.TabixFile when available for efficient remote range queries.
    Falls back to streaming the full file (slower) when tabix index is unavailable.
    """

    def __init__(
        self,
        ftp_base: str = DEFAULT_FTP_BASE,
        timeout_s: float = DEFAULT_TIMEOUT_S,
        session: requests.Session | None = None,
    ) -> None:
        self.ftp_base = ftp_base.rstrip("/")
        self.timeout_s = timeout_s
        self.session = session or requests.Session()

    def fetch_region(
        self,
        accession: str,
        chromosome: str,
        start_bp: int,
        end_bp: int,
    ) -> RegionResult:
        """Range-fetch the harmonised file for a genomic window."""
        url = harmonised_file_url(accession, ftp_base=self.ftp_base,
                                  session=self.session, timeout_s=self.timeout_s)
        notes: list[str] = []
        rows = self._tabix_fetch_rows(url, chromosome, start_bp, end_bp, notes)
        variants = [_normalise_row(r) for r in rows]
        # Drop any whose normalisation failed.
        variants = [v for v in variants if v is not None]
        release = GWASCatalogRelease(
            accession=accession,
            harmonised_url=url,
            fetched_at_utc=_now_utc(),
        )
        return RegionResult(
            accession=accession,
            chromosome=chromosome,
            region_start_bp=start_bp,
            region_end_bp=end_bp,
            n_variants=len(variants),
            variants=variants,
            release=release,
            notes=notes,
        )

    def _tabix_fetch_rows(
        self,
        url: str,
        chromosome: str,
        start_bp: int,
        end_bp: int,
        notes: list[str],
    ) -> list[dict[str, Any]]:
        """Tabix range fetch via pysam. Header columns come from the first row
        of the bgzipped file via HTTP range request (the GWAS-SSF v1.0 format
        does not prefix the column-name line with `#`, so tabix's `.header`
        attribute is empty for these files).
        """
        try:
            import pysam  # noqa: F401  imported lazily
        except ImportError as e:
            raise GWASCatalogFetchError(
                "pysam is required for tabix range fetches; install via `pip install pysam`"
            ) from e
        from pysam import TabixFile  # local import keeps top-level cheap

        try:
            tbx = TabixFile(url)
        except (OSError, ValueError) as e:
            raise GWASCatalogFetchError(
                f"could not open tabix index for {url}: {e!s}"
            ) from e

        try:
            header_cols = _header_columns(list(tbx.header))
            if not header_cols:
                header_cols = _fetch_header_via_http(url, self.session, self.timeout_s)
            chrom_q = chromosome.lstrip("chr")
            try:
                rows = list(tbx.fetch(chrom_q, max(0, start_bp - 1), end_bp))
            except ValueError as e:
                # Tabix indexes can use 'chr' prefix or not. Retry the other form.
                if not chrom_q.startswith("chr"):
                    rows = list(tbx.fetch(f"chr{chrom_q}", max(0, start_bp - 1), end_bp))
                else:
                    raise GWASCatalogFetchError(
                        f"tabix fetch failed for {chromosome}:{start_bp}-{end_bp} on {url}: {e!s}"
                    ) from e
            parsed = []
            for row_str in rows:
                fields = row_str.split("\t")
                if len(fields) != len(header_cols):
                    notes.append(
                        f"row column count {len(fields)} != header {len(header_cols)}; skipping"
                    )
                    continue
                parsed.append(dict(zip(header_cols, fields)))
            return parsed
        finally:
            tbx.close()


def _fetch_header_via_http(url: str, session: requests.Session, timeout_s: float) -> list[str]:
    """Fetch the first ~16 KB of a bgzipped harmonised file via HTTP Range,
    decompress just enough to read line 1, and return its tab-split tokens.

    Used when pysam's tabix `.header` is empty (the GWAS-SSF v1.0 format
    does not `#`-prefix the column-name line).
    """
    headers = {"Range": "bytes=0-32767"}
    resp = session.get(url, headers=headers, timeout=timeout_s)
    if resp.status_code not in (200, 206):
        raise GWASCatalogFetchError(
            f"could not fetch column-name row from {url}: HTTP {resp.status_code}"
        )
    raw = resp.content
    # bgzip is gzip-format-compatible, so streaming gzip decompression of the
    # leading bytes works.
    import gzip
    import io
    try:
        with gzip.GzipFile(fileobj=io.BytesIO(raw)) as gz:
            first_line = gz.readline().decode("utf-8", errors="replace").rstrip("\r\n")
    except (OSError, EOFError) as e:
        raise GWASCatalogFetchError(
            f"could not decompress leading bytes of {url}: {e!s}"
        ) from e
    if not first_line:
        raise GWASCatalogFetchError(f"first line of {url} is empty")
    cols = [c.lstrip("#") for c in first_line.split("\t")]
    return cols


def _header_columns(header_lines: list[str]) -> list[str]:
    """The harmonised TSV header is the last `#`-or-tab-stripped line of the
    file's preamble. Most files put it at line 1 with no `#` prefix.
    """
    if not header_lines:
        return []
    last = header_lines[-1]
    if last.startswith("#"):
        last = last[1:]
    return last.split("\t")


def _normalise_row(row: dict[str, Any]) -> RegionVariant | None:
    """Convert one harmonised TSV row to RegionVariant.

    Schema reference: GWAS Catalog harmonised columns include
        hm_variant_id, hm_chrom, hm_pos, hm_other_allele, hm_effect_allele,
        hm_beta, hm_odds_ratio, standard_error, p_value, hm_effect_allele_frequency
        (+ the optional ci_lower / ci_upper of an odds-ratio study)
    Returns None if essentials are missing.

    The parse is faithful: a binary-trait study's `odds_ratio` and its
    `ci_lower` / `ci_upper` are passed through as their own fields, and no
    `beta` is synthesised that the source did not report. Deriving log(OR) is
    left to the consumer.
    """
    chrom = (row.get("hm_chrom") or row.get("chromosome") or "").lstrip("chr")
    pos = _maybe_int(row.get("hm_pos") or row.get("base_pair_location"))
    ref = (row.get("hm_other_allele") or row.get("other_allele") or "").upper()
    alt = (row.get("hm_effect_allele") or row.get("effect_allele") or "").upper()
    if not (chrom and pos and ref and alt):
        return None
    variant_id = row.get("hm_variant_id") or f"{chrom}_{pos}_{ref}_{alt}"
    return RegionVariant(
        variant_id=str(variant_id),
        chromosome=chrom,
        position=int(pos),
        ref=ref,
        alt=alt,
        beta=_maybe_float(row.get("hm_beta") or row.get("beta")),
        se=_maybe_float(row.get("standard_error") or row.get("se")),
        p_value=_maybe_float(row.get("p_value") or row.get("pvalue")),
        odds_ratio=_maybe_float(row.get("hm_odds_ratio") or row.get("odds_ratio")),
        ci_lower=_maybe_float(row.get("hm_ci_lower") or row.get("ci_lower")),
        ci_upper=_maybe_float(row.get("hm_ci_upper") or row.get("ci_upper")),
        effect_allele_frequency=_maybe_float(
            row.get("hm_effect_allele_frequency") or row.get("effect_allele_frequency")
        ),
        raw=dict(row),
    )


def _maybe_float(v: Any) -> float | None:
    if v is None or v == "" or v == "NA":
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if (f == f) else None


def _maybe_int(v: Any) -> int | None:
    if v is None or v == "" or v == "NA":
        return None
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def _now_utc() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


# Local cache — mirrors ClawBio variant-annotation convention.
import os as _os  # noqa: E402

DEFAULT_CACHE_DIR = Path(
    _os.environ.get(
        "GWAS_CATALOG_CACHE_DIR",
        Path.home() / ".clawbio" / "gwas_catalog_region_fetch_cache",
    )
).expanduser()


def main(argv: list[str] | None = None) -> int:
    """Standard skill CLI: --input <config> --output <dir> --demo.

    Config schema (JSON or YAML):
        accession: GCST90269602
        chromosome: "1"
        start_bp: 108774968
        end_bp: 109774968

    Writes <output>/{variants.tsv, manifest.yaml, report.md}.
    """
    parser = argparse.ArgumentParser(
        prog="gwas-catalog-region-fetch",
        description="Fetch a region of GWAS summary statistics from NHGRI-EBI GWAS Catalog harmonised via tabix-on-FTP.",
    )
    parser.add_argument("--input", type=Path, help="JSON or YAML config (see docstring).")
    parser.add_argument("--output", type=Path,
                        help="Output directory; created if missing. Required unless --list-demos.")
    parser.add_argument("--demo", nargs="?", const="__default__", default=None,
                        metavar="NAME",
                        help="Run a bundled demo. Bare --demo runs the default; "
                             "pass a name to choose a specific one. See --list-demos.")
    parser.add_argument("--list-demos", action="store_true",
                        help="List bundled demo configs in this skill's examples/ directory.")
    parser.add_argument("--no-cache", action="store_true",
                        help="Bypass the local cache; always fetch fresh from FTP.")
    args = parser.parse_args(argv)

    if args.list_demos:
        _print_available_demos()
        return 0
    if args.demo is None and args.input is None:
        parser.error("either --input <config> or --demo [NAME] or --list-demos is required")
    if args.output is None:
        parser.error("--output is required")
    args.output.mkdir(parents=True, exist_ok=True)

    if args.demo is not None:
        cfg_path = _resolve_demo_path(args.demo)
        cfg = _load_config(cfg_path)
        print(f"info: using bundled demo {cfg_path.name}", file=sys.stderr)
    else:
        cfg = _load_config(args.input)

    cache_dir = None if args.no_cache else DEFAULT_CACHE_DIR
    client = GWASCatalogClient()
    result = _fetch_with_cache(client=client, cfg=cfg, cache_dir=cache_dir)

    tsv_path = args.output / "variants.tsv"
    _write_canonical_tsv(result, tsv_path)

    manifest = {
        "skill": "gwas-catalog-region-fetch",
        "version": "0.2.0",
        "accession": cfg["accession"],
        "region": {"chromosome": str(cfg["chromosome"]),
                   "start_bp": int(cfg["start_bp"]),
                   "end_bp": int(cfg["end_bp"])},
        "n_variants": result.n_variants,
        "release": {
            "accession": result.release.accession,
            "harmonised_url": result.release.harmonised_url,
            "harmoniser_version": result.release.harmoniser_version,
            "fetched_at_utc": result.release.fetched_at_utc,
        },
        "outputs": {"variants_tsv": "variants.tsv"},
    }
    try:
        import yaml as _yaml
        (args.output / "manifest.yaml").write_text(_yaml.safe_dump(manifest, sort_keys=False))
    except ImportError:
        (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str))

    report = [
        "# gwas-catalog-region-fetch report",
        "",
        f"- **Study accession:** `{cfg['accession']}`",
        f"- **Region:** chr{cfg['chromosome']}:{int(cfg['start_bp']):,}-{int(cfg['end_bp']):,}",
        f"- **Variants returned:** {result.n_variants}",
        f"- **Source URL:** {result.release.harmonised_url}",
        f"- **Output TSV:** {tsv_path.name}",
    ]
    (args.output / "report.md").write_text("\n".join(report) + "\n")

    print(f"gwas-catalog-region-fetch: {result.n_variants} variants -> {tsv_path}")
    print(f"  source: GWAS Catalog harmonised | accession {result.release.accession}")
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
        for cand in ("default.json", "default.yaml", "default.yml", "input.json"):
            p = examples / cand
            if p.is_file():
                return p
        files = _list_demos()
        if not files:
            raise FileNotFoundError(f"no bundled demo configs found in {examples}")
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


# Bumped when the cached row shape changes, so an entry written by an earlier
# version (for example one without the odds-ratio CI bounds) is never served.
_CACHE_SCHEMA = "v2"


def _cache_key(cfg: dict) -> str:
    chrom = str(cfg["chromosome"]).lstrip("chr")
    return (
        f"{cfg['accession']}__chr{chrom}_{int(cfg['start_bp'])}_{int(cfg['end_bp'])}"
        f"__{_CACHE_SCHEMA}.json"
    )


def _fetch_with_cache(*, client, cfg: dict, cache_dir: Path | None) -> "RegionResult":
    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path = cache_dir / _cache_key(cfg)
        if cache_path.is_file():
            return _region_result_from_cache(json.loads(cache_path.read_text()))
    result = client.fetch_region(
        accession=cfg["accession"],
        chromosome=str(cfg["chromosome"]),
        start_bp=int(cfg["start_bp"]),
        end_bp=int(cfg["end_bp"]),
    )
    if cache_dir is not None:
        cache_path = cache_dir / _cache_key(cfg)
        cache_path.write_text(json.dumps(result.to_dict(), default=str))
    return result


def _region_result_from_cache(d: dict) -> "RegionResult":
    rel = d.get("release") or {}
    release = GWASCatalogRelease(
        accession=rel.get("accession", d["accession"]),
        harmonised_url=rel.get("harmonised_url", ""),
        fetched_at_utc=rel.get("fetched_at_utc", ""),
        harmoniser_version=rel.get("harmoniser_version"),
    )
    variants = [
        RegionVariant(
            variant_id=v["variant_id"], chromosome=v["chromosome"], position=int(v["position"]),
            ref=v["ref"], alt=v["alt"],
            beta=v.get("beta"), se=v.get("se"), p_value=v.get("p_value"),
            odds_ratio=v.get("odds_ratio"),
            ci_lower=v.get("ci_lower"), ci_upper=v.get("ci_upper"),
            effect_allele_frequency=v.get("effect_allele_frequency"),
            raw=v.get("raw") or {},
        ) for v in d.get("variants", [])
    ]
    return RegionResult(
        accession=d["accession"], chromosome=d["chromosome"],
        region_start_bp=int(d["region_start_bp"]), region_end_bp=int(d["region_end_bp"]),
        release=release, n_variants=int(d["n_variants"]),
        variants=variants, notes=list(d.get("notes") or []),
    )


def _write_canonical_tsv(result, tsv_path: Path) -> None:
    """Emit the canonical sumstats-slice TSV (variant_id, chromosome,
    position_bp, allele_a, allele_b, beta, se, p, eaf, study_id) for
    downstream coloc / fine-mapping / regional plotting tools."""
    cols = ["variant_id", "chromosome", "position_bp",
            "allele_a", "allele_b", "beta", "se", "p", "eaf", "study_id"]
    with tsv_path.open("w") as f:
        f.write("# locuscompare-schema-version: 1.0\n")
        f.write("# source: gwas_catalog\n")
        f.write(f"# accession: {result.accession}\n")
        f.write("\t".join(cols) + "\n")
        for v in result.variants:
            row = [
                v.variant_id, v.chromosome, str(v.position),
                v.ref, v.alt,
                "" if v.beta is None else f"{v.beta:.6g}",
                "" if v.se is None else f"{v.se:.6g}",
                "" if v.p_value is None else f"{v.p_value:.6g}",
                "" if v.effect_allele_frequency is None else f"{v.effect_allele_frequency:.6g}",
                result.accession,
            ]
            f.write("\t".join(row) + "\n")


if __name__ == "__main__":
    sys.exit(main())
