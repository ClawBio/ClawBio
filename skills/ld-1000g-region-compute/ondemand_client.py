"""On-demand 1000G LD region fetch + plink 1.9 r² compute.

Per-region tabix fetch from EBI's 1000G FTP, super-pop-filtered, run through
plink 1.9 to get r² between the lead and every variant in the window.
Caches both the region VCF and the LD output to
`~/.clawbio/locuscompare_cache/1000g/`.

This means a fresh ClawBio install can render LD-coloured plots without
asking the user to download a 3 GB PLINK panel first: the EBI fetch is a
~5-50 MB byte-range request per locus, and plink 1.9 is sub-second on the
resulting region VCF despite being single-threaded.

Constructor arg `super_pop` selects which 1000G samples to keep
(EUR/AFR/AMR/EAS/SAS). The sample-to-super-pop mapping is fetched once
(small TSV) and cached.

License: 1000G data are open-access with attribution (Auton 2015, Clarke
2017). plink 1.9 binary is GPL-3 (subprocess invocation only; no GPL
contamination of MIT-licensed locuscompare code).
"""
from __future__ import annotations

import csv
import os
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import requests

# 1000G GRCh38 phased genotypes (NYGC re-imputed, 2019-03-12 release).
# These are the canonical liftover-free GRCh38 panel for LD reference work.
ONEKG_VCF_URL_TEMPLATE = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/"
    "ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
)
ONEKG_PANEL_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "integrated_call_samples_v3.20130502.ALL.panel"
)

DEFAULT_CACHE_DIR = Path(
    os.environ.get(
        "LOCUSCOMPARE_CACHE_DIR",
        Path.home() / ".clawbio" / "locuscompare_cache",
    )
) / "1000g"

# plink 1.9 is the supported binary (ubiquitous via brew / apt / conda;
# ships --ld-snp + --r2 + --ld-window-r2 natively).
DEFAULT_PLINK_BIN = os.environ.get("PLINK_BIN", "plink")

PLINK_NOT_FOUND_HINT = (
    "plink binary not found. Install via `brew install brewsci/bio/plink` "
    "(macOS), `apt-get install plink1.9` (Ubuntu/Debian), or "
    "`conda install -c bioconda plink` (any platform); direct binary "
    "downloads are available at https://www.cog-genomics.org/plink/1.9/. "
    "Then either ensure it is on PATH or set PLINK_BIN to its absolute path."
)

PANEL_ID_DEFAULT = "1000g_phase3_v5b_grch38_basic"
PANEL_VERSION_DEFAULT = "5b_remote_2019_03_12"


@dataclass
class OnDemandLDPair:
    partner_variant_id: str
    r2: float
    dprime: float | None = None


@dataclass
class OnDemandLDResult:
    panel_id: str
    panel_version: str
    super_pop: str
    plink_version: str  # binary version string returned by `plink --version`
    chromosome: str
    lead_variant_id: str
    window_bp: int
    n_partners_requested: int
    n_partners_returned: int
    pairs: list[OnDemandLDPair]
    fetched_at_utc: str
    notes: list[str] = field(default_factory=list)


def _alpha_plink_id(variant_id: str) -> str:
    """Convert a `chr_pos_ref_alt` id to plink 1.9's `chr:pos:allele1:allele2`,
    alleles in alphabetical order.

    plink 1.9's `--set-missing-var-ids '@:#:$1:$2'` writes the two alleles in
    alphabetical order, not REF/ALT order, so REF=G ALT=A becomes `chr:pos:A:G`.
    A naive `_` -> `:` swap cannot find such a variant with `--ld-snp` / `--extract`.
    """
    parts = variant_id.split("_", 3)
    if len(parts) < 4:
        return variant_id.replace("_", ":", 3)
    chrom, pos, a, b = parts
    a_sorted, b_sorted = sorted([a, b])
    return f"{chrom}:{pos}:{a_sorted}:{b_sorted}"


class OnDemandLDError(Exception):
    """Raised when on-demand LD compute cannot proceed."""


def _detect_plink_version(plink_bin: str) -> str:
    """Return the plink `--version` string for the manifest."""
    if shutil.which(plink_bin) is None:
        raise OnDemandLDError(f"{PLINK_NOT_FOUND_HINT} (looked for: {plink_bin})")
    try:
        proc = subprocess.run(
            [plink_bin, "--version"],
            capture_output=True, text=True, check=False, timeout=10,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        raise OnDemandLDError(f"could not run `{plink_bin} --version`: {e!s}") from e
    line = (proc.stdout or proc.stderr).strip().splitlines()
    return line[0] if line else "unknown"


def _resolve_super_pop_samples(super_pop: str, cache_dir: Path) -> list[str]:
    """Fetch (or read from cache) the 1000G sample to super_pop mapping;
    return the list of sample IDs in the requested super-pop.
    """
    panel_path = cache_dir / "integrated_call_samples_v3.20130502.ALL.panel"
    if not panel_path.is_file():
        cache_dir.mkdir(parents=True, exist_ok=True)
        resp = requests.get(ONEKG_PANEL_URL, timeout=60)
        resp.raise_for_status()
        panel_path.write_bytes(resp.content)
    samples: list[str] = []
    with panel_path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("super_pop") == super_pop:
                samples.append(row["sample"])
    if not samples:
        raise OnDemandLDError(
            f"no samples found for super_pop={super_pop!r} in {panel_path}. "
            f"Valid super-pop codes: EUR, AFR, AMR, EAS, SAS."
        )
    return samples


def _fetch_region_vcf(
    chromosome: str,
    start_bp: int,
    end_bp: int,
    cache_dir: Path,
) -> Path:
    """Tabix-fetch a 1000G region VCF via pysam. Cached per (chrom, start, end)."""
    chrom_bare = chromosome.removeprefix("chr") if chromosome.startswith("chr") else chromosome
    cache_dir.mkdir(parents=True, exist_ok=True)
    region_vcf = cache_dir / f"chr{chrom_bare}_{start_bp}_{end_bp}.vcf.gz"
    if region_vcf.is_file() and region_vcf.stat().st_size > 0:
        return region_vcf

    import pysam  # lazy import; pysam is optional for non-LD-colouring runs

    url = ONEKG_VCF_URL_TEMPLATE.format(chrom=chrom_bare)
    tmp_vcf = region_vcf.with_suffix(".tmp.vcf")
    try:
        # pysam.VariantFile supports remote tabix-indexed VCFs; the .tbi at
        # the source URL is loaded transparently.
        with pysam.VariantFile(url) as src:
            # 1000G files use bare contig names (`1`, `2`, ...), no `chr`.
            chrom_used = chrom_bare
            with tmp_vcf.open("w") as out:
                out.write(str(src.header))
                for rec in src.fetch(chrom_used, start_bp, end_bp):
                    out.write(str(rec))
    except Exception as e:
        if tmp_vcf.exists():
            tmp_vcf.unlink()
        raise OnDemandLDError(
            f"tabix-fetch from {url} for {chrom_bare}:{start_bp}-{end_bp} failed: {e!s}"
        ) from e

    import gzip
    with tmp_vcf.open("rb") as f_in, gzip.open(region_vcf, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    tmp_vcf.unlink()
    return region_vcf


class OnDemand1000GLDClient:
    """LD client that fetches per-region 1000G VCFs on demand from EBI FTP.

    Shells out to plink 1.9 for r² compute. Construction validates the
    plink binary is reachable; per-call work performs the tabix region
    fetch, super-pop sample filter, and r² compute.
    """

    def __init__(
        self,
        super_pop: str = "EUR",
        plink_bin: str = DEFAULT_PLINK_BIN,
        cache_dir: Path | None = None,
        panel_id: str = PANEL_ID_DEFAULT,
        panel_version: str = PANEL_VERSION_DEFAULT,
    ) -> None:
        self.super_pop = super_pop
        self.plink_bin = plink_bin
        self.cache_dir = (cache_dir or DEFAULT_CACHE_DIR).expanduser()
        self.panel_id = panel_id
        self.panel_version = panel_version
        self.plink_version = _detect_plink_version(plink_bin)
        self._super_pop_samples: list[str] | None = None

    def _samples(self) -> list[str]:
        if self._super_pop_samples is None:
            self._super_pop_samples = _resolve_super_pop_samples(self.super_pop, self.cache_dir)
        return self._super_pop_samples

    def r2_with_lead(
        self,
        lead: str,
        partners: Iterable[str],
        chromosome: str,
        window_bp: int,
    ) -> OnDemandLDResult:
        """Compute r² between `lead` and each `partner` via on-demand region fetch.

        `lead` and `partners` use OT-style chr_pos_ref_alt ids. `chromosome`
        is the chromosome name (with or without `chr` prefix). `window_bp`
        determines the region around the lead to fetch.
        """
        partner_list = list(partners)
        notes: list[str] = []
        if not partner_list:
            return OnDemandLDResult(
                panel_id=self.panel_id, panel_version=self.panel_version,
                super_pop=self.super_pop, plink_version=self.plink_version,
                chromosome=chromosome, lead_variant_id=lead, window_bp=window_bp,
                n_partners_requested=0, n_partners_returned=0, pairs=[],
                fetched_at_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                notes=["no partners requested"],
            )

        chrom_bare = chromosome.removeprefix("chr") if chromosome.startswith("chr") else chromosome
        # Parse lead position from variant_id (chr_pos_ref_alt)
        lead_parts = lead.split("_")
        if len(lead_parts) < 4:
            raise OnDemandLDError(f"cannot parse lead variant id: {lead!r}")
        lead_pos = int(lead_parts[1])
        half = max(window_bp // 2, 1)
        start_bp = max(0, lead_pos - half)
        end_bp = lead_pos + half

        region_vcf = _fetch_region_vcf(chrom_bare, start_bp, end_bp, self.cache_dir)
        notes.append(f"fetched 1000G region VCF to {region_vcf}")

        samples = self._samples()
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            keep_path = td_path / "keep.txt"
            with keep_path.open("w") as f:
                # plink 1.9 `--vcf` assigns FID=IID=sample-id (per
                # https://www.cog-genomics.org/plink/1.9/data#vcf). `--keep`
                # expects "FID\tIID" pairs, so write the sample id twice.
                for s in samples:
                    f.write(f"{s}\t{s}\n")
            out_prefix = td_path / "ld_out"
            # 1000G GRCh38 VCFs hold rsids in the ID column, not chr:pos:ref:alt.
            # plink 1.9 `--set-missing-var-ids '@:#:$1:$2'` rewrites IDs into
            # the canonical chr:pos:ref:alt form (@=chr, #=bp, $1/$2=alleles).
            # `--set-all-var-ids` (plink2's spelling) is NOT available in 1.9;
            # `--set-missing-var-ids` covers the case because every VCF row's
            # ID field is `.` in this distribution, which 1.9 treats as missing.
            lead_panel = _alpha_plink_id(lead)
            ot_form_by_plink: dict[str, str] = {}
            partner_panel: list[str] = []
            for p in partner_list:
                pf = _alpha_plink_id(p)
                ot_form_by_plink[pf] = p
                partner_panel.append(pf)

            extract_path = td_path / "extract.txt"
            with extract_path.open("w") as f:
                f.write(lead_panel + "\n")
                for p in partner_panel:
                    f.write(p + "\n")

            cmd = [
                self.plink_bin,
                "--vcf", str(region_vcf),
                "--keep", str(keep_path),
                "--set-missing-var-ids", "@:#:$1:$2",
                "--extract", str(extract_path),
                "--r2",
                "--ld-snp", lead_panel,
                "--ld-window-kb", str(max(window_bp // 1000, 1)),
                "--ld-window", "99999",
                "--ld-window-r2", "0",
                "--out", str(out_prefix),
            ]
            try:
                proc = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=300)
            except (subprocess.TimeoutExpired, FileNotFoundError) as e:
                raise OnDemandLDError(f"plink invocation failed: {e!s}") from e
            if proc.returncode != 0:
                # Lead not in the 1000G panel (rare allele, indel encoded
                # differently, variant added after 1000G): plink cannot compute
                # LD against a variant it has no genotypes for. Return no pairs
                # with a note, so the caller renders without LD colouring
                # instead of failing the whole run.
                stderr = proc.stderr or ""
                if (
                    "No valid variants specified by --ld-snp" in stderr
                    or "No variants remaining after --extract" in stderr
                ):
                    notes.append(
                        f"lead variant {lead!r} (or every partner) not in 1000G "
                        f"{self.super_pop} panel; returning empty pairs"
                    )
                    return OnDemandLDResult(
                        panel_id=self.panel_id, panel_version=self.panel_version,
                        super_pop=self.super_pop, plink_version=self.plink_version,
                        chromosome=chrom_bare, lead_variant_id=lead, window_bp=window_bp,
                        n_partners_requested=len(partner_list), n_partners_returned=0,
                        pairs=[],
                        fetched_at_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                        notes=notes,
                    )
                raise OnDemandLDError(
                    f"plink exited with code {proc.returncode}. "
                    f"stderr (truncated): {stderr[:1000]}"
                )

            ld_path = Path(f"{out_prefix}.ld")
            if not ld_path.exists():
                raise OnDemandLDError(
                    f"plink produced no .ld output at {out_prefix}; "
                    f"stdout (truncated): {(proc.stdout or '')[:1000]}"
                )
            pairs = _parse_ld(ld_path, lead_panel, notes)
            for p in pairs:
                p.partner_variant_id = ot_form_by_plink.get(
                    p.partner_variant_id, p.partner_variant_id.replace(":", "_", 3)
                )

        return OnDemandLDResult(
            panel_id=self.panel_id, panel_version=self.panel_version,
            super_pop=self.super_pop, plink_version=self.plink_version,
            chromosome=chrom_bare, lead_variant_id=lead, window_bp=window_bp,
            n_partners_requested=len(partner_list), n_partners_returned=len(pairs),
            pairs=pairs,
            fetched_at_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            notes=notes,
        )


def _parse_ld(path: Path, lead: str, notes: list[str]) -> list[OnDemandLDPair]:
    """Parse a plink 1.9 `--r2`-emitted `.ld` file.

    Columns: `CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2`. Whitespace-separated
    (multiple spaces between fields), so we tokenise with `str.split()`
    rather than csv.
    """
    out: list[OnDemandLDPair] = []
    with path.open() as f:
        header_tokens: list[str] | None = None
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            tokens = line.split()
            if header_tokens is None:
                header_tokens = tokens
                continue
            row = dict(zip(header_tokens, tokens))
            id_a = row.get("SNP_A", "")
            id_b = row.get("SNP_B", "")
            partner = id_b if id_a == lead else (id_a if id_b == lead else None)
            if partner is None:
                continue
            r2_str = row.get("R2", "")
            if not r2_str or r2_str == "NA":
                notes.append(f"missing r² for partner {partner}")
                continue
            try:
                r2 = float(r2_str)
            except ValueError:
                continue
            out.append(OnDemandLDPair(partner_variant_id=partner, r2=r2))
    return out


__all__ = [
    "DEFAULT_PLINK_BIN",
    "OnDemand1000GLDClient",
    "OnDemandLDError",
    "OnDemandLDPair",
    "OnDemandLDResult",
    "PLINK_NOT_FOUND_HINT",
]
