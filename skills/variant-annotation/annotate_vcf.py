#!/usr/bin/env python3
"""
Variant Annotation Pipeline — ClawBio skill.

Reads a VCF file, queries the Ensembl VEP REST API (GRCh38) in batches
of up to 200 variants, and produces:
  1. An annotated TSV (CHROM, POS, REF, ALT, Gene, Consequence, ClinVar, gnomAD_AF)
  2. A markdown summary highlighting pathogenic/likely pathogenic variants
     with gnomAD AF < 0.001

Coordinate system notes:
  - VCF is 1-based, inclusive.
  - VEP region notation is 1-based, inclusive.
  - Indels in VCF use a left-aligned padding base.  This script strips the
    padding base and adjusts coordinates before sending to VEP.

Rate limiting:
  - Max 15 requests/second to comply with Ensembl's fair-use policy.
  - Exponential backoff on HTTP 429 (Too Many Requests) and 503 (Service Unavailable).

Author: ClawBio community
License: MIT
"""

import csv
import hashlib
import json
import sys
import time
from urllib.error import HTTPError
from urllib.request import Request, urlopen

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

VEP_URL = "https://rest.ensembl.org/vep/homo_sapiens/region"
BATCH_SIZE = 200        # VEP accepts up to 200 variants per POST
MAX_RPS = 15            # Ensembl fair-use: 15 requests per second
RETRY_MAX = 3           # Max retries on transient errors
GNOMAD_AF_THRESHOLD = 0.001  # Flag variants rarer than this

# ---------------------------------------------------------------------------
# VCF parsing
# ---------------------------------------------------------------------------


def parse_vcf(path: str) -> list[tuple[str, int, str, str]]:
    """Parse a VCF file and return a list of (chrom, pos, ref, alt) tuples.

    Multi-allelic sites (comma-separated ALTs) are split into separate
    entries.  All other VCF columns are ignored.

    Args:
        path: Filesystem path to a VCF 4.x file.

    Returns:
        List of (chrom, pos, ref, alt) tuples.  *pos* is a 1-based integer.
    """
    variants: list[tuple[str, int, str, str]] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue  # malformed line
            chrom = fields[0]
            try:
                pos = int(fields[1])
            except ValueError:
                continue  # skip unparseable POS
            ref = fields[3].upper()
            alt_field = fields[4].upper()
            # Split multi-allelic ALTs
            for alt in alt_field.split(","):
                alt = alt.strip()
                if alt in (".", "*"):
                    continue  # skip missing / spanning-deletion markers
                variants.append((chrom, pos, ref, alt))
    return variants


# ---------------------------------------------------------------------------
# Coordinate conversion — VCF → VEP region notation
# ---------------------------------------------------------------------------


def format_vep_region(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Convert a single VCF variant to VEP POST region notation.

    VEP region format:  "chrom start end allele strand"
    All coordinates are 1-based inclusive.

    Indel rules (VCF uses a left-aligned padding base):
      SNV       — REF and ALT both length 1.
                   Region: chrom pos pos ALT 1
      Deletion  — len(REF) > len(ALT), ALT is a prefix of REF.
                   The padding base is the shared prefix.
                   Deleted bases start at (pos + len(ALT)).
                   Region: chrom (pos+len(ALT)) (pos+len(REF)-1) - 1
      Insertion — len(ALT) > len(REF), REF is a prefix of ALT.
                   Inserted bases are ALT[len(REF):].
                   The insertion point is between pos+len(REF)-1 and pos+len(REF).
                   Region: chrom (pos+len(REF)) (pos+len(REF)-1) inserted 1
                   (start > end signals an insertion to VEP)
      Complex   — REF and ALT differ in length but share no prefix.
                   Treat as MNV / block substitution.
                   Region: chrom pos (pos+len(REF)-1) ALT 1

    Raises:
        ValueError: If ref or alt is empty.
    """
    if not ref or not alt:
        raise ValueError(f"Empty REF or ALT: ref={ref!r}, alt={alt!r}")

    # Strip 'chr' prefix — VEP uses bare chromosome names for GRCh38
    c = chrom.removeprefix("chr").removeprefix("Chr").removeprefix("CHR")

    ref_len = len(ref)
    alt_len = len(alt)

    if ref_len == 1 and alt_len == 1:
        # SNV
        return f"{c} {pos} {pos} {ref}/{alt} 1"

    # Determine shared prefix length (the padding bases)
    shared = 0
    for r, a in zip(ref, alt):
        if r == a:
            shared += 1
        else:
            break

    if ref_len > alt_len and shared == alt_len:
        # Pure deletion: ALT is a prefix of REF
        del_start = pos + shared          # first deleted base
        del_end = pos + ref_len - 1       # last deleted base
        return f"{c} {del_start} {del_end} {ref[shared:]}/- 1"

    elif alt_len > ref_len and shared == ref_len:
        # Pure insertion: REF is a prefix of ALT
        inserted = alt[shared:]
        ins_after = pos + shared - 1      # base AFTER which the insertion occurs
        # VEP convention: start = ins_after + 1, end = ins_after  (start > end)
        return f"{c} {ins_after + 1} {ins_after} -/{inserted} 1"

    else:
        # Complex substitution / MNV
        # Strip shared prefix, adjust coordinates
        new_ref = ref[shared:]
        new_alt = alt[shared:]
        new_start = pos + shared
        new_end = new_start + len(new_ref) - 1
        return f"{c} {new_start} {new_end} {new_ref}/{new_alt} 1"


# ---------------------------------------------------------------------------
# VEP REST API
# ---------------------------------------------------------------------------


def query_vep(batch: list[str], attempt: int = 0) -> list[dict]:
    """POST a batch of variants to the Ensembl VEP REST API.

    Args:
        batch: List of VEP region strings.
        attempt: Current retry attempt (for exponential backoff).

    Returns:
        List of VEP result dicts.

    Raises:
        HTTPError: On non-retryable HTTP errors.
    """
    payload = json.dumps({"variants": batch, "check_existing": 1, "af_gnomade": 1, "af_gnomadg": 1, "clinical_significance": 1}).encode("utf-8")
    req = Request(
        VEP_URL,
        data=payload,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
        },
    )
    try:
        with urlopen(req, timeout=120) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except HTTPError as exc:
        if exc.code in (429, 503) and attempt < RETRY_MAX:
            wait = 2 ** (attempt + 1)  # 2s, 4s, 8s
            print(
                f"  ⏳ VEP returned {exc.code}, retrying in {wait}s "
                f"(attempt {attempt + 1}/{RETRY_MAX})…",
                file=sys.stderr,
            )
            time.sleep(wait)
            return query_vep(batch, attempt + 1)
        raise


# ---------------------------------------------------------------------------
# Annotation extraction
# ---------------------------------------------------------------------------


def extract_annotations(result: dict) -> tuple[str, str, str, str]:
    """Extract gene, consequence, ClinVar, and gnomAD AF from a VEP result.

    Uses VEP's own ``most_severe_consequence`` field rather than
    re-implementing the Ensembl severity ranking, which is non-trivial
    and version-dependent.

    Args:
        result: A single variant dict from the VEP response.

    Returns:
        (gene_symbol, most_severe_consequence, clinvar_significance, gnomad_af)
        Fields that are unavailable are returned as ".".
    """
    consequence = result.get("most_severe_consequence", ".")
    gene = "."
    clinvar = "."
    gnomad_af = "."

    # ---- Gene symbol: pick from the transcript whose consequence matches ----
    for tc in result.get("transcript_consequences", []):
        if consequence in tc.get("consequence_terms", []):
            gene = tc.get("gene_symbol", ".")
            break
    # Fallback: take the first transcript consequence with a gene symbol
    if gene == ".":
        for tc in result.get("transcript_consequences", []):
            if "gene_symbol" in tc:
                gene = tc["gene_symbol"]
                break

    # ---- ClinVar and gnomAD from colocated variants ----
    # Determine the ALT allele from the result for frequency lookup
    alt_allele = None
    allele_str = result.get("allele_string", "")
    if "/" in allele_str:
        parts = allele_str.split("/")
        if len(parts) >= 2:
            alt_allele = parts[1]  # first ALT allele

    for cv in result.get("colocated_variants", []):
        # ClinVar significance — VEP uses "clin_sig_allele" (per-allele)
        # or "clin_sig" (legacy). Try both.
        if clinvar == ".":
            clin_sig = cv.get("clin_sig_allele") or cv.get("clin_sig")
            if clin_sig:
                if isinstance(clin_sig, list):
                    clinvar = ",".join(str(s) for s in clin_sig)
                elif isinstance(clin_sig, str):
                    # clin_sig_allele format: "A:benign;A:drug_response"
                    # Extract unique significances, stripping allele prefixes
                    sigs = set()
                    for entry in clin_sig.split(";"):
                        entry = entry.strip()
                        if ":" in entry:
                            sigs.add(entry.split(":", 1)[1])
                        else:
                            sigs.add(entry)
                    clinvar = ",".join(sorted(sigs))
                else:
                    clinvar = str(clin_sig)

        # gnomAD allele frequency — VEP nests under "frequencies" -> allele -> population
        if gnomad_af == "." and alt_allele:
            freqs = cv.get("frequencies", {})
            allele_freqs = freqs.get(alt_allele, {})
            # Try gnomAD exomes first, then genomes, then generic af
            for key in ("gnomade", "gnomadg", "af"):
                val = allele_freqs.get(key)
                if val is not None:
                    gnomad_af = str(val)
                    break
            # Fallback: flat keys (older VEP versions)
            if gnomad_af == ".":
                for key in ("gnomade_af", "gnomad_af", "minor_allele_freq"):
                    val = cv.get(key)
                    if val is not None:
                        gnomad_af = str(val)
                        break

    return gene, consequence, clinvar, gnomad_af


# ---------------------------------------------------------------------------
# Result matching
# ---------------------------------------------------------------------------


def build_result_index(
    vep_results: list[dict],
) -> dict[tuple[str, int, str], dict]:
    """Index VEP results by (chrom, start, allele_string) for fast lookup.

    The *input* field echoed back by VEP is the most reliable way to match
    results to the original variants, but its format depends on the input
    notation.  We also index by (seq_region_name, start) as a fallback.
    """
    idx: dict[tuple[str, int, str], dict] = {}
    for r in vep_results:
        # Primary key: the input string we sent
        inp = r.get("input", "")
        parts = inp.split()
        if len(parts) >= 4:
            key = (parts[0], int(parts[1]), parts[3])  # (chrom, start, allele)
            idx[key] = r
        # Secondary key
        sr = str(r.get("seq_region_name", ""))
        st = r.get("start")
        als = r.get("allele_string", "")
        if sr and st:
            idx[(sr, int(st), als)] = r
    return idx


def match_variant(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    index: dict[tuple[str, int, str], dict],
) -> dict | None:
    """Try to match a VCF variant back to its VEP result."""
    c = chrom.removeprefix("chr").removeprefix("Chr").removeprefix("CHR")
    ref_len = len(ref)
    alt_len = len(alt)

    if ref_len == 1 and alt_len == 1:
        # SNV: VEP start == VCF pos, allele == ALT
        for allele in (alt, f"{ref}/{alt}"):
            hit = index.get((c, pos, allele))
            if hit:
                return hit
    elif ref_len > alt_len:
        # Deletion: start == pos + len(shared prefix)
        shared = 0
        for r, a in zip(ref, alt):
            if r == a:
                shared += 1
            else:
                break
        del_start = pos + shared
        for allele in ("-", f"{ref[shared:]}/{alt[shared:] or '-'}"):
            hit = index.get((c, del_start, allele))
            if hit:
                return hit
    elif alt_len > ref_len:
        # Insertion: start == pos + len(shared prefix)
        shared = 0
        for r, a in zip(ref, alt):
            if r == a:
                shared += 1
            else:
                break
        ins_start = pos + shared
        inserted = alt[shared:]
        for allele in (inserted, f"-/{inserted}"):
            hit = index.get((c, ins_start, allele))
            if hit:
                return hit
    else:
        # Complex / MNV
        shared = 0
        for r, a in zip(ref, alt):
            if r == a:
                shared += 1
            else:
                break
        new_start = pos + shared
        new_alt = alt[shared:]
        hit = index.get((c, new_start, new_alt))
        if hit:
            return hit

    # Brute-force fallback: search by input string containing our pos
    pos_str = str(pos)
    for key, r in index.items():
        inp = r.get("input", "")
        if pos_str in inp and c in inp:
            return r

    return None


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def main() -> None:
    if len(sys.argv) != 4:
        print(
            "Usage: python3 annotate_vcf.py <input.vcf> <output.tsv> <summary.md>",
            file=sys.stderr,
        )
        sys.exit(1)

    vcf_path = sys.argv[1]
    tsv_path = sys.argv[2]
    summary_path = sys.argv[3]

    # 1. Parse VCF
    print(f"📂 Parsing VCF: {vcf_path}", file=sys.stderr)
    variants = parse_vcf(vcf_path)
    print(f"   Found {len(variants)} variant(s)", file=sys.stderr)

    if not variants:
        print("⚠️  No variants found in VCF. Exiting.", file=sys.stderr)
        sys.exit(0)

    # 2. Batch and query VEP with rate limiting
    print(f"🌐 Querying Ensembl VEP (batches of {BATCH_SIZE})…", file=sys.stderr)
    all_results: list[dict] = []
    interval = 1.0 / MAX_RPS
    n_batches = (len(variants) + BATCH_SIZE - 1) // BATCH_SIZE

    for batch_num, i in enumerate(range(0, len(variants), BATCH_SIZE), start=1):
        batch_variants = variants[i : i + BATCH_SIZE]
        batch_regions = []
        for chrom, pos, ref, alt in batch_variants:
            try:
                region = format_vep_region(chrom, pos, ref, alt)
                batch_regions.append(region)
            except ValueError as exc:
                print(f"   ⚠️  Skipping malformed variant: {exc}", file=sys.stderr)

        if not batch_regions:
            continue

        print(
            f"   Batch {batch_num}/{n_batches} ({len(batch_regions)} variants)…",
            file=sys.stderr,
        )
        time.sleep(interval)  # rate-limit
        try:
            results = query_vep(batch_regions)
            all_results.extend(results)
        except HTTPError as exc:
            print(
                f"   ❌ VEP error (HTTP {exc.code}) on batch {batch_num}. "
                f"Skipping this batch.",
                file=sys.stderr,
            )

    # 3. Index results for matching
    result_index = build_result_index(all_results)

    # 4. Write annotated TSV
    print(f"📝 Writing annotated TSV: {tsv_path}", file=sys.stderr)
    flagged: list[tuple[str, int, str, str, str, str, str, str]] = []

    with open(tsv_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            ["CHROM", "POS", "REF", "ALT", "Gene", "Consequence", "ClinVar", "gnomAD_AF"]
        )
        for chrom, pos, ref, alt in variants:
            matched = match_variant(chrom, pos, ref, alt, result_index)
            if matched:
                gene, cons, clin, af = extract_annotations(matched)
            else:
                gene, cons, clin, af = ".", ".", ".", "."

            writer.writerow([chrom, pos, ref, alt, gene, cons, clin, af])

            # Flag: pathogenic/likely pathogenic AND rare
            if af != "." and clin != ".":
                try:
                    af_val = float(af)
                except (ValueError, TypeError):
                    af_val = None
                if af_val is not None and af_val < GNOMAD_AF_THRESHOLD:
                    clin_lower = clin.lower()
                    if "pathogenic" in clin_lower:
                        flagged.append((chrom, pos, ref, alt, gene, cons, clin, af))

    # 5. Write markdown summary
    print(f"📋 Writing summary: {summary_path}", file=sys.stderr)
    with open(tsv_path, "rb") as fh:
        tsv_sha256 = hashlib.sha256(fh.read()).hexdigest()

    with open(summary_path, "w") as fh:
        fh.write("# Variant Annotation Summary\n\n")
        fh.write(f"**Input VCF**: `{vcf_path}`\n\n")
        fh.write(f"**Total variants annotated**: {len(variants)}\n\n")
        fh.write(
            f"**Flagged** (ClinVar pathogenic/likely pathogenic + "
            f"gnomAD AF < {GNOMAD_AF_THRESHOLD}): **{len(flagged)}**\n\n"
        )

        if flagged:
            fh.write("## Flagged Variants\n\n")
            fh.write(
                "| CHROM | POS | REF | ALT | Gene | Consequence | ClinVar | gnomAD_AF |\n"
            )
            fh.write(
                "|-------|-----|-----|-----|------|-------------|---------|----------|\n"
            )
            for row in flagged:
                fh.write("| " + " | ".join(str(x) for x in row) + " |\n")
            fh.write("\n")
        else:
            fh.write("_No variants met the flagging criteria._\n\n")

        fh.write("## Reproducibility\n\n")
        fh.write(f"- **Annotated TSV**: `{tsv_path}`\n")
        fh.write(f"- **TSV SHA-256**: `{tsv_sha256}`\n")
        fh.write(f"- **Genome assembly**: GRCh38\n")
        fh.write(f"- **VEP endpoint**: `{VEP_URL}`\n")
        fh.write(f"- **Batch size**: {BATCH_SIZE}\n")
        fh.write(f"- **Flagging threshold**: gnomAD AF < {GNOMAD_AF_THRESHOLD}\n")
        fh.write(f"- **Date**: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}\n")

    print(f"✅ Done. {len(flagged)} variant(s) flagged.", file=sys.stderr)


if __name__ == "__main__":
    main()
