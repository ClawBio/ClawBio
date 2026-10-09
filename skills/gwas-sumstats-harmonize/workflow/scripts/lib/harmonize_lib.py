"""
scripts/lib/harmonize_lib.py

Plain-python helpers shared by the harmonization stage scripts: table I/O,
column-name detection, allele logic and the p/beta/SE conversions. Standard
library only, so every stage runs in any Python 3.10+ environment.
"""
import csv
import gzip
import math
import os
import re
from decimal import Decimal, InvalidOperation
from statistics import NormalDist, median

CANONICAL = ["SNP", "CHR", "BP", "EA", "NEA", "EAF", "BETA", "SE", "P", "N"]
HELPERS = ["OR", "LOG10P", "Z"]  # carried from map_columns to derive_effects only

# Lower-cased header aliases, in priority order. Deliberately absent: MAF and
# MINOR_ALLELE_FREQUENCY (not allele-specific, so never an effect-allele frequency).
ALIASES = {
    "SNP": ["snp", "rsid", "rs_id", "rsids", "markername", "variant_id", "id", "snpid", "markerid", "marker", "rs"],
    "CHR": ["chr", "chrom", "#chrom", "chromosome", "#chr", "chr_name", "hg19chrc"],
    "BP": ["bp", "pos", "position", "base_pair_location", "genpos",
           "bp_hg19", "bp_grch37", "pos_b37", "bp_hg38", "bp_grch38", "pos_b38"],
    "EA": ["ea", "effect_allele", "a1", "allele1", "tested_allele", "alt", "eff_allele"],
    "NEA": ["nea", "other_allele", "non_effect_allele", "a2", "allele2", "allele0", "oa", "ref"],
    "EAF": ["eaf", "effect_allele_frequency", "a1_freq", "a1freq", "freq1", "frq", "freq", "af_alt", "alt_freq"],
    "BETA": ["beta", "b", "effect", "log_or", "logor"],
    "OR": ["or", "odds_ratio"],
    "SE": ["se", "stderr", "standard_error", "log(or)_se", "se_beta", "sebeta"],
    "P": ["p", "pval", "p_value", "pvalue", "p-value", "p.value", "p_bolt_lmm"],
    "LOG10P": ["log10p", "mlog10p", "neg_log10_p_value", "-log10p"],
    "Z": ["z", "zscore", "z_score", "z_stat"],
    "N": ["n", "n_total", "obs_ct", "nobs", "n_samples", "totaln"],
}
NOT_EAF = {"maf", "minor_allele_frequency"}
# Build-labelled position/chromosome headers: the label is recorded, never ignored.
BUILD_HINTS = {
    "bp_hg19": "GRCh37", "bp_grch37": "GRCh37", "pos_b37": "GRCh37", "hg19chrc": "GRCh37",
    "bp_hg38": "GRCh38", "bp_grch38": "GRCh38", "pos_b38": "GRCh38",
}
# A palindromic SNV is only strand-resolvable when both the effect-allele and
# the reference frequency are clearly away from 0.5.
PALINDROMIC_BAND = (0.4, 0.6)
MISSING_TOKENS = {"", "na", "nan", "none", "null", "."}
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
CHR_ORDER = {str(i): i for i in range(1, 23)}
CHR_ORDER.update({"X": 23, "Y": 24, "XY": 25, "MT": 26})
_NORMAL = NormalDist()
_SNP_POS_RE = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|XY|MT|M)[:_](\d+)", re.IGNORECASE)


class HarmonizeError(Exception):
    """Raised when an input cannot be harmonized (unmappable columns, bad mapping)."""


# ── I/O ───────────────────────────────────────────────────────────────────────


def _open(path, mode):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode + "t", newline="")
    return open(path, mode, newline="")


def sniff_delimiter(first_line):
    if "\t" in first_line:
        return "\t"
    if "," in first_line:
        return ","
    return None  # whitespace


def read_raw(path):
    """Read a delimited summary-statistics file of unknown delimiter -> (header, rows)."""
    with _open(path, "r") as fh:
        first = fh.readline().rstrip("\r\n")
        delim = sniff_delimiter(first)
        header = first.split(delim) if delim else first.split()
        if delim:
            rows = [r for r in csv.reader(fh, delimiter=delim)]
        else:
            rows = [line.split() for line in fh if line.strip()]
    rows = [[c.strip() for c in r] for r in rows if any(c.strip() for c in r)]
    return [h.strip() for h in header], rows


def read_table(path):
    """Read a tab-separated file with a header -> list of dicts."""
    with _open(path, "r") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_table(path, columns, rows):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with _open(path, "w") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow([row.get(c, "") for c in columns])


def is_missing(value):
    return value is None or str(value).strip().lower() in MISSING_TOKENS


# ── Column detection ──────────────────────────────────────────────────────────


def detect_columns(header, explicit=None):
    """Map canonical names to header names.

    Returns {"mapping": {canonical: header_name}, "plink2_other": (ref, alt) | None,
             "notes": [...]}. `explicit` ({canonical: header_name}) always wins.
    """
    lower = {h.lower(): h for h in header}
    mapping, notes, used = {}, [], set()
    for canon, name in (explicit or {}).items():
        canon = canon.upper()
        if name not in header:
            raise HarmonizeError(
                f"explicit column mapping {canon} -> '{name}' not found in header {header}"
            )
        mapping[canon] = name
        used.add(name)

    plink2_other = None
    if "EA" not in mapping and "NEA" not in mapping and {"allele1", "allele2"} <= set(lower)             and {"af_allele2", "ac_allele2"} & set(lower):
        # SAIGE step 2: "association results are with regard to Allele2", so
        # BETA and AF_Allele2 belong to Allele2. METAL and BOLT use Allele1.
        mapping["EA"], mapping["NEA"] = lower["allele2"], lower["allele1"]
        used.update({lower["allele1"], lower["allele2"]})
        if "EAF" not in mapping and "af_allele2" in lower:
            mapping["EAF"] = lower["af_allele2"]
            used.add(lower["af_allele2"])
        notes.append("SAIGE Allele1/Allele2 detected: EA = Allele2, NEA = Allele1 (SAIGE reports BETA for Allele2)")
    elif "EA" not in mapping and "NEA" not in mapping and {"a1", "ref", "alt"} <= set(lower):
        # PLINK2 --glm: A1 is the tested (effect) allele and may be REF or ALT;
        # the other allele is whichever of REF/ALT A1 is not, row by row.
        mapping["EA"] = lower["a1"]
        plink2_other = (lower["ref"], lower["alt"])
        used.update({lower["a1"], lower["ref"], lower["alt"]})
        notes.append("PLINK2 REF/ALT/A1 detected: EA = A1, NEA = the other of REF/ALT per row")

    for canon, aliases in ALIASES.items():
        if canon in mapping or (canon == "NEA" and plink2_other):
            continue
        for alias in aliases:
            name = lower.get(alias)
            if name is not None and name not in used:
                mapping[canon] = name
                used.add(name)
                break

    for h in header:
        if h.lower() in NOT_EAF and "EAF" not in mapping:
            notes.append(f"column '{h}' is MAF, not an effect-allele frequency; EAF left empty")

    builds = {BUILD_HINTS[mapping[c].lower()] for c in ("BP", "CHR")
              if c in mapping and mapping[c].lower() in BUILD_HINTS}
    if len(builds) > 1:
        raise HarmonizeError(f"CHR and BP columns name different genome builds: {sorted(builds)}")
    build = builds.pop() if builds else None
    if build:
        notes.append(f"genome build {build} inferred from the column name '{mapping['BP']}'")
    test = lower.get("test")
    return {"mapping": mapping, "plink2_other": plink2_other, "notes": notes, "header": list(header),
            "build": build, "test_column": test}


def check_derivable(detected):
    """Raise HarmonizeError unless every required canonical field is present or derivable."""
    m = set(detected["mapping"])
    missing = []
    if not ({"CHR", "BP"} <= m or "SNP" in m):
        missing.append("chromosome + position (CHR/BP, or a chr:pos SNP id)")
    if "EA" not in m:
        missing.append("effect allele (EA)")
    if "NEA" not in m and not detected.get("plink2_other"):
        missing.append("other allele (NEA)")
    if not ({"BETA", "OR"} & m or {"Z", "EAF", "N"} <= m):
        # Z alone fixes neither the scale nor the size of an effect; with EAF
        # and N it gives a per-SD (standardised) BETA and SE.
        missing.append("effect size (BETA or OR, or Z together with EAF and N)")
    if not ({"P", "LOG10P", "Z"} & m) and not ({"BETA", "OR"} & m and "SE" in m):
        missing.append("p-value (P, LOG10P or Z)")
    if not ({"SE"} & m) and not ({"P", "LOG10P", "Z"} & m):
        missing.append("standard error (SE, or a p-value to derive it)")
    if missing:
        raise HarmonizeError(
            "cannot harmonize: missing " + "; ".join(missing)
            + f". Header seen: {detected['header']}. Pass an explicit column mapping."
        )


def split_snp_position(snp):
    """'1:12345', 'chr1_12345:A:G' -> ('1', '12345'); otherwise (None, None)."""
    m = _SNP_POS_RE.match(snp or "")
    return (m.group(1), m.group(2)) if m else (None, None)


# ── Alleles and chromosomes ───────────────────────────────────────────────────


def normalize_chr(value):
    v = str(value).strip().upper()
    if v.startswith("CHR"):
        v = v[3:]
    v = {"23": "X", "24": "Y", "25": "XY", "26": "MT", "M": "MT"}.get(v, v)
    return v.lstrip("0") or v


def chr_sort_key(chrom, bp):
    c = normalize_chr(chrom)
    try:
        pos = int(float(bp))
    except ValueError:
        pos = 0
    return (CHR_ORDER.get(c, 99), c, pos)


def is_snv(allele):
    return len(allele) == 1 and allele in COMPLEMENT


def is_palindromic(a, b):
    return is_snv(a) and is_snv(b) and COMPLEMENT[a] == b


def complement(allele):
    """Reverse complement: the same allele read off the other strand ("AC" -> "GT")."""
    return "".join(COMPLEMENT.get(base, base) for base in reversed(allele))


def classify_alignment(ea, nea, ref, alt):
    """How (EA, NEA) relates to a reference (REF, ALT), with EA expected to be ALT.

    Returns aligned | swapped | strand_flipped | strand_flipped_swapped |
    palindromic | mismatch. A palindromic pair (A/T, C/G) matching the reference
    alleles is returned as "palindromic", never as aligned or swapped: reverse-
    strand A/T is indistinguishable from a forward-strand swap, so only the
    allele frequency can decide (see resolve_palindromic).
    """
    if is_palindromic(ea, nea):
        return "palindromic" if {ea, nea} == {ref, alt} else "mismatch"
    if (ea, nea) == (alt, ref):
        return "aligned"
    if (ea, nea) == (ref, alt):
        return "swapped"
    cea, cnea = complement(ea), complement(nea)
    if (cea, cnea) == (alt, ref):
        return "strand_flipped"
    if (cea, cnea) == (ref, alt):
        return "strand_flipped_swapped"
    return "mismatch"


def resolve_palindromic(ea, nea, eaf, ref, alt, af):
    """Strand-resolve a palindromic SNV by frequency; None when it cannot be.

    The effect allele is the reference ALT when EAF and the reference ALT
    frequency fall on the same side of 0.5, and REF otherwise. Unresolvable
    when either frequency is missing or inside PALINDROMIC_BAND. The return
    value names the operation align_reference applies, as classify_alignment does.
    """
    lo, hi = PALINDROMIC_BAND
    if eaf is None or af is None or lo <= eaf <= hi or lo <= af <= hi:
        return None
    ea_is_alt = (eaf > 0.5) == (af > 0.5)
    if ea == alt:
        return "aligned" if ea_is_alt else "strand_flipped_swapped"
    return "strand_flipped" if ea_is_alt else "swapped"


# ── Statistics ────────────────────────────────────────────────────────────────


def to_float(value):
    if is_missing(value):
        return None
    try:
        f = float(value)
    except ValueError:
        return None
    return f if math.isfinite(f) else None


def valid_p(value):
    """0 < p <= 1, evaluated exactly so values below the float64 minimum survive."""
    if is_missing(value):
        return False
    try:
        d = Decimal(str(value).strip())
    except InvalidOperation:
        return False
    return d.is_finite() and Decimal(0) < d <= Decimal(1)


def p_from_log10p(log10p):
    """-log10(p) -> p as a decimal string, exact even where float64 underflows."""
    x = Decimal(str(log10p).strip())
    exponent = int(x.to_integral_value(rounding="ROUND_CEILING"))
    mantissa = 10 ** float(exponent - x)  # in (1, 10]
    if mantissa >= 10:
        mantissa /= 10
        exponent -= 1
    return f"{mantissa:.4f}e-{exponent}" if exponent > 0 else f"{10 ** float(-x):.6g}"


def z_from_p(p):
    # -inv_cdf(p/2), not inv_cdf(1 - p/2): 1 - p/2 rounds to exactly 1.0 once
    # p < ~1e-16, which loses every genome-wide-significant tail.
    return -_NORMAL.inv_cdf(p / 2)


def se_from_beta_p(beta, p):
    z = z_from_p(p)
    return abs(beta) / z if z > 0 else None


def p_from_beta_se(beta, se):
    return 2 * _NORMAL.cdf(-abs(beta / se))


def beta_se_from_z(z, eaf, n):
    """Standardised (per-SD) BETA and SE from Z, EAF and N (Zhu et al. 2016, SMR)."""
    denom = math.sqrt(2 * eaf * (1 - eaf) * (n + z * z))
    return z / denom, 1 / denom


def p_from_z(z):
    return 2 * _NORMAL.cdf(-abs(z))


def lambda_gc(pvalues):
    """Genomic inflation: median chi-square(1 df) / 0.4549364."""
    chis = []
    for p in pvalues:
        p = max(float(p), 1e-300)
        if 0 < p <= 1:
            chis.append(z_from_p(p) ** 2 if p < 1 else 0.0)
    return median(chis) / 0.4549364 if chis else float("nan")


def fmt(x):
    return f"{x:.6g}"
