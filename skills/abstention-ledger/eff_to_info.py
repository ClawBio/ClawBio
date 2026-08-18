"""Adapter: legacy SnpEff ``EFF`` -> the INFO keys existing ClawBio skills read.

Why
---
`rare-high-impact-variants` takes its consequence from ``MC``, ``Consequence`` or
``ANN`` (``rare_high_impact_variants.py:126``). Classic SnpEff writes ``EFF``, so
on legacy-annotated input that chain yields the empty string, the record is
dropped at line 160, and the run finishes with exit code 0 and every impact
metric at zero. This module is the translation layer that makes the existing
skill work on that corpus, without modifying it.

The mapping is explicit on purpose
----------------------------------
The matcher is a case-insensitive **substring** test over eight terms
(``rare_high_impact_variants.py:33-36``). That has a trap: three of the six
loss-of-function effect names SnpEff emits happen to contain a matching
substring, and three do not.

    STOP_GAINED           -> contains "stop_gained"     -> matches by accident
    START_LOST            -> contains "start_lost"      -> matches by accident
    STOP_LOST             -> contains "stop_lost"       -> matches by accident
    FRAME_SHIFT           -> "frameshift" != "frame_shift"        -> no match
    SPLICE_SITE_DONOR     -> "splice_donor" != "splice_site_donor" -> no match
    SPLICE_SITE_ACCEPTOR  -> likewise                              -> no match

So piping raw effect names through produces a *partially* correct count. On the
Berlin challenge data that is 29 of 68 — a number that looks plausible, carries
no warning, and is wrong. A silent zero is at least obviously broken; a silent
29 is not. Both failure modes are reproduced side by side in ``prove_gap.py``.

SO accessions are decorative here: the skill never parses them. They are emitted
anyway so the output is honest about which ontology term was intended.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from legacy_eff import EffAnnotation, parse_eff  # noqa: E402

# SnpEff classic effect name -> (SO accession, Sequence Ontology term).
# Terms are the ones the downstream matcher recognises where a match is intended.
SNPEFF_TO_SO: dict[str, tuple[str, str]] = {
    # HIGH — loss of function
    "STOP_GAINED": ("SO:0001587", "nonsense"),
    "FRAME_SHIFT": ("SO:0001589", "frameshift_variant"),
    "SPLICE_SITE_DONOR": ("SO:0001575", "splice_donor_variant"),
    "SPLICE_SITE_ACCEPTOR": ("SO:0001574", "splice_acceptor_variant"),
    "START_LOST": ("SO:0002012", "start_lost"),
    "STOP_LOST": ("SO:0001578", "stop_lost"),
    # MODERATE / LOW
    "NON_SYNONYMOUS_CODING": ("SO:0001583", "missense_variant"),
    "SYNONYMOUS_CODING": ("SO:0001819", "synonymous_variant"),
    "CODON_DELETION": ("SO:0001822", "inframe_deletion"),
    "CODON_INSERTION": ("SO:0001821", "inframe_insertion"),
    # MODIFIER
    "INTRON": ("SO:0001627", "intron_variant"),
    "EXON": ("SO:0001791", "exon_variant"),
    "UPSTREAM": ("SO:0001631", "upstream_gene_variant"),
    "DOWNSTREAM": ("SO:0001632", "downstream_gene_variant"),
    "UTR_5_PRIME": ("SO:0001623", "5_prime_UTR_variant"),
    "UTR_3_PRIME": ("SO:0001624", "3_prime_UTR_variant"),
    "INTERGENIC": ("SO:0001628", "intergenic_variant"),
}

MODES = ("raw", "passthrough", "mapped")


class UnmappedEffect(KeyError):
    """An effect name absent from SNPEFF_TO_SO.

    Raised rather than guessed. Emitting an unrecognised name would land in the
    passthrough failure mode this module exists to document.
    """


def most_severe(annotations: list[EffAnnotation]) -> EffAnnotation:
    """The annotation the record would be filed under.

    Reduced to exactly one, because the downstream label is the first matching
    term in the matcher's own declaration order rather than the most severe. A
    concatenated consequence string silently reassigns the reported effect.
    """
    return max(annotations, key=lambda a: a.impact_rank)


def build_info(
    annotations: list[EffAnnotation],
    *,
    mode: str,
    raw_eff: str,
    frequency: float | None = None,
    frequency_key: str = "gnomAD_AF",
) -> tuple[str, str | None]:
    """Return ``(info_string, unmapped_effect_or_None)``."""
    if mode == "raw":
        return f"EFF={raw_eff}", None

    best = most_severe(annotations)
    fields: list[str] = []
    unmapped: str | None = None

    if mode == "passthrough":
        # Deliberately wrong: the SnpEff name, unmapped. Reproduces what a
        # plausible-looking one-line "fix" does.
        fields.append(f"MC={best.effect}")
    else:
        try:
            accession, term = SNPEFF_TO_SO[best.effect]
        except KeyError:
            unmapped = best.effect
        else:
            fields.append(f"MC={accession}|{term}")

    if best.gene:
        # The skill takes .split(":")[0].split("|")[0], so a bare symbol is safe.
        fields.append(f"GENEINFO={best.gene}")

    # Frequency is emitted only when a documented value exists. Supplying a
    # placeholder would manufacture the rarity claim this project refuses to make;
    # omitting the key lands the record in the skill's own `frequency_unknown`
    # bucket, which is the honest destination.
    if frequency is not None:
        fields.append(f"{frequency_key}={frequency:.6g}")

    fields.append(f"EFF={raw_eff}")
    return ";".join(fields), unmapped


def load_frequencies(evidence: pathlib.Path | None) -> dict[str, float]:
    """Highest documented population frequency per record, from the VEP cache.

    Keyed on the response's ``input`` field, which echoes the query string we
    sent verbatim. That is an exact join. Keying on ``(chrom, start)`` is not:
    VEP normalises indels and reports ``start`` shifted by one for insertions,
    so an exact positional join silently drops them and a tolerant one risks
    matching a neighbouring record. Both modules that read this cache use the
    same key so their counts cannot drift apart.
    """
    if evidence is None or not evidence.exists():
        return {}
    blob = json.loads(evidence.read_text())
    if blob.get("status") != "ok":
        return {}
    out: dict[str, float] = {}
    for rec in blob.get("records", []):
        key = rec.get("input")
        if not key:
            continue
        best = None
        for col in rec.get("colocated_variants") or []:
            for fval in (col.get("frequencies") or {}).values():
                if isinstance(fval, dict):
                    for val in fval.values():
                        if isinstance(val, (int, float)):
                            best = float(val) if best is None else max(best, float(val))
        if best is not None:
            out[key] = best
    return out


def region_key(row: dict) -> str:
    """The exact string ``fetch_evidence.to_region_string`` posted for this row."""
    ident = row["ID"] if row.get("ID") and row["ID"] != "." else "."
    return f"{row['CHROM']} {row['POS']} {ident} {row['REF']} {row['ALT']} . . ."


VCF_HEADER = """##fileformat=VCFv4.2
##reference=human_g1k_v37.fasta
##source=abstention-ledger/eff_to_info.py mode={mode}
##INFO=<ID=EFF,Number=.,Type=String,Description="legacy SnpEff annotation, verbatim from input">
##INFO=<ID=MC,Number=.,Type=String,Description="molecular consequence (SO term)">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="gene symbol">
##INFO=<ID=gnomAD_AF,Number=1,Type=Float,Description="highest documented population frequency, Ensembl GRCh37 colocated variants (gnomAD exomes r2.1.1); absent when no record exists">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSON
"""


def convert(
    tsv: pathlib.Path,
    out_vcf: pathlib.Path,
    *,
    mode: str,
    evidence: pathlib.Path | None = None,
) -> dict:
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}")

    freqs = load_frequencies(evidence)
    with tsv.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    unmapped: dict[str, int] = {}
    written = 0
    with_freq = 0
    out_vcf.parent.mkdir(parents=True, exist_ok=True)

    with out_vcf.open("w") as out:
        out.write(VCF_HEADER.format(mode=mode))
        for r in rows:
            annotations = parse_eff(r["EFF"])
            if not annotations:
                continue
            freq = freqs.get(region_key(r))
            if freq is not None:
                with_freq += 1
            info, bad = build_info(
                annotations, mode=mode, raw_eff=r["EFF"], frequency=freq
            )
            if bad:
                unmapped[bad] = unmapped.get(bad, 0) + 1
            gt = r["SON_GT_DP_GQ"].split(":")[0]
            out.write(
                "\t".join(
                    [r["CHROM"], r["POS"], r["ID"] or ".", r["REF"], r["ALT"],
                     ".", "PASS", info, "GT", gt]
                )
                + "\n"
            )
            written += 1

    return {
        "mode": mode,
        "records_written": written,
        "records_with_documented_frequency": with_freq,
        "unmapped_effects": unmapped,
        "output": str(out_vcf),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, type=pathlib.Path, help="segregation TSV")
    ap.add_argument("--output", required=True, type=pathlib.Path, help="VCF to write")
    ap.add_argument("--mode", default="mapped", choices=MODES)
    ap.add_argument("--evidence", type=pathlib.Path, help="cached VEP JSON, for documented frequencies")
    args = ap.parse_args()

    stats = convert(args.input, args.output, mode=args.mode, evidence=args.evidence)
    print(json.dumps(stats, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
