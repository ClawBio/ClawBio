"""Tests for the abstention ledger.

Run: python3 skills/abstention-ledger/tests/test_abstention_ledger.py
(or via pytest; no plugins or fixtures required)
"""

from __future__ import annotations

import json
import pathlib
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
SKILL = HERE.parent
sys.path.insert(0, str(SKILL))

from abstention_ledger import analyse, infer_origin, load_tsv, reproduce_segregation  # noqa: E402
from gates import evaluate_variant, non_proband_carriers  # noqa: E402
from legacy_eff import parse_eff, summarise  # noqa: E402

FAILURES: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        print(f"  ok   {name}")
    else:
        print(f"  FAIL {name} {detail}")
        FAILURES.append(name)


# --------------------------------------------------------------------------
# legacy EFF parsing
# --------------------------------------------------------------------------

def test_parse_single() -> None:
    eff = "STOP_GAINED(HIGH|NONSENSE|Cga/Tga|R36*|150|PLA2G2C||CODING|NM_001105572.1|1|1)"
    (a,) = parse_eff(eff)
    check("effect parsed", a.effect == "STOP_GAINED", a.effect)
    check("impact parsed", a.impact == "HIGH", a.impact)
    check("gene parsed", a.gene == "PLA2G2C", a.gene)
    check("transcript parsed", a.transcript_id == "NM_001105572.1", a.transcript_id)
    check("aa change parsed", a.aa_change == "R36*", a.aa_change)


def test_parse_multi_transcript() -> None:
    """The real FBXO7 record: HIGH on one transcript, MODERATE on two others."""
    eff = (
        "NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|atG/atA|M115I|522|FBXO7||CODING|NM_012179.3|2|1),"
        "NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|atG/atA|M36I|443|FBXO7||CODING|NM_001033024.1|2|1),"
        "START_LOST(HIGH|MISSENSE|atG/atA|M1I|408|FBXO7||CODING|NM_001257990.1|2|1)"
    )
    s = summarise(eff)
    check("three annotations", len(s.annotations) == 3, str(len(s.annotations)))
    check("max impact HIGH", s.max_impact == "HIGH", s.max_impact)
    check("single gene", s.genes == ["FBXO7"], str(s.genes))
    check("transcript-dependent detected", s.is_impact_transcript_dependent)


def test_empty_fields_survive() -> None:
    eff = "DOWNSTREAM(MODIFIER||||847|CLCN6||CODING|NM_001256959.1||1)"
    (a,) = parse_eff(eff)
    check("modifier impact", a.impact == "MODIFIER", a.impact)
    check("empty codon tolerated", a.codon_change == "", repr(a.codon_change))
    check("aa length kept", a.aa_length == "847", a.aa_length)


def test_two_genes_is_not_disagreement() -> None:
    """A HIGH call in one gene next to a MODIFIER in another is not a conflict."""
    eff = (
        "DOWNSTREAM(MODIFIER||||847|CLCN6||CODING|NM_001256959.1||1),"
        "STOP_LOST(HIGH|MISSENSE|Tga/Cga|*152R|151|NPPA||CODING|NM_006172.3|3|1)"
    )
    s = summarise(eff)
    check("two genes seen", set(s.genes) == {"CLCN6", "NPPA"}, str(s.genes))
    check("no false transcript conflict", not s.is_impact_transcript_dependent)


def test_unparseable_token_skipped() -> None:
    check("garbage skipped", parse_eff("not an annotation") == [])
    check("mixed input keeps good token", len(parse_eff("junk,STOP_GAINED(HIGH|||||G||CODING|NM_1|1|1)")) == 1)


# --------------------------------------------------------------------------
# gates
# --------------------------------------------------------------------------

def test_low_complexity_by_family() -> None:
    rec = {"variant_id": "x", "chrom": "1", "pos": "1000", "genotypes": {}}
    s = summarise("FRAME_SHIFT(HIGH||-/A|-77?|310|OR2DEMO1||CODING|NM_1.1|1|1)")
    v = evaluate_variant(rec, s)
    check("OR family gated", "LOW_COMPLEXITY_LOCUS" in v.codes, str(v.codes))


def test_low_complexity_by_position() -> None:
    rec = {"variant_id": "x", "chrom": "6", "pos": "31000500", "genotypes": {}}
    s = summarise("STOP_GAINED(HIGH|NONSENSE|Cag/Aag|Q90*|265|CLEANGENE||CODING|NM_1.1|2|1)")
    v = evaluate_variant(rec, s)
    check("MHC position gated", "LOW_COMPLEXITY_LOCUS" in v.codes, str(v.codes))


def test_clean_record_passes() -> None:
    rec = {"variant_id": "x", "chrom": "1", "pos": "1000", "genotypes": {}}
    s = summarise("STOP_GAINED(HIGH|NONSENSE|Cga/Tga|R41*|300|CLEANGENE||CODING|NM_1.1|2|1)")
    v = evaluate_variant(rec, s)
    check("clean record reviewable", v.verdict == "REVIEWABLE", str(v.codes))


def test_non_proband_carriers_is_not_a_gate() -> None:
    """Regression: this fired on 100% of records when it was a gate."""
    gts = {"SON": "0/1", "FATHER": "0/1", "SISTER": "0/0", "MOTHER": "0/0"}
    check("father listed", non_proband_carriers(gts) == ["father"], str(non_proband_carriers(gts)))
    rec = {"variant_id": "x", "chrom": "1", "pos": "1000", "genotypes": gts}
    s = summarise("STOP_GAINED(HIGH|NONSENSE|Cga/Tga|R41*|300|CLEANGENE||CODING|NM_1.1|2|1)")
    check("carrier status does not withhold", evaluate_variant(rec, s).verdict == "REVIEWABLE")


# --------------------------------------------------------------------------
# segregation
# --------------------------------------------------------------------------

def test_origin_inference() -> None:
    def rec(f, m):
        return {"genotypes": {"FATHER": f, "MOTHER": m, "SON": "0/1", "SISTER": "0/0"}}

    check("paternal", infer_origin(rec("0/1", "0/0"), father="FATHER", mother="MOTHER") == "paternal")
    check("maternal", infer_origin(rec("0/0", "0/1"), father="FATHER", mother="MOTHER") == "maternal")
    check("ambiguous", infer_origin(rec("0/1", "0/1"), father="FATHER", mother="MOTHER") == "ambiguous")
    check("no carrier", infer_origin(rec("0/0", "0/0"), father="FATHER", mother="MOTHER") == "no_carrier_parent")
    check("hom parent still carries", infer_origin(rec("1/1", "0/0"), father="FATHER", mother="MOTHER") == "paternal")


def test_demo_end_to_end() -> None:
    tsv = SKILL / "examples" / "demo_segregation.tsv"
    records = load_tsv(tsv)
    check("demo has 8 records", len(records) == 8, str(len(records)))

    seg = reproduce_segregation(records)
    check("demo labels reproduce exactly", seg["as_labelled"]["mismatches_vs_supplied"] == 0,
          str(seg["as_labelled"]))

    with tempfile.TemporaryDirectory() as td:
        out = pathlib.Path(td)
        result = analyse(tsv, None, None, out)
        check("report written", (out / "report.md").exists())
        check("result written", (out / "result.json").exists())
        check("ledger table written", (out / "tables" / "abstention_ledger.tsv").exists())
        check("checksums written", (out / "reproducibility" / "checksums.sha256").exists())
        check("4 withheld in demo", result["withheld"] == 4, str(result["withheld"]))
        check("4 reviewable in demo", result["reviewable"] == 4, str(result["reviewable"]))
        check("json round-trips", isinstance(json.loads((out / "result.json").read_text()), dict))
        # Without an evidence layer the report must say "not checked", never
        # "no record found" — the latter claims we looked.
        review_section = (out / "report.md").read_text().split("## 5.")[1].split("## 6.")[0]
        check("unchecked is not reported as absent", "no record found" not in review_section,
              review_section[:200])
        check("unchecked is labelled not checked", "not checked" in review_section)


def main() -> int:
    for fn in [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]:
        print(f"\n{fn.__name__}")
        fn()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILED: {', '.join(FAILURES)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
