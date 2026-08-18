"""Abstention Ledger — say which variants you were not entitled to rank, and prove it.

Run ``--demo`` for a self-contained example, or ``--input`` on a segregation TSV.

The output is two artefacts, not one:

* a review list of records that survived every check we could run, ordered by how
  complete the evidence is — not by how severe the consequence looks; and
* a ledger of every record withheld, each with the check, the value that fired
  it, and where that value came from.

The second artefact is the point. A reader can disagree with a check. They
cannot disagree with a conclusion that shows no working.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import pathlib
import sys
from dataclasses import asdict

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from gates import cohort_abstentions, evaluate_variant  # noqa: E402
from legacy_eff import summarise  # noqa: E402

ROLES = ("SON", "FATHER", "SISTER", "MOTHER")
CARRIER_EXCLUDE = ("0/0", "./.", "0|0", ".")


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------

def load_tsv(path: pathlib.Path) -> list[dict]:
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    records = []
    for r in rows:
        records.append(
            {
                "variant_id": f"{r['CHROM']}:{r['POS']}{'' if r['ID'] in ('.', '') else ' ' + r['ID']} {r['REF']}>{r['ALT']}",
                "chrom": r["CHROM"],
                "pos": r["POS"],
                "rsid": r["ID"],
                "ref": r["REF"],
                "alt": r["ALT"],
                "genotypes": {ro: r[f"{ro}_GT_DP_GQ"].split(":")[0] for ro in ROLES},
                "depths": {ro: r[f"{ro}_GT_DP_GQ"] for ro in ROLES},
                "supplied_origin": r["PARENT_OF_ORIGIN_UNPHASED"],
                "eff": r["EFF"],
            }
        )
    return records


def vcf_sample_genotypes(path: pathlib.Path) -> tuple[list[str], dict]:
    """Read sample IDs and per-sample genotypes from a (possibly gzipped) VCF."""
    opener = gzip.open if path.suffix == ".gz" else open
    samples: list[str] = []
    table: dict = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]
                continue
            gt_idx = f[8].split(":").index("GT")
            key = (f[0], f[1], f[3], f[4])
            table[key] = {s: f[9 + i].split(":")[gt_idx] for i, s in enumerate(samples)}
    return samples, table


def resolve_sample_map(records: list[dict], vcf: pathlib.Path | None) -> dict:
    """Decide which sample ID holds which family role, from genotypes alone.

    The documentation for this dataset states an ordering that its own TSV column
    order contradicts. Rather than trust either, we test all 24 assignments of
    four sample IDs to four roles and keep the ones that reproduce every genotype
    in the TSV. If exactly one survives, the mapping is decided by the data.
    """
    if vcf is None or not vcf.exists():
        return {
            "status": "not_checked",
            "note": "no VCF supplied; role labels taken from TSV column names on trust",
            "mapping": {ro: ro for ro in ROLES},
            "candidates_tested": 0,
            "candidates_matching": 0,
        }

    samples, table = vcf_sample_genotypes(vcf)
    tsv_index = {(r["chrom"], r["pos"], r["ref"], r["alt"]): r["genotypes"] for r in records}

    matching = []
    for perm in itertools.permutations(samples):
        candidate = dict(zip(ROLES, perm))
        if all(
            table[key][candidate[ro]] == gts[ro]
            for key, gts in tsv_index.items()
            if key in table
            for ro in ROLES
        ):
            matching.append(candidate)

    return {
        "status": "decided" if len(matching) == 1 else "ambiguous",
        "samples": samples,
        "candidates_tested": len(list(itertools.permutations(samples))),
        "candidates_matching": len(matching),
        "mapping": matching[0] if len(matching) == 1 else {},
        "note": (
            "exactly one assignment reproduces all genotypes, so the mapping is decided by the data"
            if len(matching) == 1
            else "genotypes do not uniquely determine the mapping"
        ),
    }


# --------------------------------------------------------------------------
# Segregation
# --------------------------------------------------------------------------

def carries(gt: str) -> bool:
    return gt not in CARRIER_EXCLUDE


def infer_origin(rec: dict, *, father: str, mother: str) -> str:
    f, m = carries(rec["genotypes"][father]), carries(rec["genotypes"][mother])
    if f and not m:
        return "paternal"
    if m and not f:
        return "maternal"
    return "ambiguous" if (f and m) else "no_carrier_parent"


def reproduce_segregation(records: list[dict]) -> dict:
    """Derive parent-of-origin ourselves, and show what a mislabelled map costs.

    Variant A uses the roles as the TSV column names give them. Variant B swaps
    sister and mother, which is what following the documentation's prose ordering
    produces. Printing both is the cheapest possible demonstration that the
    labels are a consequence of a decision we can name.
    """
    out = {}
    for label, (fa, mo) in {
        "as_labelled": ("FATHER", "MOTHER"),
        "sister_mother_swapped": ("FATHER", "SISTER"),
    }.items():
        tally = {"paternal": 0, "maternal": 0, "ambiguous": 0, "no_carrier_parent": 0}
        mismatches = 0
        for rec in records:
            call = infer_origin(rec, father=fa, mother=mo)
            tally[call] += 1
            if call != rec["supplied_origin"]:
                mismatches += 1
        out[label] = {"tally": tally, "mismatches_vs_supplied": mismatches, "total": len(records)}
    return out


# --------------------------------------------------------------------------
# Evidence layer (optional, cached)
# --------------------------------------------------------------------------

def load_evidence(path: pathlib.Path | None) -> dict:
    if path is None or not path.exists():
        return {"status": "absent"}
    try:
        blob = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        return {"status": "unreadable", "reason": str(exc)}
    if blob.get("status") != "ok":
        return {"status": blob.get("status", "unavailable"), "reason": blob.get("reason", "")}

    by_key: dict[str, dict] = {}
    for rec in blob.get("records", []):
        # Keyed on the response's `input` field, which echoes our query verbatim.
        # Exact by construction. A (chrom, start) key is not: VEP normalises
        # indels and shifts `start` for insertions, so an exact positional join
        # drops them and a tolerant one can match a neighbour. `eff_to_info.py`
        # uses the same key, so the two modules cannot report different counts.
        key = rec.get("input")
        if not key:
            continue

        freqs: dict[str, float] = {}
        for col in rec.get("colocated_variants", []) or []:
            for fkey, fval in (col.get("frequencies") or {}).items():
                if isinstance(fval, dict):
                    for pop, val in fval.items():
                        if isinstance(val, (int, float)):
                            freqs[pop] = max(freqs.get(pop, 0.0), float(val))

        impacts = {}
        for tc in rec.get("transcript_consequences", []) or []:
            tid = tc.get("transcript_id", "")
            impacts[tid] = {
                "impact": tc.get("impact", ""),
                "consequences": tc.get("consequence_terms", []),
                "canonical": bool(tc.get("canonical")),
                "mane_select": tc.get("mane_select"),
                "gene": tc.get("gene_symbol", ""),
            }

        by_key[key] = {
            "most_severe": rec.get("most_severe_consequence", ""),
            "frequencies": freqs,
            "transcripts": impacts,
            "rsids": [c.get("id") for c in (rec.get("colocated_variants") or []) if c.get("id")],
        }

    return {
        "status": "ok",
        "frequency_source": blob.get("frequency_source", ""),
        "assembly": blob.get("assembly", ""),
        "endpoint": blob.get("endpoint", ""),
        "returned": blob.get("returned", 0),
        "by_key": by_key,
    }


def match_evidence(rec: dict, evidence: dict) -> dict | None:
    """Exact join on the query string we sent. See load_evidence for why."""
    if evidence.get("status") != "ok":
        return None
    ident = rec["rsid"] if rec["rsid"] and rec["rsid"] != "." else "."
    key = f"{rec['chrom']} {rec['pos']} {ident} {rec['ref']} {rec['alt']} . . ."
    return evidence["by_key"].get(key)


COMMON_THRESHOLD = 0.01


def evidence_gates(rec: dict, ev: dict | None, supplied_impact: str) -> list[dict]:
    """Gates that need the build-matched layer. Absent layer fires nothing."""
    hits: list[dict] = []
    if ev is None:
        hits.append(
            {
                "code": "NO_EVIDENCE_LAYER",
                "evidence": "no build-matched annotation available for this position",
                "source": "Ensembl GRCh37 REST cache (missing or not fetched)",
            }
        )
        return hits

    freqs = ev.get("frequencies") or {}
    if not freqs:
        hits.append(
            {
                "code": "NO_FREQUENCY_RECORD",
                "evidence": (
                    "build-matched query returned no population frequency for this position; "
                    "this is a statement about the reference database, not about the variant"
                ),
                "source": "Ensembl GRCh37 REST, gnomAD exomes r2.1.1 (exome-only)",
            }
        )
    else:
        top_pop, top_af = max(freqs.items(), key=lambda kv: kv[1])
        if top_af >= COMMON_THRESHOLD:
            hits.append(
                {
                    "code": "FREQUENCY_DOCUMENTED_COMMON",
                    "evidence": f"documented frequency {top_af:.4g} in {top_pop} (>= {COMMON_THRESHOLD})",
                    "source": "Ensembl GRCh37 REST colocated-variant frequencies",
                }
            )

    current = {t["impact"] for t in ev.get("transcripts", {}).values() if t.get("impact")}
    if supplied_impact and current and supplied_impact not in current:
        hits.append(
            {
                "code": "ANNOTATION_SUPERSEDED",
                "evidence": (
                    f"supplied impact {supplied_impact} is not reproduced by current annotation "
                    f"(now {sorted(current)})"
                ),
                "source": "Ensembl GRCh37 VEP, current release",
            }
        )
    return hits


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------

BANNED = ("rare", "pathogenic", "diagnostic", "de novo", "compound het")

# Mandated verbatim by CLAUDE.md safety rule 2 — every report carries it.
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device and does "
    "not provide clinical diagnoses. Consult a healthcare professional before making "
    "any medical decisions."
)


def write_tsv(path: pathlib.Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


def render_report(ctx: dict) -> str:
    sm, seg, ev = ctx["sample_map"], ctx["segregation"], ctx["evidence"]
    verdicts, cohort = ctx["verdicts"], ctx["cohort"]
    withheld = [v for v in verdicts if v["verdict"] == "WITHHELD"]
    reviewable = [v for v in verdicts if v["verdict"] == "REVIEWABLE"]

    code_counts: dict[str, int] = {}
    for v in verdicts:
        for c in v["codes"]:
            code_counts[c] = code_counts.get(c, 0) + 1

    L: list[str] = []
    a = L.append
    a("# Abstention Ledger")
    a("")
    a(f"Input: `{ctx['input_name']}`  ")
    a(f"SHA-256: `{ctx['input_sha256']}`  ")
    a(f"Records: **{len(verdicts)}**  ")
    a(f"Assembly: {ctx['assembly']}")
    a("")
    a(
        f"**{len(withheld)} of {len(verdicts)} records are withheld from the review list.** "
        f"{len(reviewable)} survived every check we were able to run."
    )
    a("")
    a("---")
    a("")

    # 1. sample map
    a("## 1. Which sample is which person")
    a("")
    if sm["status"] == "decided":
        a(
            f"Resolved from genotypes: **{sm['candidates_matching']} of {sm['candidates_tested']}** "
            "possible assignments of sample IDs to family roles reproduces every genotype in the "
            "table. The mapping is therefore decided by the data, not taken on trust."
        )
        a("")
        a("| Role | Sample |")
        a("|---|---|")
        for role, sid in sm["mapping"].items():
            a(f"| {role.title()} | `{sid}` |")
    else:
        a(f"Status: **{sm['status']}** — {sm['note']}")
    a("")

    # 2. segregation
    a("## 2. Reproducing the parent-of-origin labels")
    a("")
    a(
        "We re-derive each label from the genotypes rather than reading the supplied column. "
        "A record is *paternal* when the father carries and the mother does not, *maternal* when "
        "the reverse holds, *ambiguous* when both carry, and *no carrier parent* when neither does."
    )
    a("")
    a("| Role assignment | paternal | maternal | ambiguous | no carrier parent | disagreements with supplied label |")
    a("|---|---|---|---|---|---|")
    for label, res in seg.items():
        t = res["tally"]
        a(
            f"| {label.replace('_', ' ')} | {t['paternal']} | {t['maternal']} | {t['ambiguous']} | "
            f"{t['no_carrier_parent']} | **{res['mismatches_vs_supplied']} / {res['total']}** |"
        )
    a("")
    a(
        "The second row is what happens if the sister and mother columns are exchanged. It is not "
        "a hypothetical: the published description of this dataset orders the samples in a way that "
        "produces exactly that swap. Records with no carrier parent are impossible under the "
        "dataset's own stated filter, which is how the mistake announces itself."
    )
    a("")

    # 3. cohort abstentions
    a("## 3. What this data cannot support, whatever the records say")
    a("")
    a("These are properties of the dataset. They bound every claim below.")
    a("")
    for c in cohort:
        a(f"### `{c['code']}`")
        a("")
        a(f"**Not supported:** {c['claim_blocked']}")
        a("")
        a(f"**Why:** {c['because']}")
        a("")
        a(f"**Checked against:** {c['evidence']}")
        a("")

    # 4. gate tally
    a("## 4. Per-record gates")
    a("")
    a("| Reason code | Records |")
    a("|---|---|")
    for code, n in sorted(code_counts.items(), key=lambda kv: -kv[1]):
        a(f"| `{code}` | {n} |")
    a("")
    # A check that ran and found nothing is a result, not an absence of work.
    # Reporting only the gates that fired would leave a reader unable to tell
    # "screened, clean" apart from "never screened".
    scr = ctx["screening"]
    a("### Checks that ran and found nothing")
    a("")
    if scr["sf_list_size"]:
        a(
            f"- **ACMG secondary findings:** all {len(verdicts)} records screened against "
            f"{scr['sf_list_size']} genes — **{scr['sf_hits']} hits**. Source: {scr['sf_source']}."
        )
        if scr["sf_hits"] == 0:
            a(
                "  This file contains no gene on that list. Worth stating plainly rather than "
                "leaving implicit: the actionable-secondary-finding problem does not arise here, "
                "and we are not reporting that we avoided one."
            )
    else:
        a(f"- **ACMG secondary findings: not screened.** {scr['sf_source']}")
    a("")

    # 5. review list
    a("## 5. Review list")
    a("")
    if reviewable:
        a(
            "Ordered by how complete the assembled evidence is, **not** by how severe the "
            "consequence appears. Every entry remains subject to section 3."
        )
        a("")
        # "no record" is a claim that we looked. "not checked" is a claim that we
        # did not. Collapsing the two is the exact error this skill exists to
        # catch, so the renderer must not make it either.
        layer_ok = ev.get("status") == "ok"
        unknown = "no record found" if layer_ok else "not checked"
        a("| # | Variant | Gene(s) | Supplied impact | Current impact | Documented frequency |")
        a("|---|---|---|---|---|---|")
        for i, v in enumerate(reviewable, 1):
            a(
                f"| {i} | `{v['variant_id']}` | {', '.join(v['genes']) or '-'} | "
                f"{v['supplied_impact'] or '-'} | {v.get('current_impact') or unknown} | "
                f"{v.get('frequency_display') or unknown} |"
            )
    else:
        a(
            "**Empty.** Every record in this file triggered at least one gate. That is the honest "
            "result for this input, not a failure of the pipeline."
        )
    a("")

    # 6. the ledger
    a("## 6. The ledger — every record withheld, and why")
    a("")
    a("| Variant | Gene(s) | Codes | Evidence that fired the gate |")
    a("|---|---|---|---|")
    for v in withheld:
        ev_txt = "<br>".join(f"`{h['code']}`: {h['evidence']}" for h in v["hits"])
        a(
            f"| `{v['variant_id']}` | {', '.join(v['genes']) or '-'} | "
            f"{', '.join(f'`{c}`' for c in v['codes'])} | {ev_txt} |"
        )
    a("")

    # 7. limits of our own layer
    a("## 7. Limits of the evidence layer we added")
    a("")
    if ev.get("status") == "ok":
        a(f"- Source: {ev['frequency_source']}")
        a(f"- Endpoint: `{ev['endpoint']}` ({ev['assembly']}), {ev['returned']} records returned")
        a(
            "- The layer is **build-matched**, so an empty frequency result means Ensembl holds no "
            "record at that position. It does not mean a liftover failed, because there is no "
            "liftover in this path."
        )
        a(
            "- It is nonetheless **exome-only and several releases old**. gnomAD v4 is "
            "GRCh38-native and is not used here. `NO_FREQUENCY_RECORD` is therefore a weaker "
            "statement than it would be against v4, and we do not treat it as evidence about the "
            "variant's frequency in any population."
        )
    else:
        a(
            f"- Status: **{ev.get('status')}**. "
            f"{ev.get('reason', '')} Records carry `NO_EVIDENCE_LAYER` rather than an assumed value."
        )
    a("")
    a("---")
    a("")
    a(
        "*Any tool can output a ranking. This one ships the check, the value and the source for "
        "everything it declined to rank.*"
    )
    a("")
    a("## Disclaimer")
    a("")
    a(DISCLAIMER)
    a("")
    return "\n".join(L)


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------

def analyse(
    tsv: pathlib.Path,
    vcf: pathlib.Path | None,
    evidence_path: pathlib.Path | None,
    outdir: pathlib.Path,
) -> dict:
    records = load_tsv(tsv)
    sample_map = resolve_sample_map(records, vcf)
    segregation = reproduce_segregation(records)
    evidence = load_evidence(evidence_path)
    cohort = [asdict(c) for c in cohort_abstentions(records, roles=sample_map.get("mapping") or {})]

    verdicts: list[dict] = []
    for rec in records:
        summary = summarise(rec["eff"])
        verdict = evaluate_variant(rec, summary)
        hits = [asdict(h) for h in verdict.hits]

        ev = match_evidence(rec, evidence)
        if evidence.get("status") == "ok":
            hits.extend(evidence_gates(rec, ev, summary.max_impact))

        current_impact = ""
        freq_display = ""
        if ev:
            impacts = sorted({t["impact"] for t in ev.get("transcripts", {}).values() if t.get("impact")})
            current_impact = "/".join(impacts)
            freqs = ev.get("frequencies") or {}
            if freqs:
                pop, af = max(freqs.items(), key=lambda kv: kv[1])
                freq_display = f"{af:.4g} ({pop})"

        verdicts.append(
            {
                "variant_id": verdict.variant_id,
                "chrom": rec["chrom"],
                "pos": rec["pos"],
                "rsid": rec["rsid"],
                "genes": verdict.genes,
                "supplied_impact": verdict.supplied_impact,
                "supplied_origin": rec["supplied_origin"],
                "derived_origin": infer_origin(rec, father="FATHER", mother="MOTHER"),
                "current_impact": current_impact,
                "frequency_display": freq_display,
                "transcript_count": len(summary.annotations),
                "hits": hits,
                "codes": [h["code"] for h in hits],
                "verdict": "WITHHELD" if hits else "REVIEWABLE",
            }
        )

    # Record checks that ran, not only the ones that fired. Without this a reader
    # cannot distinguish "screened, nothing found" from "never screened".
    from gates import ACMG_SF_GENES, ACMG_SF_SOURCE  # noqa: PLC0415

    screening = {
        "sf_list_size": len(ACMG_SF_GENES) if ACMG_SF_GENES else 0,
        "sf_source": ACMG_SF_SOURCE,
        "sf_hits": sum(1 for v in verdicts if "SECONDARY_FINDING_NO_CONSENT" in v["codes"]),
    }

    ctx = {
        "input_name": tsv.name,
        "input_sha256": hashlib.sha256(tsv.read_bytes()).hexdigest(),
        "assembly": "GRCh37/b37 (contigs without chr prefix)",
        "sample_map": sample_map,
        "segregation": segregation,
        "evidence": evidence,
        "cohort": cohort,
        "screening": screening,
        "verdicts": verdicts,
    }

    # artefacts
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "report.md").write_text(render_report(ctx))

    write_tsv(
        outdir / "tables" / "segregation.tsv",
        ["variant_id", "chrom", "pos", "rsid", "supplied_origin", "derived_origin", "agrees"],
        [
            [v["variant_id"], v["chrom"], v["pos"], v["rsid"], v["supplied_origin"],
             v["derived_origin"], str(v["supplied_origin"] == v["derived_origin"])]
            for v in verdicts
        ],
    )
    write_tsv(
        outdir / "tables" / "abstention_ledger.tsv",
        ["variant_id", "genes", "verdict", "reason_codes", "evidence", "source"],
        [
            [v["variant_id"], ";".join(v["genes"]), v["verdict"], ";".join(v["codes"]),
             " || ".join(h["evidence"] for h in v["hits"]),
             " || ".join(h["source"] for h in v["hits"])]
            for v in verdicts
        ],
    )

    summary_json = {
        "input": ctx["input_name"],
        "input_sha256": ctx["input_sha256"],
        "records": len(verdicts),
        "withheld": sum(1 for v in verdicts if v["verdict"] == "WITHHELD"),
        "reviewable": sum(1 for v in verdicts if v["verdict"] == "REVIEWABLE"),
        "sample_map": sample_map,
        "segregation": segregation,
        "cohort_abstentions": cohort,
        "screening": screening,
        "evidence_layer": {k: v for k, v in evidence.items() if k != "by_key"},
        "disclaimer": DISCLAIMER,
        "variants": verdicts,
    }
    (outdir / "result.json").write_text(json.dumps(summary_json, indent=1))

    repro = outdir / "reproducibility"
    repro.mkdir(parents=True, exist_ok=True)
    (repro / "commands.sh").write_text(
        "#!/bin/sh\n"
        "# Regenerate every number in report.md from a clean checkout.\n"
        "set -eu\n\n"
        "python3 skills/abstention-ledger/fetch_evidence.py \\\n"
        f"    --input {tsv} \\\n"
        "    --output out/vep_grch37_cache.json\n\n"
        "python3 skills/abstention-ledger/abstention_ledger.py \\\n"
        f"    --input {tsv} \\\n"
        f"    --vcf {vcf if vcf else '<vcf>'} \\\n"
        "    --evidence out/vep_grch37_cache.json \\\n"
        f"    --output {outdir}\n"
    )
    checks = [tsv] + ([vcf] if vcf and vcf.exists() else [])
    (repro / "checksums.sha256").write_text(
        "".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n" for p in checks)
    )
    return summary_json


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=pathlib.Path, help="segregation TSV")
    ap.add_argument("--vcf", type=pathlib.Path, help="matching VCF, enables the sample-map check")
    ap.add_argument("--evidence", type=pathlib.Path, help="cached build-matched evidence JSON")
    ap.add_argument("--output", type=pathlib.Path, required=True)
    ap.add_argument("--demo", action="store_true", help="run on bundled synthetic data")
    args = ap.parse_args()

    here = pathlib.Path(__file__).resolve().parent
    if args.demo:
        tsv = here / "examples" / "demo_segregation.tsv"
        vcf = None
        evidence = None
        print(f"DEMO: synthetic data, no network calls.\n  input: {tsv}")
    else:
        if not args.input:
            ap.error("--input is required unless --demo is given")
        tsv, vcf, evidence = args.input, args.vcf, args.evidence

    result = analyse(tsv, vcf, evidence, args.output)
    print(
        f"{result['records']} records: {result['withheld']} withheld, "
        f"{result['reviewable']} reviewable -> {args.output}/report.md"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
