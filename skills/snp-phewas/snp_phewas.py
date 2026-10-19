#!/usr/bin/env python3
"""snp-phewas: one variant looked up in every outcome of Pan-UK Biobank and FinnGen R12.

    python skills/snp-phewas/snp_phewas.py --variant rs12740374 --output /tmp/snp_phewas
    python skills/snp-phewas/snp_phewas.py --demo --output /tmp/snp_phewas_demo
    python skills/snp-phewas/snp_phewas.py --variant rs12740374 --all-outcomes --output DIR
    python skills/snp-phewas/snp_phewas.py --variant rs12740374 --quick --output DIR
    python skills/snp-phewas/snp_phewas.py --list-datasets

Writes report.md, result.json, tables/, figures/ and reproducibility/ under --output and
nowhere else (downloaded tabix indexes go to the index cache, `--index-cache`,
`SNP_PHEWAS_INDEX_CACHE`, default ~/.clawbio/snp_phewas/index).
"""

from __future__ import annotations

import argparse
import csv
import logging
import math
import shlex
import sys
import time
from collections import Counter
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = SKILL_DIR.parents[1]
for _p in (SKILL_DIR, _PROJECT_ROOT):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from clawbio.common.report import (  # noqa: E402
    generate_report_footer,
    generate_report_header,
    write_result_json,
)
from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

import snp_phewas_schema as S  # noqa: E402
from snp_phewas_catalog import (  # noqa: E402
    get_dataset,
    load_datasets,
    load_outcome_catalog,
    not_served,
    resolve_dataset_ids,
)
from snp_phewas_harmonize import basis_label  # noqa: E402
from snp_phewas_lookup import SnpPhewasResult, lookup_variant  # noqa: E402
from snp_phewas_replication import VERDICT_LABELS, ReplicationTable, build_replication  # noqa: E402
from snp_phewas_sources import DEFAULT_WORKERS, FINNGEN_ROUTES  # noqa: E402
from snp_phewas_tabix import INDEX_STATS, SourceUnreachable, reset_index_stats  # noqa: E402
from snp_phewas_variant import VariantResolutionError  # noqa: E402

SKILL_NAME = "snp-phewas"
VERSION = "0.1.0"
SCRIPT_REL = "skills/snp-phewas/snp_phewas.py"
PIP_DEPS = ["requests>=2.28", "matplotlib>=3.7", "numpy>=1.24", "pyyaml>=6.0"]
TOP_N = 12

log = logging.getLogger("snp_phewas")


# --- argument handling ---------------------------------------------------------------

INDEX_SIZE_NOTE = (
    "Index files (downloaded once per release into the index cache; measured 2026-10-18 from\nthe sources' Content-Length): default scope 2.08 GB (589 Pan-UK Biobank + 622 FinnGen);\n--all-outcomes 17.9 GB (7,160 + 2,469); --quick none."
)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="snp-phewas",
        description="Variant PheWAS across Pan-UK Biobank and FinnGen R12: every outcome, "
                    "nulls included, on one effect allele, with replication verdicts.",
        epilog=INDEX_SIZE_NOTE, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variant", help="rsID (rs12740374) or chr_pos_ref_alt (1_109274968_G_T)")
    ap.add_argument("--input", help="a text file whose first non-comment line is the variant")
    ap.add_argument("--output", help="output directory")
    ap.add_argument("--demo", action="store_true",
                    help="run offline on the bundled rs12740374 slices (3 + 3 outcomes)")
    ap.add_argument("--build", default="GRCh38", choices=("GRCh38", "GRCh37"),
                    help="genome build of a chr_pos_ref_alt variant (default GRCh38)")
    ap.add_argument("--dataset", default="default",
                    help="comma-separated dataset ids, 'default' (pan_ukbb_eur,finngen_r12) "
                         "or 'all'; see --list-datasets")
    ap.add_argument("--all-outcomes", action="store_true",
                    help="read every outcome in both biobanks, not only those matched across "
                         "the two (the default)")
    ap.add_argument("--quick", action="store_true",
                    help="FinnGen web API only: one request, values only where p < 0.05")
    ap.add_argument("--outcomes",
                    help="restrict reads: dataset:outcome_id,... (other outcomes are listed "
                         "as not part of the lookup)")
    ap.add_argument("--finngen-route", default="hybrid", choices=FINNGEN_ROUTES,
                    help="hybrid (default): API then files for censored endpoints; api; tabix")
    ap.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                    help=f"parallel file reads (default {DEFAULT_WORKERS})")
    ap.add_argument("--index-cache", help="directory for downloaded tabix indexes")
    ap.add_argument("--list-datasets", action="store_true",
                    help="list the curated datasets with license and attribution, then exit")
    return ap


def _variant_from_input(path: str) -> str:
    lines = [ln.strip() for ln in Path(path).read_text(encoding="utf-8").splitlines()]
    lines = [ln for ln in lines if ln and not ln.startswith("#")]
    if not lines:
        raise SystemExit(f"{path}: no variant found (one rsID or chr_pos_ref_alt per file)")
    if len(lines) > 1:
        raise SystemExit(f"{path}: {len(lines)} variants found; this skill looks up one "
                         f"variant per run")
    return lines[0].split()[0]


def _parse_outcomes(spec: str | None) -> dict[str, list[str]] | None:
    if not spec:
        return None
    out: dict[str, list[str]] = {}
    for tok in spec.split(","):
        if ":" not in tok:
            raise SystemExit(f"--outcomes: {tok!r} is not dataset_id:outcome_id")
        ds, oid = tok.split(":", 1)
        out.setdefault(ds.strip(), []).append(oid.strip())
    return out


# --- writers -------------------------------------------------------------------------

def _write_records(recs: list[dict], path: Path) -> None:
    cols = list(recs[0]) if recs else ["harmonized_outcome_id"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in recs:
            w.writerow({k: ("" if v is None else v) for k, v in r.items()})


def format_p(neglog10_p: float | None) -> str:
    """p from -log10 p, without underflow (Pan-UKBB LDL at rs12740374 is p = 1.6e-465)."""
    if neglog10_p is None:
        return ""
    if neglog10_p < 300:
        return f"{10.0 ** -neglog10_p:.2e}"
    e = math.ceil(neglog10_p)
    return f"{10.0 ** (e - neglog10_p):.2f}e-{e}"


def _fmt_p(row: S.OutcomeRow) -> str:
    return format_p(row.neglog10_p)


def _fmt(x: float | None, nd: int = 4) -> str:
    return "" if x is None else f"{x:.{nd}g}"


def _top_rows(result: SnpPhewasResult, dsid: str, n: int) -> list[S.OutcomeRow]:
    rows = [r for r in result.dataset_rows(dsid)
            if r.status in S.TESTED_STATUSES and r.neglog10_p is not None]
    return sorted(rows, key=lambda r: -r.neglog10_p)[:n]


def _scope_sentences(result: SnpPhewasResult, args: argparse.Namespace) -> list[str]:
    """The plain statement of what was read, when it is less than the user may assume."""
    out = []
    if args.demo:
        out.append("**Demo data.** This run reads bundled slices of 3 Pan-UK Biobank outcome "
                   "files and 3 FinnGen endpoint files around rs12740374, offline; the rsID and "
                   "the GRCh37 position use answers recorded from Open Targets and Ensembl. "
                   "Every other outcome is listed as not part of this lookup.")
    if args.quick:
        out.append("**Quick mode: FinnGen web API only.** The API reports values only where "
                   "p < 0.05; every other FinnGen endpoint is labelled 'above the API reporting "
                   "threshold' (tested, value not shown). Pan-UK Biobank was not read and there "
                   "is no replication view. Run without --quick for full values.")
    for dsid in result.dataset_ids:
        prov = result.provenance[dsid]
        if prov.get("scope") == "harmonized" and prov["n_queried"] < prov["n_catalog"]:
            out.append(f"**{get_dataset(dsid).label}: {prov['n_queried']:,} of "
                       f"{prov['n_catalog']:,} outcomes were read**, the ones matched to the other "
                       f"biobank (the outcomes the replication view compares). The other "
                       f"{prov['n_catalog'] - prov['n_queried']:,} are rows with status "
                       f"'not part of this lookup', not nulls. Run with --all-outcomes to read "
                       f"every outcome in both biobanks.")
    return out


def write_report(result: SnpPhewasResult, table: ReplicationTable | None, out: Path,
                 args: argparse.Namespace, wall_s: float) -> Path:
    v = result.variant
    title = f"Variant PheWAS across biobanks: {v.rs_id or v.variant_id_grch38}"
    lines = [generate_report_header(
        title, SKILL_NAME, VERSION,
        extra_metadata={"Variant (GRCh38)": v.variant_id_grch38,
                        "Effect allele": f"{v.grch38.alt} (other allele {v.grch38.ref})",
                        "Datasets": ", ".join(get_dataset(d).label for d in result.dataset_ids)})]
    scope = _scope_sentences(result, args)
    if scope:
        lines += ["## What this answer covers", "", *[f"- {s}" for s in scope], ""]
    lines += ["## Summary", ""]
    for dsid in result.dataset_ids:
        c = result.status_counts(dsid)
        tested = sum(c[s] for s in S.TESTED_STATUSES)
        prov = result.provenance[dsid]
        line = (f"- **{get_dataset(dsid).label}**: tested in {tested:,} of "
                f"{prov['n_queried']:,} outcomes read ({prov['n_catalog']:,} in the catalog)")
        rows = result.dataset_rows(dsid)
        n_bonf = S.phenome_wide_n(rows)
        if n_bonf:
            alpha = 0.05 / n_bonf
            # Low-confidence values (too few carriers) are shown but never counted as a pass.
            n_sig = sum(1 for r in rows
                        if r.status == S.ST_OK and r.p_value is not None and r.p_value < alpha)
            n_low = c[S.ST_LOW_CONFIDENCE]
            # The threshold is phenome-wide only when no outcome was left out of the read.
            where = ("phenome-wide" if not c[S.ST_NOT_QUERIED]
                     else f"over the {n_bonf:,} outcomes looked up")
            line += f"; {n_sig:,} below the threshold 0.05 / {n_bonf:,} ({where})"
            if n_low:
                line += (" (1 low-confidence value is not counted)" if n_low == 1 else
                         f" ({n_low:,} low-confidence values are not counted)")
        lines.append(line + ".")
    if table is not None and table.pairs:
        pair = table.pairs[0]
        c = table.counts(pair)
        lines.append(f"- **Replication** ({get_dataset(pair[0]).label} vs "
                     f"{get_dataset(pair[1]).label}, threshold 0.05 / {table.n_shared[pair]:,} "
                     f"shared outcomes): {c['replicates']:,} replicate, "
                     f"{c['concordant_one_significant']:,} same direction significant in one "
                     f"only, {c['discordant']:,} opposite directions, {c['both_null']:,} null in "
                     f"both.")
    lines += ["", "## Notes", "", *[f"- {n}" for n in result.notes], ""]

    for dsid in result.dataset_ids:
        top = _top_rows(result, dsid, TOP_N)
        if not top:
            continue
        lines += [f"## Strongest associations: {get_dataset(dsid).label}", "",
                  "| Outcome | Outcome id | beta | SE | p | EAF | Status |",
                  "|---|---|---|---|---|---|---|"]
        for r in top:
            beta = _fmt(r.beta) if r.beta is not None else "withheld"
            lines.append(f"| {S.display_label(r.outcome_label)} | `{r.outcome_id}` | {beta} | {_fmt(r.se)} | "
                         f"{_fmt_p(r)} | {_fmt(r.eaf, 3)} | {r.status_label} |")
        lines.append("")

    if table is not None and table.pairs:
        pair = table.pairs[0]
        lines += [f"## Replication across biobanks ({len(table.rows):,} harmonized outcomes)",
                  "", f"Significance threshold: 0.05 / {table.n_shared[pair]:,} outcomes both "
                  f"datasets share. Full table: `tables/replication.tsv`.", ""]
        lines += [f"- {VERDICT_LABELS[k]}: {n:,}" for k, n in table.counts(pair).items() if n]
        shown = [r for r in table.rows if r.verdicts[pair] in
                 ("replicates", "concordant_one_significant", "discordant",
                  "single_dataset_signal", "direction_unknown")][:TOP_N]
        if shown:
            a, b = pair
            lines += ["", f"| Harmonized outcome | {get_dataset(a).label} p | "
                      f"{get_dataset(b).label} p | Verdict |", "|---|---|---|---|"]
            for r in shown:
                ca, cb = r.cells.get(a), r.cells.get(b)
                lines.append(f"| {r.harmonized_label} | {_fmt_p(ca) if ca else ''} | "
                             f"{_fmt_p(cb) if cb else ''} | {VERDICT_LABELS[r.verdicts[pair]]} |")
        basis = Counter(r.harmonization_basis for r in result.rows)
        lines += ["", "Harmonization basis across all rows:", ""]
        lines += [f"- {basis_label(k)}: {n:,}" for k, n in basis.most_common()]
        lines.append("")

    lines += ["## What was read", ""]
    for dsid in result.dataset_ids:
        prov = result.provenance[dsid]
        route = prov.get("route", "none (no position on this build)")
        files = prov.get("n_files_read")
        bits = [f"release `{prov['release']}` ({prov['build']})", f"route {route}"]
        if files is not None:
            bits.append(f"{files:,} files read")
        if prov.get("n_requests"):
            bits.append(f"{prov['n_requests']} API request(s)")
        if prov.get("elapsed_s") is not None:
            bits.append(f"{prov['elapsed_s']:.0f} s")
        lines.append(f"- {get_dataset(dsid).label}: " + "; ".join(bits) + ".")
        hist = ", ".join(f"{S.STATUS_LABELS.get(k, k)}: {n:,}"
                         for k, n in Counter(prov.get("status_counts", {})).most_common())
        lines.append(f"  - Status of every catalog outcome: {hist}.")
    lines.append(f"- Tabix indexes: {INDEX_STATS['downloaded']:,} downloaded "
                 f"({INDEX_STATS['bytes_downloaded'] / 1e6:,.1f} MB), "
                 f"{INDEX_STATS['cache_hits']:,} already cached.")
    lines.append(f"- Wall time: {wall_s:.1f} s.")
    lines += ["", "## Data sources", "", *[f"- {a}" for a in result.attributions()]]
    lines += ["", "## Files", "",
              "- `tables/snp_phewas_rows.tsv`: one row per catalog outcome per dataset, every status.",
              "- `tables/excluded.tsv`: the rows with no value, each with its reason.",
              "- `figures/manhattan.png`: outcome Manhattan, one panel per dataset."]
    if table is not None:
        lines += ["- `tables/replication.tsv`, `figures/replication_matrix.png`: the replication view."]
    lines.append(generate_report_footer())
    path = out / "report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def run(args: argparse.Namespace, argv: list[str]) -> int:
    t0 = time.time()
    out = Path(args.output).expanduser().resolve()
    (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    reset_index_stats()

    kwargs: dict = dict(build=args.build, workers=args.workers, index_cache=args.index_cache,
                        finngen_route=args.finngen_route,
                        scope="all" if args.all_outcomes else "harmonized",
                        outcomes=_parse_outcomes(args.outcomes))
    if args.quick:
        kwargs["finngen_route"] = "api"
    if args.demo:
        from snp_phewas_demo import DEMO_VARIANT, demo_kwargs, demo_readers

        variant = args.variant or DEMO_VARIANT
        kwargs.update(demo_kwargs(readers=demo_readers(finngen_route=kwargs["finngen_route"])))
        datasets = "default"
    else:
        variant = args.variant or _variant_from_input(args.input)
        datasets = args.dataset
        if not args.quick:
            # Before any byte moves, so nobody is surprised by the first run's download.
            size = "17.9 GB" if args.all_outcomes else "2.08 GB"
            print(f"snp-phewas: a first run downloads up to {size} of index files into "
                  f"{args.index_cache or 'the index cache'} (reused afterwards); "
                  f"--quick downloads none.", file=sys.stderr)
    if args.quick:
        datasets = "finngen_r12"
    result = lookup_variant(variant, datasets=datasets, **kwargs)

    produced: list[Path] = []
    rows_tsv = out / "tables" / "snp_phewas_rows.tsv"
    S.write_rows_tsv(result.rows, str(rows_tsv))
    excl_tsv = out / "tables" / "excluded.tsv"
    S.write_rows_tsv(result.excluded, str(excl_tsv))
    produced += [rows_tsv, excl_tsv]
    if result.uncatalogued:
        p = out / "tables" / "uncataloged.tsv"
        S.write_rows_tsv(result.uncatalogued, str(p))
        produced.append(p)

    import snp_phewas_plots as plots

    produced.append(Path(plots.manhattan(result, str(out / "figures" / "manhattan.png"))))
    table = None
    if len(result.dataset_ids) >= 2:
        table = build_replication(result)
        rep = out / "tables" / "replication.tsv"
        _write_records(table.to_records(), rep)
        produced.append(rep)
        produced.append(Path(plots.replication_matrix(
            table, str(out / "figures" / "replication_matrix.png"))))
        if table.pairs and any(r.verdicts[table.pairs[0]] not in ("not_tested_in_one",
                                                                  "not_read_in_one",
                                                                  "low_confidence_in_one")
                               for r in table.rows):
            produced.append(Path(plots.paired_effects(
                table, str(out / "figures" / "paired_effects.png"))))

    wall = time.time() - t0
    produced.append(write_report(result, table, out, args, wall))
    v = result.variant
    summary = {
        "variant_id_grch38": v.variant_id_grch38, "rs_id": v.rs_id,
        "effect_allele": v.grch38.alt, "other_allele": v.grch38.ref,
        "scope": result.scope, "quick": bool(args.quick), "demo": bool(args.demo),
        "wall_time_s": round(wall, 1),
        "datasets": {d: {"n_catalog": p["n_catalog"], "n_queried": p["n_queried"],
                         "status_counts": p.get("status_counts", {})}
                     for d, p in result.provenance.items()},
    }
    if table is not None and table.pairs:
        pair = table.pairs[0]
        summary["replication"] = {"pair": list(pair), "n_shared": table.n_shared[pair],
                                  "alpha": table.alpha[pair], "verdicts": table.counts(pair)}
    # "ok" only when some outcome carries a value: an all-outage run is not a success.
    any_value = any(r.status in S.TESTED_STATUSES for r in result.rows)
    data = {"query": v.query, "resolution_note": v.resolution_note, "notes": result.notes,
            "provenance": result.provenance,
            "index_cache": dict(INDEX_STATS),
            "uncataloged": [r.outcome_id for r in result.uncatalogued]}
    produced.append(write_result_json(
        out, SKILL_NAME, VERSION, summary, data,
        datasets={d: get_dataset(d).release for d in result.dataset_ids},
        status="ok" if any_value else "error", ok=any_value))

    cmd = shlex.join(["python", SCRIPT_REL, *argv])
    produced.append(write_commands_sh(out, cmd))
    produced.append(write_environment_yml(out, "clawbio-snp-phewas", PIP_DEPS))
    write_checksums(produced, out, anchor=out)

    print((out / "report.md").read_text(encoding="utf-8"))
    print(f"wrote {out}")
    return 0


def list_datasets() -> int:
    for dsid, d in load_datasets().items():
        n = len(load_outcome_catalog(dsid))
        default = " (default)" if d.default else ""
        print(f"{dsid}{default}: {d.label}; release {d.release}; {d.build}; {n:,} outcomes; "
              f"license {d.license} ({d.license_verdict}; checked {d.license_verified})")
        print(f"    cite: {d.attribution}")
    for k, why in not_served().items():
        print(f"{k}: not served. {why}")
    return 0


def main(argv: list[str] | None = None) -> int:
    raw = list(sys.argv[1:] if argv is None else argv)
    ap = build_parser()
    args = ap.parse_args(raw)
    logging.basicConfig(level=logging.INFO, format="%(message)s", stream=sys.stderr)
    if args.list_datasets:
        return list_datasets()
    if not args.output:
        ap.error("--output is required")
    if not (args.demo or args.variant or args.input):
        ap.error("give --variant, --input or --demo")
    if args.quick and args.all_outcomes:
        ap.error("--quick reads FinnGen only; it cannot be combined with --all-outcomes")
    if args.quick and args.dataset not in ("default", "finngen_r12"):
        ap.error("--quick reads FinnGen only; drop --dataset or use --dataset finngen_r12")
    try:
        if not args.demo:
            resolve_dataset_ids(args.dataset)
        return run(args, raw)
    except VariantResolutionError as e:
        print(f"snp-phewas: {e}", file=sys.stderr)
        return 2
    except SourceUnreachable as e:
        print(f"snp-phewas: a source could not be reached ({e}); this is an outage, "
              f"not an absence. Try again later.", file=sys.stderr)
        return 3


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
