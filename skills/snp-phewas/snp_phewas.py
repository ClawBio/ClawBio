#!/usr/bin/env python3
"""snp-phewas: one variant looked up in every outcome of Pan-UK Biobank, FinnGen R12 and
the Million Veteran Program.

    python skills/snp-phewas/snp_phewas.py --variant rs12740374 --output /tmp/snp_phewas
    python skills/snp-phewas/snp_phewas.py --demo --output /tmp/snp_phewas_demo
    python skills/snp-phewas/snp_phewas.py --variant rs12740374 --all-outcomes --output DIR
    python skills/snp-phewas/snp_phewas.py --variant rs12740374 --ancestry AFR --output DIR
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
    ANCESTRY_WORDS,
    datasets_for_ancestry,
    get_dataset,
    label_in_list,
    labels_listed,
    load_datasets,
    load_outcome_catalog,
    not_served,
    resolve_dataset_ids,
)
from snp_phewas_harmonize import basis_label  # noqa: E402
from snp_phewas_lookup import (  # noqa: E402
    SnpPhewasResult,
    lookup_variant,
    matched_outcome_ids,
    matched_scope_sentence,
)
from snp_phewas_replication import (  # noqa: E402
    THREE_WAY_LABELS,
    V_NO_COUNTERPART,
    VERDICT_LABELS,
    ReplicationTable,
    build_replication,
    replication_lines,
)
from snp_phewas_sources import DEFAULT_WORKERS, FINNGEN_ROUTES  # noqa: E402
from snp_phewas_tabix import INDEX_STATS, SourceUnreachable, reset_index_stats  # noqa: E402
from snp_phewas_variant import VariantResolutionError  # noqa: E402

SKILL_NAME = "snp-phewas"
VERSION = "0.2.0"
SCRIPT_REL = "skills/snp-phewas/snp_phewas.py"
PIP_DEPS = ["requests>=2.28", "matplotlib>=3.7", "numpy>=1.24", "pyyaml>=6.0"]
TOP_N = 12

log = logging.getLogger("snp_phewas")


# --- argument handling ---------------------------------------------------------------

INDEX_SIZE_NOTE = (
    "Index files (downloaded once per release into the index cache; summed from the sources'\n"
    "Content-Length on 2026-10-19 for the default three biobanks): default matched scope\n"
    "4.53 GB (1,156 Pan-UK Biobank + 622 FinnGen + 727 Million Veteran Program files);\n"
    "--all-outcomes 20.0 GB (7,160 + 2,469 + the 1,221 readable MVP files); --quick none."
)
# The same numbers, for the message printed before a live run starts.
INDEX_SIZE_DEFAULT = "4.53 GB"
INDEX_SIZE_ALL = "20.0 GB"


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="snp-phewas",
        description="Variant PheWAS across Pan-UK Biobank, FinnGen R12 and the Million Veteran "
                    "Program: every outcome, nulls included, on one effect allele, with "
                    "replication verdicts.",
        epilog=INDEX_SIZE_NOTE, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variant", help="rsID (rs12740374) or chr_pos_ref_alt (1_109274968_G_T)")
    ap.add_argument("--input", help="a text file whose first non-comment line is the variant")
    ap.add_argument("--output", help="output directory")
    ap.add_argument("--demo", action="store_true",
                    help="run offline on the bundled rs12740374 slices (3 + 3 + 2 outcomes)")
    ap.add_argument("--build", default="GRCh38", choices=("GRCh38", "GRCh37"),
                    help="genome build of a chr_pos_ref_alt variant (default GRCh38)")
    ap.add_argument("--dataset", default="default",
                    help="comma-separated dataset ids, 'default' (pan_ukbb_eur,finngen_r12,"
                         "mvp_eur) or 'all'; see --list-datasets")
    ap.add_argument("--ancestry", type=str.upper, choices=sorted(ANCESTRY_WORDS),
                    help="one genetic-ancestry group in every biobank that has it (AFR, AMR, "
                         "CSA, EAS, EUR, MID, or META for MVP's multi-population "
                         "meta-analysis); a biobank without the group is named and left out")
    ap.add_argument("--all-outcomes", action="store_true",
                    help="read every outcome in every biobank, not only those matched across "
                         "them (the default)")
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


def _scope_sentences(result: SnpPhewasResult, args: argparse.Namespace,
                     ancestry_notes: list[str] | None = None) -> list[str]:
    """The plain statement of what was read, when it is less than the user may assume."""
    out = []
    if args.demo:
        out.append("**Demo data.** This run reads bundled slices of 3 Pan-UK Biobank outcome "
                   "files, 3 FinnGen endpoint files and 2 Million Veteran Program study files "
                   "around rs12740374, offline; the rsID and the GRCh37 position use answers "
                   "recorded from Open Targets and Ensembl. Every other outcome is listed as "
                   "not part of this lookup.")
    if args.quick:
        out.append("**Quick mode: FinnGen web API only.** The API reports values only where "
                   "p < 0.05; every other FinnGen endpoint is labeled 'above the API reporting "
                   "threshold' (tested, value not shown). No other biobank was read and there "
                   "is no replication view. Run without --quick for full values.")
    out += [f"**{s}**" for s in (ancestry_notes or [])]
    narrowed = [d for d in result.dataset_ids
                if result.provenance[d].get("scope") == "harmonized"]
    if narrowed:
        out.append(f"**This was the quicker read: {labels_listed(narrowed)} "
                   f"{'was' if len(narrowed) == 1 else 'were'} read only over the outcomes "
                   f"matched across the biobanks, the same rule in each.** Run with "
                   f"--all-outcomes to read every outcome in every biobank.")
        for dsid in narrowed:
            out.append(matched_scope_sentence(label_in_list(dsid), result.dataset_rows(dsid),
                                              result.provenance[dsid])
                       .split(" Run with --all-outcomes")[0])
    return out


def _other_groups_sentence(result: SnpPhewasResult) -> str | None:
    """On a default lookup, the ancestry groups that can be asked for instead, read from the
    dataset catalog, so a user can find them without knowing they exist."""
    ds = load_datasets()
    if not all(ds[d].default for d in result.dataset_ids if d in ds):
        return None
    if not any(ds[d].ancestry == "EUR" for d in result.dataset_ids if d in ds):
        return None  # a FinnGen-only read (--quick) read no European-ancestry group
    served = {s.ancestry for s in ds.values() if not s.default}
    codes = [a for a in ANCESTRY_WORDS if a in served and a not in ("EUR", "META")]
    if not codes:
        return None
    words = [f"{ANCESTRY_WORDS[a]} ({a})" for a in codes]
    listed = words[0] if len(words) == 1 else f"{', '.join(words[:-1])} or {words[-1]}"
    meta = (", or META for the Million Veteran Program's multi-population meta-analysis"
            if "META" in served else "")
    label = result.variant.rs_id or result.variant.variant_id_grch38
    return (f"This lookup read the European-ancestry datasets. Other ancestry groups can be "
            f"read with --ancestry: {listed}{meta}. For example: `python {SCRIPT_REL} "
            f"--variant {label} --ancestry {codes[0]} --output DIR`.")


def _read_counts(result: SnpPhewasResult, dsid: str) -> tuple[int, int, int]:
    """(tested, read, unreadable) counted from the rows: a study the source holds with no
    readable file is never "read", whatever the scope asked for."""
    rows = result.dataset_rows(dsid)
    tested = sum(1 for r in rows if r.status in S.TESTED_STATUSES)
    unreadable = sum(1 for r in rows if r.status == S.ST_NOT_READABLE)
    read = sum(1 for r in rows if r.status not in (S.ST_NOT_QUERIED, S.ST_NOT_READABLE))
    return tested, read, unreadable


def _short(dsid: str) -> str:
    return get_dataset(dsid).short_label


def _replication_rows_table(table: ReplicationTable) -> list[str]:
    """The harmonized outcomes significant in at least one dataset: p per dataset and the
    verdict of each pair (and of all three)."""
    shown = [r for r in table.rows if r.significant_anywhere(table.alpha)][:TOP_N]
    if not shown:
        return []
    ids = table.dataset_ids
    head = ["Harmonized outcome", *[f"{_short(d)} p" for d in ids],
            *[f"{_short(a)} × {_short(b)}" for a, b in table.pairs]]
    if table.n_shared_3 is not None:
        head.append("All three")
    lines = ["", "| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for r in shown:
        cells = [S.display_label(r.harmonized_label)]
        cells += [_fmt_p(c) if (c := r.cells.get(d)) is not None and c.neglog10_p is not None
                  else ("" if c is None else c.status_label) for d in ids]
        cells += ["" if r.verdicts[p] == V_NO_COUNTERPART else VERDICT_LABELS[r.verdicts[p]]
                  for p in table.pairs]
        if table.n_shared_3 is not None:
            cells.append(THREE_WAY_LABELS[r.three_way] if r.three_way else "")
        lines.append("| " + " | ".join(cells) + " |")
    return lines


def write_report(result: SnpPhewasResult, table: ReplicationTable | None, out: Path,
                 args: argparse.Namespace, wall_s: float,
                 ancestry_notes: list[str] | None = None) -> Path:
    v = result.variant
    title = f"Variant PheWAS across biobanks: {v.rs_id or v.variant_id_grch38}"
    lines = [generate_report_header(
        title, SKILL_NAME, VERSION,
        extra_metadata={"Variant (GRCh38)": v.variant_id_grch38,
                        "Effect allele": f"{v.grch38.alt} (other allele {v.grch38.ref})",
                        "Datasets": labels_listed(result.dataset_ids)})]
    scope = _scope_sentences(result, args, ancestry_notes)
    if scope:
        lines += ["## What this answer covers", "", *[f"- {s}" for s in scope], ""]
    lines += ["## Summary", ""]
    for dsid in result.dataset_ids:
        c = result.status_counts(dsid)
        tested, read, unreadable = _read_counts(result, dsid)
        prov = result.provenance[dsid]
        cat = f"{prov['n_catalog']:,} in the catalog"
        if unreadable:
            cat += f", {unreadable:,} of them with no readable file"
        line = (f"- **{get_dataset(dsid).label}**: tested in {tested:,} of "
                f"{read:,} outcomes read ({cat})")
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
        for pair in table.pairs:
            n = table.n_shared[pair]
            c = table.counts(pair)
            if not n:
                lines.append(f"- **Replication, {labels_listed(list(pair))}**: no shared "
                             f"outcomes.")
                continue
            lines.append(f"- **Replication, {labels_listed(list(pair))}** (threshold 0.05 / "
                         f"{n:,} shared outcomes): {c['replicates']:,} replicate, "
                         f"{c['concordant_one_significant']:,} same direction significant in "
                         f"one only, {c['discordant']:,} opposite directions, "
                         f"{c['both_null']:,} null in both.")
        if table.n_shared_3:
            c3 = table.counts_three_way()
            lines.append(f"- **Replication, all three** (threshold 0.05 / {table.n_shared_3:,} "
                         f"outcomes all three share): {c3['replicates_in_all_three']:,} "
                         f"replicate in all three, {c3['replicates_in_two']:,} in two.")
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
        lines += replication_lines(table, heading="##")
        lines += _replication_rows_table(table)
        lines += ["", "Full table: `tables/replication.tsv`."]
        basis = Counter(r.harmonization_basis for r in result.rows)
        lines += ["", "Harmonization basis across all rows:", ""]
        lines += [f"- {basis_label(k)}: {n:,}" for k, n in basis.most_common()]
        lines.append("")

    other = _other_groups_sentence(result)
    if other:
        lines += ["## Other ancestry groups", "", other, ""]

    lines += ["## What was read", ""]
    for dsid in result.dataset_ids:
        prov = result.provenance[dsid]
        route = prov.get("route", "none (no position on this build)")
        files = prov.get("n_files_read")
        bits = [f"release `{prov['release']}` ({prov['build']})", f"route {route}"]
        if files is not None:
            bits.append(f"{files:,} files read")
        if prov.get("n_not_readable"):
            bits.append(f"{prov['n_not_readable']:,} studies with no readable file (not read)")
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
        lines += ["- `tables/replication.tsv`, `figures/replication_matrix.png`: the replication "
                  "view (every dataset pair, and all three when three datasets are read).",
                  "- `figures/paired_effects.png`: effect against effect, one row per dataset "
                  "pair with shared outcomes."]
    lines.append(generate_report_footer())
    path = out / "report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _resolve_datasets(args: argparse.Namespace) -> tuple[str | list[str], list[str]]:
    """(the datasets to read, the sentences that say what an ancestry choice substituted)."""
    if args.ancestry:
        if args.dataset != "default":
            raise SystemExit("snp-phewas: give --ancestry or --dataset, not both")
        return datasets_for_ancestry(args.ancestry)
    return args.dataset, []


def _start_message(datasets: str | list[str], all_outcomes: bool, index_cache: str | None) -> str:
    """What a live run will read and download, said before any byte moves."""
    ids = resolve_dataset_ids(datasets)
    where = index_cache or "the index cache"
    if all_outcomes or len({get_dataset(d).family for d in ids}) < 2:
        n = [sum(1 for o in load_outcome_catalog(d) if o.access_status in (None, "harmonised"))
             for d in ids]
        size = (f"up to {INDEX_SIZE_ALL} of index files (the default three biobanks)"
                if set(ids) == {"pan_ukbb_eur", "finngen_r12", "mvp_eur"}
                else "one index file per outcome file")
        return (f"snp-phewas: reading every readable outcome ({', '.join(f'{k:,}' for k in n)} "
                f"files in {labels_listed(ids)}). A first run downloads {size} into {where} "
                f"(reused afterwards); the time of a full read is not measured. --quick "
                f"downloads none.")
    n = [len(matched_outcome_ids(d, ids)) for d in ids]
    default3 = set(ids) == {"pan_ukbb_eur", "finngen_r12", "mvp_eur"}
    size = f"{INDEX_SIZE_DEFAULT} of index files" if default3 else "their index files"
    timing = (" Measured on a Mac mini for the default three biobanks: 159 s and 613 s with an "
              "empty index cache (the spread is the GWAS Catalog's download speed), 57 s and "
              "66 s once the indexes are cached." if default3 else "")
    return (f"snp-phewas: reading the outcomes matched across the biobanks "
            f"({', '.join(f'{k:,}' for k in n)} files in {labels_listed(ids)}). A first run "
            f"downloads {size} into {where} (reused afterwards).{timing} "
            f"--all-outcomes reads everything; --quick downloads none.")


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
    datasets, ancestry_notes = _resolve_datasets(args)
    if args.demo:
        from snp_phewas_demo import DEMO_VARIANT, demo_kwargs, demo_readers

        variant = args.variant or DEMO_VARIANT
        kwargs.update(demo_kwargs(readers=demo_readers(finngen_route=kwargs["finngen_route"])))
        datasets, ancestry_notes = "default", []
    else:
        variant = args.variant or _variant_from_input(args.input)
        if not args.quick:
            # Before any byte moves, so nobody is surprised by the first run's download.
            print(_start_message(datasets, args.all_outcomes, args.index_cache),
                  file=sys.stderr)
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
        if table.pairs and any(table.n_shared.values()):
            produced.append(Path(plots.paired_effects(
                table, str(out / "figures" / "paired_effects.png"))))

    wall = time.time() - t0
    produced.append(write_report(result, table, out, args, wall, ancestry_notes))
    v = result.variant
    summary = {
        "variant_id_grch38": v.variant_id_grch38, "rs_id": v.rs_id,
        "effect_allele": v.grch38.alt, "other_allele": v.grch38.ref,
        "scope": result.scope, "quick": bool(args.quick), "demo": bool(args.demo),
        "ancestry": args.ancestry or None,
        "wall_time_s": round(wall, 1),
        "datasets": {d: {"n_catalog": p["n_catalog"], "n_queried": p["n_queried"],
                         "n_read": _read_counts(result, d)[1],
                         "status_counts": p.get("status_counts", {})}
                     for d, p in result.provenance.items()},
    }
    if table is not None and table.pairs:
        summary["replication"] = {
            "pairs": [{"pair": list(pair), "n_shared": table.n_shared[pair],
                       "alpha": table.alpha[pair], "verdicts": table.counts(pair)}
                      for pair in table.pairs]}
        if table.n_shared_3 is not None:
            summary["replication"]["three_way"] = {
                "n_shared": table.n_shared_3, "alpha": table.alpha_3,
                "classes": table.counts_three_way()}
    # "ok" only when some outcome carries a value: an all-outage run is not a success.
    any_value = any(r.status in S.TESTED_STATUSES for r in result.rows)
    data = {"query": v.query, "resolution_note": v.resolution_note, "notes": result.notes,
            "ancestry_notes": ancestry_notes,
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
    if args.demo and (args.ancestry or args.dataset != "default"):
        ap.error("--demo reads the bundled European-ancestry slices only; drop --ancestry "
                 "and --dataset")
    if args.quick and args.ancestry:
        ap.error("--quick reads FinnGen only (one population); drop --ancestry")
    if args.ancestry and args.dataset != "default":
        ap.error("give --ancestry or --dataset, not both")
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
