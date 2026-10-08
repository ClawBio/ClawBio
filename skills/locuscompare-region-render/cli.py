"""locuscompare-region-render — two-trait regional colocalization visualization.

Reads a JSON or YAML config, validates against the schema, builds a
LocusCompareSpec, picks the right fetcher / LD client / gene-track source,
and runs the render pipeline.

CLI follows the standard skill convention: `--input <file> --output <dir>
--demo`.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import yaml


# Resolve sibling-skill imports for the fork's flat layout. Inject this
# skill's own dir + each sibling-skill dir onto sys.path.
_SKILL_DIR = Path(__file__).resolve().parent
_SKILLS_ROOT = _SKILL_DIR.parent
for _p in (
    _SKILL_DIR,
    _SKILLS_ROOT / "eqtl-catalogue-region-fetch",
    _SKILLS_ROOT / "gwas-catalog-region-fetch",
    _SKILLS_ROOT / "ld-1000g-region-compute",
):
    if _p.exists() and str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="locuscompare",
        description=(
            "Two-trait regional colocalization visualization in the Liu 2019 "
            "LocusCompare convention. Renders a 4-panel publication-grade plot "
            "from a config describing the lead variant, exposure + outcome "
            "sumstats, and optional LD + gene-track sources."
        ),
    )
    parser.add_argument(
        "--input",
        type=Path,
        help="Path to the run config (JSON or YAML). See INPUT_SCHEMA.md.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output directory; created if missing. Required unless --list-demos.",
    )
    parser.add_argument(
        "--demo",
        nargs="?",
        const="__default__",
        default=None,
        metavar="NAME",
        help="Run a bundled demo. Bare --demo runs the default; pass a name "
             "(e.g. --demo 02_eqtl_catalogue_x_gwas_catalog) to choose a "
             "specific one. See --list-demos.",
    )
    parser.add_argument(
        "--list-demos",
        action="store_true",
        help="List bundled demos in this skill's examples/ directory.",
    )
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
        try:
            config_path = _resolve_demo_path(args.demo)
        except FileNotFoundError as e:
            print(f"error: {e}", file=sys.stderr)
            return 2
        config = _load_config(config_path)
        config_dir = config_path.parent
        print(f"info: using bundled demo {config_path.parent.name}/{config_path.name}", file=sys.stderr)
    else:
        config = _load_config(args.input)
        config_dir = args.input.parent

    return _run(config=config, config_dir=config_dir, output=args.output)


def _examples_dir() -> Path:
    return Path(__file__).resolve().parent / "examples"


def _list_demos() -> list[Path]:
    """Each demo for the orchestrator is a subdirectory containing
    `config.{json,yaml,yml}`. Return the sorted list of those configs."""
    out: list[Path] = []
    base = _examples_dir()
    if not base.is_dir():
        return out
    for d in sorted(base.iterdir()):
        if not d.is_dir() or d.name in {"recipes", "chains"}:
            continue
        for cand in ("config.yaml", "config.yml", "config.json"):
            p = d / cand
            if p.is_file():
                out.append(p)
                break
    return out


def _resolve_demo_path(name: str) -> Path:
    """Map --demo NAME to a bundled config path. Special name `__default__`
    (bare --demo) picks `02_eqtl_catalogue_x_gwas_catalog/`'s config (or the
    first listed demo)."""
    base = _examples_dir()
    if name == "__default__":
        for cand in ("02_eqtl_catalogue_x_gwas_catalog",):
            for ext in ("yaml", "yml", "json"):
                p = base / cand / f"config.{ext}"
                if p.is_file():
                    return p
        files = _list_demos()
        if not files:
            raise FileNotFoundError(f"no bundled demo configs found in {base}")
        return files[0]
    # Named demo: try base/<name>/config.{yaml,yml,json}
    for ext in ("yaml", "yml", "json"):
        p = base / name / f"config.{ext}"
        if p.is_file():
            return p
    # Or try short prefix match (e.g. --demo 02 picks 02_eqtl_catalogue_x_gwas_catalog)
    for d in _examples_dir().iterdir():
        if d.is_dir() and d.name.startswith(f"{name}_"):
            for ext in ("yaml", "yml", "json"):
                p = d / f"config.{ext}"
                if p.is_file():
                    return p
    available = ", ".join(p.parent.name for p in _list_demos())
    raise FileNotFoundError(
        f"no bundled demo named {name!r} under {base}. Available: {available}"
    )


def _print_available_demos() -> None:
    paths = _list_demos()
    if not paths:
        print(f"no bundled demos in {_examples_dir()}")
        return
    print(f"Bundled demos in {_examples_dir()}:")
    for p in paths:
        is_default = p.parent.name in ("02_eqtl_catalogue_x_gwas_catalog",)
        marker = " (default)" if is_default else ""
        print(f"  {p.parent.name}{marker}    [{p.name}]")


def _run(*, config: dict, config_dir: Path, output: Path) -> int:
    """Translate a parsed config into a LocusCompareSpec + clients and render.

    Returns: 0 on success, non-zero on user-visible error.
    """
    # Lazy imports (kept lazy so `--help` doesn't pull in heavy deps).
    from locuscompare_region_render import (
        EXPOSURE_KIND_EQTL_CATALOGUE,
        EXPOSURE_KIND_UKB_PPP,
        LocusCompareSpec,
        Tier2NotAvailable,
        render_locuscompare_for_lead,
    )
    from eqtl_catalogue_region_fetch import EQTLCatalogueClient
    from gwas_catalog_region_fetch import GWASCatalogClient
    from ld_1000g_region_compute import (
        SuperPop,
    )
    from ondemand_client import (
        DEFAULT_PLINK_BIN,
        OnDemand1000GLDClient,
        OnDemandLDError,
    )
    from ukb_ppp_region_fetch import UKBPPPClient

    # ----- Required blocks
    if "lead" not in config:
        print("config missing required 'lead' block", file=sys.stderr)
        return 2
    lead = config["lead"]
    lead_variant_id = lead["variant_id"]
    chromosome = str(lead["chromosome"])
    lead_position_bp = int(lead["position_bp"])
    window_bp = int(lead["window_bp"])
    lead_rs_id = lead.get("rs_id") or None

    exposure = config.get("exposure") or {}
    outcome = config.get("outcome") or {}
    if not exposure or not outcome:
        print("config requires both 'exposure' and 'outcome' blocks", file=sys.stderr)
        return 2

    exposure_label = exposure.get("trait_label", "")
    outcome_label = outcome.get("trait_label", "")

    # ----- Resolve exposure + outcome inputs
    # Two entry vectors per side, mix-and-match supported (per INPUT_SCHEMA.md):
    #   1. `fetch:` block -> live tabix call via the bundled eQTL Cat / UKB-PPP
    #      / GWAS Cat fetchers.
    #   2. `sumstats_path:` -> pre-fetched TSV in the canonical INPUT_SCHEMA
    #      format; resolved via a duck-typed stub client so the core
    #      orchestrator's fetcher dispatch is unchanged.
    if "sumstats_path" not in exposure and "fetch" not in exposure:
        print(
            "exposure requires either a `fetch:` block (live tabix) or a "
            "`sumstats_path:` field (pre-fetched canonical TSV; see "
            "INPUT_SCHEMA.md).",
            file=sys.stderr,
        )
        return 2
    if "sumstats_path" not in outcome and "fetch" not in outcome:
        print(
            "outcome requires either a `fetch:` block (live tabix) or a "
            "`sumstats_path:` field (pre-fetched canonical TSV; see "
            "INPUT_SCHEMA.md).",
            file=sys.stderr,
        )
        return 2

    exposure_prefetched = "sumstats_path" in exposure
    outcome_prefetched = "sumstats_path" in outcome

    # Defaults; overridden below per branch.
    exposure_source = "prefetched"
    eqtl_dataset_id = ""
    molecular_trait_id: str | None = None
    pqtl_protein_label = ""
    pqtl_ancestry = "EUR"
    gwas_accession = "prefetched"

    if exposure_prefetched:
        # Force the eqtl-catalogue dispatch branch -- shape-compatible with
        # the loaded variants; the renderer is fetcher-agnostic downstream.
        exposure_source = "eqtl_catalogue"
        eqtl_dataset_id = exposure.get("study_id") or "prefetched"
        molecular_trait_id = exposure.get("molecular_trait_id")
    else:
        exposure_source = exposure["fetch"]["source"]
        if exposure_source not in ("eqtl_catalogue", "ukb_ppp"):
            print(
                f"unsupported exposure.fetch.source: {exposure_source} "
                f"(expected one of: eqtl_catalogue, ukb_ppp)",
                file=sys.stderr,
            )
            return 2
        eqtl_dataset_id = exposure["fetch"].get("dataset_id", "")
        molecular_trait_id = exposure["fetch"].get("molecular_trait_id")
        pqtl_protein_label = exposure["fetch"].get("protein_label", "")
        pqtl_ancestry = exposure["fetch"].get("ancestry", "EUR")
        if exposure_source == "ukb_ppp" and not pqtl_protein_label:
            print(
                "exposure.fetch.protein_label required when source=ukb_ppp",
                file=sys.stderr,
            )
            return 2
        if exposure_source == "eqtl_catalogue" and not eqtl_dataset_id:
            print(
                "exposure.fetch.dataset_id required when source=eqtl_catalogue",
                file=sys.stderr,
            )
            return 2

    if outcome_prefetched:
        gwas_accession = outcome.get("study_id") or "prefetched"
    else:
        if outcome["fetch"]["source"] != "gwas_catalog":
            print(
                f"unsupported outcome.fetch.source: {outcome['fetch']['source']}",
                file=sys.stderr,
            )
            return 2
        gwas_accession = outcome["fetch"]["accession"]

    # ----- LD client (optional; gracefully degrades if unset / unavailable)
    # Three modes, in priority order:
    # 1. `ld.source: synthetic` (+ `ld_matrix_path`): pre-loaded r2 dict from
    #    a two-column TSV. Used by the offline synthetic demo and any other
    #    case where the caller already has an LD matrix in hand.
    # 2. On-demand region fetch from EBI 1000G FTP via plink 1.9 (default;
    #    matches ClawBio "local-first install, no multi-GB pre-download" UX).
    # 3. None: render without LD coloring (graceful fallback).
    ld_client = None
    ld_block = config.get("ld") or {}
    super_pop_str = ld_block.get("super_pop", "EUR")
    super_pop = SuperPop[super_pop_str]
    ld_source = ld_block.get("source", "1000g_phase3_grch38")
    if ld_source == "synthetic":
        from _prefetched import PrefetchedLDClient, load_synthetic_ld
        ld_matrix_path = ld_block.get("ld_matrix_path")
        if not ld_matrix_path:
            print(
                "ld.source=synthetic requires ld.ld_matrix_path pointing at "
                "a 2-column TSV (partner_variant_id, r2).",
                file=sys.stderr,
            )
            return 2
        ld_matrix_full = _resolve_path(ld_matrix_path, config_dir)
        r2_by_partner = load_synthetic_ld(ld_matrix_full)
        ld_client = PrefetchedLDClient(r2_by_partner=r2_by_partner)
        print(
            f"info: using synthetic LD matrix ({len(r2_by_partner)} partners; "
            f"no plink call): {ld_matrix_full}",
            file=sys.stderr,
        )
    elif ld_source == "1000g_phase3_grch38":
        try:
            plink_bin = ld_block.get("plink_bin") or DEFAULT_PLINK_BIN
            ld_client = OnDemand1000GLDClient(
                super_pop=super_pop_str,
                plink_bin=plink_bin,
            )
            print(
                f"info: using on-demand 1000G LD client (super_pop={super_pop_str}); "
                "first run will fetch a region VCF (~5-50 MB) and cache it.",
                file=sys.stderr,
            )
        except (OnDemandLDError, ImportError) as e:
            print(
                f"warning: on-demand LD client unavailable ({e!s}). "
                "Rendering without LD coloring; variants will appear grey. "
                "Install plink 1.9 (brew / apt / conda) and pysam, or set "
                "PLINK_BIN to a plink binary path.",
                file=sys.stderr,
            )
            ld_client = None

    # ----- Gene track. Four sources in priority order:
    # 1. `gene_track.source: synthetic` (+ `genes_path`): load from a TSV
    #    (gene_symbol, start, end, strand, [biotype]). Used by the offline
    #    synthetic demo and any case where the caller already has a gene
    #    track in hand. No Ensembl REST call.
    # 2. Caller-supplied local GTF path (legacy; fastest).
    # 3. On-demand fetch from Ensembl REST (default).
    # 4. None: render without gene track + caveat (graceful fallback).
    gene_track_block = config.get("gene_track") or {}
    gencode_gtf_path = gene_track_block.get("gtf_path") or os.environ.get("GENCODE_GTF")
    if gencode_gtf_path:
        gencode_gtf_path = Path(gencode_gtf_path)
    gene_biotypes = tuple(gene_track_block.get("biotypes") or ("protein_coding",))
    prefetched_gene_track = None
    gt_source = gene_track_block.get("source", "gencode_v39")
    if gt_source == "synthetic":
        from _prefetched import load_synthetic_gene_track
        genes_path = gene_track_block.get("genes_path")
        if not genes_path:
            print(
                "gene_track.source=synthetic requires gene_track.genes_path "
                "pointing at a TSV (gene_symbol, start, end, strand, [biotype]).",
                file=sys.stderr,
            )
            return 2
        genes_full = _resolve_path(genes_path, config_dir)
        prefetched_gene_track = load_synthetic_gene_track(genes_full)
        print(
            f"info: gene track from synthetic TSV ({len(prefetched_gene_track)} "
            f"genes): {genes_full}",
            file=sys.stderr,
        )
    elif gencode_gtf_path is None and gt_source in ("gencode_v39", "ensembl_rest"):
        try:
            from _fetchers.gencode_ondemand import (
                fetch_region_genes_remote,
            )
            from locuscompare_region_render import GeneTrackEntry
            half = max(window_bp // 2, 1)
            start_bp = max(0, lead_position_bp - half)
            end_bp = lead_position_bp + half
            genes, release_meta, gt_notes = fetch_region_genes_remote(
                chromosome=chromosome,
                start_bp=start_bp,
                end_bp=end_bp,
                biotypes=gene_biotypes,
            )
            prefetched_gene_track = [
                GeneTrackEntry(
                    gene_symbol=g.gene_symbol,
                    start=g.start, end=g.end, strand=g.strand,
                    exons=[(e.start, e.end) for e in g.exons],
                    biotype=g.biotype,
                )
                for g in genes
            ]
            print(
                f"info: gene track from Ensembl REST ({len(prefetched_gene_track)} genes, "
                f"{release_meta['release_label']})",
                file=sys.stderr,
            )
        except Exception as e:
            print(
                f"warning: on-demand gene-track fetch failed ({e!s}); "
                "rendering without a gene track.",
                file=sys.stderr,
            )
            prefetched_gene_track = None

    # ----- Build the spec (exposure-kind dispatch)
    if exposure_source == "ukb_ppp":
        exposure_kind = EXPOSURE_KIND_UKB_PPP
    else:
        exposure_kind = EXPOSURE_KIND_EQTL_CATALOGUE

    spec = LocusCompareSpec(
        lead_variant_id=lead_variant_id,
        chromosome=chromosome,
        lead_position_bp=lead_position_bp,
        window_bp=window_bp,
        eqtl_dataset_id=eqtl_dataset_id,
        molecular_trait_id=molecular_trait_id,
        gwas_accession=gwas_accession,
        exposure_kind=exposure_kind,
        pqtl_protein_label=pqtl_protein_label,
        pqtl_ancestry=pqtl_ancestry,
        exposure_gene_symbol=exposure.get("gene_symbol", "") or pqtl_protein_label,
        outcome_trait_label=outcome_label,
        exposure_id_extra="",
        outcome_id_extra="",
        provenance_prefix=_format_provenance_prefix(config),
        release_tag=str((config.get("provenance") or {}).get("ot_release") or ""),
        notes=list((config.get("caveats") or [])),
        super_pop=super_pop,
        gencode_gtf_path=gencode_gtf_path,
        gene_biotypes=gene_biotypes,
        prefetched_gene_track=prefetched_gene_track,
        lead_rs_id=lead_rs_id,
    )

    # ----- Build clients
    # Exposure side: stub from a pre-fetched TSV when the config block uses
    # `sumstats_path:`, else the live fetcher matching the dispatched kind.
    if exposure_prefetched:
        from _prefetched import PrefetchedEQTLClient, load_sumstats_tsv
        exposure_tsv = _resolve_path(exposure["sumstats_path"], config_dir)
        exposure_variants = load_sumstats_tsv(exposure_tsv)
        eqtl_client = PrefetchedEQTLClient(
            variants=exposure_variants,
            dataset_id=eqtl_dataset_id,
            study_label=exposure.get("trait_label") or None,
        )
        ukb_ppp_client = None
        print(
            f"info: exposure sumstats from pre-fetched TSV "
            f"({len(exposure_variants)} variants): {exposure_tsv}",
            file=sys.stderr,
        )
    else:
        eqtl_client = EQTLCatalogueClient()
        ukb_ppp_client = UKBPPPClient() if exposure_kind == EXPOSURE_KIND_UKB_PPP else None

    # Outcome side: stub from a pre-fetched TSV when the config block uses
    # `sumstats_path:`, else the live GWAS Catalog fetcher.
    if outcome_prefetched:
        from _prefetched import (
            PrefetchedGWASClient,
            gwas_variants_from_eqtl,
            load_sumstats_tsv,
        )
        outcome_tsv = _resolve_path(outcome["sumstats_path"], config_dir)
        outcome_variants = gwas_variants_from_eqtl(load_sumstats_tsv(outcome_tsv))
        gwas_client = PrefetchedGWASClient(
            variants=outcome_variants,
            accession=gwas_accession,
        )
        print(
            f"info: outcome sumstats from pre-fetched TSV "
            f"({len(outcome_variants)} variants): {outcome_tsv}",
            file=sys.stderr,
        )
    else:
        gwas_client = GWASCatalogClient()

    # ----- Render
    plot_path = output / f"{lead_variant_id}_full_locuscompare.png"
    try:
        result = render_locuscompare_for_lead(
            spec,
            eqtl_client=eqtl_client,
            gwas_client=gwas_client,
            ld_client=ld_client,
            out_path=plot_path,
            ukb_ppp_client=ukb_ppp_client,
        )
    except Tier2NotAvailable as e:
        print(f"render failed: {e!s}", file=sys.stderr)
        return 1

    # ----- Manifest + report
    manifest = {
        "skill": "locuscompare",
        "version": "0.1.0",
        "lead_variant_id": lead_variant_id,
        "lead_rs_id": lead_rs_id,
        "n_pairs": result.n_pairs,
        "n_palindromic_excluded": result.n_palindromic_excluded,
        "plot_path": str(plot_path.relative_to(output)),
        "notes": result.notes,
        "render_block": result.manifest_block,
    }
    (output / "manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False))

    lead_line = (
        f"- **Lead variant:** `{lead_variant_id}`"
        + (f" ({lead_rs_id}; " if lead_rs_id else " (")
        + f"chr{chromosome}:{lead_position_bp}, ±{window_bp//2000} kb)"
    )
    report_lines = [
        "# locuscompare report",
        "",
        lead_line,
        f"- **Exposure:** {exposure_label}",
        f"- **Outcome:** {outcome_label}",
        f"- **n_pairs:** {result.n_pairs}",
        f"- **n_palindromic_excluded:** {result.n_palindromic_excluded}",
        f"- **Plot:** {plot_path.name}",
    ]
    if result.notes:
        report_lines.append("")
        report_lines.append("## Notes")
        for n in result.notes:
            report_lines.append(f"- {n}")
    (output / "report.md").write_text("\n".join(report_lines) + "\n")

    print(f"locuscompare: rendered {plot_path}")
    print(f"  n_pairs = {result.n_pairs}, n_palindromic_excluded = {result.n_palindromic_excluded}")
    return 0


def _format_provenance_prefix(config: dict) -> str:
    """Build a provenance label prefix from optional config blocks."""
    prov = config.get("provenance") or {}
    bits: list[str] = []
    if prov.get("ot_release"):
        bits.append(f"OT release: {prov['ot_release']}")
    if prov.get("gwas_lookup_run_dir"):
        bits.append(f"gwas-lookup chain: {prov['gwas_lookup_run_dir']}")
    return " | ".join(bits) + (" | " if bits else "")


def _load_config(path: Path) -> dict:
    """Load JSON or YAML by extension. Returns the parsed dict."""
    text = path.read_text()
    if path.suffix.lower() in (".yaml", ".yml"):
        return yaml.safe_load(text) or {}
    if path.suffix.lower() == ".json":
        return json.loads(text)
    raise ValueError(f"unsupported config extension: {path.suffix}")


def _resolve_path(p: str | Path, config_dir: Path) -> Path:
    """Resolve a possibly-relative path against the config file's directory."""
    pp = Path(p)
    if pp.is_absolute():
        return pp
    return (config_dir / pp).resolve()


if __name__ == "__main__":
    sys.exit(main())
