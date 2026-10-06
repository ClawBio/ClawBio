"""Write report.md, result.json and the reproducibility bundle.

Two rules govern the prose:

* Every number is read from a table the pipeline wrote or from the effective
  config. Nothing is typed into a template, so changing
  `mag.checkm2.min_completeness` changes the report.
* An empty upstream cell renders as `not measured`. Upstream writes one where a
  measurement could not be made, and a report that printed 0 there would be
  asserting something the data does not support.

The command for a `--check` run arrives through `pipeline_source["command"]`,
because the documented `write_report_md` signature has no other channel for it
and a planned run must still show what it would have executed.
"""

from __future__ import annotations

import json
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

SKILL_DIR = Path(__file__).resolve().parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

# The repo root holds the `clawbio` package; the documented direct invocation
# (`python skills/metagenome-mag-pipeline/metagenome_mag_pipeline.py`) puts only
# this directory on sys.path.
_PROJECT_ROOT = SKILL_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from magpipe_outputs import thresholds_from_config
from magpipe_preflight import network_requirements
from magpipe_schemas import (
    DISCLAIMER,
    SKILL_NAME,
    SKILL_VERSION,
    SNAKEMAKE_PIN,
)

from clawbio.common.reproducibility import (
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

NOT_MEASURED = "not measured"

# Metrics worth showing from QUAST, in report order. Anything the run produced is
# in the JSON; the table stays readable.
_QUAST_METRICS = (
    "# contigs",
    "Largest contig",
    "Total length",
    "GC (%)",
    "# contigs (>= 1000 bp)",
    "N50 (scaffolds)",
)


def _cell(value: Any) -> str:
    if value is None:
        return NOT_MEASURED
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:g}"
    return str(value)


def _table(columns: list[str], rows: list[list[Any]]) -> str:
    if not rows:
        return "_No rows._\n"
    header = "| " + " | ".join(columns) + " |"
    divider = "|" + "|".join(["---"] * len(columns)) + "|"
    body = ["| " + " | ".join(_cell(value) for value in row) + " |" for row in rows]
    return "\n".join([header, divider, *body]) + "\n"


def stage_plan(config: dict) -> list[dict[str, Any]]:
    """Which stages this config enables, and why any are off."""
    assembly = bool(config.get("assembly", {}).get("enabled", True))
    binning = bool(config.get("binning", {}).get("enabled", True))
    profiler = config.get("profile", {}).get("profiler", "none")
    mag = bool(config.get("mag", {}).get("enabled", True))
    taxonomy = bool(config.get("taxonomy", {}).get("enabled", True))
    annotation = bool(config.get("annotation", {}).get("enabled", True))
    heterogeneity = bool(config.get("heterogeneity", {}).get("enabled", True))

    def off(reason: str) -> str:
        return reason

    return [
        {"stage": "read QC", "enabled": True, "reason": ""},
        {
            "stage": "taxonomic profiling",
            "enabled": profiler != "none",
            "reason": "" if profiler != "none" else off("profile.profiler is none"),
        },
        {
            "stage": "assembly (metaSPAdes, QUAST)",
            "enabled": assembly,
            "reason": "" if assembly else off("assembly.enabled is false"),
        },
        {
            "stage": "binning (MetaBAT2/MaxBin2/VAMB, DAS Tool)",
            "enabled": binning,
            "reason": "" if binning else off("binning.enabled is false"),
        },
        {
            "stage": "MAG assessment (CheckM2, dRep, GTDB-Tk)",
            "enabled": mag,
            "reason": ""
            if mag
            else off(
                "needs --checkm2-db and --gtdbtk-db"
                if binning
                else "binning is off, and MAG assessment requires binning"
            ),
        },
        {
            "stage": "contig taxonomy (DIAMOND, LCA)",
            "enabled": taxonomy and assembly,
            "reason": ""
            if taxonomy and assembly
            else off(
                "needs --taxonomy-db and --taxdump"
                if assembly
                else "requires assembly, which is off"
            ),
        },
        {
            "stage": "functional annotation (DIAMOND, hmmsearch)",
            "enabled": annotation and assembly,
            "reason": ""
            if annotation and assembly
            else off(
                "needs --kegg-db, --cog-db and --pfam-db"
                if assembly
                else "requires assembly, which is off"
            ),
        },
        {
            "stage": "strain heterogeneity",
            "enabled": heterogeneity and binning,
            "reason": ""
            if heterogeneity and binning
            else off(
                "heterogeneity.enabled is false"
                if not heterogeneity
                else "requires binning, which is off"
            ),
        },
        {
            "stage": "read-sketch similarity (sourmash)",
            "enabled": bool(config.get("qc", {}).get("sourmash", {}).get("enabled", False)),
            "reason": ""
            if config.get("qc", {}).get("sourmash", {}).get("enabled", False)
            else off("qc.sourmash.enabled is false"),
        },
        {
            "stage": "consolidated QC report (MultiQC)",
            "enabled": bool(config.get("report", {}).get("multiqc", True)),
            "reason": ""
            if config.get("report", {}).get("multiqc", True)
            else off("report.multiqc is false"),
        },
        {"stage": "provenance and timings", "enabled": True, "reason": ""},
    ]


def _stage_lines(plan: list[dict[str, Any]]) -> str:
    ran = [entry["stage"] for entry in plan if entry["enabled"]]
    skipped = [entry for entry in plan if not entry["enabled"]]
    lines = ["**Stages run**: " + (", ".join(ran) if ran else "none") + "\n"]
    lines.append("**Stages skipped**: " + (", ".join(entry["stage"] for entry in skipped) if skipped else "none"))
    if skipped:
        lines.append("")
        lines += [f"- {entry['stage']} — {entry['reason']}" for entry in skipped]
    return "\n".join(lines) + "\n"


def _section_qc(findings: dict) -> str:
    rows = [
        [
            row.get("sample"),
            row.get("layout"),
            row.get("retained_reads"),
            row.get("retained_bases"),
            row.get("retained_pairs"),
        ]
        for row in findings["qc"]["read_counts"]
    ]
    text = "## Read QC\n\nRetained reads and bases, counted from the post-QC FASTQ "\
        "files. `retained_reads` counts both mates; `retained_pairs` counts fragments.\n\n"
    text += _table(["sample", "layout", "retained_reads", "retained_bases", "retained_pairs"], rows)
    sourmash = findings["qc"].get("sourmash")
    if sourmash and sourmash.get("samples"):
        text += "\n### Read-sketch similarity (sourmash)\n\n"
        text += _table(
            ["sample"] + sourmash["samples"],
            [
                [sample] + [sourmash["values"][sample].get(other) for other in sourmash["samples"]]
                for sample in sourmash["samples"]
                if sample in sourmash["values"]
            ],
        )
    return text


def _section_profile(findings: dict) -> str:
    profile = findings["profile"]
    text = (
        "## Taxonomic profiling (MetaPhlAn)\n\n"
        f"Profile database: `{profile['version']}`. Every value is a **percentage "
        "relative abundance**, not a read count.\n\n"
    )
    top = profile.get("top_taxa", {})
    rows = []
    for sample in profile["samples"]:
        for entry in top.get(sample, [])[:10]:
            rows.append([sample, entry["clade_name"], entry["relative_abundance_pct"]])
    text += _table(["sample", "clade", "relative_abundance (%)"], rows)
    text += (
        "\nA clade at the top of one sample's profile and absent from another's is "
        "a compositional difference. It is not evidence of a different organism "
        "unless the two profiles are compared on the same basis.\n"
    )
    return text


def _section_assembly(findings: dict) -> str:
    assembly = findings["assembly"]
    metrics = assembly["metrics"]
    names = [name for name in _QUAST_METRICS if name in metrics]
    extra = sorted(set(metrics) - set(names))
    rows = [
        [metric] + [assembly["quast"].get(unit, {}).get(metric) for unit in assembly["units"]]
        for metric in [*names, *extra]
    ]
    text = "## Assembly (QUAST)\n\n"
    text += _table(["metric"] + assembly["units"], rows)
    text += (
        "\nContig counts describe the assembly after `assembly.min_contig_length` "
        "filtering, not a genome size.\n"
    )
    return text


def _section_bins(findings: dict, config: dict) -> str:
    bins = findings["bins"]
    # Read from the config this report was called with, not from the parse: the
    # report must echo the thresholds the run used even when it is re-rendered.
    thresholds = thresholds_from_config(config)
    text = (
        "## Bins and MAG assessment\n\n"
        f"{bins['count']} bin(s) recovered; {bins['passing_qa']} pass the configured "
        "quality gate (CheckM2 completeness >= "
        f"{_cell(thresholds.get('checkm2_min_completeness'))}, contamination <= "
        f"{_cell(thresholds.get('checkm2_max_contamination'))}). dRep "
        "dereplication uses completeness >= "
        f"{_cell(thresholds.get('drep_completeness'))} and contamination <= "
        f"{_cell(thresholds.get('drep_contamination'))}.\n\n"
    )
    text += "### Genome quality (CheckM2)\n\n"
    text += _table(
        ["genome", "completeness (%)", "contamination (%)", "passes_qa"],
        [
            [row.get("genome"), row.get("completeness"), row.get("contamination"), row.get("passes_qa")]
            for row in bins["checkm2"]
        ],
    )
    text += "\n### Taxonomy (GTDB-Tk)\n\n"
    text += _table(
        ["genome", "classification"],
        [[row.get("genome"), row.get("classification")] for row in bins["gtdb"]],
    )
    classified = sum(1 for row in bins["gtdb"] if row.get("classification"))
    text += (
        f"\n{classified} of {len(bins['gtdb'])} bin(s) received a classification. "
        "Passing the completeness and contamination gate is a filter, not a "
        "finding: it does not by itself establish a high-quality MAG, and no "
        "coverage, marker-gene or contamination estimate beyond CheckM2's is "
        "reported here.\n"
    )
    if bins["drep"]:
        text += "\n### Dereplicated representatives (dRep)\n\n"
        text += _table(
            [key for key in bins["drep"][0] if key],
            [[row.get(key) for key in bins["drep"][0] if key] for row in bins["drep"]],
        )
    return text


def _section_diagnostic(findings: dict) -> str:
    rows = [
        [
            row.get("sample"),
            row.get("layout"),
            row.get("reads"),
            row.get("assembly_capture_pct"),
            row.get("catalogue_capture_pct"),
        ]
        for row in findings["diagnostic"]["bottleneck"]
    ]
    text = (
        "## Where the reads were lost\n\n"
        "Each column is the percentage of post-QC reads placed exactly once on the "
        "sample's own assembly, and on the concatenated bin catalogue, over the "
        "same denominator.\n\n"
    )
    text += _table(
        ["sample", "layout", "reads", "assembly_capture_pct (%)", "catalogue_capture_pct (%)"],
        rows,
    )
    text += (
        f"\n`{NOT_MEASURED}` means the pipeline could not compute that axis — for "
        "example a single-end sample has no assembly to place reads against. It is "
        "not zero: reporting 0% would place the sample on the axis and imply its "
        "reads were definitively not in any bin. A high assembly figure with a low "
        "catalogue figure means reads were assembled and then lost at binning; a "
        "low assembly figure means the loss happened at or before assembly.\n"
    )
    return text


def _section_taxonomy(findings: dict) -> str:
    rows = [
        [
            row.get("contig"),
            row.get("taxonomy"),
            row.get("assigned_genes"),
            row.get("total_genes"),
            row.get("disparity"),
        ]
        for row in findings["taxonomy"]["contigs"]
    ]
    text = (
        "## Contig taxonomy\n\n"
        "Each contig gets the most specific label a quorum of its genes agrees on; "
        "a contig whose genes disagree is left `unclassified` rather than given a "
        "majority label. High disparity suggests a chimera or a mixed bin.\n\n"
    )
    return text + _table(
        ["contig", "taxonomy", "assigned_genes", "total_genes", "disparity"], rows
    )


def _section_annotation(findings: dict) -> str:
    rows = [
        [row.get("bin"), row.get("source"), row.get("function"), row.get("genes")]
        for row in findings["annotation"]["functions"]
    ]
    text = (
        "## Functional annotation\n\n"
        "Gene counts per function, taken from the assembly's gene calls. This is a "
        "presence count, not a per-sample abundance: `unassigned` collects functions "
        "on contigs no binner claimed.\n\n"
    )
    return text + _table(["bin", "source", "function", "genes"], rows)


def _section_heterogeneity(findings: dict) -> str:
    rows = [
        [
            row.get("sample"),
            row.get("bin"),
            row.get("genes"),
            row.get("covered_positions"),
            row.get("heterogeneity_pct"),
        ]
        for row in findings["heterogeneity"]["rows"]
    ]
    text = (
        "## Strain heterogeneity\n\n"
        "The percentage of covered coding positions whose codon would change the "
        "protein between the two commonest alleles. A low value is evidence the bin "
        "is one strain; it is not a completeness estimate.\n\n"
    )
    text += _table(
        ["sample", "bin", "genes", "covered_positions", "heterogeneity (%)"], rows
    )
    text += (
        f"\n`{NOT_MEASURED}` means no coding position passed the configured depth "
        "filter, which is a different claim from a measured zero.\n"
    )
    return text


def _section_provenance(findings: dict) -> str:
    provenance = findings["provenance"]
    pipeline = provenance["pipeline"]
    text = "## Provenance\n\n"
    text += _table(
        ["field", "value"],
        [[key, value] for key, value in pipeline.items()],
    )
    text += "\n### Tool versions\n\n"
    text += _table(
        [key for key in ("tool", "version", "evidence") if key in (provenance["tools"][0] if provenance["tools"] else {})],
        [
            [row.get(key) for key in ("tool", "version", "evidence")]
            for row in provenance["tools"]
        ],
    )
    if provenance["references"]:
        text += "\n### Reference checksums\n\n"
        text += _table(
            ["resource", "md5"],
            [[row.get("resource"), row.get("md5")] for row in provenance["references"]],
        )
    manifest = provenance["manifest"]
    if manifest:
        text += (
            f"\nRun manifest: samplesheet SHA-256 `{manifest.get('samplesheet_sha256')}`, "
            f"{manifest.get('source_file_count')} pipeline source file(s) hashed, "
            f"{manifest.get('conda_environment_count')} conda environment(s) recorded.\n"
        )
    if provenance["runtime"]:
        text += "\n### Slowest rules\n\n"
        rows = sorted(
            provenance["runtime"],
            key=lambda row: row.get("s") if isinstance(row.get("s"), (int, float)) else 0,
            reverse=True,
        )[:10]
        text += _table(
            ["rule", "seconds", "max_rss"],
            [[row.get("rule"), row.get("s"), row.get("max_rss")] for row in rows],
        )
    return text


def _caveats(mode: str, findings: dict | None, config: dict) -> str:
    items: list[str] = []
    if mode == "demo":
        items.append(
            "**Synthetic data.** The demo's reads are random fragments generated at "
            "run time. They have no biological truth: no taxonomy and no marker "
            "genes, so the numbers above prove that the pipeline executes, not that "
            "it is accurate. "
            "Consensus binning cannot complete on such input at all, because "
            "MaxBin2 and DAS Tool need single-copy marker genes."
        )
    if findings and findings.get("profile"):
        items.append(
            "**Relative abundance, not counts.** MetaPhlAn values are percentages of "
            "the sample, so they cannot be compared across samples as if they were "
            "read counts."
        )
    if findings and findings.get("bins"):
        items.append(
            "**Threshold passing is not a quality verdict.** A bin that clears the "
            "configured completeness and contamination gate has passed a filter. "
            "That alone does not establish a high-quality MAG."
        )
    if findings:
        has_blank = any(
            row.get("heterogeneity_pct") is None
            for row in (findings.get("heterogeneity") or {}).get("rows", [])
        ) or any(
            row.get("assembly_capture_pct") is None
            for row in findings.get("diagnostic", {}).get("bottleneck", [])
        )
        if has_blank:
            items.append(
                f"**Blank is not zero.** `{NOT_MEASURED}` marks an axis the pipeline "
                "could not measure. Never read it as 0% and never plot it as 0."
            )
    if findings and findings.get("taxonomy"):
        items.append(
            "**Consensus labels hide disagreement.** A contig whose genes disagree is "
            "reported `unclassified` rather than given a majority label; treat a high "
            "disparity value as a signal of a chimera or a mixed bin."
        )
    if not items:
        return ""
    return "\n## Caveats\n\n" + "\n".join(f"- {item}" for item in items) + "\n"


def write_report_md(
    output_dir: Path,
    *,
    status: str,
    mode: str,
    findings: dict | None,
    config: dict,
    pipeline_source: dict,
    warnings: list[str],
) -> Path:
    """Write report.md. `pipeline_source` may carry a `command` key for a plan."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    plan = stage_plan(config)
    commit = pipeline_source.get("commit", "unknown")
    ref = pipeline_source.get("ref", "")

    parts = ["# Shotgun metagenomics: MAG pipeline report\n"]
    label = "demo (synthetic reads)" if mode == "demo" else mode
    parts.append(f"**Mode**: {label} · **Status**: {status}\n")
    if status == "error":
        parts.append(
            "**The run failed.** No findings below come from a completed pipeline run; "
            "the cause is under Warnings.\n"
        )
    if ref and ref != commit:
        parts.append(f"**Pipeline**: metagenomics-workflow @ `{ref}` ({commit})\n")
    else:
        parts.append(f"**Pipeline**: metagenomics-workflow @ `{commit}`\n")
    kinds = {"local_checkout": "local checkout", "cached_clone": "cached clone"}
    kind = pipeline_source.get("source_kind", "unknown")
    parts.append(f"**Source**: {kinds.get(kind, kind)}\n")
    if base_config := pipeline_source.get("base_config"):
        parts.append(f"**Effective config derived from**: `{base_config}`\n")
    if mode == "check":
        parts.append(
            "**Nothing was executed.** `--check` validated the inputs and rendered a "
            "dry run only, so no results below come from a pipeline run.\n"
        )
    parts.append(_stage_lines(plan))

    if mode == "check":
        command = pipeline_source.get("command", "")
        parts.append("## Planned stages\n\n" + _table(
            ["stage", "enabled", "reason"],
            [[entry["stage"], "yes" if entry["enabled"] else "no", entry["reason"]] for entry in plan],
        ))
        parts.append("\n## Command that would run\n\n```bash\n" + (command or "(not recorded)") + "\n```\n")
        parts.append("\n## Network use\n\n")
        parts += [f"- {line}" for line in network_requirements(config, [])] + ["\n"]
    elif findings is not None:
        parts.append(_section_qc(findings))
        if findings.get("profile"):
            parts.append(_section_profile(findings))
        if findings.get("assembly"):
            parts.append(_section_assembly(findings))
        if findings.get("bins"):
            parts.append(_section_bins(findings, config))
        parts.append(_section_diagnostic(findings))
        if findings.get("taxonomy"):
            parts.append(_section_taxonomy(findings))
        if findings.get("annotation"):
            parts.append(_section_annotation(findings))
        if findings.get("heterogeneity"):
            parts.append(_section_heterogeneity(findings))
        parts.append(_section_provenance(findings))

    links = [line for line in (
        "- `results/reports/multiqc/multiqc_report.html` — MultiQC read-QC report"
        if _present(output_dir, "results/reports/multiqc/multiqc_report.html")
        else "",
        "- `results/reports/assembly/quast.html` — QUAST assembly report"
        if _present(output_dir, "results/reports/assembly/quast.html")
        else "",
    ) if line]
    if links:
        parts.append("\n## Pipeline reports\n\n" + "\n".join(links) + "\n")

    if warnings:
        parts.append("\n## Warnings\n\n" + "\n".join(f"- {line}" for line in warnings) + "\n")

    parts.append(_caveats(mode, findings, config))

    parts.append(
        "\n## Reproducibility\n\n"
        "The exact command, the environment and the SHA-256 of every reported file "
        "are in `reproducibility/` (`commands.sh`, `environment.yml`, "
        "`checksums.sha256`). The effective config is `pipeline/config.yaml` and "
        "the inputs are `pipeline/samplesheet.tsv`.\n"
    )
    parts.append(f"\n*{DISCLAIMER}*\n")

    path = output_dir / "report.md"
    path.write_text("\n".join(parts), encoding="utf-8")
    return path


def _present(output_dir: Path, relative: str) -> bool:
    return (output_dir / relative).exists()


def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def write_result_json(
    output_dir: Path,
    *,
    status: str,
    mode: str,
    findings: dict | None,
    config: dict,
    pipeline_source: dict,
    command: str,
    warnings: list[str],
    config_source: dict | None = None,
) -> Path:
    """Write result.json: what ran, what it found, and how to reproduce it.

    `config_source` records which upstream configuration the effective config was
    derived from, so a demo run and a real run are distinguishable from the file
    alone. It is optional because a caller that failed before rendering a config
    has none to record.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "status": status,
        "mode": mode,
        "skill": SKILL_NAME,
        "version": SKILL_VERSION,
        "generated_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "pipeline": _jsonable(pipeline_source),
        "command": command,
        "stages": stage_plan(config),
        "config": _jsonable(config),
        "config_source": _jsonable(config_source) if config_source else None,
        "findings": _jsonable(findings),
        "warnings": list(warnings),
        "disclaimer": DISCLAIMER,
        "tool_versions": {"snakemake_pin": SNAKEMAKE_PIN},
    }
    path = output_dir / "result.json"
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path


def write_reproducibility(output_dir: Path, *, command: str, files: list[Path]) -> None:
    """Write the ClawBio reproducibility bundle: commands, environment, checksums."""
    output_dir = Path(output_dir)
    write_commands_sh(output_dir, command)
    write_environment_yml(
        output_dir,
        "clawbio-metagenome-mag-pipeline",
        pip_deps=[],
        conda_deps=[f"snakemake={SNAKEMAKE_PIN}", "pandas", "pyyaml"],
        channels=["conda-forge", "bioconda"],
    )
    write_checksums([Path(item) for item in files], output_dir, anchor=output_dir)
