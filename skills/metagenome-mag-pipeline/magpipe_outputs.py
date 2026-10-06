"""Read the pipeline's final tables, and decide which of them must exist.

`expected_outputs` mirrors upstream's `analysis_targets()` in
`workflow/rules/helpers.smk` at the pinned commit: the always-present tables plus
the ones the enabled stages produce. A path may contain a `*` where upstream
expands over samples; such a path counts as present when at least one file
matches.

`parse_outputs` keeps upstream's column names verbatim, converts numeric strings
to numbers, and maps an empty cell to None. Upstream writes an empty cell where a
measurement could not be made, and reports that as "not measured"; turning it
into 0 would assert something the data does not support.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path
from typing import Any

SKILL_DIR = Path(__file__).resolve().parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

from magpipe_errors import ErrorCode, SkillError
from magpipe_schemas import ALWAYS_EXPECTED_OUTPUTS, FINAL_OUTPUTS

STAGE = "outputs"

# A plain integer or decimal, without a leading zero that would turn an
# identifier such as "007" into a number.
_NUMERIC = re.compile(r"-?(?:0|[1-9]\d*)(?:\.\d+)?\Z")
_TOP_TAXA = 10


def _fail(message: str, fix: str, **details: object) -> SkillError:
    return SkillError(
        stage=STAGE,
        error_code=ErrorCode.EXPECTED_OUTPUTS_NOT_FOUND,
        message=message,
        fix=fix,
        details=dict(details),
    )


def expected_outputs(config: dict) -> list[str]:
    """Relative paths the run must have produced, given the effective config."""
    assembly = bool(config.get("assembly", {}).get("enabled", True))
    binning = bool(config.get("binning", {}).get("enabled", True))
    paths = [FINAL_OUTPUTS["read_counts"]]

    if config.get("profile", {}).get("profiler", "none") != "none":
        paths.append(FINAL_OUTPUTS["profile"])
    if assembly:
        paths.append(FINAL_OUTPUTS["quast"])
    if binning:
        paths.append("intermediate/binning/dastool/*/contig2bin.tsv")
    if config.get("mag", {}).get("enabled", True):
        paths += [
            FINAL_OUTPUTS["bins"],
            FINAL_OUTPUTS["checkm2_quality"],
            FINAL_OUTPUTS["drep_genome_info"],
            FINAL_OUTPUTS["gtdb_summary"],
        ]
    paths.append(FINAL_OUTPUTS["bottleneck"])
    if config.get("taxonomy", {}).get("enabled", True) and assembly:
        paths.append(FINAL_OUTPUTS["contig_taxonomy"])
    if config.get("annotation", {}).get("enabled", True) and assembly:
        paths.append(FINAL_OUTPUTS["function_abundance"])
        if config.get("annotation", {}).get("minpath", {}).get("data_dir"):
            paths.append("intermediate/annotation/minpath/*.minpath")
    if config.get("heterogeneity", {}).get("enabled", True) and binning:
        paths.append(FINAL_OUTPUTS["strain_heterogeneity"])
    if config.get("qc", {}).get("sourmash", {}).get("enabled", False):
        paths.append(FINAL_OUTPUTS["sourmash_similarities"])
    if config.get("report", {}).get("multiqc", True):
        paths.append(FINAL_OUTPUTS["multiqc"])
    paths += list(ALWAYS_EXPECTED_OUTPUTS)
    return paths


def _resolve(results_dir: Path, relative: str) -> list[Path]:
    return sorted(results_dir.glob(relative)) if "*" in relative else [results_dir / relative]


def _check_present(results_dir: Path, config: dict) -> dict[str, Path]:
    missing: list[str] = []
    found: dict[str, Path] = {}
    for relative in expected_outputs(config):
        matches = [path for path in _resolve(results_dir, relative) if path.exists()]
        if not matches:
            missing.append(relative)
        else:
            found[relative] = matches[0]
    if missing:
        raise _fail(
            "the pipeline did not write "
            + str(len(missing))
            + " expected output(s): "
            + ", ".join(missing),
            "Read <output>/pipeline/snakemake.log for the failing rule. Snakemake "
            "keeps its state under <output>/pipeline/workdir, so re-running the "
            "same command resumes rather than starting again.",
            missing=missing,
            results_dir=str(results_dir),
        )
    return found


def _value(cell: str | None) -> Any:
    """Upstream's empty cell is "not measured"; a number is a number."""
    if cell is None:
        return None
    text = cell.strip()
    if not text:
        return None
    if text in ("True", "False"):
        return text == "True"
    if _NUMERIC.fullmatch(text):
        return float(text) if "." in text else int(text)
    return text


def _read_table(path: Path, delimiter: str | None = None) -> list[dict[str, Any]]:
    if delimiter is None:
        delimiter = "," if path.suffix == ".csv" else "\t"
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows = []
        for row in reader:
            rows.append(
                {
                    (key or "").strip(): _value(value)
                    for key, value in row.items()
                    if key is not None
                }
            )
    return rows


def _matrix(path: Path) -> dict[str, Any]:
    """A sourmash comparison matrix: an index column, then one column per sample."""
    rows = _read_table(path)
    if not rows:
        return {"samples": [], "values": {}}
    columns = [key for key in rows[0] if key]
    samples = columns[1:]
    return {
        "samples": samples,
        "values": {
            row[columns[0]]: {name: row.get(name) for name in samples}
            for row in rows
            if row.get(columns[0]) is not None
        },
    }


def _profile(path: Path) -> dict[str, Any]:
    """MetaPhlAn output: a version line, a header line, then percentage rows."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise _fail(f"{path} is too short to be a MetaPhlAn table", "Check the profiler rule's log.")
    version = lines[0].strip()
    header = lines[1].split("\t")
    samples = [name for name in header[1:] if name]
    clades = []
    for line in lines[2:]:
        if not line.strip():
            continue
        cells = line.split("\t")
        clades.append(
            {
                "clade_name": cells[0],
                **{
                    name: _value(cells[index + 1] or None)
                    for index, name in enumerate(samples)
                },
            }
        )
    top = {
        sample: [
            {"clade_name": row["clade_name"], "relative_abundance_pct": row.get(sample)}
            for row in sorted(
                (row for row in clades if isinstance(row.get(sample), (int, float))),
                key=lambda row: row[sample],
                reverse=True,
            )[:_TOP_TAXA]
        ]
        for sample in samples
    }
    return {"version": version, "samples": samples, "clades": clades, "top_taxa": top}


def _quast(path: Path) -> tuple[list[str], dict[str, dict[str, Any]]]:
    """QUAST writes metrics down the first column and one column per assembly."""
    rows = _read_table(path)
    if not rows:
        return [], {}
    header = [key for key in rows[0] if key]
    units = header[1:]
    metrics: dict[str, dict[str, Any]] = {}
    for row in rows:
        metric = row.get(header[0])
        if metric:
            metrics[metric] = {unit: row.get(unit) for unit in units}
    return units, metrics


def thresholds_from_config(config: dict) -> dict[str, Any]:
    mag = config.get("mag", {})
    return {
        "checkm2_min_completeness": mag.get("checkm2", {}).get("min_completeness"),
        "checkm2_max_contamination": mag.get("checkm2", {}).get("max_contamination"),
        "drep_completeness": mag.get("drep", {}).get("completeness"),
        "drep_contamination": mag.get("drep", {}).get("contamination"),
    }


def _pipeline_revision(path: Path) -> dict[str, Any]:
    return {row["field"]: row["value"] for row in _read_table(path) if row.get("field")}


def _references(path: Path) -> list[dict[str, Any]]:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    lines = [line for line in lines if not line.startswith("#")]
    if not lines:
        return []
    rows = list(csv.DictReader(lines, delimiter="\t"))
    return [
        {key.strip(): _value(value) for key, value in row.items() if key}
        for row in rows
    ]


def _manifest(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}
    return {
        "samplesheet_sha256": data.get("samplesheet_sha256"),
        "source_file_count": len(data.get("source_sha256") or {}),
        "conda_environment_count": len(data.get("resolved_conda_packages") or {}),
    }


def parse_outputs(results_dir: Path, config: dict) -> dict:
    """Read every table the run was expected to produce.

    A stage that did not run maps to None, never to an empty table that reads
    like "nothing was found".
    """
    results_dir = Path(results_dir)
    found = _check_present(results_dir, config)

    read_counts = _read_table(found[FINAL_OUTPUTS["read_counts"]])
    sourmash = (
        _matrix(found[FINAL_OUTPUTS["sourmash_similarities"]])
        if FINAL_OUTPUTS["sourmash_similarities"] in found
        else None
    )

    profile = None
    if FINAL_OUTPUTS["profile"] in found:
        profile = _profile(found[FINAL_OUTPUTS["profile"]])

    assembly = None
    if FINAL_OUTPUTS["quast"] in found:
        units, metrics = _quast(found[FINAL_OUTPUTS["quast"]])
        assembly = {
            "units": units,
            "quast": {
                unit: {metric: values.get(unit) for metric, values in metrics.items()}
                for unit in units
            },
            "metrics": metrics,
        }

    bins = None
    if FINAL_OUTPUTS["bins"] in found:
        mags = _read_table(found[FINAL_OUTPUTS["bins"]])
        bins = {
            "count": len(mags),
            "genomes": [row.get("genome") for row in mags],
            "mags": mags,
            "checkm2": _read_table(found[FINAL_OUTPUTS["checkm2_quality"]]),
            "drep": _read_table(found[FINAL_OUTPUTS["drep_genome_info"]]),
            "gtdb": _read_table(found[FINAL_OUTPUTS["gtdb_summary"]]),
            "passing_qa": sum(
                1 for row in mags if str(row.get("passes_qa")).lower() == "true"
            ),
            "thresholds": thresholds_from_config(config),
        }

    taxonomy = None
    if FINAL_OUTPUTS["contig_taxonomy"] in found:
        taxonomy = {"contigs": _read_table(found[FINAL_OUTPUTS["contig_taxonomy"]])}

    annotation = None
    if FINAL_OUTPUTS["function_abundance"] in found:
        annotation = {"functions": _read_table(found[FINAL_OUTPUTS["function_abundance"]])}

    heterogeneity = None
    if FINAL_OUTPUTS["strain_heterogeneity"] in found:
        heterogeneity = {
            "rows": _read_table(found[FINAL_OUTPUTS["strain_heterogeneity"]]),
        }

    provenance = {
        "pipeline": _pipeline_revision(found[FINAL_OUTPUTS["pipeline_revision"]]),
        "tools": _read_table(found[FINAL_OUTPUTS["software_versions"]]),
        "references": _references(found[FINAL_OUTPUTS["reference_db_md5"]]),
        "runtime": _read_table(found[FINAL_OUTPUTS["runtime"]]),
        "manifest": _manifest(found[FINAL_OUTPUTS["run_manifest"]]),
    }

    return {
        "qc": {"read_counts": read_counts, "sourmash": sourmash},
        "profile": profile,
        "assembly": assembly,
        "bins": bins,
        "diagnostic": {"bottleneck": _read_table(found[FINAL_OUTPUTS["bottleneck"]])},
        "taxonomy": taxonomy,
        "annotation": annotation,
        "heterogeneity": heterogeneity,
        "provenance": provenance,
    }
