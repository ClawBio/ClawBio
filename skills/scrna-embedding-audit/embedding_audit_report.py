"""Output rendering for scrna-embedding-audit.

No scanpy dependency at import time. Takes the audit result produced by
``scrna_embedding_audit.run_audit`` and writes strictly valid JSON, CSV tables,
one figure and the markdown report. Every table row carries a ``provenance``
column so that a detached CSV still says where each number came from.
"""

from __future__ import annotations

import csv
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

from clawbio.common.report import DISCLAIMER

GROUP_TITLES = {
    "bio": "biological conservation",
    "batch": "batch mixing",
    "transfer": "cross-batch label transfer",
}
GROUP_ORDER = ("bio", "batch", "transfer")

SCOPE_SENTENCE = (
    "Resampling varies only which cells are scored, with all embeddings held fixed; it does not "
    "cover embedding training randomness, donor-to-donor variation, or the choice of k. Verdicts "
    "describe this cell sample and these fixed embeddings."
)
BASELINE_SENTENCE = (
    "The baseline is an unintegrated PCA of the same cells, so any batch-aware method is expected "
    "to score higher on batch mixing; the question this run answers is by how much and at what "
    "cost to biological conservation."
)
CIRCULARITY_SENTENCE = (
    "If an embedding was trained on these labels (scANVI, a fine-tuned model) or the labels were "
    "derived from one of the scored embeddings, its biological-conservation scores are circular; "
    "this skill cannot detect that."
)
NO_BATCH_SENTENCE = "Batch mixing was not compared: no `--batch-key` was given."


def _as_json_value(value: Any) -> Any:
    """Convert NumPy/Pandas values to JSON values, representing non-finite as null."""
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, np.ndarray):
        return _as_json_value(value.tolist())
    if isinstance(value, dict):
        return {str(key): _as_json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_as_json_value(item) for item in value]
    return value


def _num(value: Any, digits: int = 3) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "null"
    return f"{number:.{digits}f}" if np.isfinite(number) else "null"


def _csv_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return "" if not np.isfinite(value) else repr(value)
    return str(value)


def _write_tables(result: dict, output_dir: Path, *, baseline_key: str) -> list[Path]:
    tables = output_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    groups = result["metric_groups"]
    written: list[Path] = []

    path = tables / "metrics.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            ["embedding", "metric", "group", "value", "reason", "provenance"]
        )
        for key, entry in result["embeddings"].items():
            for metric, value in entry["metrics"].items():
                writer.writerow(
                    [
                        key,
                        metric,
                        groups[metric],
                        _csv_value(value),
                        entry["reasons"].get(metric, ""),
                        entry["provenance"],
                    ]
                )
    written.append(path)

    path = tables / "verdicts.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "embedding",
                "metric",
                "group",
                "verdict",
                "full_value",
                "baseline_value",
                "median_diff",
                "sd_diff",
                "frac_positive",
                "n_resamples",
                "reason",
                "provenance",
            ]
        )
        for key, entry in result["verdicts"].items():
            for metric, v in entry["per_metric"].items():
                writer.writerow(
                    [
                        key,
                        metric,
                        v["group"],
                        v["verdict"],
                        _csv_value(v["full_value"]),
                        _csv_value(v["baseline_value"]),
                        _csv_value(v["median_diff"]),
                        _csv_value(v["sd_diff"]),
                        _csv_value(v["frac_positive"]),
                        v["n_resamples"],
                        v["reason"],
                        f"{result['embeddings'][key]['provenance']} vs {baseline_key}",
                    ]
                )
    written.append(path)

    path = tables / "resamples.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["resample", "embedding", "metric", "value", "provenance"])
        for key, per_metric in result["resamples"].items():
            prov = result["embeddings"][key]["provenance"]
            for metric, values in per_metric.items():
                for r, value in enumerate(values):
                    writer.writerow([r, key, metric, _csv_value(value), prov])
    written.append(path)
    return written


def _write_figure(
    result: dict, output_dir: Path, *, baseline_key: str, control_key: str
) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures = output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    metrics = list(result["metric_definitions"].keys())
    keys = list(result["embedding_keys"]) + [baseline_key, control_key]
    n_keys = len(keys)
    width = 0.8 / n_keys
    fig, ax = plt.subplots(figsize=(max(7.0, 1.1 * len(metrics)), 4.2))
    x = np.arange(len(metrics))
    for i, key in enumerate(keys):
        values = [result["embeddings"][key]["metrics"][m] for m in metrics]
        heights = [v if v is not None else 0.0 for v in values]
        colour = (
            "0.55" if key == baseline_key else ("0.85" if key == control_key else None)
        )
        bars = ax.bar(
            x + (i - (n_keys - 1) / 2) * width,
            heights,
            width,
            label=key,
            color=colour,
            edgecolor="black",
            linewidth=0.4,
        )
        for bar, v in zip(bars, values):
            if v is None:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    0.02,
                    "null",
                    ha="center",
                    va="bottom",
                    fontsize=6,
                    rotation=90,
                )
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("metric value (higher is better)")
    ax.set_ylim(
        min(
            -0.05,
            min(
                (
                    v
                    for k in keys
                    for v in result["embeddings"][k]["metrics"].values()
                    if v is not None
                ),
                default=0.0,
            )
            - 0.05,
        ),
        1.05,
    )
    ax.axhline(0.0, color="black", linewidth=0.5)
    ax.set_title("Embeddings vs in-run PCA baseline (full data)")
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    path = figures / "metrics_vs_baseline.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return [path]


def _counts_sentence(key: str, counts_by_group: dict) -> str:
    parts = []
    for g in GROUP_ORDER:
        c = counts_by_group[g]
        parts.append(
            f"{GROUP_TITLES[g]}: higher than the PCA baseline on {c['higher']} of {c['n_metrics']} metrics, "
            f"lower on {c['lower']}, not separable on {c['not separable']}, tied on {c['tied']}, "
            f"not evaluated on {c['not evaluated']}"
        )
    return f"`{key}` — " + "; ".join(parts) + "."


def _render_report(
    result: dict,
    *,
    source_label: str,
    demo: bool,
    skill_name: str,
    skill_version: str,
    timestamp: str,
    baseline_key: str,
    control_key: str,
) -> str:
    emb = result["embeddings"]
    batch_line = (
        f"`{result['batch']['key']}` ({result['batch']['n_unique']} unique)"
        if result["batch"]["key"]
        else "none given"
    )
    cells_line = f"{result['n_cells']} in file, {result['cells_scored']} scored"
    if result["subsampled"]:
        cells_line += " (stratified subsample)"
    if result["cells_dropped"]:
        cells_line += f", dropped {result['cells_dropped']}"
    lines: list[str] = [
        f"# {skill_name} report",
        "",
        f"- **Generated**: {timestamp}",
        f"- **Skill**: {skill_name} v{skill_version}",
        f"- **Input**: {source_label}" + (" (no measured data)" if demo else ""),
        f"- **Cells**: {cells_line}",
        f"- **Labels**: `{result['labels']['key']}` ({result['labels']['n_unique']} unique)",
        f"- **Batch**: {batch_line}",
        (
            f"- **Baseline**: {baseline_key}, {result['baseline']['recipe']} "
            f"({result['baseline']['n_hvg_used']} genes, {result['baseline']['n_pcs_used']} PCs); "
            f"provenance: {result['baseline']['provenance']}"
        ),
        (
            f"- **Metric definitions**: scib-metrics {result['references']['scib_metrics_fixture_version']} "
            "(fixture-checked re-implementation; see table below)"
        ),
        "",
        "## Summary",
        "",
    ]
    if demo:
        lines += [
            "This run used synthetic data generated from `examples/demo_spec.json`; it demonstrates the report, not a biological result.",
            "",
        ]
    for key in result["embedding_keys"]:
        lines.append(
            "- " + _counts_sentence(key, result["verdicts"][key]["counts_by_group"])
        )
    lines += [
        "",
        "NMI and ARI share one KMeans clustering. Counts are not comparable across embeddings; read the table.",
        "",
        SCOPE_SENTENCE,
        "",
        BASELINE_SENTENCE,
        "",
        CIRCULARITY_SENTENCE,
        "",
    ]
    if not result["batch"]["key"]:
        lines += [NO_BATCH_SENTENCE, ""]

    lines += [
        "## Verdicts",
        "",
        "Rule: "
        + f"{result['parameters']['n_resamples']} stratified half-samples; a metric is `higher`/`lower` when every "
        f"half-sample difference shares the sign and |median| > {result['effective']['sd_factor']:g} x SD of the differences "
        f"(differences within {result['effective']['diff_tolerance']:g} count as zero); `tied` when both values are at or above "
        f"{result['effective']['ceiling']:g} or every difference is zero; otherwise `not separable`. "
        "The 2 x SD factor, the tolerance and the ceiling are this skill's choices.",
        "",
        "| embedding | metric | group | verdict | value | baseline | median diff | SD diff | frac > 0 | note |",
        "|---|---|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for key in result["embedding_keys"]:
        for metric, v in result["verdicts"][key]["per_metric"].items():
            lines.append(
                f"| `{key}` | {metric} | {GROUP_TITLES[v['group']]} | **{v['verdict']}** | {_num(v['full_value'])} | "
                f"{_num(v['baseline_value'])} | {_num(v['median_diff'], 4)} | {_num(v['sd_diff'], 4)} | "
                f"{_num(v['frac_positive'], 2)} | {v['reason']} |"
            )

    lines += [
        "",
        "## Full-data values",
        "",
        "| embedding | provenance | "
        + " | ".join(result["metric_definitions"].keys())
        + " |",
        "|---|---|" + "---:|" * len(result["metric_definitions"]),
    ]
    for key, entry in emb.items():
        values = " | ".join(
            _num(entry["metrics"][m]) for m in result["metric_definitions"]
        )
        lines.append(f"| `{key}` | {entry['provenance']} | {values} |")
    lines += [
        "",
        f"`{control_key}` is a seeded random embedding scored for scale only; it never enters a verdict.",
        "",
    ]

    null_rows = [(k, m, r) for k, e in emb.items() for m, r in e["reasons"].items()]
    if null_rows:
        lines += [
            "### Metrics not computed",
            "",
            "| embedding | metric | reason |",
            "|---|---|---|",
        ]
        lines += [f"| `{k}` | {m} | {r} |" for k, m, r in null_rows]
        lines.append("")

    lines += [
        "## Metric definitions",
        "",
        "| metric | group | definition | source |",
        "|---|---|---|---|",
    ]
    for m, d in result["metric_definitions"].items():
        lines.append(
            f"| {m} | {GROUP_TITLES[result['metric_groups'][m]]} | {d['definition']} | {d['source']} |"
        )
    lines += [
        "",
        "Not included: "
        + "; ".join(
            f"`{k}` ({v})" for k, v in result["metric_definitions_omitted"].items()
        ),
        "",
    ]

    if result["warnings"]:
        lines += ["## Warnings", ""] + [f"- {w}" for w in result["warnings"]] + [""]

    lines += [
        "## Parameters",
        "",
        "| parameter | value |",
        "|---|---|",
    ]
    for k, v in result["parameters"].items():
        lines.append(f"| `{k}` | {v if v is not None else 'default'} |")
    for k, v in result["effective"].items():
        lines.append(f"| effective `{k}` | {v} |")

    lines += [
        "",
        "## Files",
        "",
        "| file | contents |",
        "|---|---|",
        "| `result.json` | every number in this report, plus definitions and provenance |",
        "| `tables/metrics.csv` | full-data metric values per embedding |",
        "| `tables/verdicts.csv` | verdict, values and resampling statistics per embedding and metric |",
        "| `tables/resamples.csv` | metric value in every half-sample |",
        "| `figures/metrics_vs_baseline.png` | full-data values, baseline and random control |",
        "| `reproducibility/` | `commands.sh`, `environment.yml`, `checksums.sha256`, `run_manifest.json` |",
        "",
        "## References",
        "",
        f"- Luecken et al., Nat Methods 2022, doi:{result['references']['luecken_2022']} (scIB metric definitions)",
        f"- Korsunsky et al., Nat Methods 2019, doi:{result['references']['korsunsky_2019_harmony_lisi']} (LISI)",
        f"- scib-metrics {result['references']['scib_metrics_fixture_version']} (reference implementation used for the fixture)",
        f"- Kedzierska et al., Genome Biology 2025, doi:{result['references']['kedzierska_2025_zero_shot_fm_benchmark']} (why a baseline is run alongside)",
        "",
        "---",
        "",
        f"*{DISCLAIMER}*",
        "",
    ]
    return "\n".join(lines)


def write_outputs(
    result: dict,
    output_dir: Path,
    *,
    source_label: str,
    demo: bool,
    skill_name: str,
    skill_version: str,
    baseline_key: str,
    control_key: str,
) -> tuple[list[Path], Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC")
    tables = _write_tables(result, output_dir, baseline_key=baseline_key)
    figures = _write_figure(
        result, output_dir, baseline_key=baseline_key, control_key=control_key
    )
    report_path = output_dir / "report.md"
    report_path.write_text(
        _render_report(
            result,
            source_label=source_label,
            demo=demo,
            skill_name=skill_name,
            skill_version=skill_version,
            timestamp=timestamp,
            baseline_key=baseline_key,
            control_key=control_key,
        ),
        encoding="utf-8",
    )
    payload = _as_json_value(
        {
            "skill": skill_name,
            "skill_version": skill_version,
            "demo": demo,
            "input": source_label,
            "timestamp": timestamp,
            **{k: v for k, v in result.items() if k != "resamples"},
            "disclaimer": DISCLAIMER,
        }
    )
    result_path = output_dir / "result.json"
    result_path.write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return tables + figures + [report_path, result_path], report_path, result_path
