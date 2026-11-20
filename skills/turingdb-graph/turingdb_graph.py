#!/usr/bin/env python3
"""turingdb-graph — Build, query, and analyse biomedical graphs in TuringDB.

ClawBio skill entry point. Subcommands:
    --build         Ingest CSV/TSV/GML/JSONL into a named TuringDB graph
    --query         Run an arbitrary Cypher query against a named graph
    --analyse-cohort  Run a fixed set of clinical-cohort analyses
    --demo          Run an end-to-end example on one of the three shipped datasets
    --stop-server   Stop a daemon started by this skill

Global flags: --host (default localhost:6666), --data-dir (default ~/.turing),
--no-auto-start, --output (output directory, default ./output; --out also works).
"""
from __future__ import annotations

import argparse
import json
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit

import pandas as pd

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    ReproCommand,
    ReproPath,
    write_checksums,
    write_environment_yml,
    write_portable_commands_sh,
)

DISCLAIMER = (
    "ClawBio is a research and educational tool. Not a medical device. "
    "This output must not be used for clinical decision-making."
)

DEFAULT_HOST = "localhost:6666"
DEFAULT_PORT = 6666
DEFAULT_DATA_DIR = Path.home() / ".turing"


# ---------------------------------------------------------------------------
# Connection / server lifecycle
# ---------------------------------------------------------------------------

def split_host(host: str) -> tuple[str, int]:
    """Split ``--host`` into the host name and port the turingdb 3.x client takes.

    Accepts ``localhost:6666``, ``localhost`` and, as earlier releases of this skill
    documented, ``http://localhost:6666``; turingdb 3.x rejects a URL as ``host``.
    """
    parts = urlsplit(host if "://" in host else f"//{host}")
    if not parts.hostname:
        raise ValueError(f"--host must be a host name or host:port, got {host!r}")
    return parts.hostname, parts.port or DEFAULT_PORT


@dataclass
class Connection:
    """Wraps a TuringDB client. Lazy-imports turingdb so --help works without it."""
    host: str
    data_dir: Path
    client: Any = None
    # In-process engine on data_dir instead of a daemon. turingdb allows one
    # embedded engine per process, so this suits a single CLI run such as --demo.
    embedded: bool = False

    def _client(self):
        from turingdb import TuringDB  # lazy import — keeps --help cheap
        if self.embedded:
            return TuringDB(type="embedded", data_dir=str(self.data_dir))
        hostname, port = split_host(self.host)
        return TuringDB(host=hostname, port=port)

    def connect(self, auto_start: bool = True, retries: int = 30) -> None:
        if self.embedded:
            self.client = self._client()
            return
        try:
            self.client = self._client()
            self.client.list_loaded_graphs()  # probe
            return
        except Exception:
            if not auto_start:
                raise RuntimeError(
                    f"Could not reach TuringDB at {self.host}. Start it manually or re-run "
                    "without --no-auto-start."
                )

        # Start the daemon
        binary = _find_turingdb_binary()
        self.data_dir.mkdir(parents=True, exist_ok=True)
        print(f"→ starting turingdb -demon (data dir: {self.data_dir})", file=sys.stderr)
        subprocess.run(
            [binary, "start", "-turing-dir", str(self.data_dir), "-demon"],
            check=True,
        )

        # Poll for the server to come up
        for _ in range(retries):
            try:
                self.client = self._client()
                self.client.list_loaded_graphs()
                return
            except Exception:
                time.sleep(0.5)
        raise RuntimeError(
            f"TuringDB daemon did not become reachable at {self.host} within "
            f"{retries * 0.5:.0f}s."
        )

    def stop(self) -> None:
        binary = _find_turingdb_binary()
        print(f"→ stopping turingdb (data dir: {self.data_dir})", file=sys.stderr)
        subprocess.run(
            [binary, "stop", "-turing-dir", str(self.data_dir)],
            check=False,
        )


def _find_turingdb_binary() -> str:
    """Locate the turingdb CLI — PATH, then common venv paths."""
    for candidate in ("turingdb", ".venv/bin/turingdb", "venv/bin/turingdb"):
        if shutil.which(candidate) or Path(candidate).exists():
            resolved = shutil.which(candidate) or str(Path(candidate).resolve())
            return resolved
    raise FileNotFoundError(
        "turingdb binary not found. Install it with `pip install turingdb` or `uv add turingdb`."
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _ensure_graph_absent(conn: Connection, graph: str) -> None:
    """Refuse to overwrite an existing graph (Safety Rule 6)."""
    if graph in conn.client.list_available_graphs():
        raise RuntimeError(
            f"Graph {graph!r} already exists. Safety Rule 6: --build is additive, "
            f"not destructive. Choose a new --graph name or drop it manually."
        )


def _write_report(out_dir: Path, md: str, summary: dict) -> tuple[Path, Path]:
    """Write the two canonical outputs: report.md + summary.json."""
    out_dir.mkdir(parents=True, exist_ok=True)
    md_path = out_dir / "report.md"
    json_path = out_dir / "summary.json"
    md_path.write_text(md + f"\n\n---\n\n*{DISCLAIMER}*\n")
    json_path.write_text(json.dumps(summary, indent=2, default=str))
    return md_path, json_path


def _infer_numeric_wrap(df_sample: pd.DataFrame, col: str) -> str:
    """Given a pandas DataFrame sample, decide whether to wrap a column in toInteger/toFloat."""
    dtype = df_sample[col].dtype
    if pd.api.types.is_integer_dtype(dtype):
        return f"toInteger(row.{col})"
    if pd.api.types.is_float_dtype(dtype):
        return f"toFloat(row.{col})"
    return f"row.{col}"


def _cypher_safe(s: str) -> str:
    """Minimal escape for string literals we interpolate into Cypher. Not a full injection guard."""
    return s.replace("\\", "\\\\").replace("'", "\\'")


# ---------------------------------------------------------------------------
# --build
# ---------------------------------------------------------------------------

def build_graph(
    conn: Connection,
    input_path: Path,
    graph: str,
    fmt: str = "auto",
    node_label: str | None = None,
    out_dir: Path = Path("output"),
) -> dict:
    """Ingest input_path into a new TuringDB graph named `graph`."""
    if fmt == "auto":
        fmt = input_path.suffix.lstrip(".").lower()
        if fmt == "tsv":
            fmt = "csv"  # treated the same after read

    _ensure_graph_absent(conn, graph)

    # Stage the file under the TuringDB data directory if it isn't already there
    data_dir_path = conn.data_dir / "data"
    data_dir_path.mkdir(parents=True, exist_ok=True)
    staged = data_dir_path / input_path.name
    if input_path.resolve() != staged.resolve():
        shutil.copy(input_path, staged)

    summary: dict = {"graph": graph, "format": fmt, "source": input_path.name}

    if fmt in ("gml",):
        conn.client.query(f"LOAD GML '{staged.name}' AS {graph}")
        summary["loader"] = "LOAD GML"
        summary["note"] = (
            "GML flattens labels/types to GMLNode/GMLEdge and all properties to strings. "
            "For typed biomedical data prefer JSONL."
        )

    elif fmt in ("jsonl",):
        conn.client.query(f"LOAD JSONL '{staged.name}' AS {graph}")
        summary["loader"] = "LOAD JSONL"

    elif fmt in ("csv",):
        if not node_label:
            raise RuntimeError(
                "CSV ingest requires --node-label. One node will be created per row with "
                "that label. Columns become properties (strings, or numeric where pandas infers)."
            )
        # Create the graph first, then LOAD CSV + CREATE inside a change
        conn.client.create_graph(graph)
        conn.client.set_graph(graph)

        df_sample = pd.read_csv(staged, nrows=200, sep="\t" if staged.suffix == ".tsv" else ",")
        column_spec = ", ".join(
            f"{col}: {_infer_numeric_wrap(df_sample, col)}" for col in df_sample.columns
        )

        change = conn.client.new_change()
        conn.client.checkout(change=change)
        conn.client.query(
            f"LOAD CSV '{staged.name}' WITH HEADERS AS row\n"
            f"CREATE (:{node_label} {{{column_spec}}})"
        )
        conn.client.query("CHANGE SUBMIT")
        conn.client.checkout()

        summary["loader"] = "LOAD CSV + CREATE"
        summary["node_label"] = node_label
        summary["columns"] = list(df_sample.columns)
        summary["row_count"] = len(pd.read_csv(
            staged, sep="\t" if staged.suffix == ".tsv" else ","
        ))

    else:
        raise RuntimeError(f"Unsupported format: {fmt!r}. Use csv, tsv, gml, or jsonl.")

    # Ensure the graph is loaded and current
    loaded = conn.client.list_loaded_graphs()
    if graph not in loaded:
        conn.client.load_graph(graph)
    conn.client.set_graph(graph)

    # Capture schema + counts
    summary.update(_summarise_graph(conn.client))

    # Resolve the latest commit for audit (surface ID so callers can time-travel)
    try:
        hist = conn.client.query("CALL db.history()")
        if not hist.empty:
            summary["commit"] = str(hist.iloc[0, 0])
    except Exception:
        pass  # history is best-effort

    # Report
    md = _render_build_report(summary)
    _write_report(out_dir, md, summary)
    return summary


def _summarise_graph(client) -> dict:
    """Return counts per label / edge type + schema introspection."""
    labels = client.query("CALL db.labels()")
    edge_types = client.query("CALL db.edgeTypes()")
    props = client.query("CALL db.propertyTypes()")

    n_nodes = int(client.query("MATCH (n) RETURN count(n)").iloc[0, 0])
    n_edges = int(client.query("MATCH ()-->() RETURN count(*)").iloc[0, 0])

    by_label = {}
    for label in labels["label"].tolist():
        c = int(client.query(f"MATCH (n:{label}) RETURN count(n)").iloc[0, 0])
        by_label[label] = c

    return {
        "nodes": n_nodes,
        "edges": n_edges,
        "labels": labels["label"].tolist(),
        "edge_types": edge_types["edgeType"].tolist() if "edgeType" in edge_types.columns else [],
        "properties": props.to_dict(orient="records"),
        "nodes_by_label": by_label,
    }


def _render_build_report(summary: dict) -> str:
    lines = [f"# Graph build: `{summary['graph']}`", ""]
    lines.append(f"- **Source**: `{summary['source']}` ({summary['format']})")
    lines.append(f"- **Loader**: `{summary['loader']}`")
    if "commit" in summary:
        lines.append(f"- **Commit**: `{summary['commit']}`")
    lines.append(f"- **Nodes**: {summary.get('nodes', 0):,}")
    lines.append(f"- **Edges**: {summary.get('edges', 0):,}")
    lines.append("")

    lines.append("## Nodes by label")
    lines.append("")
    if summary.get("nodes_by_label"):
        lines.append("| Label | Count |")
        lines.append("|-------|------:|")
        for label, c in sorted(summary["nodes_by_label"].items(), key=lambda x: -x[1]):
            lines.append(f"| `{label}` | {c:,} |")
    else:
        lines.append("_No labels found._")
    lines.append("")

    if summary.get("edge_types"):
        lines.append("## Edge types")
        lines.append("")
        lines.append(", ".join(f"`{e}`" for e in summary["edge_types"]))
        lines.append("")

    if "note" in summary:
        lines.append(f"> **Note**: {summary['note']}")
        lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# --query
# ---------------------------------------------------------------------------

def run_query(
    conn: Connection,
    graph: str,
    cypher: str,
    fmt: str = "md",
    out_dir: Path = Path("output"),
) -> dict:
    """Run a Cypher query and serialise the result."""
    if graph not in conn.client.list_available_graphs():
        raise RuntimeError(f"Graph {graph!r} not found. Available: {conn.client.list_available_graphs()}")
    if graph not in conn.client.list_loaded_graphs():
        conn.client.load_graph(graph)
    conn.client.set_graph(graph)

    df = conn.client.query(cypher)

    summary = {
        "graph": graph,
        "cypher": cypher,
        "rows": len(df),
        "columns": df.columns.tolist(),
    }

    # Write outputs in the requested format(s)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "result.json").write_text(df.to_json(orient="records", indent=2))
    (out_dir / "result.tsv").write_text(df.to_csv(sep="\t", index=False))

    if fmt == "md":
        md = "# Query result\n\n"
        md += f"```cypher\n{cypher}\n```\n\n"
        md += f"**{len(df)} row(s), {len(df.columns)} column(s)**\n\n"
        md += df.to_markdown(index=False) if len(df) else "_(empty result)_"
        md_path, json_path = _write_report(out_dir, md, summary)
        summary["outputs"] = {"report": str(md_path), "summary": str(json_path)}
    else:
        summary["outputs"] = {
            "json": str(out_dir / "result.json"),
            "tsv": str(out_dir / "result.tsv"),
        }

    return summary


# ---------------------------------------------------------------------------
# --analyse-cohort
# ---------------------------------------------------------------------------

# Cypher recipes for the fixed cohort analyses. Explicit rather than templated so an
# auditor can read exactly what aggregation was performed.
COHORT_QUERIES = {
    "patient_count":
        "MATCH (p:Patient) RETURN count(p) AS n",
    "age_stats":
        "MATCH (p:Patient) WHERE p.age IS NOT NULL RETURN p.age AS age",
    "junior_count":
        "MATCH (p:Patient) WHERE p.age < 18 RETURN count(p) AS n",
    "senior_count":
        "MATCH (p:Patient) WHERE p.age > 65 RETURN count(p) AS n",
    "top_conditions":
        "MATCH (p:Patient)-[:HAS]->(c:MedicalCondition) "
        "RETURN p.id AS patient, c.displayName AS condition",
    "top_medications":
        "MATCH (p:Patient)-[:TOOK_MEDICATION]->(m:Medication) "
        "RETURN p.id AS patient, m.displayName AS medication",
    "top_comorbidities":
        "MATCH (p:Patient)-[:HAS]->(c1:MedicalCondition), "
        "(p)-[:HAS]->(c2:MedicalCondition) "
        "RETURN p.id AS patient, c1.displayName AS condition_a, c2.displayName AS condition_b",
    "top_comedications":
        "MATCH (p:Patient)-[:TOOK_MEDICATION]->(m1:Medication), "
        "(p)-[:TOOK_MEDICATION]->(m2:Medication) "
        "RETURN p.id AS patient, m1.displayName AS medication_a, m2.displayName AS medication_b",
}


def analyse_cohort(
    conn: Connection,
    graph: str,
    top_n: int = 10,
    out_dir: Path = Path("output"),
) -> dict:
    """Run the fixed cohort analyses on a patient-centric graph."""
    if graph not in conn.client.list_loaded_graphs():
        conn.client.load_graph(graph)
    conn.client.set_graph(graph)

    # Schema precondition — fail fast with a helpful message
    labels = conn.client.query("CALL db.labels()")["label"].tolist()
    if "Patient" not in labels:
        raise RuntimeError(
            f"Graph {graph!r} has no `Patient` nodes. `--analyse-cohort` expects a "
            f"patient-centric clinical schema — see reference/biomedical.md § 2a. "
            f"Labels found: {labels}"
        )

    patient_count = int(conn.client.query(COHORT_QUERIES["patient_count"]).iloc[0, 0])
    ages = conn.client.query(COHORT_QUERIES["age_stats"])["age"]

    age_stats = None
    if len(ages) > 0:
        ages_num = pd.to_numeric(ages, errors="coerce").dropna()
        if len(ages_num) > 0:
            age_stats = {
                "n": int(len(ages_num)),
                "min": int(ages_num.min()),
                "max": int(ages_num.max()),
                "mean": float(round(ages_num.mean(), 1)),
                "median": float(ages_num.median()),
            }

    # Only run analyses whose edge type exists in the graph
    edge_types = set(conn.client.query("CALL db.edgeTypes()")["edgeType"].tolist())

    top_conditions = []
    if "HAS" in edge_types and "MedicalCondition" in labels:
        df = conn.client.query(COHORT_QUERIES["top_conditions"])
        top_conditions = (df.drop_duplicates().groupby("condition")["patient"].nunique()
                          .rename("patients").sort_values(ascending=False)
                          .head(top_n).reset_index().to_dict(orient="records"))

    top_medications = []
    if "TOOK_MEDICATION" in edge_types and "Medication" in labels:
        df = conn.client.query(COHORT_QUERIES["top_medications"])
        top_medications = (df.drop_duplicates().groupby("medication")["patient"].nunique()
                           .rename("patients").sort_values(ascending=False)
                           .head(top_n).reset_index().to_dict(orient="records"))

    top_comorbidities = []
    if "HAS" in edge_types and "MedicalCondition" in labels:
        df = conn.client.query(COHORT_QUERIES["top_comorbidities"])
        df = df[df["condition_a"] != df["condition_b"]].copy()
        lo = df[["condition_a", "condition_b"]].min(axis=1)
        hi = df[["condition_a", "condition_b"]].max(axis=1)
        df["condition_a"] = lo.values
        df["condition_b"] = hi.values
        df_pair = (df.drop_duplicates()
                   .groupby(["condition_a", "condition_b"])["patient"].nunique()
                   .rename("patients").reset_index())
        df_pair = df_pair[df_pair["patients"] >= 2].sort_values("patients", ascending=False).head(top_n)
        top_comorbidities = df_pair.to_dict(orient="records")

    top_comedications = []
    if "TOOK_MEDICATION" in edge_types and "Medication" in labels:
        df = conn.client.query(COHORT_QUERIES["top_comedications"])
        df = df[df["medication_a"] != df["medication_b"]].copy()
        lo = df[["medication_a", "medication_b"]].min(axis=1)
        hi = df[["medication_a", "medication_b"]].max(axis=1)
        df["medication_a"] = lo.values
        df["medication_b"] = hi.values
        df_pair = (df.drop_duplicates()
                   .groupby(["medication_a", "medication_b"])["patient"].nunique()
                   .rename("patients").reset_index())
        df_pair = df_pair[df_pair["patients"] >= 2].sort_values("patients", ascending=False).head(top_n)
        top_comedications = df_pair.to_dict(orient="records")

    junior = int(conn.client.query(COHORT_QUERIES["junior_count"]).iloc[0, 0])
    senior = int(conn.client.query(COHORT_QUERIES["senior_count"]).iloc[0, 0])

    summary = {
        "graph": graph,
        "patient_count": patient_count,
        "age_stats": age_stats,
        "junior_count": junior,
        "senior_count": senior,
        "top_conditions": top_conditions,
        "top_medications": top_medications,
        "top_comorbidities": top_comorbidities,
        "top_comedications": top_comedications,
    }

    md = _render_cohort_report(summary, top_n)
    _write_report(out_dir, md, summary)
    return summary


def _render_cohort_report(s: dict, top_n: int) -> str:
    lines = [f"# Cohort analysis: `{s['graph']}`", ""]
    lines.append(f"- **Patients**: {s['patient_count']:,}")
    if s["age_stats"]:
        a = s["age_stats"]
        lines.append(
            f"- **Ages** (n={a['n']:,}): min {a['min']}, max {a['max']}, "
            f"mean {a['mean']}, median {a['median']}"
        )
    lines.append(f"- **Under 18**: {s['junior_count']:,}")
    lines.append(f"- **Over 65**: {s['senior_count']:,}")
    lines.append("")

    def _table(title: str, rows: list, cols: list):
        if not rows:
            return
        lines.append(f"## {title}")
        lines.append("")
        lines.append("| " + " | ".join(cols) + " |")
        lines.append("|" + "|".join(["---"] * len(cols)) + "|")
        for r in rows:
            lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
        lines.append("")

    _table(f"Top {top_n} conditions", s["top_conditions"], ["condition", "patients"])
    _table(f"Top {top_n} medications", s["top_medications"], ["medication", "patients"])
    _table(f"Top {top_n} comorbidities (co-occurring in ≥ 2 patients)",
           s["top_comorbidities"], ["condition_a", "condition_b", "patients"])
    _table(f"Top {top_n} comedications (co-prescribed in ≥ 2 patients)",
           s["top_comedications"], ["medication_a", "medication_b", "patients"])

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# --demo
# ---------------------------------------------------------------------------

DEMO_SPECS = {
    "cohort": {
        "file": "demo/cohort.csv",
        "fmt": "csv",
        "node_label": "PatientRow",  # intermediate — rewritten below via CREATE-from-rows
        "graph": "demo_cohort",
        "analyse": True,
    },
    "pathway": {
        "file": "demo/pathway.gml",
        "fmt": "gml",
        "graph": "demo_pathway",
        "analyse": False,
    },
    "antibody": {
        "file": "demo/antibody.csv",
        "fmt": "csv",
        "node_label": "AntibodyRow",
        "graph": "demo_antibody",
        "analyse": False,
    },
}


def run_demo(
    conn: Connection,
    which: str,
    skill_root: Path,
    out_dir: Path = Path("output"),
) -> dict:
    if which not in DEMO_SPECS:
        raise RuntimeError(f"Unknown demo {which!r}. Choose from {list(DEMO_SPECS)}.")
    spec = DEMO_SPECS[which]
    input_path = skill_root / spec["file"]
    if not input_path.exists():
        raise FileNotFoundError(f"Demo file missing: {input_path}")

    # Pick a fresh graph name so repeated demo runs don't collide
    base = spec["graph"]
    existing = set(conn.client.list_available_graphs()) if conn.client else set()
    graph = base
    i = 1
    while graph in existing:
        i += 1
        graph = f"{base}_v{i}"

    demo_out = out_dir / f"demo-{which}"

    if which == "cohort":
        # Cohort demo needs a richer graph (Patient + Medication + MedicalCondition + edges).
        # Build it with explicit CREATE from the DataFrame so --analyse-cohort has something
        # to chew on. This is the one demo where --build's generic CSV path is too flat.
        summary = _build_cohort_demo(conn, input_path, graph, demo_out)
        if spec["analyse"]:
            cohort_summary = analyse_cohort(conn, graph, out_dir=demo_out / "analysis")
            summary["cohort_analysis"] = cohort_summary
    else:
        summary = build_graph(
            conn, input_path, graph, fmt=spec["fmt"],
            node_label=spec.get("node_label"), out_dir=demo_out,
        )
        # Run one illustrative query
        sample_q = _demo_sample_query(which, graph)
        if sample_q:
            qres = run_query(conn, graph, sample_q, fmt="md",
                             out_dir=demo_out / "sample_query")
            summary["sample_query"] = qres

    print(f"✓ demo '{which}' complete → {demo_out}")
    return summary


def _demo_sample_query(which: str, graph: str) -> str | None:
    if which == "pathway":
        return (
            "MATCH (a)-->(b)-->(c) "
            "RETURN a.`displayName (String)` AS start, "
            "b.`displayName (String)` AS via, "
            "c.`displayName (String)` AS end LIMIT 10"
        )
    if which == "antibody":
        return (
            "MATCH (a:AntibodyRow) "
            "RETURN a.antibody_name AS antibody, a.target_protein AS protein, "
            "a.gene_name AS gene, a.cited_by AS citations "
            "ORDER BY a.cited_by DESC LIMIT 10"
        )
    return None


def _build_cohort_demo(conn: Connection, csv_path: Path, graph: str, out_dir: Path) -> dict:
    """Build a patient-centric clinical graph from the demo cohort CSV."""
    _ensure_graph_absent(conn, graph)

    df = pd.read_csv(csv_path)
    conn.client.create_graph(graph)
    conn.client.set_graph(graph)

    change = conn.client.new_change()
    conn.client.checkout(change=change)

    # Secondary nodes first (dedup'd by Python)
    for col, label in [
        ("condition", "MedicalCondition"),
        ("medication", "Medication"),
        ("doctor", "Doctor"),
        ("hospital", "Hospital"),
        ("blood_type", "BloodType"),
        ("gender", "Gender"),
    ]:
        for val in df[col].dropna().unique():
            conn.client.query(
                f"CREATE (:{label} {{displayName: '{_cypher_safe(str(val))}'}})"
            )
    conn.client.query("COMMIT")

    # Primary nodes + edges
    for _, row in df.iterrows():
        pid = _cypher_safe(str(row["patient_id"]))
        name = _cypher_safe(str(row["name"]))
        age = int(row["age"])
        conn.client.query(
            f"CREATE (:Patient {{id: '{pid}', displayName: '{name}', age: {age}}})"
        )
        conn.client.query("COMMIT")

        conn.client.query(f"""
            MATCH (p:Patient {{id: '{pid}'}}),
                  (c:MedicalCondition {{displayName: '{_cypher_safe(str(row['condition']))}'}}),
                  (m:Medication {{displayName: '{_cypher_safe(str(row['medication']))}'}}),
                  (d:Doctor {{displayName: '{_cypher_safe(str(row['doctor']))}'}}),
                  (h:Hospital {{displayName: '{_cypher_safe(str(row['hospital']))}'}}),
                  (b:BloodType {{displayName: '{_cypher_safe(str(row['blood_type']))}'}}),
                  (g:Gender {{displayName: '{_cypher_safe(str(row['gender']))}'}})
            CREATE (p)-[:HAS]->(c),
                   (p)-[:TOOK_MEDICATION]->(m),
                   (p)-[:IS_TREATED_BY]->(d),
                   (p)-[:IS_TREATED_IN]->(h),
                   (p)-[:IS]->(b),
                   (p)-[:IS]->(g)
        """)

    conn.client.query("CHANGE SUBMIT")
    conn.client.checkout()

    summary = {"graph": graph, "source": csv_path.name, "format": "csv (cohort demo builder)", "loader": "CREATE (cohort demo builder)"}
    summary.update(_summarise_graph(conn.client))
    md = _render_build_report(summary)
    _write_report(out_dir, md, summary)
    return summary


# ---------------------------------------------------------------------------
# Reproducibility bundle
# ---------------------------------------------------------------------------

# Mirrors the turingdb-graph extra in pyproject.toml.
REPRO_PIP_DEPS = ["turingdb>=3.1,<4", "pandas>=2.0", "tabulate>=0.9"]


def _shell_quote(value: str) -> str:
    """Quote value for bash so commands.sh hands it back unchanged.

    commands.sh is written with LF line endings, which turns a carriage return
    inside an ordinary quoted string into a newline; bash's $'...' quoting
    spells it as an escape instead.
    """
    if "\r" not in value:
        return shlex.quote(value)
    escaped = value.replace("\\", "\\\\").replace("'", "\\'").replace("\r", "\\r")
    return f"$'{escaped}'"


def _replay_args(args: argparse.Namespace, out_dir: Path) -> list[str | ReproPath]:
    """The CLI arguments that rerun this invocation into the bundle's own output directory."""
    q = _shell_quote  # free-text values (graph names, Cypher) reach commands.sh verbatim
    replay: list[str | ReproPath] = []
    if not args.demo:  # the demo runs in-process and never reads these
        if args.host != DEFAULT_HOST:
            replay += ["--host", q(args.host)]
        if args.no_auto_start:
            replay.append("--no-auto-start")
    if args.demo:
        replay += ["--demo", args.demo]
    elif args.build:
        replay += ["--build", "--input", ReproPath(args.input.resolve()), "--graph", q(args.graph),
                   "--format", args.format]
        if args.node_label:
            replay += ["--node-label", q(args.node_label)]
    elif args.query:
        replay += ["--query", "--graph", q(args.graph), "--cypher", q(args.cypher),
                   "--result-format", args.result_format]
    elif args.analyse_cohort:
        replay += ["--analyse-cohort", "--graph", q(args.graph), "--top-n", str(args.top_n)]
    replay += ["--output", ReproPath(out_dir, anchor="output_dir")]
    return replay


def write_reproducibility(args: argparse.Namespace, out_dir: Path) -> None:
    """Write reproducibility/ (commands.sh, environment.yml, checksums.sha256) under out_dir."""
    out_dir = out_dir.resolve()
    repro_dir = out_dir / "reproducibility"
    outputs = sorted(
        p for p in out_dir.rglob("*") if p.is_file() and repro_dir not in p.parents
    )
    commands = write_portable_commands_sh(
        out_dir,
        ReproCommand(
            script_path=Path("skills/turingdb-graph/turingdb_graph.py"),
            args=_replay_args(args, out_dir),
            comment="Replay this turingdb-graph run (needs the turingdb-graph extra).",
        ),
        repo_root=_PROJECT_ROOT,
    )
    environment = write_environment_yml(
        out_dir, env_name="clawbio-turingdb-graph", pip_deps=REPRO_PIP_DEPS,
    )
    write_checksums([*outputs, commands, environment], out_dir, anchor=out_dir)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="turingdb-graph",
        description="Build, query, and analyse biomedical graphs in TuringDB.",
    )
    # Global flags
    parser.add_argument("--host", default=DEFAULT_HOST, help=f"TuringDB host (default {DEFAULT_HOST})")
    parser.add_argument("--data-dir", type=Path, default=None,
                        help=f"TuringDB data directory (default {DEFAULT_DATA_DIR}); "
                             "--demo ignores it and uses a temporary one")
    parser.add_argument("--no-auto-start", action="store_true",
                        help="Do not launch the daemon if the server is unreachable")
    parser.add_argument("--output", "--out", dest="out", type=Path, default=Path("output"),
                        help="Output directory (default ./output)")

    sub = parser.add_mutually_exclusive_group(required=True)
    sub.add_argument("--build", action="store_true")
    sub.add_argument("--query", action="store_true")
    sub.add_argument("--analyse-cohort", action="store_true")
    sub.add_argument("--demo", choices=list(DEMO_SPECS))
    sub.add_argument("--stop-server", action="store_true")

    # --build args
    parser.add_argument("--input", type=Path, help="Input file for --build")
    parser.add_argument("--graph", help="Target graph name")
    parser.add_argument("--format", default="auto", choices=["auto", "csv", "tsv", "gml", "jsonl"])
    parser.add_argument("--node-label", help="Node label for CSV ingest (one node per row)")

    # --query args
    parser.add_argument("--cypher", help="Cypher query string for --query")
    parser.add_argument("--result-format", default="md", choices=["md", "json", "tsv"])

    # --analyse-cohort args
    parser.add_argument("--top-n", type=int, default=10)

    args = parser.parse_args(argv)

    if args.stop_server:
        Connection(host=args.host, data_dir=args.data_dir or DEFAULT_DATA_DIR).stop()
        return 0

    if args.demo:
        # The demo runs on an in-process engine in a throwaway data directory: no
        # daemon to start or stop, and nothing written outside --output.
        with tempfile.TemporaryDirectory(prefix="turingdb-graph-demo-") as data_dir:
            conn = Connection(host=args.host, data_dir=Path(data_dir), embedded=True)
            conn.connect()
            _run_subcommand(parser, args, conn)
    else:
        conn = Connection(host=args.host, data_dir=args.data_dir or DEFAULT_DATA_DIR)
        conn.connect(auto_start=not args.no_auto_start)
        _run_subcommand(parser, args, conn)

    write_reproducibility(args, args.out)
    print(f"✓ done — outputs: {args.out}")
    return 0


def _run_subcommand(parser: argparse.ArgumentParser, args: argparse.Namespace,
                    conn: Connection) -> None:
    if args.build:
        if not args.input or not args.graph:
            parser.error("--build requires --input and --graph")
        build_graph(
            conn, args.input, args.graph,
            fmt=args.format, node_label=args.node_label, out_dir=args.out,
        )
    elif args.query:
        if not args.graph or not args.cypher:
            parser.error("--query requires --graph and --cypher")
        run_query(conn, args.graph, args.cypher,
                  fmt=args.result_format, out_dir=args.out)
    elif args.analyse_cohort:
        if not args.graph:
            parser.error("--analyse-cohort requires --graph")
        analyse_cohort(conn, args.graph, top_n=args.top_n, out_dir=args.out)
    elif args.demo:
        skill_root = Path(__file__).resolve().parent
        run_demo(conn, args.demo, skill_root, out_dir=args.out)
    else:
        parser.error("No subcommand selected")


if __name__ == "__main__":
    sys.exit(main())
