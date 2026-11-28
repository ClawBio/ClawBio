#!/usr/bin/env python3
"""
polars_bio_runner.py: Genomic interval operations, multi-format I/O, DataFusion
SQL, and pileup on Polars DataFrames via polars-bio.

Wraps polars-bio (https://github.com/biodatageeks/polars-bio).
Requires: uv pip install polars-bio   (Python 3.12-3.14)

Citation: Wiewiorka M, Khamutou P, Zbysinski M, Gambin T. polars-bio - fast,
scalable, and out-of-core operations on large genomic interval datasets.
Bioinformatics, 2025, 41(12):btaf640. https://doi.org/10.1093/bioinformatics/btaf640

Subcommands:
  overlap nearest merge coverage cluster complement subtract count-overlaps
                - genomic interval operations
  io            - read/scan a format (--describe for schema-only inspection)
  sql           - register file(s) as tables and run DataFusion SQL (--query)
  pileup        - per-base read depth (mosdepth-compatible) from BAM/CRAM

Usage:
  python polars_bio_runner.py overlap --a a.bed --b b.bed --output /tmp/out
  python polars_bio_runner.py io --input s.vcf --format vcf --describe --output /tmp/out
  python polars_bio_runner.py sql --input s.vcf --query "SELECT * FROM t" --output /tmp/out
  python polars_bio_runner.py --demo --output /tmp/polars_bio_demo
"""
from __future__ import annotations

__version__ = "0.1.0"
__author__ = "ClawBio (adapted from K-Dense scientific-agent-skills/polars-bio)"
__license__ = "Apache-2.0"

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

# No hardcoded polars-bio version anywhere: the version that actually ran is
# detected at runtime via polars_bio.__version__ (see _detect_pb_version) and
# recorded in result.json/report.md. Dependencies are declared in SKILL.md.

DISCLAIMER = (
    "ClawBio is a research and educational tool. "
    "It is not a medical device and does not provide clinical diagnoses. "
    "Consult a healthcare professional before making any medical decisions."
)

CITATION = (
    "Wiewiorka M, Khamutou P, Zbysinski M, Gambin T. polars-bio - fast, scalable, "
    "and out-of-core operations on large genomic interval datasets. "
    "Bioinformatics, 2025, 41(12):btaf640. https://doi.org/10.1093/bioinformatics/btaf640"
)

# expand / sort_bedframe are documented conceptually upstream but are NOT exposed
# as polars_bio.* functions in the current API, so they are not subcommands.
INTERVAL_OPS = {
    "overlap", "nearest", "merge", "coverage", "cluster",
    "complement", "subtract", "count-overlaps",
}
ALL_SUBCOMMANDS = INTERVAL_OPS | {"io", "sql", "pileup"}


def _import_pb():
    """Lazily import polars-bio; emit install hint and exit(2) if unavailable."""
    if os.environ.get("POLARS_BIO_FORCE_IMPORT_ERROR") == "1":
        _dependency_error()
    try:
        import polars_bio as pb  # noqa: F401
        return pb
    except ImportError:
        _dependency_error()


def _dependency_error():
    sys.stderr.write(
        "ERROR: polars-bio is not installed.\n"
        "Install it with:  uv pip install polars-bio\n"
        "Requires Python 3.12-3.14.\n"
    )
    sys.exit(2)


def _detect_pb_version() -> str:
    """Actual installed polars-bio version (provenance), or 'unknown'."""
    try:
        import polars_bio
        return getattr(polars_bio, "__version__", "unknown")
    except ImportError:
        return "unknown"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="polars_bio_runner.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="polars-bio: genomic interval ops, I/O, SQL, pileup.",
    )
    p.add_argument("subcommand", nargs="?", choices=sorted(ALL_SUBCOMMANDS),
                   help="Operation to run (omit with --demo).")
    p.add_argument("--a", help="Primary interval file (BED, >=4 columns).")
    p.add_argument("--b", help="Secondary interval file (BED), for two-input ops.")
    p.add_argument("--input", help="Input file for io/sql/pileup.")
    p.add_argument("--format", help="Format for io (bed,vcf,vcf_zarr,gff,gtf,fasta,fastq,bam,cram,sam,pairs,bigwig,bigbed).")
    p.add_argument("--reference", help="Reference FASTA path (required for CRAM).")
    p.add_argument("--describe", action="store_true", help="io: schema-only inspection (no full read).")
    p.add_argument("--query", help="DataFusion SQL query for the sql subcommand.")
    p.add_argument("--k", type=int, default=1, help="nearest: number of nearest intervals (default 1).")
    p.add_argument("--min-mapping-quality", dest="min_mapping_quality", type=int, default=0,
                   help="pileup: minimum mapping quality (default 0).")
    p.add_argument("--one-based", action="store_true",
                   help="Output 1-based closed coordinates. Default is 0-based half-open "
                        "(matches BED on disk; controls output representation, not input parsing).")
    p.add_argument("--genome",
                   help="complement: contig-sizes file (chrom<TAB>size per line) used to bound gaps.")
    p.add_argument("--output", help="Output directory.")
    p.add_argument("--demo", action="store_true", help="Run overlap on bundled synthetic BED data.")
    p.add_argument("--verbose", action="store_true", help="Verbose output.")
    return p


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def build_result(subcommand, params, inputs, output_rows, output_schema, figure,
                 polars_bio_version=None, caveats=None):
    return {
        "skill": "polars-bio",
        "subcommand": subcommand,
        "params": params,
        "polars_bio_version": polars_bio_version or _detect_pb_version(),
        "inputs": list(inputs),
        "output_rows": output_rows,
        "output_schema": output_schema,
        "figure": figure,
        "caveats": list(caveats or []),
        "report": "report.md",
        "timestamp": _now_iso(),
    }


def write_result_json(result: dict, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "result.json").write_text(json.dumps(result, indent=2))


def write_report(result: dict, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    schema = result.get("output_schema") or {}
    schema_rows = "\n".join(f"| `{c}` | {t} |" for c, t in schema.items()) or "| (n/a) | |"
    fig = f"\n![figure]({result['figure']})\n" if result.get("figure") else ""
    params_json = json.dumps(result["params"], indent=2)
    inputs = ", ".join(result["inputs"]) or "(none)"
    caveats = result.get("caveats") or []
    caveat_block = (
        "## Caveats\n\n" + "".join(f"- {c}\n" for c in caveats) + "\n" if caveats else ""
    )
    report = (
        f"# polars-bio — {result['subcommand']}\n\n"
        f"**polars-bio version:** {result['polars_bio_version']}\n"
        f"**Inputs:** {inputs}\n"
        f"**Output rows:** {result['output_rows']}\n"
        f"**Generated:** {result['timestamp']}\n\n"
        f"{caveat_block}"
        f"## Parameters\n\n"
        f"```json\n{params_json}\n```\n\n"
        f"## Output schema\n\n"
        f"| Column | Type |\n|--------|------|\n{schema_rows}\n"
        f"{fig}\n"
        f"## Reproducibility\n\n"
        f"See `reproducibility/commands.sh`.\n\n"
        f"## Citation\n\n{CITATION}\n\n"
        f"---\n\n_{DISCLAIMER}_\n"
    )
    (out_dir / "report.md").write_text(report)


def write_figure(subcommand, df, out_dir: Path):
    """Render an interval/coverage figure. Returns filename or None (io/sql/empty)."""
    if subcommand in {"io", "sql"} or df is None or df.height == 0:
        return None
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    cols = df.columns
    start_col = next((c for c in cols if "start" in c), None)
    end_col = next((c for c in cols if "end" in c), None)
    if not start_col or not end_col:
        return None
    rows = df.head(40)
    fig, ax = plt.subplots(figsize=(8, max(2, rows.height * 0.25)))
    for i, (s, e) in enumerate(zip(rows[start_col], rows[end_col])):
        ax.plot([s, e], [i, i], lw=4)
    ax.set_xlabel("position")
    ax.set_ylabel("interval")
    ax.set_title(f"polars-bio {subcommand} (first {rows.height} rows)")
    fig.tight_layout()
    fig.savefig(out_dir / "figure.png", dpi=100)
    plt.close(fig)
    return "figure.png"


def write_reproducibility(args, out_dir: Path) -> None:
    repro = out_dir / "reproducibility"
    repro.mkdir(parents=True, exist_ok=True)
    argv = " ".join(sys.argv[1:])
    (repro / "commands.sh").write_text(
        "#!/usr/bin/env bash\n"
        f"# polars-bio runner v{__version__}, polars-bio {_detect_pb_version()}\n"
        f"python polars_bio_runner.py {argv}\n"
    )


def main() -> None:
    args = build_parser().parse_args()
    if not args.output:
        sys.stderr.write("ERROR: --output is required.\n")
        sys.exit(1)
    if not args.demo and not args.subcommand:
        sys.stderr.write("ERROR: a subcommand is required (or use --demo).\n")
        sys.exit(1)
    raise SystemExit(_dispatch(args))


DEMO_A = Path(__file__).parent / "examples" / "demo_a.bed"
DEMO_B = Path(__file__).parent / "examples" / "demo_b.bed"

TWO_INPUT_OPS = {"overlap", "nearest", "coverage", "subtract", "count-overlaps"}
SINGLE_INPUT_OPS = {"merge", "cluster", "complement"}
_DF = {"output_type": "polars.DataFrame"}


def _zero_based(args) -> bool:
    """Coordinate output representation: 0-based half-open by default (BED-native),
    1-based closed when --one-based is passed. Controls output display only; both
    inputs share the setting so interval logic is identical across modes."""
    return not args.one_based


def read_bed_lf(pb, path, zero_based):
    return pb.scan_bed(str(path), use_zero_based=bool(zero_based))


def _load_genome_view(path):
    """Build a contig-bounds view_df (chrom,start=0,end=size) from a chrom-sizes
    file (chrom<TAB>size per line; '#' comments and blank lines ignored)."""
    import polars as pl
    chroms, ends = [], []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t") if "\t" in line else line.split()
        if len(parts) < 2:
            raise ValueError(f"Bad chrom-sizes line (need chrom<TAB>size): {line!r}")
        chroms.append(parts[0])
        ends.append(int(parts[1]))
    return pl.DataFrame({"chrom": chroms, "start": [0] * len(chroms), "end": ends})


def run_interval_op(pb, op, args):
    zb = _zero_based(args)
    a = read_bed_lf(pb, args.a, zb)
    if op in TWO_INPUT_OPS:
        b = read_bed_lf(pb, args.b, zb)
        if op == "overlap":
            return pb.overlap(a, b, **_DF)
        if op == "nearest":
            return pb.nearest(a, b, k=args.k, **_DF)
        if op == "coverage":
            return pb.coverage(a, b, **_DF)
        if op == "subtract":
            return pb.subtract(a, b, **_DF)
        if op == "count-overlaps":
            return pb.count_overlaps(a, b, **_DF)
    if op == "merge":
        return pb.merge(a, **_DF)
    if op == "cluster":
        return pb.cluster(a, **_DF)
    if op == "complement":
        view = _load_genome_view(args.genome) if args.genome else None
        return pb.complement(a, view_df=view, **_DF)
    raise ValueError(f"Unsupported interval op: {op}")


def _schema_of(df):
    return {c: str(t) for c, t in zip(df.columns, df.dtypes)}


def _write_result_table(df, out_dir: Path) -> str:
    """Write the result table, falling back to NDJSON for nested columns."""
    try:
        df.write_csv(out_dir / "result.csv")
        return "result.csv"
    except Exception:
        df.write_ndjson(out_dir / "result.ndjson")
        return "result.ndjson"


_COMPLEMENT_UNBOUNDED = (
    "complement was run without contig bounds: trailing gaps are unbounded "
    "(span to i64::MAX) and are not genomically meaningful. Pass "
    "--genome <chrom-sizes> (chrom<TAB>size per line) to bound the gaps."
)


def _run_interval_and_write(pb, op, args, out_dir, inputs, params):
    df = run_interval_op(pb, op, args)
    caveats = []
    if op == "complement" and not args.genome:
        caveats.append(_COMPLEMENT_UNBOUNDED)
        sys.stderr.write("WARNING: " + _COMPLEMENT_UNBOUNDED + "\n")
    figure = write_figure(op, df, out_dir)
    result = build_result(op, params, inputs, df.height, _schema_of(df), figure,
                          caveats=caveats)
    write_result_json(result, out_dir)
    write_report(result, out_dir)
    write_reproducibility(args, out_dir)
    _write_result_table(df, out_dir)
    return 0


_SCAN = {
    "bed": "scan_bed", "vcf": "scan_vcf", "vcf_zarr": "scan_vcf_zarr",
    "gff": "scan_gff", "gtf": "scan_gtf", "fasta": "scan_fasta",
    "fastq": "scan_fastq", "bam": "scan_bam", "cram": "scan_cram",
    "sam": "scan_sam", "pairs": "scan_pairs",
    "bigwig": "scan_bigwig", "bigbed": "scan_bigbed",
}
_DESCRIBE = {
    "vcf": "describe_vcf", "vcf_zarr": "describe_vcf_zarr", "bam": "describe_bam",
    "cram": "describe_cram", "sam": "describe_sam",
}


def run_io(pb, args):
    """Return (df_or_None, schema_dict). df is None when --describe (schema-only)."""
    fmt = (args.format or "").lower()
    if args.describe:
        if fmt not in _DESCRIBE:
            raise ValueError(f"--describe not supported for format {fmt!r}")
        desc = getattr(pb, _DESCRIBE[fmt])(args.input)
        ddf = desc.collect() if hasattr(desc, "collect") else desc
        return None, _schema_of(ddf)
    if fmt not in _SCAN:
        raise ValueError(f"Unsupported io format: {fmt!r}")
    scan = getattr(pb, _SCAN[fmt])
    if fmt == "cram":
        if not args.reference:
            raise ValueError("CRAM requires --reference (path to the reference FASTA).")
        lf = scan(args.input, reference_path=args.reference)
    elif fmt == "bed":
        # Honor the same coordinate default as the interval-op path.
        lf = scan(args.input, use_zero_based=_zero_based(args))
    else:
        lf = scan(args.input)
    df = lf.head(50).collect()
    return df, _schema_of(df)


_REGISTER = {
    ".bed": "register_bed", ".vcf": "register_vcf", ".gff": "register_gff",
    ".gff3": "register_gff", ".gtf": "register_gtf", ".fasta": "register_fasta",
    ".fa": "register_fasta", ".fastq": "register_fastq", ".fq": "register_fastq",
    ".bam": "register_bam", ".sam": "register_sam", ".pairs": "register_pairs",
    ".bw": "register_bigwig", ".bigwig": "register_bigwig",
    ".bb": "register_bigbed", ".bigbed": "register_bigbed",
}


def _register_fn(pb, path):
    name = Path(path).name.lower()
    for ext, fn in _REGISTER.items():
        if name.endswith(ext):
            return getattr(pb, fn)
    raise ValueError(f"Cannot infer register function for {path!r}")


def run_sql(pb, args):
    # register_bed has no use_zero_based param, so reconcile the coordinate system
    # via the global option before registration (keeps sql consistent with io/interval).
    try:
        pb.set_option("datafusion.bio.coordinate_system_zero_based",
                      "true" if _zero_based(args) else "false")
    except Exception:
        pass
    _register_fn(pb, args.input)(args.input, name="t")
    res = pb.sql(args.query)
    return res.collect() if hasattr(res, "collect") else res


def run_pileup(pb, args):
    bam = Path(args.input)
    if bam.suffix.lower() == ".bam" and not (
        bam.with_suffix(".bam.bai").exists()
        or bam.with_suffix(".bai").exists()
        or Path(str(bam) + ".bai").exists()
    ):
        raise FileNotFoundError(
            f"BAM index not found for {bam}. Create one with: samtools index {bam}"
        )
    res = pb.depth(str(bam), min_mapping_quality=args.min_mapping_quality)
    return res.collect() if hasattr(res, "collect") else res


def _dispatch(args) -> int:
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    pb = _import_pb()

    if args.demo:
        args.subcommand, args.a, args.b = "overlap", str(DEMO_A), str(DEMO_B)

    op = args.subcommand
    if op in INTERVAL_OPS:
        if not args.a:
            sys.stderr.write(f"ERROR: '{op}' needs --a.\n")
            return 1
        if op in TWO_INPUT_OPS and not args.b:
            sys.stderr.write(f"ERROR: '{op}' needs both --a and --b.\n")
            return 1
        params = {"k": args.k, "zero_based": _zero_based(args)}
        if op == "complement":
            params["genome"] = args.genome
        inputs = [args.a] + ([args.b] if args.b else [])
        return _run_interval_and_write(pb, op, args, out_dir, inputs, params)

    if op == "io":
        if not args.input or not args.format:
            sys.stderr.write("ERROR: 'io' needs --input and --format.\n")
            return 1
        df, schema = run_io(pb, args)
        rows = df.height if df is not None else None
        params = {"format": args.format, "describe": bool(args.describe)}
        result = build_result("io", params, [args.input], rows, schema, None)
        write_result_json(result, out_dir)
        write_report(result, out_dir)
        write_reproducibility(args, out_dir)
        if df is not None:
            _write_result_table(df, out_dir)
        return 0

    if op == "sql":
        if not args.input or not args.query:
            sys.stderr.write("ERROR: 'sql' needs --input and --query.\n")
            return 1
        df = run_sql(pb, args)
        params = {"query": args.query}
        result = build_result("sql", params, [args.input], df.height, _schema_of(df), None)
        write_result_json(result, out_dir)
        write_report(result, out_dir)
        write_reproducibility(args, out_dir)
        _write_result_table(df, out_dir)
        return 0

    if op == "pileup":
        if not args.input:
            sys.stderr.write("ERROR: 'pileup' needs --input (BAM/CRAM).\n")
            return 1
        try:
            df = run_pileup(pb, args)
        except FileNotFoundError as exc:
            sys.stderr.write(f"ERROR: {exc}\n")
            return 1
        figure = write_figure("coverage", df, out_dir)
        params = {"min_mapping_quality": args.min_mapping_quality}
        result = build_result("pileup", params, [args.input], df.height, _schema_of(df), figure)
        write_result_json(result, out_dir)
        write_report(result, out_dir)
        write_reproducibility(args, out_dir)
        _write_result_table(df, out_dir)
        return 0

    sys.stderr.write(f"ERROR: subcommand '{op}' not recognized.\n")
    return 1


if __name__ == "__main__":
    main()
