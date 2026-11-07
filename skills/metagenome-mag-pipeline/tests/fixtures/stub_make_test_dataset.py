#!/usr/bin/env python3
"""A stand-in for upstream `test/make_test_dataset.py` with the same CLI.

The real generator draws fragments from simulated genomes; this one writes a few
records of fixed sequence so the tests exercise samplesheet discovery and config
rendering without a genome simulator in the way. It accepts the same flags and
writes the same artifacts:

    <outdir>/<SAMPLE>/<SAMPLE>_1.fastq.gz
    <outdir>/<SAMPLE>/<SAMPLE>_2.fastq.gz
    <outdir>/samplesheet.tsv          columns: sample fastq_1 fastq_2 group assembly
    <dbdir>/decontamination/test_host.fa
"""

import argparse
import gzip
from pathlib import Path

SEQUENCE = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
QUALITY = "I" * len(SEQUENCE)


def write_fastq(path: Path, records: int, mate: int, sample: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        for index in range(records):
            handle.write(f"@{sample}:{index + 1}/{mate}\n{SEQUENCE}\n+\n{QUALITY}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--dbdir", default=None)
    parser.add_argument("--samples", type=int, default=2)
    parser.add_argument("--genomes", type=int, default=4)
    parser.add_argument("--genome-length", type=int, default=12000)
    parser.add_argument("--pairs", type=int, default=2000)
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--insert", type=int, default=400)
    parser.add_argument("--host-fraction", type=float, default=0.05)
    parser.add_argument("--error-rate", type=float, default=0.002)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    dbdir = Path(args.dbdir).resolve() if args.dbdir else outdir / "dbs"

    rows = []
    for index in range(args.samples):
        sample = f"TEST_{index + 1:03d}"
        r1 = outdir / sample / f"{sample}_1.fastq.gz"
        r2 = outdir / sample / f"{sample}_2.fastq.gz"
        write_fastq(r1, min(args.pairs, 8), 1, sample)
        write_fastq(r2, min(args.pairs, 8), 2, sample)
        rows.append((sample, str(r1), str(r2), sample, "individual"))

    sheet = outdir / "samplesheet.tsv"
    with sheet.open("w") as handle:
        handle.write("sample\tfastq_1\tfastq_2\tgroup\tassembly\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")

    host = dbdir / "decontamination"
    host.mkdir(parents=True, exist_ok=True)
    (host / "test_host.fa").write_text(">test_host_contig_1\n" + SEQUENCE + "\n")

    print(f"Wrote samplesheet: {sheet}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
