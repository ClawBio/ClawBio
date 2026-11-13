#!/usr/bin/env python3
"""Rebuild the deterministic, fully simulated methylation demo fixture.

The bundled probe identifiers come from a public Illumina 450K-style example.
No donor methylation measurements are read by this generator.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
from pathlib import Path


DATA_DIR = Path(__file__).resolve().parent
SEED = "20261019"
SAMPLES = (
    ("synthetic_01", "blood", 35, "F"),
    ("synthetic_02", "blood", 65, "M"),
)


def generate(output: Path) -> None:
    with gzip.open(DATA_DIR / "probe_ids.txt.gz", "rt", encoding="utf-8") as handle:
        probe_ids = [line.strip() for line in handle if line.strip()]
    if not probe_ids or len(probe_ids) != len(set(probe_ids)):
        raise ValueError("Probe identifier list must be nonempty and unique")

    # Integer-only transforms of hashes give a stable bimodal beta matrix
    # across Python and NumPy versions. They carry no biological signal.
    values_by_sample: list[list[str]] = [[], []]
    scale = 1 << 32
    for probe_id in probe_ids:
        digest = hashlib.sha256(f"{SEED}:{probe_id}".encode("ascii")).digest()
        u = int.from_bytes(digest[0:4], "big")
        spread = u * u * u * 350_000 // (scale * scale * scale)
        baseline = spread if (digest[12] & 1) == 0 else 1_000_000 - spread
        for sample_number in range(2):
            start = 4 + 4 * sample_number
            v = int.from_bytes(digest[start:start + 4], "big")
            noise = v * 80_001 // scale - 40_000
            value = max(0, min(1_000_000, baseline + noise))
            values_by_sample[sample_number].append(f"{value // 1_000_000}.{value % 1_000_000:06d}")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.writer(text, lineterminator="\n")
                writer.writerow(["sample_id", "tissue_type", "age", "gender", *probe_ids])
                for sample, values in zip(SAMPLES, values_by_sample):
                    writer.writerow([*sample, *values])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    generate(args.output)


if __name__ == "__main__":
    main()
