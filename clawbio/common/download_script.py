"""Emit a bash download script for FASTQ files listed in a samplesheet.

Shared by the archive-fetch skills (`ena-fetch`, `geo-fetch`,
`arrayexpress-fetch`, `pride-fetch`), whose `samplesheet` command produces the
input this consumes. Adapted from UKDRI/informatics_data_skills @ 7cc3e6e
(`fastq-download-script/scripts/make_download_script.py`).
Copyright (c) 2026 UK Dementia Research Institute. Licensed MIT.

This module only *writes* a script. Running or submitting it is a separate,
explicitly requested step — see `clawbio.common.job_exec`.

On the SLURM partition: Slurm defines no standard default partition name. The
partition marked `Default=YES` in slurm.conf is a per-site choice, so emitting a
guess such as `batch` makes sbatch fail outright with "invalid partition
specified" wherever that name does not exist. When no partition is given the
directive is emitted commented out, and the controller picks the site default.
"""

from __future__ import annotations

import csv
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path

# Clean output fields: control chars (incl. CR/LF/tab) -> _, so no value can
# break the line it is written into.
_CTRL_RE = re.compile(r"[\x00-\x1f\x7f-\x9f]")

Groups = list[tuple[str, list[str]]]


def _safe_field(value: str | None) -> str:
    """Value safe for one line: control chars -> _, runs of spaces collapsed."""
    return re.sub(r" +", " ", _CTRL_RE.sub("_", value or "")).strip()


def _is_safe_url(url: str) -> bool:
    """Reject anything that would break out of the double quotes in bash."""
    return not (_CTRL_RE.search(url) or '"' in url)


@dataclass
class SlurmOptions:
    """SLURM header settings. Every unset value is emitted commented out."""

    job_name: str = "fastq_download"
    partition: str | None = None
    account: str | None = None
    cpus: str = "1"
    mem: str = "4G"
    time: str = "24:00:00"
    email: str | None = None


def urls_from_samplesheet(path: Path | str) -> Groups:
    """Return [(sample, [urls])] in row order, from the fastq_* columns."""
    path = Path(path)
    out: Groups = []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        cols = reader.fieldnames or []
        fastq_cols = [c for c in cols if c and c.lower().startswith("fastq")]
        if not fastq_cols:
            raise SystemExit(
                f"No 'fastq_*' columns found in {path}. Is this a samplesheet? "
                "For a plain URL list, pass a URL list instead.")
        sample_col = next((c for c in cols if c and c.lower() == "sample"), None)
        for i, row in enumerate(reader):
            sample = (row.get(sample_col) if sample_col else None) or f"row{i + 1}"
            urls = []
            for col in fastq_cols:
                value = (row.get(col) or "").strip()
                if not value:
                    continue
                if not (value.startswith("http") or value.startswith("ftp")):
                    print(f"warning: skipping non-URL value in {col}: {value}",
                          file=sys.stderr)
                    continue
                if not _is_safe_url(value):
                    print(f"warning: skipping unsafe URL in {col}: {value!r}",
                          file=sys.stderr)
                    continue
                urls.append(value)
            if urls:
                out.append((sample, urls))
    return out


def urls_from_url_list(path: Path | str) -> Groups:
    """Return [(label, [url])] from a plain one-URL-per-line file."""
    out: Groups = []
    for i, line in enumerate(Path(path).read_text().splitlines()):
        url = line.strip()
        if not url or url.startswith("#"):
            continue
        if not (url.startswith("http") or url.startswith("ftp")):
            continue
        if not _is_safe_url(url):
            print(f"warning: skipping unsafe URL: {url!r}", file=sys.stderr)
            continue
        out.append((f"file{i + 1}", [url]))
    return out


def _download_cmd(tool: str, url: str, outdir: str) -> str:
    name = _safe_field(url.rstrip("/").split("/")[-1])
    dest = f'"{outdir}/{name}"'
    if tool == "curl":
        # -f fail on HTTP errors, -s silent (no progress bar), -S still show
        # errors, -L follow redirects, --retry for transient failures.
        return f'curl -fsSL --retry 3 --create-dirs -o {dest} "{url}"'
    # wget: -q fully quiet, retries with a wait, -O explicit output path.
    return f'wget -q --tries=3 --waitretry=5 -O {dest} "{url}"'


def _slurm_header(opts: SlurmOptions) -> str:
    lines = ["#!/bin/bash", "#", f"#SBATCH --job-name={opts.job_name}"]
    if opts.partition:
        lines.append(f"#SBATCH --partition={opts.partition}")
    else:
        lines.append("# #SBATCH --partition=<your_partition>   # set for your cluster")
    if opts.account:
        lines.append(f"#SBATCH --account={opts.account}")
    else:
        lines.append("# #SBATCH --account=<your_account>        # set if required")
    lines += [
        f"#SBATCH --cpus-per-task={opts.cpus}",
        f"#SBATCH --mem={opts.mem}",
        f"#SBATCH --time={opts.time}",
        "#SBATCH --output=slurm-%j.out",
        "#SBATCH --error=slurm-%j.err",
    ]
    if opts.email:
        lines += ["#SBATCH --mail-type=END,FAIL", f"#SBATCH --mail-user={opts.email}"]
    else:
        lines += ["# #SBATCH --mail-type=END,FAIL",
                  "# #SBATCH --mail-user=<you@example.org>"]
    return "\n".join(lines)


def build_download_script(
    groups: Groups,
    *,
    tool: str = "wget",
    outdir: str = "fastq",
    slurm: SlurmOptions | None = None,
) -> tuple[str, int]:
    """Render the script body. `slurm=None` omits the SLURM header entirely."""
    parts = [_slurm_header(slurm), ""] if slurm is not None else ["#!/bin/bash"]
    parts += ["set -euo pipefail", "", f'OUTDIR="{outdir}"', 'mkdir -p "$OUTDIR"', ""]
    n_files = 0
    for sample, urls in groups:
        parts.append(f"# sample: {_safe_field(sample)}")
        for url in urls:
            parts.append(_download_cmd(tool, url, "$OUTDIR"))
            n_files += 1
        parts.append("")
    parts += [f'echo "Downloaded {n_files} file(s) to $OUTDIR"', ""]
    return "\n".join(parts), n_files


def write_download_script(
    groups: Groups,
    out_path: Path | str,
    *,
    tool: str = "wget",
    outdir: str = "fastq",
    slurm: SlurmOptions | None = None,
) -> tuple[Path, int]:
    """Write the script and mark it executable. Returns (path, file count)."""
    if not groups:
        raise SystemExit("No FASTQ URLs found in the input.")
    body, n_files = build_download_script(groups, tool=tool, outdir=outdir, slurm=slurm)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(body, encoding="utf-8")
    mode = out_path.stat().st_mode
    out_path.chmod(mode | ((mode & 0o444) >> 2))  # mirror read bits into execute
    return out_path, n_files
