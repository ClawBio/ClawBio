"""Emit a bash download script for FASTQ files listed in a samplesheet.

Shared by the archive-fetch skills (`ena-fetch`, `geo-fetch`,
`arrayexpress-fetch`, `pride-fetch`), whose `samplesheet` command produces the
input this consumes. Adapted from UKDRI/informatics_data_skills @ 7cc3e6e
(`fastq-download-script/scripts/make_download_script.py`).
Copyright (c) 2026 UK Dementia Research Institute. Licensed MIT.

This module only *writes* a script. Running or submitting it is a separate step
the caller takes; nothing here executes anything.

On the SLURM partition: Slurm defines no standard default partition name. The
partition marked `Default=YES` in slurm.conf is a per-site choice, so emitting a
guess such as `batch` makes sbatch fail outright with "invalid partition
specified" wherever that name does not exist. When no partition is given the
directive is emitted commented out, and the controller picks the site default.
"""

from __future__ import annotations

import csv
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


# --- transfer robustness defaults ---
# These scripts fetch multi-gigabyte archive files, usually unattended on a
# compute node. The failure that matters is not "the URL was wrong" but "the
# transfer dropped at 90% after two hours". So: retry, wait between attempts,
# refuse to hang forever, and resume rather than restart.
RETRIES = 5           # attempts per file
RETRY_DELAY = 10      # seconds between attempts
CONNECT_TIMEOUT = 30  # seconds to establish a connection
STALL_BYTES = 1024    # a transfer below this rate...
STALL_SECONDS = 120   # ...for this long is dead; abort so the retry can fire

# Plain --retry already covers 429 and 5xx, which is most archive throttling.
# It does NOT cover 403, which is what EBI returned under a burst of requests
# when this was tested (2026-09-22). --retry-all-errors closes that gap.
#
# The cost: a genuinely missing URL is now retried too, so it takes about
# RETRIES x RETRY_DELAY before failing. That is bounded, not multiplied --
# `set -e` aborts the script at the first failure, so it is ~50 s once, not
# once per file. Worth it to survive a rate limit mid-job.
#
# The flag is curl >= 7.71. Hardcoding it breaks the script outright on an
# older curl (RHEL 7 ships 7.29), so the script probes for it at run time.
CURL_PROBE = """# curl >= 7.71 retries HTTP errors (403/429/503) too, not just network
# failures. Older curl rejects the flag outright, so probe before using it.
RETRY_ALL=""
if curl --help all 2>/dev/null | grep -q -- --retry-all-errors; then
  RETRY_ALL="--retry-all-errors"
fi
"""


def _download_cmd(tool: str, url: str, outdir: str) -> str:
    name = _safe_field(url.rstrip("/").split("/")[-1])
    dest = f'"{outdir}/{name}"'
    if tool == "curl":
        # -f fail on HTTP errors, -s silent (no progress bar), -S still show
        # errors, -L follow redirects (BioStudies' file route 302s), --retry
        # with a delay for transient failures, --connect-timeout and the
        # --speed-limit/--speed-time pair so a stalled transfer aborts and
        # retries instead of hanging until the job's walltime, and -C - to
        # resume a partial file rather than restart it. -C - is safe on a
        # missing or already-complete file: both exit 0.
        return (f'curl -fsSL --retry {RETRIES} --retry-delay {RETRY_DELAY} '
                f'$RETRY_ALL --connect-timeout {CONNECT_TIMEOUT} '
                f'--speed-limit {STALL_BYTES} --speed-time {STALL_SECONDS} '
                f'-C - --create-dirs -o {dest} "{url}"')
    # wget: -q fully quiet, --tries/--waitretry for transient failures,
    # --timeout bounds both the connect and the read, -O the explicit output
    # path, and -c to resume.
    #
    # -c with -O resumes correctly -- verified against GNU Wget 1.25.0 on
    # 2026-09-22 and re-verified 2026-09-23: a truncated file produced `206
    # Partial Content` with only the remainder transferred, and the result was
    # byte-identical to a fresh download and a valid gzip. An earlier comment
    # here claimed the man page documents -c with -O as unsupported; it does
    # not. It documents that restriction for -N and for -nc, not for -c.
    #
    # NOT -nc. That is --no-clobber, not "no continue": it SKIPS a file that
    # already exists, which would strand a partial download truncated forever
    # -- the opposite of what -c is for.
    #
    # -c only resumes transfers left by a *previous* invocation; mid-transfer
    # retry within one run is already wget's default. That is exactly the
    # resubmitted-job case these scripts are written for.
    return (f'wget -q --tries={RETRIES} --waitretry={RETRY_DELAY} '
            f'--timeout={CONNECT_TIMEOUT} -c -O {dest} "{url}"')


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
    tool: str = "curl",
    outdir: str = "fastq",
    slurm: SlurmOptions | None = None,
) -> tuple[str, int]:
    """Render the script body. `slurm=None` omits the SLURM header entirely."""
    parts = [_slurm_header(slurm), ""] if slurm is not None else ["#!/bin/bash"]
    parts += ["set -euo pipefail", "", f'OUTDIR="{outdir}"', 'mkdir -p "$OUTDIR"', ""]
    if tool == "curl":
        parts += [CURL_PROBE]
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
    tool: str = "curl",
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
