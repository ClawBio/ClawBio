#!/usr/bin/env python3
"""Look at the machine you are on and say how to run BJ-Expression here.

Read-only. It runs a few version commands and reads /proc, and it changes nothing.
It reports what it found, picks a container profile, and prints the next commands with
YOUR paths. Nothing is assumed about the host: a laptop with Docker, a workstation with
Apptainer and a locked-down container each get a different answer.

Profiles it can recommend:
  docker             Docker is installed and its daemon answers.
  singularity-sif    Apptainer/Singularity is present and the usual .sif images should work.
  singularity-sandbox
                     Apptainer/Singularity is present, but this host looks like one where
                     .sif images cannot be mounted (no FUSE, a gVisor kernel, or an
                     overlay setting turned off). Pre-built sandbox directories are used.
  none               No container runtime was found.

The sandbox recommendation is a judgement from the signals listed in the report. It is
confirmed or ruled out by --probe-image, which runs one tiny command inside an image you name.

Usage:
  bj_expression_wrapper.py --detect-environment [--output DIR] [--probe-image PATH]
  python environment_check.py --path ~/runs --json-out env.json
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

# Free space worth having before a run. These figures come from the author's full-read runs
# (about 100 GB of kept results per 20 pairs, plus working space) and are a starting point.
# Change them with --min-free-gb. Subsampled runs need far less.
DEFAULT_MIN_FREE_GB = 300

# Filesystems that commonly cause trouble for Nextflow work directories or container mounts.
SLOW_OR_LIMITED_FS = {
    "9p": "a 9p mount (common in gVisor sandboxes): slow metadata, overlay mounts may fail",
    "fuse": "a FUSE mount: slow for many small files",
    "fuse.s3fs": "an s3fs mount: object storage, unsuitable for a Nextflow work directory",
    "nfs": "an NFS mount: file timestamps and locking can differ from local disk",
    "nfs4": "an NFS mount: file timestamps and locking can differ from local disk",
    "cifs": "a CIFS mount: file timestamps and locking can differ from local disk",
}


@dataclass
class Finding:
    """One thing learned about the host."""
    name: str
    value: str
    status: str = "ok"          # ok | warn | problem | info
    note: str = ""


@dataclass
class Report:
    findings: list[Finding] = field(default_factory=list)
    profile: str = "none"
    reasons: list[str] = field(default_factory=list)
    next_steps: list[str] = field(default_factory=list)

    def add(self, name, value, status="ok", note=""):
        self.findings.append(Finding(name, str(value), status, note))

    def worst(self) -> str:
        order = ["ok", "info", "warn", "problem"]
        return max((f.status for f in self.findings), key=order.index, default="ok")


# ──────────────────────────────────────────────
# Probes. Each takes what it needs, so tests can feed fixed inputs.
# ──────────────────────────────────────────────
def run_version(cmd: list[str], pattern: str, timeout: int = 30) -> str:
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired):
        return ""
    m = re.search(pattern, p.stdout + p.stderr)
    return m.group(1) if m else ""


def read_text(path: str) -> str:
    try:
        return Path(path).read_text(errors="replace")
    except OSError:
        return ""


def is_gvisor(kernel_release: str, version_text: str) -> bool:
    return "gvisor" in kernel_release.lower() or "gvisor" in version_text.lower()


def mem_total_gb(meminfo: str) -> float | None:
    m = re.search(r"^MemTotal:\s+(\d+)\s+kB", meminfo, re.M)
    return round(int(m.group(1)) / 1024 / 1024, 1) if m else None


def filesystem_of(path: Path, mounts_text: str) -> str:
    """Filesystem type of the mount that holds `path` (longest matching mount point)."""
    best, best_len = "unknown", -1
    target = str(path)
    for line in mounts_text.splitlines():
        parts = line.split()
        if len(parts) < 3:
            continue
        mount_point, fstype = parts[1].replace("\\040", " "), parts[2]
        if (target == mount_point or target.startswith(mount_point.rstrip("/") + "/")
                or mount_point == "/") and len(mount_point) > best_len:
            best, best_len = fstype, len(mount_point)
    return best


def apptainer_overlay_setting(conf_text: str) -> str | None:
    """Value of 'enable overlay' in apptainer.conf (yes, no, try, driver), if set."""
    m = re.search(r"^\s*enable overlay\s*=\s*(\S+)", conf_text, re.M)
    return m.group(1) if m else None


def find_apptainer_conf() -> str:
    exe = shutil.which("apptainer") or shutil.which("singularity")
    if not exe:
        return ""
    root = Path(exe).resolve().parent.parent
    for rel in ("etc/apptainer/apptainer.conf", "etc/singularity/singularity.conf"):
        for base in (root, Path("/")):
            p = base / rel
            if p.is_file():
                return read_text(str(p))
    return ""


def probe_image(runtime: str, image: Path, timeout: int = 120) -> tuple[bool, str]:
    """Run `true` inside an image or sandbox directory. Returns (worked, message)."""
    cmd = [runtime, "exec"]
    if image.is_dir():
        cmd.append("--writable")
    cmd += [str(image), "true"]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, f"could not run: {exc}"
    if p.returncode == 0:
        return True, "ran a command inside the image"
    return False, (p.stderr.strip().splitlines() or ["exit " + str(p.returncode)])[-1]


# ──────────────────────────────────────────────
# Detection
# ──────────────────────────────────────────────
def detect(path: Path, min_free_gb: int = DEFAULT_MIN_FREE_GB,
           image_to_probe: Path | None = None, probes: dict | None = None) -> Report:
    """Collect findings and choose a profile. `probes` overrides the readers in tests."""
    pr = {
        "kernel": platform.release(),
        "version_text": read_text("/proc/version"),
        "meminfo": read_text("/proc/meminfo"),
        "mounts": read_text("/proc/mounts"),
        "cpus": os.cpu_count() or 0,
        "which": shutil.which,
        "run_version": run_version,
        "apptainer_conf": None,       # lazily read
        "dev_fuse": os.path.exists("/dev/fuse"),
        "disk_free_gb": None,
        "docker_ok": None,
        "probe_image": probe_image,
    }
    pr.update(probes or {})
    which = pr["which"]
    rep = Report()

    # Host
    gv = is_gvisor(pr["kernel"], pr["version_text"])
    rep.add("operating system", f"{platform.system()} {pr['kernel']}")
    rep.add("sandboxed kernel (gVisor)", "yes" if gv else "no",
            "warn" if gv else "ok",
            "gVisor does not support FUSE or overlay mounts, so .sif images usually fail"
            if gv else "")
    rep.add("CPUs", pr["cpus"], "warn" if pr["cpus"] < 8 else "ok",
            "fewer than 8 CPUs will make STAR and HTSeq slow" if pr["cpus"] < 8 else "")
    mem = mem_total_gb(pr["meminfo"])
    if mem is None:
        rep.add("memory", "unknown", "info")
    else:
        rep.add("memory", f"{mem} GB", "warn" if mem < 64 else "ok",
                "STAR needs about 32 GB for the human index; 64 GB or more is comfortable"
                if mem < 64 else "")

    # Storage
    target = path.expanduser().resolve()
    probe_at = target
    while not probe_at.exists() and probe_at != probe_at.parent:
        probe_at = probe_at.parent
    free = pr["disk_free_gb"]
    if free is None:
        free = round(shutil.disk_usage(probe_at).free / 1e9, 1)
    rep.add("free disk at output", f"{free} GB ({probe_at})",
            "warn" if free < min_free_gb else "ok",
            f"under {min_free_gb} GB; full-read runs keep about 100 GB of results per 20 pairs, "
            "so lower this only for subsampled runs" if free < min_free_gb else "")
    fs = filesystem_of(probe_at, pr["mounts"])
    note = SLOW_OR_LIMITED_FS.get(fs, "")
    rep.add("filesystem at output", fs, "warn" if note else "ok", note)
    tmp_fs = filesystem_of(Path("/tmp"), pr["mounts"])
    rep.add("filesystem at /tmp", tmp_fs, "info",
            "keep work and scratch on the output disk, not /tmp, when /tmp is a small overlay"
            if tmp_fs in ("overlay", "tmpfs") else "")

    # Python and uv
    rep.add("python", platform.python_version(),
            "problem" if sys.version_info < (3, 10) else "ok",
            "the wrapper needs Python 3.10 or newer" if sys.version_info < (3, 10) else "")
    uv = which("uv")
    rep.add("uv", uv or "not found", "ok" if uv else "info",
            "" if uv else "optional; used to run the ClawBio tests and the skill linter")

    # Java and Nextflow
    java = which("java")
    rep.add("java", pr["run_version"](["java", "-version"], r'version "([^"]+)"') or "found"
            if java else "not found", "ok" if java else "problem",
            "" if java else "Nextflow needs a JDK, for example: conda install -c conda-forge openjdk")
    nf = which("nextflow")
    rep.add("nextflow", pr["run_version"](["nextflow", "-version"], r"version\s+(\S+)") or "found"
            if nf else "not found", "ok" if nf else "problem",
            "" if nf else "install from https://www.nextflow.io/docs/latest/install.html")

    # Container runtimes
    docker = which("docker")
    docker_ok = pr["docker_ok"]
    if docker and docker_ok is None:
        try:
            docker_ok = subprocess.run(["docker", "info"], capture_output=True,
                                       timeout=30).returncode == 0
        except (OSError, subprocess.TimeoutExpired):
            docker_ok = False
    if docker:
        rep.add("docker", "daemon answers" if docker_ok else "installed, daemon not reachable",
                "ok" if docker_ok else "warn")
    else:
        rep.add("docker", "not found", "info")

    sif_rt = which("apptainer") or which("singularity")
    sif_name = Path(sif_rt).name if sif_rt else ""
    if sif_rt:
        ver = pr["run_version"]([sif_name, "--version"], r"(\d+\.\d+(?:\.\d+)?\S*)")
        rep.add(sif_name, ver or "found", "ok")
    else:
        rep.add("apptainer/singularity", "not found", "info")

    overlay = fuse_helper = None
    if sif_rt:
        conf = pr["apptainer_conf"] if pr["apptainer_conf"] is not None else find_apptainer_conf()
        overlay = apptainer_overlay_setting(conf)
        fuse_helper = which("fusermount3") or which("fusermount")
        rep.add("/dev/fuse", "present" if pr["dev_fuse"] else "absent",
                "ok" if pr["dev_fuse"] else "warn")
        rep.add("fusermount", fuse_helper or "not found", "ok" if fuse_helper else "warn")
        rep.add("apptainer 'enable overlay'", overlay or "default", "info")

    # Optional functional probe
    probe_result = None
    if image_to_probe is not None:
        if not sif_rt:
            rep.add("image probe", "skipped", "info", "no apptainer/singularity to run it with")
        else:
            worked, msg = pr["probe_image"](sif_name, image_to_probe)
            probe_result = worked
            rep.add(f"image probe ({image_to_probe.name})", "works" if worked else "failed",
                    "ok" if worked else "problem", msg)

    # Choose a profile
    if docker and docker_ok:
        rep.profile = "docker"
        rep.reasons.append("Docker is installed and its daemon answered.")
    elif sif_rt:
        sandbox_signals = []
        if gv:
            sandbox_signals.append("the kernel is gVisor")
        if not pr["dev_fuse"]:
            sandbox_signals.append("/dev/fuse is absent")
        if not fuse_helper:
            sandbox_signals.append("no fusermount helper")
        if overlay == "no":
            sandbox_signals.append("apptainer overlay is turned off")
        if probe_result is True and image_to_probe is not None and not image_to_probe.is_dir():
            rep.profile = "singularity-sif"
            rep.reasons.append("a .sif image ran, so the standard Singularity profile works.")
        elif probe_result is False or sandbox_signals:
            rep.profile = "singularity-sandbox"
            why = ", ".join(sandbox_signals) or "the image probe failed"
            rep.reasons.append(
                f"{sif_name} is present but .sif images are unlikely to mount ({why}). "
                "Pre-built sandbox directories avoid the mount.")
        else:
            rep.profile = "singularity-sif"
            rep.reasons.append(f"{sif_name} is present and no signal suggests .sif mounts fail. "
                               "Confirm with --probe-image on any .sif or sandbox you have.")
    else:
        rep.profile = "none"
        rep.reasons.append("No container runtime was found.")

    rep.next_steps = next_steps(rep.profile, target, sif_name)
    return rep


def next_steps(profile: str, out_dir: Path, sif_name: str) -> list[str]:
    home = Path.home()
    if profile == "docker":
        return ["Run the wrapper with --profile docker (the default)."]
    if profile == "singularity-sif":
        return ["Run the wrapper with --profile singularity.",
                "If a run fails with a mount, overlay or FUSE error, rerun this check with "
                "--probe-image <a .sif file> and follow the sandbox steps."]
    if profile == "singularity-sandbox":
        cache = home / "bj-expression-sandboxes"
        cfg = cache / "sandbox.config"
        return [
            f"Build the sandboxes once (large download): python apptainer_sandboxes.py "
            f"--cache-dir {cache} --bind {home} --config-out {cfg}",
            f"Run the wrapper with --profile singularity --extra-config {cfg}",
            "Folders outside your home (FASTQs, references) must be bound too: add --bind <path> to the "
            "apptainer_sandboxes.py command for each one, then rebuild the config.",
        ]
    return ["Install Docker (or Apptainer), then rerun this check.",
            "Docker needs a host where you can run a daemon, such as a VM or workstation."]


# ──────────────────────────────────────────────
# Output
# ──────────────────────────────────────────────
SYMBOL = {"ok": "ok  ", "info": "info", "warn": "WARN", "problem": "FAIL"}


def format_report(rep: Report) -> str:
    width = max((len(f.name) for f in rep.findings), default=10)
    lines = ["Environment check (read-only)", ""]
    for f in rep.findings:
        line = f"  [{SYMBOL[f.status]}] {f.name.ljust(width)}  {f.value}"
        if f.note:
            line += f"\n         {' ' * width}  {f.note}"
        lines.append(line)
    lines += ["", f"Recommended profile: {rep.profile}"]
    lines += [f"  why: {r}" for r in rep.reasons]
    lines += ["", "Next steps:"] + [f"  {i}. {s}" for i, s in enumerate(rep.next_steps, 1)]
    return "\n".join(lines)


def to_json(rep: Report) -> dict:
    return {"profile": rep.profile, "reasons": rep.reasons, "next_steps": rep.next_steps,
            "worst_status": rep.worst(), "findings": [asdict(f) for f in rep.findings]}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--path", type=Path, default=Path.cwd(),
                    help="where results and work files will go (default: current directory)")
    ap.add_argument("--min-free-gb", type=int, default=DEFAULT_MIN_FREE_GB)
    ap.add_argument("--probe-image", type=Path,
                    help="a .sif file or sandbox directory to run one tiny command in")
    ap.add_argument("--json-out", type=Path)
    args = ap.parse_args(argv)
    rep = detect(args.path, args.min_free_gb, args.probe_image)
    print(format_report(rep))
    if args.json_out:
        args.json_out.write_text(json.dumps(to_json(rep), indent=2) + "\n")
    return 1 if rep.worst() == "problem" else 0


if __name__ == "__main__":
    sys.exit(main())
