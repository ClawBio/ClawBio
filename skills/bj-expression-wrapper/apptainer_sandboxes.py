#!/usr/bin/env python3
"""Build Apptainer sandbox directories so Nextflow can run containers on hosts
that cannot mount .sif images.

Why: in a gVisor container (no FUSE, no overlay) `apptainer exec image.sif`
fails. A sandbox *directory* run with `--writable` works. Nextflow looks in
`singularity.cacheDir` for a file named after the image and uses it without
pulling, and it will use a directory there. So: build one sandbox per image,
named the way Nextflow names it, and add the config this module writes.

Without an overlay, Apptainer cannot create bind-mount points inside the image,
so each bound path must already exist inside every sandbox; this module makes
them. Nextflow's own `autoMounts` is turned off for the same reason.

Tested end to end with Nextflow 25.04.6 and one real image (seqtk 1.3) only.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
DEFAULT_IMAGES_FILE = SKILL_DIR / "bj_expression_2_0_0_rna_images.txt"


class SandboxError(Exception):
    """A user-facing failure."""


def cache_name(image: str) -> str:
    """The file name Nextflow looks for in singularity.cacheDir."""
    name = image.split("://", 1)[-1]
    return name.replace(":", "-").replace("/", "-") + ".img"


def read_images(path: Path) -> list[str]:
    lines = (ln.strip() for ln in Path(path).read_text().splitlines())
    return [ln for ln in lines if ln and not ln.startswith("#")]


def sandbox_config(cache_dir: Path, bind_paths: list[Path]) -> str:
    binds = " ".join(f"--bind {p}:{p}" for p in bind_paths)
    return (
        "// Written by apptainer_sandboxes.py: run pre-built sandbox directories\n"
        "// where .sif images cannot be mounted (no FUSE / overlay).\n"
        "singularity {\n"
        "  enabled = true\n"
        "  autoMounts = false\n"
        f"  cacheDir = '{cache_dir}'\n"
        f"  runOptions = '--writable {binds}'\n"
        "}\n"
    )


def _find_apptainer() -> str:
    exe = shutil.which("apptainer") or shutil.which("singularity")
    if not exe:
        raise SandboxError("apptainer not found on PATH (conda install -c conda-forge apptainer)")
    return exe


def build_all(images: list[str], cache_dir: Path, bind_paths: list[Path],
              errors: dict | None = None) -> dict:
    """Build any sandbox that is missing. One failure does not stop the rest."""
    apptainer = _find_apptainer()
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    scratch = cache_dir / ".apptainer-tmp"
    env = {**os.environ,
           "APPTAINER_CACHEDIR": str(cache_dir / ".apptainer-cache"),
           "APPTAINER_TMPDIR": str(scratch)}  # keep big temp files off /tmp
    scratch.mkdir(exist_ok=True)

    status: dict = {"built": [], "skipped": [], "failed": []}
    for image in images:
        name = cache_name(image)
        final, partial = cache_dir / name, cache_dir / f"{name}.partial"
        if final.is_dir():
            status["skipped"].append(name)
            continue
        shutil.rmtree(partial, ignore_errors=True)
        proc = subprocess.run(
            [apptainer, "build", "--sandbox", "--force", str(partial), f"docker://{image}"],
            capture_output=True, text=True, env=env)
        if proc.returncode != 0:
            shutil.rmtree(partial, ignore_errors=True)
            status["failed"].append(name)
            if errors is not None:
                errors[image] = proc.stderr.strip()[-500:]
            continue
        for p in bind_paths:
            (partial / str(p).lstrip("/")).mkdir(parents=True, exist_ok=True)
        partial.rename(final)  # only a finished sandbox ever has the final name
        status["built"].append(name)
    return status


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="apptainer-sandboxes",
        description="Build Apptainer sandboxes for Nextflow on hosts without FUSE/overlay.")
    ap.add_argument("--cache-dir", type=Path,
                    default=Path.home() / "bj-expression-sandboxes",
                    help="singularity.cacheDir to fill (default ~/bj-expression-sandboxes)")
    ap.add_argument("--images-file", type=Path, default=DEFAULT_IMAGES_FILE)
    ap.add_argument("--bind", type=Path, action="append",
                    help="Host path to bind at the same path in containers "
                         "(repeatable; default: your home directory)")
    ap.add_argument("--config-out", type=Path,
                    help="Write the Nextflow config here (use it with --extra-config)")
    ap.add_argument("--dry-run", action="store_true", help="List what would be built")
    args = ap.parse_args(argv)

    try:
        images = read_images(args.images_file)
        binds = [p.resolve() for p in (args.bind or [Path.home()])]
        cache = args.cache_dir.expanduser().resolve()
        if args.dry_run:
            for im in images:
                state = "exists" if (cache / cache_name(im)).is_dir() else "to build"
                print(f"{cache_name(im)}  <-  {im}  [{state}]")
            return 0
        errors: dict = {}
        status = build_all(images, cache, binds, errors)
    except SandboxError as exc:
        print(f"apptainer-sandboxes: error: {exc}", file=sys.stderr)
        return 1

    if args.config_out:
        args.config_out.parent.mkdir(parents=True, exist_ok=True)
        args.config_out.write_text(sandbox_config(cache, binds))
        print(f"Config: {args.config_out}")
    print(f"built {len(status['built'])}, skipped {len(status['skipped'])}, "
          f"failed {len(status['failed'])}")
    for image, err in errors.items():
        print(f"apptainer-sandboxes: error: {image}: {err}", file=sys.stderr)
    return 1 if status["failed"] else 0


if __name__ == "__main__":
    sys.exit(main())
