"""Tests for apptainer_sandboxes.py (red/green TDD).

Background: on hosts where Apptainer cannot mount .sif images (no FUSE, no
overlay; e.g. a gVisor container), Nextflow can still use a pre-built sandbox
*directory* found in `singularity.cacheDir` under the name it would have given
the pulled image. This module builds those directories and writes the config.
No test pulls a real image: they use a fake `apptainer`.
"""

import importlib.util
import json
import os
import re
import stat
import sys
from pathlib import Path

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

SKILL_DIR = Path(__file__).resolve().parent.parent
MODULE = SKILL_DIR / "apptainer_sandboxes.py"
IMAGES_FILE = SKILL_DIR / "bj_expression_2_0_0_rna_images.txt"


def load():
    spec = importlib.util.spec_from_file_location("apptainer_sandboxes", MODULE)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["apptainer_sandboxes"] = mod
    spec.loader.exec_module(mod)
    return mod


FAKE_APPTAINER = f"""#!{sys.executable}
import json, os, sys
args = sys.argv[1:]
with open(os.environ["FAKE_APPTAINER_LOG"], "a") as fh:
    fh.write(json.dumps({{"args": args,
        "cachedir": os.environ.get("APPTAINER_CACHEDIR"),
        "tmpdir": os.environ.get("APPTAINER_TMPDIR")}}) + "\\n")
if args[:1] == ["--version"]:
    print("apptainer version 1.5.4"); sys.exit(0)
if args[:2] == ["build", "--sandbox"]:
    dest = args[args.index("--force") + 1] if "--force" in args else args[2]
    if os.environ.get("FAKE_APPTAINER_FAIL") and os.environ["FAKE_APPTAINER_FAIL"] in " ".join(args):
        print("FATAL: pull failed", file=sys.stderr); sys.exit(1)
    os.makedirs(dest + "/bin", exist_ok=True)
    sys.exit(0)
sys.exit(2)
"""


@pytest.fixture
def fake(tmp_path, monkeypatch):
    bindir = tmp_path / "bin"
    bindir.mkdir()
    exe = bindir / "apptainer"
    exe.write_text(FAKE_APPTAINER)
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    log = tmp_path / "apptainer.log"
    monkeypatch.setenv("PATH", str(bindir))
    monkeypatch.setenv("FAKE_APPTAINER_LOG", str(log))
    monkeypatch.delenv("FAKE_APPTAINER_FAIL", raising=False)
    return {"log": log, "tmp": tmp_path}


def calls(log: Path) -> list[dict]:
    return [json.loads(ln) for ln in log.read_text().splitlines()] if log.exists() else []


# ──────────────────────────────────────────────
# Cache names
# ──────────────────────────────────────────────

def test_cache_name_matches_what_nextflow_looked_up_for_a_real_image():
    # Verified against Nextflow 25.04.6: this is the file name it used.
    mod = load()
    assert (mod.cache_name("quay.io/biocontainers/seqtk:1.3--h5bf99c6_3")
            == "quay.io-biocontainers-seqtk-1.3--h5bf99c6_3.img")


def test_docker_scheme_is_ignored():
    mod = load()
    a = mod.cache_name("docker://quay.io/biocontainers/star:2.7.6a--0")
    assert a == mod.cache_name("quay.io/biocontainers/star:2.7.6a--0")


IMAGE = st.from_regex(r"[a-z0-9.]{1,12}(/[a-z0-9_]{1,10}){1,2}:[A-Za-z0-9_.-]{1,12}",
                      fullmatch=True)


@settings(deadline=None, max_examples=80)
@given(image=IMAGE)
def test_property_cache_names_are_single_path_components(image):
    mod = load()
    name = mod.cache_name(image)
    assert name.endswith(".img")
    assert "/" not in name and ":" not in name
    assert mod.cache_name(image) == name  # deterministic


# ──────────────────────────────────────────────
# The committed image list
# ──────────────────────────────────────────────

def test_image_list_is_the_17_public_images_with_distinct_cache_names():
    mod = load()
    images = mod.read_images(IMAGES_FILE)
    assert len(images) == 17
    assert len(set(images)) == 17
    for im in images:
        assert re.fullmatch(r"(public\.ecr\.aws/v2u5x5f1/basejumper|quay\.io/biocontainers/[a-z0-9_.-]+):\S+", im), im
    assert len({mod.cache_name(i) for i in images}) == 17


def test_read_images_skips_comments_and_blanks(tmp_path):
    mod = load()
    p = tmp_path / "l.txt"
    p.write_text("# header\n\nquay.io/biocontainers/a:1\n  # note\nquay.io/biocontainers/b:2  \n")
    assert mod.read_images(p) == ["quay.io/biocontainers/a:1", "quay.io/biocontainers/b:2"]


# ──────────────────────────────────────────────
# Config text
# ──────────────────────────────────────────────

def test_config_turns_off_automounts_and_binds_explicitly(tmp_path):
    mod = load()
    text = mod.sandbox_config(tmp_path / "cache", [Path("/home/x"), Path("/data")])
    assert "autoMounts = false" in text
    assert f"cacheDir = '{tmp_path / 'cache'}'" in text
    assert "--writable" in text
    assert "--bind /home/x:/home/x" in text and "--bind /data:/data" in text
    assert "enabled = true" in text


# ──────────────────────────────────────────────
# Building
# ──────────────────────────────────────────────

def test_build_creates_sandbox_with_mount_points_and_no_partial_left(fake):
    mod = load()
    cache = fake["tmp"] / "cache"
    status = mod.build_all(["quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"], cache,
                           [Path("/home/x")])
    final = cache / "quay.io-biocontainers-seqtk-1.3--h5bf99c6_3.img"
    assert final.is_dir() and (final / "bin").is_dir()
    assert (final / "home" / "x").is_dir(), "bind destination must exist in the sandbox"
    assert not list(cache.glob("*.partial"))
    assert status == {"built": [str(final.name)], "skipped": [], "failed": []}


def test_existing_sandboxes_are_skipped(fake):
    mod = load()
    cache = fake["tmp"] / "cache"
    img = "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"
    mod.build_all([img], cache, [Path("/home/x")])
    n = len(calls(fake["log"]))
    status = mod.build_all([img], cache, [Path("/home/x")])
    assert status["skipped"] and not status["built"]
    assert len(calls(fake["log"])) == n, "must not rebuild"


def test_a_failed_build_leaves_nothing_that_looks_complete(fake, monkeypatch):
    mod = load()
    cache = fake["tmp"] / "cache"
    monkeypatch.setenv("FAKE_APPTAINER_FAIL", "star")
    status = mod.build_all(["quay.io/biocontainers/star:2.7.6a--0",
                            "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"],
                           cache, [Path("/home/x")])
    assert [Path(f).name for f in status["failed"]] == ["quay.io-biocontainers-star-2.7.6a--0.img"] \
        or any("star" in f for f in status["failed"])
    assert not (cache / "quay.io-biocontainers-star-2.7.6a--0.img").exists()
    assert (cache / "quay.io-biocontainers-seqtk-1.3--h5bf99c6_3.img").is_dir(), \
        "one failure must not stop the others"


def test_apptainer_scratch_stays_under_the_cache_dir_not_tmp(fake):
    mod = load()
    cache = fake["tmp"] / "cache"
    mod.build_all(["quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"], cache, [Path("/home/x")])
    build = [c for c in calls(fake["log"]) if c["args"][:1] == ["build"]][0]
    assert build["cachedir"].startswith(str(cache))
    assert build["tmpdir"].startswith(str(cache))


def test_missing_apptainer_is_a_clear_error(tmp_path, monkeypatch):
    mod = load()
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(mod.SandboxError, match="apptainer"):
        mod.build_all(["quay.io/biocontainers/a:1"], tmp_path / "c", [Path("/home/x")])


# ──────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────

def test_cli_dry_run_lists_work_and_runs_nothing(fake, capsys):
    mod = load()
    cache = fake["tmp"] / "cache"
    code = mod.main(["--cache-dir", str(cache), "--dry-run"])
    out = capsys.readouterr().out
    assert code == 0
    assert out.count(".img") >= 17
    assert not any(c["args"][:1] == ["build"] for c in calls(fake["log"]))
    assert not cache.exists() or not list(cache.glob("*.img"))


def test_cli_writes_the_config_file(fake):
    mod = load()
    cache, cfg = fake["tmp"] / "cache", fake["tmp"] / "out" / "sandbox.config"
    code = mod.main(["--cache-dir", str(cache), "--images-file", str(IMAGES_FILE),
                     "--config-out", str(cfg), "--bind", str(fake["tmp"])])
    assert code == 0
    assert "--writable" in cfg.read_text()
    assert f"--bind {fake['tmp']}:{fake['tmp']}" in cfg.read_text()


def test_cli_exit_code_is_nonzero_when_any_build_fails(fake, monkeypatch, capsys):
    mod = load()
    monkeypatch.setenv("FAKE_APPTAINER_FAIL", "multiqc")
    code = mod.main(["--cache-dir", str(fake["tmp"] / "cache")])
    assert code == 1
    assert "multiqc" in capsys.readouterr().err
