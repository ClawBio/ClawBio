"""Tests for the bj-expression-wrapper skill (red/green TDD).

BioSkryb's bj-expression pipeline needs Java, Nextflow and Docker/Singularity,
none of which CI has. These tests therefore use a fake `nextflow` (plus fake
`java` and `docker`) on a restricted PATH. The fake writes the output tree that
`workflow/rnaseq.nf` publishes. Nothing here runs the real pipeline.
"""

import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from hypothesis import example, given, settings
from hypothesis import strategies as st

SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT = SKILL_DIR / "bj_expression_wrapper.py"
LICENSE_URL = "https://www.bioskryb.com/basejumper/local-license/"
PIPELINE_TAG = "2.0.0"


def _load_module():
    spec = importlib.util.spec_from_file_location("bj_expression_wrapper", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["bj_expression_wrapper"] = mod
    spec.loader.exec_module(mod)
    return mod


# ──────────────────────────────────────────────
# Fixtures: fake tools, fake pipeline checkout, fake FASTQs
# ──────────────────────────────────────────────

FAKE_NEXTFLOW = f"""#!{sys.executable}
import json, os, sys
args = sys.argv[1:]
log = os.environ.get("FAKE_NF_LOG")
if log:
    with open(log, "a") as fh:
        fh.write(json.dumps(args) + "\\n")
envlog = os.environ.get("FAKE_NF_ENVLOG")
if envlog:
    with open(envlog, "a") as fh:
        fh.write(os.environ.get("NXF_VER", "<unset>") + "\\n")
if args[:1] == ["-version"]:
    print("nextflow version 25.04.0 build 5000")
    sys.exit(0)
mode = os.environ.get("FAKE_NF_MODE", "ok")
if mode == "fail":
    if "-work-dir" in args:  # a real failed run still leaves its work/ behind
        wdf = args[args.index("-work-dir") + 1]
        os.makedirs(os.path.join(wdf, "ab", "cd"), exist_ok=True)
        open(os.path.join(wdf, "ab", "cd", "task.txt"), "w").write("scratch\\n")
    print("boom: process STARALIGN failed", file=sys.stderr)
    sys.exit(3)
if mode == "empty":
    sys.exit(0)
params = json.load(open(args[args.index("-params-file") + 1]))
res = params["publish_dir"] + "_" + params["timestamp"]
sa = res + "/secondary_analyses"
drop = os.environ.get("FAKE_NF_DROP", "")
import csv
samples = [r["biosampleName"] for r in csv.DictReader(open(params["input_csv"]))]
files = [
    res + "/multiqc/multiqc_report.html",
    sa + "/secondary_metrics/qualimap_stats_mqc.csv",
    sa + "/secondary_metrics/pipeline_metrics_summary_percents.csv",
    sa + "/secondary_metrics/df_dynamicrange_expression.tsv",
    sa + "/quantification_htseq/df_gene_counts_starhtseq.tsv",
    res + "/tertiary_analyses/classification_cell_typing/df_cell_typing_summary_singler_hpca_gtex_tcga.tsv",
]
# exon duplicates (cell typing and MultiQC built from exon counts) go to <publish_dir>_exon_<ts>
rex = params["publish_dir"] + "_exon_" + params["timestamp"]
sx = rex + "/secondary_analyses"
exon_files = [
    rex + "/multiqc/multiqc_report.html",
    sx + "/quantification_htseq/df_gene_counts_starhtseq.tsv",
    rex + "/tertiary_analyses/classification_cell_typing/df_cell_typing_summary_singler_hpca_gtex_tcga.tsv",
]
if os.environ.get("FAKE_NF_DROP_EXON"):
    exon_files = [f for f in exon_files if os.environ["FAKE_NF_DROP_EXON"] not in f]
files += exon_files
for s in samples:
    if s != drop:
        files += [sa + "/alignment_htseq/" + s + ".bam",
                  sa + "/alignment_htseq/" + s + ".bam.bai"]
wd = args[args.index("-work-dir") + 1] if "-work-dir" in args else None
if wd:
    os.makedirs(os.path.join(wd, "ab", "cd"), exist_ok=True)
    open(os.path.join(wd, "ab", "cd", "task.txt"), "w").write("scratch\\n")
pm = os.environ.get("FAKE_NF_PUBLISH", "copy")
for f in files:
    os.makedirs(os.path.dirname(f), exist_ok=True)
    if pm == "symlink":
        w = os.path.join(os.getcwd(), "fake_work", f.lstrip("/").replace("/", "__"))
        os.makedirs(os.path.dirname(w), exist_ok=True)
        open(w, "w").write("fake\\n")
        os.symlink(w, f)
    else:
        open(f, "w").write("fake\\n")
if pm == "symlink":
    d = os.path.join(os.getcwd(), "fake_work", "multiqc_data")
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "x.txt"), "w").write("data\\n")
    os.symlink(d, res + "/multiqc/multiqc_data")
tr = os.environ.get("FAKE_NF_TRACE", "clean")
if "-with-trace" in args and tr != "none" and not os.path.exists(args[args.index("-with-trace") + 1]):
    tp = args[args.index("-with-trace") + 1]
    rows = ["1\\tab\\t-\\tRNASEQ_WF:STARALIGN (s1)\\tCOMPLETED\\t0"]
    if tr == "ignored":
        rows += ["2\\tcd\\t-\\tRNASEQ_WF:CELL_TYPING (CELL_TYPING)\\tFAILED\\t1",
                 "2\\tcd\\t-\\tRNASEQ_WF:CELL_TYPING (CELL_TYPING)\\tFAILED\\t1"]
    if tr == "retry_ok":
        rows += ["3\\tef\\t-\\tRNASEQ_WF:SALMONQUANT (s2)\\tFAILED\\t137",
                 "3\\tef\\t-\\tRNASEQ_WF:SALMONQUANT (s2)\\tCOMPLETED\\t0"]
    with open(tp, "w") as fh:
        fh.write("task_id\\thash\\tnative_id\\tname\\tstatus\\texit\\n" + "\\n".join(rows) + "\\n")
"""

FAKE_JAVA = '#!/bin/sh\necho \'openjdk version "17.0.9" 2023-10-17\' >&2\n'
FAKE_DOCKER = "#!/bin/sh\nexit 0\n"
FAKE_SIF = "#!/bin/sh\necho 'apptainer version 1.3.0'\n"


def make_tools(tmp_path: Path, *, java=True, docker=True, singularity=False,
               nextflow=True) -> dict:
    """A restricted PATH holding only the fakes requested plus real git."""
    bindir = tmp_path / "tools"
    bindir.mkdir(exist_ok=True)

    def put(name, body):
        p = bindir / name
        p.write_text(body)
        p.chmod(0o755)

    if nextflow:
        put("nextflow", FAKE_NEXTFLOW)
    if java:
        put("java", FAKE_JAVA)
    if docker:
        put("docker", FAKE_DOCKER)
    if singularity:
        put("apptainer", FAKE_SIF)
    git = shutil.which("git")
    if git and not (bindir / "git").exists():
        (bindir / "git").symlink_to(git)
    return {"PATH": str(bindir)}


GIT = ["git", "-c", "user.name=t", "-c", "user.email=t@t",
       "-c", "commit.gpgsign=false", "-c", "tag.gpgsign=false"]


def make_pipeline(tmp_path: Path, tag: str | None = PIPELINE_TAG,
                  dirty: bool = False) -> Path:
    d = tmp_path / "bj-expression"
    d.mkdir()
    (d / "main.nf").write_text("// fake bj-expression\n")
    (d / "nextflow.config").write_text("manifest { version = '2.0.0' }\n")
    subprocess.run([*GIT, "init", "-q", "-b", "main"], cwd=d, check=True)
    subprocess.run([*GIT, "add", "."], cwd=d, check=True)
    subprocess.run([*GIT, "commit", "-q", "-m", "init"], cwd=d, check=True)
    if tag:
        subprocess.run([*GIT, "tag", tag], cwd=d, check=True)
    if dirty:
        (d / "main.nf").write_text("// edited\n")
    return d


def make_fastqs(tmp_path: Path, names=("s1", "s2")) -> dict:
    d = tmp_path / "fq"
    d.mkdir(exist_ok=True)
    out = {}
    for n in names:
        r1, r2 = d / f"{n}_R1_001.fastq.gz", d / f"{n}_R2_001.fastq.gz"
        gz_fastq(r1, 10, 1)
        gz_fastq(r2, 10, 2)
        out[n] = (r1, r2)
    return out


def write_csv(path: Path, rows, groups=True) -> Path:
    head = "biosampleName,read1,read2" + (",groups" if groups else "")
    lines = [head]
    for name, r1, r2 in rows:
        lines.append(f"{name},{r1},{r2}" + (",g1" if groups else ""))
    path.write_text("\n".join(lines) + "\n")
    return path


def run_cli(args, env=None, expect=0, cwd=None):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), *map(str, args)],
        capture_output=True, text=True, env=env, cwd=cwd,
    )
    if expect is not None:
        assert result.returncode == expect, (
            f"expected exit {expect}, got {result.returncode}\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}")
    return result


@pytest.fixture
def setup(tmp_path):
    fq = make_fastqs(tmp_path)
    csv_path = write_csv(tmp_path / "in.csv",
                         [(n, *fq[n]) for n in ("s1", "s2")])
    return {
        "tmp": tmp_path,
        "csv": csv_path,
        "pipeline": make_pipeline(tmp_path),
        "env": make_tools(tmp_path),
        "out": tmp_path / "out",
    }


def base_args(s, *extra):
    choice = [] if any(a in extra for a in ("--reads-per-file", "--skip-subsampling")) \
        else ["--reads-per-file", "5"]
    return ["--input", s["csv"], "--output", s["out"],
            "--pipeline-dir", s["pipeline"], *choice, *extra]


def result_json(out: Path) -> dict:
    return json.loads((out / "result.json").read_text())


# ──────────────────────────────────────────────
# Input CSV validation
# ──────────────────────────────────────────────

def test_valid_csv_is_normalised_to_absolute_paths(tmp_path, monkeypatch):
    mod = _load_module()
    fq = make_fastqs(tmp_path)
    monkeypatch.chdir(tmp_path)
    csv_path = write_csv(tmp_path / "in.csv",
                         [("s1", "fq/s1_R1_001.fastq.gz", "fq/s1_R2_001.fastq.gz")])
    rows = mod.validate_input_csv(csv_path, require_groups=True)
    assert [r["biosampleName"] for r in rows] == ["s1"]
    assert Path(rows[0]["read1"]).is_absolute()
    assert Path(rows[0]["read1"]) == fq["s1"][0].resolve()


def test_missing_required_column_is_named(tmp_path):
    mod = _load_module()
    p = tmp_path / "in.csv"
    p.write_text("biosampleName,read1\ns1,a.fastq.gz\n")
    with pytest.raises(mod.WrapperError, match="read2"):
        mod.validate_input_csv(p, require_groups=False)


def test_header_only_csv_is_rejected(tmp_path):
    mod = _load_module()
    p = tmp_path / "in.csv"
    p.write_text("biosampleName,read1,read2\n")
    with pytest.raises(mod.WrapperError, match="no samples"):
        mod.validate_input_csv(p, require_groups=False)


def test_duplicate_sample_names_are_rejected(tmp_path):
    mod = _load_module()
    fq = make_fastqs(tmp_path, ("a", "b"))
    p = write_csv(tmp_path / "in.csv",
                  [("dup", *fq["a"]), ("dup", *fq["b"])])
    with pytest.raises(mod.WrapperError, match="dup"):
        mod.validate_input_csv(p, require_groups=True)


@pytest.mark.parametrize("bad", ["has space", "a/b", "..", ".hidden", "-lead",
                                 "semi;colon", ""])
def test_unsafe_sample_names_are_rejected(tmp_path, bad):
    mod = _load_module()
    fq = make_fastqs(tmp_path, ("a",))
    p = write_csv(tmp_path / "in.csv", [(bad, *fq["a"])])
    with pytest.raises(mod.WrapperError):
        mod.validate_input_csv(p, require_groups=True)


def test_missing_local_fastq_is_named(tmp_path):
    mod = _load_module()
    fq = make_fastqs(tmp_path, ("a",))
    p = write_csv(tmp_path / "in.csv",
                  [("a", fq["a"][0], tmp_path / "nope_R2.fastq.gz")])
    with pytest.raises(mod.WrapperError, match="nope_R2"):
        mod.validate_input_csv(p, require_groups=True)


def test_same_file_for_read1_and_read2_is_rejected(tmp_path):
    mod = _load_module()
    fq = make_fastqs(tmp_path, ("a",))
    p = write_csv(tmp_path / "in.csv", [("a", fq["a"][0], fq["a"][0])])
    with pytest.raises(mod.WrapperError, match="same file"):
        mod.validate_input_csv(p, require_groups=True)


def test_s3_paths_are_accepted_and_flagged_remote(tmp_path):
    mod = _load_module()
    p = write_csv(tmp_path / "in.csv",
                  [("a", "s3://bucket/a_R1.fastq.gz", "s3://bucket/a_R2.fastq.gz")])
    rows = mod.validate_input_csv(p, require_groups=True)
    assert rows[0]["read1"] == "s3://bucket/a_R1.fastq.gz"
    assert rows[0]["_remote"] is True


def test_groups_column_is_required_unless_qc_plots_are_skipped(tmp_path):
    mod = _load_module()
    fq = make_fastqs(tmp_path, ("a",))
    p = write_csv(tmp_path / "in.csv", [("a", *fq["a"])], groups=False)
    with pytest.raises(mod.WrapperError, match="groups"):
        mod.validate_input_csv(p, require_groups=True)
    assert mod.validate_input_csv(p, require_groups=False)


# ──────────────────────────────────────────────
# Properties
# ──────────────────────────────────────────────

NAME = st.from_regex(r"[A-Za-z0-9][A-Za-z0-9_-]{0,12}", fullmatch=True)


@settings(deadline=None, max_examples=40)
@given(names=st.lists(NAME, min_size=1, max_size=6, unique=True),
       data=st.data())
@example(names=["a"], data=None)
def test_property_validated_rows_roundtrip_and_ignore_row_order(
        tmp_path_factory, names, data):
    mod = _load_module()
    tmp = tmp_path_factory.mktemp("prop")
    fq = make_fastqs(tmp, tuple(names))
    rows = [(n, *fq[n]) for n in names]
    shuffled = list(rows)
    if data is not None:
        shuffled = data.draw(st.permutations(rows))
    a = mod.validate_input_csv(write_csv(tmp / "a.csv", rows), require_groups=True)
    b = mod.validate_input_csv(write_csv(tmp / "b.csv", shuffled), require_groups=True)
    assert sorted(r["biosampleName"] for r in a) == sorted(names)
    key = lambda r: r["biosampleName"]  # noqa: E731
    assert sorted(a, key=key) == sorted(b, key=key)
    # writing the normalised rows back out and re-reading is stable
    out = tmp / "norm.csv"
    mod.write_normalised_csv(a, out)
    c = mod.validate_input_csv(out, require_groups=True)
    assert sorted(c, key=key) == sorted(a, key=key)


@settings(deadline=None, max_examples=60)
@given(name=st.text(min_size=1, max_size=15))
def test_property_accepted_sample_names_are_path_and_shell_safe(name):
    mod = _load_module()
    if mod.is_safe_sample_name(name):
        assert re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", name)
        assert ".." not in name and "/" not in name and " " not in name


# ──────────────────────────────────────────────
# Params and command
# ──────────────────────────────────────────────

def test_params_only_carry_required_values_and_explicit_overrides(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    params = json.loads((setup["out"] / "reproducibility" / "params.json").read_text())
    assert set(params) == {"input_csv", "publish_dir", "timestamp", "genome", "n_reads"}
    assert params["genome"] == "GRCh38"
    assert params["n_reads"] == 10, "--reads-per-file 5 x 2: the pipeline counts R1+R2 together"
    for k in ("read_length", "min_reads", "seqtk_sample_seed"):
        assert k not in params, f"{k} would silently override the pipeline default"


def test_explicit_overrides_are_passed_through(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms", "--reads-per-file", "8",
                      "--read-length", "75", "--min-reads", "100", "--seed", "123",
                      "--skip-cutadapt", "--genome", "GRCm39",
                      "--max-cpus", "8", "--max-memory", "50.GB"),
            env=setup["env"])
    p = json.loads((setup["out"] / "reproducibility" / "params.json").read_text())
    assert p["n_reads"] == 16 and p["read_length"] == 75
    assert p["min_reads"] == 100 and p["seqtk_sample_seed"] == 123
    assert p["skip_cutadapt"] is True and p["genome"] == "GRCm39"
    assert p["max_cpus"] == 8 and p["max_memory"] == "50.GB"


@pytest.mark.parametrize("bad", [
    ["--reads-per-file", "0"], ["--reads-per-file", "-3"], ["--reads-per-file", "lots"],
    ["--n-reads", "0"], ["--read-length", "0"],
    ["--min-reads", "-1"], ["--genome", "hg19"], ["--max-cpus", "0"],
    ["--max-memory", "lots"], ["--profile", "podman"],
])
def test_invalid_options_are_rejected_before_anything_runs(setup, bad):
    res = run_cli(base_args(setup, "--accept-bioskryb-terms", *bad),
                  env=setup["env"], expect=None)
    assert res.returncode in (1, 2)
    # the wrapper's / argparse's own message, not a missing-script traceback
    assert "error:" in res.stderr
    assert not (setup["out"] / "result.json").exists()


def test_command_uses_params_file_work_dir_and_pinned_source(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    cmd = (setup["out"] / "reproducibility" / "commands.sh").read_text()
    assert "nextflow run" in cmd and "-params-file" in cmd
    assert str(setup["pipeline"] / "main.nf") in cmd
    assert "-work-dir" in cmd and str(setup["out"]) in cmd
    assert "-resume" not in cmd
    assert "-profile" not in cmd, "docker is the pipeline default"


def test_singularity_profile_and_resume_are_forwarded(tmp_path):
    fq = make_fastqs(tmp_path)
    s = {"tmp": tmp_path, "pipeline": make_pipeline(tmp_path),
         "env": make_tools(tmp_path, docker=False, singularity=True),
         "out": tmp_path / "out",
         "csv": write_csv(tmp_path / "in.csv", [("s1", *fq["s1"])])}
    run_cli(base_args(s, "--accept-bioskryb-terms", "--profile", "singularity",
                      "--resume"), env=s["env"])
    cmd = (s["out"] / "reproducibility" / "commands.sh").read_text()
    assert "-profile singularity" in cmd and "-resume" in cmd


# ──────────────────────────────────────────────
# Pipeline source: never fetched, must be the pinned tag
# ──────────────────────────────────────────────

def test_missing_pipeline_dir_explains_how_to_get_it_and_the_licence(setup):
    res = run_cli(["--input", setup["csv"], "--output", setup["out"],
                   "--pipeline-dir", setup["tmp"] / "nowhere", "--reads-per-file", "5",
                   "--accept-bioskryb-terms"], env=setup["env"], expect=None)
    assert res.returncode != 0
    assert "git clone" in res.stderr and LICENSE_URL in res.stderr


def test_checkout_not_at_the_pinned_tag_is_refused(tmp_path):
    fq = make_fastqs(tmp_path)
    s = {"tmp": tmp_path, "pipeline": make_pipeline(tmp_path, tag=None),
         "env": make_tools(tmp_path), "out": tmp_path / "out",
         "csv": write_csv(tmp_path / "in.csv", [("s1", *fq["s1"])])}
    res = run_cli(base_args(s, "--accept-bioskryb-terms"), env=s["env"], expect=None)
    assert res.returncode != 0
    assert PIPELINE_TAG in res.stderr


def test_dirty_checkout_is_refused_unless_allowed(tmp_path):
    fq = make_fastqs(tmp_path)
    s = {"tmp": tmp_path, "pipeline": make_pipeline(tmp_path, dirty=True),
         "env": make_tools(tmp_path), "out": tmp_path / "out",
         "csv": write_csv(tmp_path / "in.csv", [("s1", *fq["s1"])])}
    res = run_cli(base_args(s, "--accept-bioskryb-terms"), env=s["env"], expect=None)
    assert res.returncode != 0 and "uncommitted" in res.stderr.lower()
    s["out"] = tmp_path / "out2"
    run_cli(base_args(s, "--accept-bioskryb-terms", "--allow-dirty-pipeline"),
            env=s["env"])
    assert result_json(s["out"])["pipeline"]["dirty"] is True


# ──────────────────────────────────────────────
# Licence acknowledgement
# ──────────────────────────────────────────────

def test_running_requires_explicit_acceptance_of_bioskryb_terms(setup):
    log = setup["tmp"] / "nf.log"
    env = {**setup["env"], "FAKE_NF_LOG": str(log)}
    res = run_cli(base_args(setup), env=env, expect=None)
    assert res.returncode != 0
    assert LICENSE_URL in res.stderr
    assert not log.exists(), "nextflow must not be launched without acceptance"


def test_check_mode_does_not_need_acceptance(setup):
    run_cli(base_args(setup, "--check"), env=setup["env"])


# ──────────────────────────────────────────────
# Running (fake nextflow)
# ──────────────────────────────────────────────

def test_successful_run_records_status_outputs_and_samples(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    r = result_json(setup["out"])
    assert r["status"] == "success"
    assert r["terms_accepted"] is True
    assert sorted(r["samples"]) == ["s1", "s2"]
    assert r["samples_missing_from_outputs"] == []
    assert r["outputs"]["multiqc_report"].endswith("multiqc_report.html")
    assert "salmon_gene_counts" not in r["outputs"], "Salmon is not run in the single-nucleus version"
    assert r["outputs"]["htseq_gene_counts"].endswith("df_gene_counts_starhtseq.tsv")
    assert r["exon_duplicate"]["results_dir"].endswith("results_exon_clawbio")
    assert r["exon_duplicate"]["missing_outputs"] == []
    assert r["exon_duplicate"]["outputs"]["multiqc_report"].endswith("multiqc/multiqc_report.html")
    assert r["outputs"]["htseq_gene_counts"].endswith("df_gene_counts_starhtseq.tsv")
    assert r["outputs"]["qualimap_stats"].endswith("qualimap_stats_mqc.csv")
    assert r["outputs"]["pipeline_metrics_summary_percents"].endswith("pipeline_metrics_summary_percents.csv")
    # BioSkryb's README names pipeline_metrics_summary.csv; 2.0.0 does not write it
    assert not any("pipeline_metrics_summary.csv" in m for m in r["missing_outputs"])
    assert Path(r["outputs"]["results_dir"]).is_dir()
    assert r["pipeline"]["version"] == PIPELINE_TAG


def test_report_has_disclaimer_and_pipeline_licence_notice(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    text = (setup["out"] / "report.md").read_text()
    assert ("ClawBio is a research and educational tool. It is not a medical "
            "device and does not provide clinical diagnoses.") in text
    assert LICENSE_URL in text


def test_samples_dropped_by_the_pipeline_are_called_out(setup):
    env = {**setup["env"], "FAKE_NF_DROP": "s2"}
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=env)
    r = result_json(setup["out"])
    assert r["samples_missing_from_outputs"] == ["s2"]
    text = (setup["out"] / "report.md").read_text().lower()
    assert "s2" in text and "min_reads" in text


def test_nextflow_failure_is_reported_with_stderr_tail(setup):
    env = {**setup["env"], "FAKE_NF_MODE": "fail"}
    res = run_cli(base_args(setup, "--accept-bioskryb-terms"), env=env, expect=1)
    r = result_json(setup["out"])
    assert r["status"] == "failed" and r["nextflow_exit_code"] == 3
    assert "STARALIGN" in r["stderr_tail"]
    assert "STARALIGN" in res.stderr


def test_exit_zero_without_key_outputs_is_not_success(setup):
    env = {**setup["env"], "FAKE_NF_MODE": "empty"}
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=env, expect=1)
    r = result_json(setup["out"])
    assert r["status"] == "incomplete"
    assert r["missing_outputs"], "must list what was expected but not produced"


# ──────────────────────────────────────────────
# Preflight and --check
# ──────────────────────────────────────────────

def test_check_validates_without_launching_nextflow(setup):
    log = setup["tmp"] / "nf.log"
    env = {**setup["env"], "FAKE_NF_LOG": str(log)}
    run_cli(base_args(setup, "--check"), env=env)
    r = result_json(setup["out"])
    assert r["status"] == "check_passed"
    assert not (setup["out"] / "results_clawbio").exists()
    # only `nextflow -version` may have been called
    calls = [json.loads(ln) for ln in log.read_text().splitlines()] if log.exists() else []
    assert all(c == ["-version"] for c in calls)


@pytest.mark.parametrize("missing,needle", [
    ({"java": False}, "java"),
    ({"nextflow": False}, "nextflow"),
    ({"docker": False}, "singularity"),
])
def test_missing_tools_are_named_with_a_fix(tmp_path, missing, needle):
    fq = make_fastqs(tmp_path)
    s = {"tmp": tmp_path, "pipeline": make_pipeline(tmp_path),
         "env": make_tools(tmp_path, **missing), "out": tmp_path / "out",
         "csv": write_csv(tmp_path / "in.csv", [("s1", *fq["s1"])])}
    res = run_cli(base_args(s, "--check"), env=s["env"], expect=None)
    assert res.returncode != 0
    assert needle in res.stderr.lower()


def test_singularity_profile_needs_apptainer_or_singularity(tmp_path):
    fq = make_fastqs(tmp_path)
    s = {"tmp": tmp_path, "pipeline": make_pipeline(tmp_path),
         "env": make_tools(tmp_path, docker=True, singularity=False),
         "out": tmp_path / "out",
         "csv": write_csv(tmp_path / "in.csv", [("s1", *fq["s1"])])}
    res = run_cli(base_args(s, "--check", "--profile", "singularity"),
                  env=s["env"], expect=None)
    assert res.returncode != 0
    assert "apptainer" in res.stderr.lower() or "singularity" in res.stderr.lower()


# ──────────────────────────────────────────────
# Output hygiene, reproducibility, demo
# ──────────────────────────────────────────────

def test_existing_results_are_not_overwritten(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    before = (setup["out"] / "result.json").read_bytes()
    res = run_cli(base_args(setup, "--accept-bioskryb-terms"),
                  env=setup["env"], expect=None)
    assert res.returncode != 0
    assert (setup["out"] / "result.json").read_bytes() == before


def test_reproducibility_bundle_is_complete_and_checksums_verify(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    repro = setup["out"] / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256",
                 "params.json", "samplesheet.normalised.csv"):
        assert (repro / name).exists(), name
    for line in (repro / "checksums.sha256").read_text().splitlines():
        digest, label = line.split("  ", 1)
        assert hashlib.sha256((setup["out"] / label).read_bytes()).hexdigest() == digest


def test_writes_only_under_output_dir(setup):
    work = setup["tmp"] / "cwd"
    work.mkdir()
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"], cwd=work)
    assert list(work.iterdir()) == []


def test_demo_needs_no_tools_and_says_nothing_was_run(tmp_path):
    out = tmp_path / "demo"
    run_cli(["--demo", "--output", out], env={"PATH": str(tmp_path)})
    r = result_json(out)
    assert r["status"] == "demo_plan_only"
    assert "plan only" in (out / "report.md").read_text().lower()


def _parse_output_contract(skill_md: Path) -> list[str]:
    text = skill_md.read_text()
    m = re.search(r"## Output Structure\s+```[^\n]*\n(.*?)```", text, re.S)
    if not m:
        return []
    files, stack = [], []
    for line in m.group(1).splitlines()[1:]:
        if not line.strip() or "(optional)" in line:
            continue
        depth = (len(line) - len(line.lstrip("│├└─ "))) // 4 - 1
        name = re.split(r"\s+#", line.lstrip("│├└─ "))[0].strip()
        stack = stack[:depth] + [name.rstrip("/")]
        if not name.endswith("/"):
            files.append("/".join(stack))
    return files


class TestOutputContract:
    def test_documented_outputs_are_produced_by_the_demo(self, tmp_path):
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        assert promised, "SKILL.md has no parseable Output Structure"
        run_cli(["--demo", "--output", tmp_path], env={"PATH": str(tmp_path / "x")})
        missing = [p for p in promised if not (tmp_path / p).exists()]
        assert not missing, f"promised but not produced: {missing}"


# ──────────────────────────────────────────────
# Nextflow version pin and extra config
# ──────────────────────────────────────────────
# Nextflow 26.x cannot parse bj-expression 2.0.0's config (it includes profile
# files that are not in the public repo); 25.04.6 and 24.10.5 can. So the
# wrapper pins NXF_VER for every nextflow call it makes.

def _env_lines(path: Path) -> list[str]:
    return path.read_text().splitlines() if path.exists() else []


def test_nextflow_version_is_pinned_by_default_for_every_call(setup):
    envlog = setup["tmp"] / "env.log"
    env = {**setup["env"], "FAKE_NF_ENVLOG": str(envlog)}
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=env)
    seen = _env_lines(envlog)
    assert len(seen) >= 2, "expected the preflight -version call and the run"
    assert set(seen) == {"25.04.6"}
    r = result_json(setup["out"])
    assert r["nextflow_version_requested"] == "25.04.6"
    cmd = (setup["out"] / "reproducibility" / "commands.sh").read_text()
    assert "NXF_VER=25.04.6" in cmd


def test_nextflow_version_can_be_changed(setup):
    envlog = setup["tmp"] / "env.log"
    env = {**setup["env"], "FAKE_NF_ENVLOG": str(envlog)}
    run_cli(base_args(setup, "--accept-bioskryb-terms",
                      "--nextflow-version", "24.10.5"), env=env)
    assert set(_env_lines(envlog)) == {"24.10.5"}


def test_nextflow_version_system_leaves_the_installed_one_alone(setup):
    envlog = setup["tmp"] / "env.log"
    env = {**setup["env"], "FAKE_NF_ENVLOG": str(envlog)}
    run_cli(base_args(setup, "--accept-bioskryb-terms",
                      "--nextflow-version", "system"), env=env)
    assert set(_env_lines(envlog)) == {"<unset>"}
    cmd = (setup["out"] / "reproducibility" / "commands.sh").read_text()
    assert "NXF_VER" not in cmd


@pytest.mark.parametrize("bad", ["abc", "25.04", "25.04.6; rm -rf /", "", "v25.04.6"])
def test_invalid_nextflow_version_is_rejected(setup, bad):
    res = run_cli(base_args(setup, "--accept-bioskryb-terms",
                            "--nextflow-version", bad), env=setup["env"], expect=None)
    assert res.returncode in (1, 2) and "error:" in res.stderr
    # the wrapper's own message, not argparse's "unrecognized arguments"
    assert "must look like 25.04.6" in res.stderr
    assert not (setup["out"] / "result.json").exists()


def test_extra_config_is_passed_with_dash_c_and_kept_in_the_bundle(setup):
    cfg = setup["tmp"] / "site.config"
    cfg.write_text("singularity { autoMounts = false }\n")
    run_cli(base_args(setup, "--accept-bioskryb-terms",
                      "--extra-config", cfg), env=setup["env"])
    repro = setup["out"] / "reproducibility"
    kept = repro / "extra.config"
    assert kept.read_text() == cfg.read_text()
    cmd = (repro / "commands.sh").read_text()
    assert f"-c {kept}" in cmd
    assert result_json(setup["out"])["extra_config"] == str(kept)
    checksums = (repro / "checksums.sha256").read_text()
    assert "reproducibility/extra.config" in checksums


def test_missing_extra_config_is_named(setup):
    res = run_cli(base_args(setup, "--accept-bioskryb-terms",
                            "--extra-config", setup["tmp"] / "nope.config"),
                  env=setup["env"], expect=None)
    assert res.returncode != 0
    assert "nope.config" in res.stderr and "not found" in res.stderr.lower()
    assert "unrecognized" not in res.stderr


def test_no_extra_config_means_no_dash_c(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    cmd = (setup["out"] / "reproducibility" / "commands.sh").read_text()
    assert " -c " not in cmd
    assert result_json(setup["out"])["extra_config"] is None


def test_relative_output_dir_works_because_nextflow_runs_elsewhere(setup):
    """Nextflow runs in <output>/nextflow_run, so every path handed to it must
    be absolute or a relative --output would lose its own params file."""
    res = run_cli(["--input", setup["csv"], "--output", "rel_out", "--reads-per-file", "5",
                   "--pipeline-dir", setup["pipeline"], "--accept-bioskryb-terms"],
                  env=setup["env"], cwd=setup["tmp"])
    r = json.loads((setup["tmp"] / "rel_out" / "result.json").read_text())
    assert r["status"] == "success", res.stderr
    assert Path(r["outputs"]["results_dir"]).is_absolute()


# ──────────────────────────────────────────────
# Ignored task failures (the pipeline retries twice, then ignores)
# ──────────────────────────────────────────────

def test_command_asks_nextflow_for_a_trace_under_the_run_dir(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    cmd = (setup["out"] / "reproducibility" / "commands.sh").read_text()
    assert re.search(r"-with-trace \S*/nextflow_run/trace_\d{8}T\d{6}\.\d{6}Z(-\d+)?\.tsv", cmd), cmd


def test_clean_trace_means_no_ignored_failures(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    r = result_json(setup["out"])
    assert r["status"] == "success"
    assert r["ignored_failures"] == [] and r["trace"]["available"] is True
    assert r["trace"]["tasks"] >= 1


def test_ignored_failure_is_reported_but_the_run_still_exits_zero(setup):
    env = {**setup["env"], "FAKE_NF_TRACE": "ignored"}
    res = run_cli(base_args(setup, "--accept-bioskryb-terms"), env=env)
    r = result_json(setup["out"])
    assert r["status"] == "success_with_ignored_failures"
    assert r["ignored_failures"] == ["RNASEQ_WF:CELL_TYPING (CELL_TYPING)"]
    text = (setup["out"] / "report.md").read_text()
    assert "CELL_TYPING" in text and "ignored" in text.lower()
    assert "ignored" in res.stdout.lower()


def test_strict_turns_ignored_failures_into_a_nonzero_exit(setup):
    env = {**setup["env"], "FAKE_NF_TRACE": "ignored"}
    res = run_cli(base_args(setup, "--accept-bioskryb-terms", "--strict"), env=env, expect=1)
    assert "CELL_TYPING" in res.stderr
    assert result_json(setup["out"])["status"] == "success_with_ignored_failures"


def test_a_task_that_failed_then_succeeded_is_a_retry_not_a_failure(setup):
    env = {**setup["env"], "FAKE_NF_TRACE": "retry_ok"}
    run_cli(base_args(setup, "--accept-bioskryb-terms", "--strict"), env=env)
    r = result_json(setup["out"])
    assert r["status"] == "success" and r["ignored_failures"] == []
    assert r["trace"]["retried_tasks"] == ["RNASEQ_WF:SALMONQUANT (s2)"]


def test_missing_trace_is_said_not_assumed_clean(setup):
    env = {**setup["env"], "FAKE_NF_TRACE": "none"}
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=env)
    r = result_json(setup["out"])
    assert r["trace"]["available"] is False and r["ignored_failures"] is None
    assert "trace" in (setup["out"] / "report.md").read_text().lower()


def test_parse_trace_uses_the_final_state_of_each_task(tmp_path):
    mod = _load_module()
    p = tmp_path / "t.tsv"
    p.write_text("task_id\thash\tnative_id\tname\tstatus\texit\n"
                 "1\ta\t-\tA (x)\tFAILED\t1\n1\ta\t-\tA (x)\tCOMPLETED\t0\n"
                 "2\tb\t-\tB (y)\tCOMPLETED\t0\n2\tb\t-\tB (y)\tFAILED\t1\n"
                 "3\tc\t-\tC (z)\tCACHED\t0\n")
    t = mod.parse_trace(p)
    assert t["ignored_failures"] == ["B (y)"] and t["retried_tasks"] == ["A (x)"]
    assert t["tasks"] == 3


def test_parse_trace_on_a_missing_or_headerless_file_is_none(tmp_path):
    mod = _load_module()
    assert mod.parse_trace(tmp_path / "nope.tsv") is None
    p = tmp_path / "bad.tsv"
    p.write_text("not a trace\n")
    assert mod.parse_trace(p) is None


@settings(deadline=None, max_examples=80)
@given(seqs=st.dictionaries(
    st.from_regex(r"[A-Z]{1,6} \([a-z0-9]{1,4}\)", fullmatch=True),
    st.lists(st.sampled_from(["COMPLETED", "FAILED", "CACHED", "ABORTED"]), min_size=1, max_size=4),
    min_size=1, max_size=6))
def test_property_a_task_is_ignored_iff_its_last_state_failed(tmp_path_factory, seqs):
    mod = _load_module()
    p = tmp_path_factory.mktemp("tr") / "t.tsv"
    rows = ["task_id\thash\tnative_id\tname\tstatus\texit"]
    for name, states in seqs.items():
        rows += [f"1\tx\t-\t{name}\t{st_}\t0" for st_ in states]
    p.write_text("\n".join(rows) + "\n")
    t = mod.parse_trace(p)
    expect = sorted(n for n, s_ in seqs.items() if s_[-1] in ("FAILED", "ABORTED"))
    assert t["ignored_failures"] == expect
    assert t["tasks"] == len(seqs)


# ──────────────────────────────────────────────
# Through the real ClawBio runner (python clawbio.py run bj-expression ...)
# ──────────────────────────────────────────────
# The runner has its own `--profile` (a patient-profile path) and only hands it
# to a few pipeline skills. Without special handling, `--profile singularity`
# is swallowed and the wrapper silently falls back to Docker.

REPO_ROOT = Path(__file__).resolve().parents[3]


def run_via_runner(args, env):
    return subprocess.run(
        [sys.executable, str(REPO_ROOT / "clawbio.py"), "run", "bj-expression",
         *map(str, args)],
        capture_output=True, text=True, env=env, cwd=REPO_ROOT)


def test_runner_forwards_profile_and_every_wrapper_flag(tmp_path):
    fq = make_fastqs(tmp_path)
    sheet = write_csv(tmp_path / "in.csv", [("s1", *fq["s1"])])
    env = make_tools(tmp_path, docker=False, singularity=True)  # no docker: needs the profile
    cfg = tmp_path / "x.config"
    cfg.write_text("// extra\n")
    out = tmp_path / "out"
    res = run_via_runner(
        ["--input", sheet, "--output", out, "--check", "--profile", "singularity",
         "--pipeline-dir", make_pipeline(tmp_path), "--genomes-base", tmp_path / "refs",
         "--extra-config", cfg, "--max-cpus", "8", "--max-memory", "60.GB",
         "--nextflow-version", "24.10.5", "--accept-bioskryb-terms", "--reads-per-file", "5"], env)
    assert (out / "result.json").exists(), res.stdout + res.stderr
    r = result_json(out)
    assert r["status"] == "check_passed"
    assert r["profile"] == "singularity", "the runner swallowed --profile"
    assert r["terms_accepted"] is True
    assert r["nextflow_version_requested"] == "24.10.5"
    assert r["extra_config"].endswith("extra.config")
    p = json.loads((out / "reproducibility" / "params.json").read_text())
    assert p["max_cpus"] == 8 and p["max_memory"] == "60.GB"
    assert p["genomes_base"] == str(tmp_path / "refs")


def test_a_resumed_run_reports_its_own_trace_not_the_previous_runs(setup):
    """Real Nextflow refuses to overwrite an existing trace file, so a fixed trace
    name made a resume report the FIRST run's failures (and miss its own)."""
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])        # clean first run
    assert result_json(setup["out"])["ignored_failures"] == []
    env = {**setup["env"], "FAKE_NF_TRACE": "ignored"}
    run_cli(base_args(setup, "--accept-bioskryb-terms", "--resume"), env=env)      # resume with a failure
    r = result_json(setup["out"])
    assert r["status"] == "success_with_ignored_failures"
    assert r["ignored_failures"] == ["RNASEQ_WF:CELL_TYPING (CELL_TYPING)"]


def test_each_invocation_gets_its_own_trace_file_and_records_the_path(setup):
    import time
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    first = result_json(setup["out"])["trace"]["path"]
    time.sleep(1.1)  # the file name carries a one-second timestamp
    run_cli(base_args(setup, "--accept-bioskryb-terms", "--resume"), env=setup["env"])
    second = result_json(setup["out"])["trace"]["path"]
    assert first != second and Path(first).exists() and Path(second).exists()


def test_a_resume_whose_trace_was_not_written_says_unknown_not_the_old_result(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])        # trace written, clean
    env = {**setup["env"], "FAKE_NF_TRACE": "none"}
    run_cli(base_args(setup, "--accept-bioskryb-terms", "--resume"), env=env)      # no trace this time
    r = result_json(setup["out"])
    assert r["trace"]["available"] is False and r["ignored_failures"] is None


# ──────────────────────────────────────────────
# Real files instead of symlinks in the results folder
# ──────────────────────────────────────────────
# BioSkryb's publishDir directives set no mode, so Nextflow publishes symlinks into
# work/ (a config override does not change that). After a successful run the wrapper
# replaces each link with a verified real copy; the targets in work/ are never touched.

import os as _os


def links_under(root: Path) -> list[Path]:
    return [Path(r) / n for r, d, f in _os.walk(root) for n in d + f if (Path(r) / n).is_symlink()]


def _results(out: Path) -> Path:
    return out / "results_clawbio"


def test_published_symlinks_become_real_files_after_a_run(setup):
    env = {**setup["env"], "FAKE_NF_PUBLISH": "symlink"}
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=env)
    res = _results(setup["out"])
    assert (res / "multiqc" / "multiqc_report.html").is_file()
    assert links_under(res) == [], "symlinks were left in the results folder"
    assert not (res / "multiqc" / "multiqc_data").is_symlink()
    assert (res / "multiqc" / "multiqc_data" / "x.txt").read_text() == "data\n"
    m = result_json(setup["out"])["results_materialised"]
    assert m["files"] >= 8 and m["dirs"] == 1 and m["broken"] == [] and m["failed"] == []
    assert result_json(setup["out"])["status"] == "success"


def test_keep_symlinks_leaves_the_links_and_says_so(setup):
    env = {**setup["env"], "FAKE_NF_PUBLISH": "symlink"}
    run_cli(base_args(setup, "--accept-bioskryb-terms", "--keep-symlinks"), env=env)
    assert len(links_under(_results(setup["out"]))) > 5
    assert result_json(setup["out"])["results_materialised"] is None


def test_the_work_directory_files_are_never_modified(setup):
    env = {**setup["env"], "FAKE_NF_PUBLISH": "symlink"}
    run_cli(base_args(setup, "--accept-bioskryb-terms", "--keep-symlinks"), env=env)
    work = setup["out"] / "nextflow_run" / "fake_work"
    before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in work.rglob("*") if p.is_file()}
    run_cli(["--output", setup["out"], "--materialise-only"], env=setup["env"])
    after = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in work.rglob("*") if p.is_file()}
    assert before == after and before


def test_a_failed_nextflow_run_is_not_materialised(setup):
    env = {**setup["env"], "FAKE_NF_MODE": "fail", "FAKE_NF_PUBLISH": "symlink"}
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=env, expect=1)
    assert result_json(setup["out"])["results_materialised"] is None


def test_materialise_only_fixes_an_existing_run_and_leaves_result_json_alone(setup):
    env = {**setup["env"], "FAKE_NF_PUBLISH": "symlink"}
    run_cli(base_args(setup, "--accept-bioskryb-terms", "--keep-symlinks"), env=env)
    before = (setup["out"] / "result.json").read_bytes()
    res = run_cli(["--output", setup["out"], "--materialise-only"], env=setup["env"])
    assert links_under(_results(setup["out"])) == []
    assert (setup["out"] / "result.json").read_bytes() == before
    rep = json.loads((setup["out"] / "materialise_report.json").read_text())
    assert rep["failed"] == [] and rep["broken"] == [] and rep["files"] >= 8
    assert "materialised" in res.stdout.lower()


def test_materialise_only_on_a_missing_results_dir_is_a_clear_error(tmp_path):
    res = run_cli(["--output", tmp_path / "nothing", "--materialise-only"], expect=None)
    assert res.returncode == 1 and "results" in res.stderr.lower()


def test_materialise_is_idempotent(tmp_path):
    mod = _load_module()
    src = tmp_path / "work" / "a.txt"
    src.parent.mkdir()
    src.write_text("hello")
    res = tmp_path / "res"
    res.mkdir()
    (res / "a.txt").symlink_to(src)
    first = mod.materialise_results(res)
    second = mod.materialise_results(res)
    assert first["files"] == 1 and second["files"] == 0 and second["already_real"] == 1
    assert (res / "a.txt").read_text() == "hello" and not (res / "a.txt").is_symlink()


def test_a_broken_link_is_reported_and_does_not_stop_the_rest(tmp_path):
    mod = _load_module()
    good = tmp_path / "work" / "good.txt"
    good.parent.mkdir()
    good.write_text("ok")
    res = tmp_path / "res"
    res.mkdir()
    (res / "good.txt").symlink_to(good)
    (res / "gone.txt").symlink_to(tmp_path / "work" / "missing.txt")
    rep = mod.materialise_results(res)
    assert rep["broken"] == ["gone.txt"] and rep["files"] == 1
    assert (res / "good.txt").read_text() == "ok" and (res / "gone.txt").is_symlink()


def test_a_copy_that_does_not_match_leaves_the_link_in_place(tmp_path, monkeypatch):
    mod = _load_module()
    src = tmp_path / "work" / "a.txt"
    src.parent.mkdir()
    src.write_text("real contents")
    res = tmp_path / "res"
    res.mkdir()
    (res / "a.txt").symlink_to(src)
    monkeypatch.setattr(mod.shutil, "copy2",
                        lambda s_, d_, **k: Path(d_).write_text("CORRUPTED"))
    rep = mod.materialise_results(res)
    assert rep["failed"] == ["a.txt"] and (res / "a.txt").is_symlink()
    assert not list(res.glob("*.materialise-tmp")), "temp copy was not cleaned up"
    assert src.read_text() == "real contents"


@settings(deadline=None, max_examples=40)
@given(tree=st.dictionaries(
    st.from_regex(r"[a-z]{1,4}(/[a-z]{1,4}){0,2}", fullmatch=True),
    st.binary(max_size=30), min_size=1, max_size=6),
    link_every=st.integers(1, 3))
def test_property_after_materialising_the_tree_reads_the_same_with_no_links(
        tmp_path_factory, tree, link_every):
    """Whatever mix of file and directory links, content is unchanged and none remain."""
    mod = _load_module()
    base = tmp_path_factory.mktemp("prop")
    work, res = base / "work", base / "res"
    work.mkdir()
    res.mkdir()
    # keep only paths that do not collide as file-vs-directory
    kept = {}
    for rel in sorted(tree):
        parts = rel.split("/")
        if any(rel == k or k.startswith(rel + "/") or rel.startswith(k + "/") for k in kept):
            continue
        kept[rel] = tree[rel]
    for i, (rel, data) in enumerate(kept.items()):
        w = work / rel
        w.parent.mkdir(parents=True, exist_ok=True)
        w.write_bytes(data)
        dest = res / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        if i % link_every == 0:
            dest.symlink_to(w)
        else:
            dest.write_bytes(data)
    expected = {rel: data for rel, data in kept.items()}
    mod.materialise_results(res)
    assert links_under(res) == []
    got = {p.relative_to(res).as_posix(): p.read_bytes() for p in res.rglob("*") if p.is_file()}
    assert got == expected


def test_trace_path_is_unique_even_within_the_same_instant(tmp_path):
    """Two invocations in the same second used to collide, and Nextflow will not
    overwrite an existing trace file (the original stale-trace bug, again)."""
    from datetime import datetime, timezone
    mod = _load_module()
    now = datetime(2026, 11, 5, 18, 5, 50, 123456, tzinfo=timezone.utc)
    first = mod.new_trace_path(tmp_path, now)
    first.write_text("x")
    second = mod.new_trace_path(tmp_path, now)
    assert first != second and not second.exists()
    assert first.parent == second.parent == tmp_path
    assert re.fullmatch(r"trace_20261105T180550\.123456Z\.tsv", first.name)


def test_trace_path_includes_microseconds_so_back_to_back_runs_differ(tmp_path):
    mod = _load_module()
    a, b = mod.new_trace_path(tmp_path), mod.new_trace_path(tmp_path)
    assert a != b or not a.exists()  # never hands out a path that already exists


# ──────────────────────────────────────────────
# Reads per file: counted up front, chosen explicitly, translated to the pipeline's units
# ──────────────────────────────────────────────
# bj-expression's SEQTK_SAMPLE divides n_reads by the reads of R1 AND R2 together, so the
# default n_reads = 100000 gave ~50,000 reads PER FILE. It only trims to read_length when it
# actually samples (otherwise it copies the files), and it hard-codes -s123 for read 2, so a
# different seed would desynchronise the mates. The wrapper therefore counts reads first,
# requires an explicit per-file choice, and refuses any seed but 123.

import gzip as _gzip


def gz_fastq(path: Path, n: int, mate: int = 1) -> Path:
    text = "".join(f"@r{i}/{mate}\nACGTACGTAC\n+\nIIIIIIIIII\n" for i in range(n))
    with _gzip.open(path, "wt") as fh:
        fh.write(text)
    return path


def counted_setup(tmp_path, counts):
    """counts = {sample: (r1_reads, r2_reads)} -> sheet + working fakes, like `setup`."""
    d = tmp_path / "cfq"
    d.mkdir()
    rows = []
    for name, (n1, n2) in counts.items():
        r1 = gz_fastq(d / f"{name}_R1.fastq.gz", n1, 1)
        r2 = gz_fastq(d / f"{name}_R2.fastq.gz", n2, 2)
        rows.append((name, r1, r2))
    return {"tmp": tmp_path, "csv": write_csv(tmp_path / "c.csv", rows),
            "pipeline": make_pipeline(tmp_path), "env": make_tools(tmp_path),
            "out": tmp_path / "out"}


def test_count_fastq_reads_on_gz_and_plain_files(tmp_path):
    mod = _load_module()
    assert mod.count_fastq_reads(gz_fastq(tmp_path / "a.fastq.gz", 7)) == 7
    plain = tmp_path / "b.fastq"
    plain.write_text("@r\nAC\n+\nII\n" * 3)
    assert mod.count_fastq_reads(plain) == 3
    assert mod.count_fastq_reads(gz_fastq(tmp_path / "empty.fastq.gz", 0)) == 0


def test_count_fastq_reads_rejects_files_that_are_not_fastq(tmp_path):
    mod = _load_module()
    bad = tmp_path / "bad.fastq.gz"
    bad.write_bytes(b"not gzip at all")
    with pytest.raises(mod.WrapperError, match="bad.fastq.gz"):
        mod.count_fastq_reads(bad)
    trunc = tmp_path / "trunc.fastq"
    trunc.write_text("@r\nAC\n+\n")  # 3 lines: not a multiple of 4
    with pytest.raises(mod.WrapperError, match="trunc.fastq"):
        mod.count_fastq_reads(trunc)


@settings(deadline=None, max_examples=40)
@given(n=st.integers(0, 60), gz=st.booleans(), final_newline=st.booleans())
def test_property_read_count_matches_the_records_written(tmp_path_factory, n, gz, final_newline):
    mod = _load_module()
    d = tmp_path_factory.mktemp("cnt")
    text = "".join(f"@r{i}\nACGT\n+\nIIII" + ("\n" if (i < n - 1 or final_newline) else "") for i in range(n))
    p = d / ("x.fastq.gz" if gz else "x.fastq")
    (_gzip.open(p, "wt") if gz else open(p, "w")).write(text)
    assert mod.count_fastq_reads(p) == n


def test_count_reads_mode_prints_the_table_and_runs_nothing(tmp_path):
    s = counted_setup(tmp_path, {"s1": (100, 100), "s2": (250, 250)})
    res = run_cli(["--input", s["csv"], "--count-reads"], env=s["env"])
    assert "s1" in res.stdout and "100" in res.stdout and "s2" in res.stdout and "250" in res.stdout
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("mode", [[], ["--check"]])
def test_a_run_without_a_read_choice_stops_and_shows_what_is_available(tmp_path, mode):
    s = counted_setup(tmp_path, {"s1": (120, 120), "s2": (80, 80)})
    # built by hand: base_args() would add a default read choice and defeat the point
    args = ["--input", s["csv"], "--output", s["out"], "--pipeline-dir", s["pipeline"],
            "--accept-bioskryb-terms", *mode]
    res = run_cli(args, env=s["env"], expect=None)
    assert res.returncode == 1
    assert "--reads-per-file" in res.stderr and "all" in res.stderr
    assert "s1" in res.stderr and "120" in res.stderr and "80" in res.stderr
    assert not (s["out"] / "result.json").exists()


def test_the_demo_does_not_need_a_read_choice(tmp_path):
    run_cli(["--demo", "--output", tmp_path / "d"], env={"PATH": str(tmp_path)})


def test_reads_per_file_is_translated_to_the_pipelines_combined_count(tmp_path):
    s = counted_setup(tmp_path, {"s1": (100, 100), "s2": (100, 100)})
    run_cli(base_args(s, "--accept-bioskryb-terms", "--reads-per-file", "40"), env=s["env"])
    p = json.loads((s["out"] / "reproducibility" / "params.json").read_text())
    assert p["n_reads"] == 80, "the pipeline divides n_reads by R1+R2 together, so per-file x 2"
    assert "skip_subsampling" not in p
    sub = result_json(s["out"])["subsampling"]
    assert sub["mode"] == "reads_per_file" and sub["reads_per_file"] == 40
    assert sub["pipeline_n_reads"] == 80


def test_all_means_no_subsampling_and_says_reads_are_not_trimmed(tmp_path):
    s = counted_setup(tmp_path, {"s1": (100, 100)})
    run_cli(base_args(s, "--accept-bioskryb-terms", "--reads-per-file", "all"), env=s["env"])
    p = json.loads((s["out"] / "reproducibility" / "params.json").read_text())
    assert p["skip_subsampling"] is True and "n_reads" not in p
    sub = result_json(s["out"])["subsampling"]
    assert sub["mode"] == "all" and sub["trim_to_bp"] is None
    assert "not trimmed" in (s["out"] / "report.md").read_text().lower()


def test_skip_subsampling_alone_counts_as_choosing_all(tmp_path):
    s = counted_setup(tmp_path, {"s1": (50, 50)})
    run_cli(base_args(s, "--accept-bioskryb-terms", "--skip-subsampling"), env=s["env"])
    assert result_json(s["out"])["subsampling"]["mode"] == "all"


def test_a_number_together_with_skip_subsampling_is_a_conflict(tmp_path):
    s = counted_setup(tmp_path, {"s1": (50, 50)})
    res = run_cli(base_args(s, "--accept-bioskryb-terms", "--reads-per-file", "10",
                            "--skip-subsampling"), env=s["env"], expect=None)
    assert res.returncode in (1, 2) and "error:" in res.stderr
    assert "unrecognized" not in res.stderr, "must be the wrapper's own conflict message"
    assert "--skip-subsampling" in res.stderr and "--reads-per-file" in res.stderr


def test_the_old_n_reads_flag_is_refused_with_the_reason(tmp_path):
    s = counted_setup(tmp_path, {"s1": (50, 50)})
    res = run_cli(base_args(s, "--accept-bioskryb-terms", "--n-reads", "100"),
                  env=s["env"], expect=None)
    assert res.returncode in (1, 2) and "--reads-per-file" in res.stderr
    assert "R1" in res.stderr  # explains that n_reads counts both files together


def test_only_seed_123_is_allowed_because_read_2_is_hardcoded(tmp_path):
    s = counted_setup(tmp_path, {"s1": (50, 50)})
    bad = run_cli(base_args(s, "--accept-bioskryb-terms", "--reads-per-file", "10", "--seed", "7"),
                  env=s["env"], expect=None)
    assert bad.returncode in (1, 2) and "123" in bad.stderr and "read 2" in bad.stderr.lower()
    assert not (s["out"] / "result.json").exists()
    run_cli(base_args(s, "--accept-bioskryb-terms", "--reads-per-file", "10", "--seed", "123"),
            env=s["env"])


def test_sampling_state_and_expected_reads_are_recorded_per_sample(tmp_path):
    s = counted_setup(tmp_path, {"big": (100, 100), "small": (30, 30)})
    run_cli(base_args(s, "--accept-bioskryb-terms", "--reads-per-file", "50"), env=s["env"])
    sub = result_json(s["out"])["subsampling"]
    assert sub["samples"]["big"] == {"r1": 100, "r2": 100, "pairs": 100,
                                     "sampling_active": True, "expected_reads_per_file": 50}
    assert sub["samples"]["small"]["sampling_active"] is False
    assert sub["samples"]["small"]["expected_reads_per_file"] == 30
    assert any("small" in w and "not trimmed" in w.lower() for w in sub["warnings"])


def test_a_sample_the_pipeline_would_drop_is_warned_about(tmp_path):
    # min_reads (5000) is compared with R1+R2 together: 2 x 2000 = 4000 < 5000
    s = counted_setup(tmp_path, {"tiny": (2000, 2000), "fine": (4000, 4000)})
    run_cli(base_args(s, "--accept-bioskryb-terms", "--reads-per-file", "all", "--min-reads", "5000"),
            env=s["env"])
    warnings = " ".join(result_json(s["out"])["subsampling"]["warnings"])
    assert "tiny" in warnings and "min_reads" in warnings
    assert "fine" not in warnings.replace("fine-", "")


def test_mates_with_different_read_counts_are_refused(tmp_path):
    s = counted_setup(tmp_path, {"s1": (100, 90)})
    res = run_cli(base_args(s, "--accept-bioskryb-terms", "--reads-per-file", "10"),
                  env=s["env"], expect=None)
    assert res.returncode == 1 and "s1" in res.stderr and "100" in res.stderr and "90" in res.stderr


def test_the_trim_length_comes_from_the_pinned_pipeline_config(tmp_path):
    s = counted_setup(tmp_path, {"s1": (100, 100)})
    (s["pipeline"] / "nextflow.config").write_text("manifest { version = '2.0.0' }\nparams {\n  read_length = 50\n}\n")
    subprocess.run([*GIT, "-C", str(s["pipeline"]), "commit", "-q", "-am", "cfg"], check=True)
    subprocess.run([*GIT, "-C", str(s["pipeline"]), "tag", "-f", PIPELINE_TAG], check=True)
    run_cli(base_args(s, "--accept-bioskryb-terms", "--reads-per-file", "10"), env=s["env"])
    assert result_json(s["out"])["subsampling"]["trim_to_bp"] == 50
    s["out"] = tmp_path / "out2"
    run_cli(base_args(s, "--accept-bioskryb-terms", "--reads-per-file", "10", "--read-length", "40"),
            env=s["env"])
    assert result_json(s["out"])["subsampling"]["trim_to_bp"] == 40


# ──────────────────────────────────────────────
# Single-nucleus version: intron primary, exon duplicates in <results>_exon
# ──────────────────────────────────────────────

def test_a_missing_exon_duplicate_makes_the_run_incomplete(setup):
    env = {**setup["env"], "FAKE_NF_DROP_EXON": "multiqc"}
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=env, expect=1)
    r = result_json(setup["out"])
    assert r["status"] == "incomplete"
    assert any("exon" in m and "multiqc" in m for m in r["missing_outputs"])
    assert r["exon_duplicate"]["missing_outputs"]


def test_exon_duplicate_folder_is_the_results_folder_with_an_exon_suffix(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    r = result_json(setup["out"])
    assert r["exon_duplicate"]["results_dir"] == str(setup["out"] / "results_exon_clawbio")
    assert (setup["out"] / "results_exon_clawbio" / "multiqc" / "multiqc_report.html").exists()


def test_the_report_says_which_counts_drive_multiqc_and_cell_typing(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    text = (setup["out"] / "report.md").read_text().lower()
    assert "intron-inclusive" in text and "_exon" in text


# ──────────────────────────────────────────────
# Remove work/ after a verified run (default); --keep-work turns it off
# ──────────────────────────────────────────────

def test_a_successful_run_removes_work_by_default_and_keeps_real_results(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    r = result_json(setup["out"])
    assert r["status"] == "success"
    assert not (setup["out"] / "work").exists(), "work/ should be removed by default"
    assert r["work_removed"] is True
    assert (setup["out"] / "results_clawbio" / "multiqc" / "multiqc_report.html").is_file()
    assert links_under(setup["out"] / "results_clawbio") == []


def test_keep_work_leaves_work_in_place(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms", "--keep-work"), env=setup["env"])
    r = result_json(setup["out"])
    assert (setup["out"] / "work").is_dir()
    assert r["work_removed"] is False and "keep-work" in r["work_cleanup_reason"]


def test_a_failed_run_keeps_work_so_it_can_be_resumed(setup):
    env = {**setup["env"], "FAKE_NF_MODE": "fail"}
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=env, expect=1)
    r = result_json(setup["out"])
    assert (setup["out"] / "work").is_dir()
    assert r["work_removed"] is False and "not successful" in r["work_cleanup_reason"]


def test_keep_symlinks_keeps_work_because_the_links_still_point_into_it(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms", "--keep-symlinks"), env=setup["env"])
    r = result_json(setup["out"])
    assert (setup["out"] / "work").is_dir()
    assert r["work_removed"] is False and "not materialised" in r["work_cleanup_reason"]


def test_cleanup_refuses_when_any_result_is_not_a_real_file(tmp_path):
    mod = _load_module()
    work = tmp_path / "work"
    work.mkdir()
    (work / "x").write_text("x")
    ok = {"files": 3, "broken": [], "failed": []}
    bad = {"files": 3, "broken": ["a.html"], "failed": []}
    assert mod.cleanup_work(work, ok, ok)[0] is True and not work.exists()
    work.mkdir(); (work / "x").write_text("x")
    removed, reason = mod.cleanup_work(work, ok, bad)
    assert removed is False and "not real files" in reason and work.exists()
    assert mod.cleanup_work(work, None, ok)[0] is False


def test_the_report_says_when_work_was_removed_and_that_resume_is_gone(setup):
    run_cli(base_args(setup, "--accept-bioskryb-terms"), env=setup["env"])
    text = (setup["out"] / "report.md").read_text().lower()
    assert "work directory was removed" in text and "--resume" in text
