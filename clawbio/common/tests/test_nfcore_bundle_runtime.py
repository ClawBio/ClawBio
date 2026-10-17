"""Relocation mechanics and emitted single-file standard-library bundles."""

import csv
import hashlib
import importlib
import re
import subprocess
import sys
from pathlib import Path

import pytest
import yaml
from hypothesis import HealthCheck, example, given, settings
from hypothesis import strategies as st


def test_csv_preview_backup_and_prefix_boundary(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    path = tmp_path / "samplesheet.csv"
    original = "sample,fastq_1\nA,/data/a.fq\nB,/data2/b.fq\nC,s3://bucket/c.fq\n"
    path.write_text(original)
    changes = api.remap_csv(path, "/data", "/moved", columns=("fastq_1",), dry_run=True)
    assert changes == [("fastq_1", "/data/a.fq", "/moved/a.fq")]
    assert path.read_text() == original
    assert not (tmp_path / "samplesheet.csv.bak").exists()
    api.remap_csv(path, "/data", "/moved", columns=("fastq_1",), dry_run=False)
    assert (tmp_path / "samplesheet.csv.bak").read_text() == original
    assert path.read_text() == original.replace("/data/a.fq", "/moved/a.fq")


@settings(deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(st.text(alphabet="abc é-_012", min_size=1, max_size=20))
def test_remap_round_trip(tmp_path, suffix):
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    path = tmp_path / "input.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["fastq_1"])
        writer.writerow(["/old/" + suffix])
    original = path.read_bytes()
    api.remap_csv(path, "/old", "/new", columns=["fastq_1"], dry_run=False)
    api.remap_csv(path, "/new", "/old", columns=["fastq_1"], dry_run=False)
    assert path.read_bytes() == original


def test_all_actions_run_and_first_failure_wins():
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    seen = []

    def action(code):
        def run():
            seen.append(code)
            return code

        return run

    assert api.run_actions([action(0), action(2), action(1), action(0)]) == 2
    assert seen == [0, 2, 1, 0]


@pytest.mark.parametrize("pipeline", ["rnaseq", "sarek", "scrnaseq"])
def test_emitted_remapper_runs_without_site_packages(pipeline, tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.bundle_emit")
    root = Path(__file__).resolve().parents[3]
    adapter = root / "skills" / f"nfcore-{pipeline}-wrapper/remap_paths.py"
    dest = tmp_path / "remap_paths.py"
    api.write_standalone_remapper(dest, adapter_path=adapter)
    result = subprocess.run(
        [sys.executable, "-I", "-S", str(dest), "--help"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "--repair-bundle" in result.stdout
    path = tmp_path / "samplesheet.valid.csv"
    path.write_text("sample,fastq_1\nA,/old/a.fq.gz\n")
    result = subprocess.run(
        [
            sys.executable,
            "-I",
            "-S",
            str(dest),
            "--old",
            "/old",
            "--new",
            "/moved",
            "--dry-run",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert path.read_text() == "sample,fastq_1\nA,/old/a.fq.gz\n"


_BLOCK = "# >>> SHARED NFCORE RUNTIME >>>\npass\n# <<< SHARED NFCORE RUNTIME <<<\n"


def test_emitter_rejects_missing_and_duplicate_markers():
    api = importlib.import_module("clawbio.common.nfcore.bundle_emit")
    with pytest.raises(ValueError):
        api.render_standalone_remapper("print(1)")
    with pytest.raises(ValueError):
        api.render_standalone_remapper(_BLOCK + _BLOCK)


@pytest.mark.parametrize(
    "adapter",
    [
        # After the block: the runtime's remap_csv would call this, not its own helper.
        _BLOCK + "def prefix_matches(path, prefix):\n    return False\n",
        _BLOCK + "escape_double_quoted = str\n",
        # Before the block: the pasted runtime would silently replace the adapter's.
        "def repair_bundle():\n    return 0\n" + _BLOCK,
    ],
)
def test_emitter_refuses_names_the_pasted_runtime_would_shadow(adapter):
    """The runtime is pasted into the adapter's own namespace, unlike the in-repo
    import. A colliding top-level name would make the shipped script behave
    differently from the tested one, so it must fail at generation time."""
    api = importlib.import_module("clawbio.common.nfcore.bundle_emit")
    with pytest.raises(ValueError, match="shadow"):
        api.render_standalone_remapper(adapter)


def test_emitter_allows_identical_stdlib_imports_and_an_adapter_entry_point():
    api = importlib.import_module("clawbio.common.nfcore.bundle_emit")
    adapter = (
        "import csv\nfrom pathlib import Path\n" + _BLOCK
        + "import re\n\ndef remap_csv(path):\n    return _nfcore_remap_csv(path)\n"
    )
    assert "def remap_csv(path):" in api.render_standalone_remapper(adapter)


@pytest.mark.parametrize("pipeline", ["rnaseq", "sarek", "scrnaseq"])
def test_standalone_remapper_mutates_verifies_and_repairs_bundle(pipeline, tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.bundle_emit")
    root = Path(__file__).resolve().parents[3]
    adapter = root / "skills" / f"nfcore-{pipeline}-wrapper/remap_paths.py"
    bundle = tmp_path / "relocated é bundle" / "reproducibility"
    dest = bundle / "remap_paths.py"
    api.write_standalone_remapper(dest, adapter_path=adapter)
    reads = tmp_path / "moved é reads"
    reads.mkdir()
    (reads / "a.fq.gz").write_bytes(b"public test read")
    samplesheet = bundle / "samplesheet.valid.csv"
    original = b"sample,fastq_1\nA,/old/a.fq.gz\n"
    samplesheet.write_bytes(original)
    (bundle / "commands.sh").write_text("#!/bin/bash\ntrue\n")
    (bundle / "params.yaml").write_text(
        "input: reproducibility/samplesheet.valid.csv\n"
    )
    (bundle.parent / "report.md").write_text("Synthetic test report\n")
    (bundle.parent / "upstream" / "results").mkdir(parents=True)
    (bundle.parent / "upstream" / "results" / "calls.vcf").write_text("synthetic\n")

    def run(*args):
        result = subprocess.run(
            [sys.executable, "-I", "-S", str(dest), *args],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, result.stdout + result.stderr

    run("--old", "/old", "--new", str(reads))
    backup = (
        samplesheet.with_suffix(".bak")
        if pipeline == "scrnaseq"
        else samplesheet.with_name(samplesheet.name + ".bak")
    )
    assert backup.read_bytes() == original
    assert list(csv.DictReader(samplesheet.read_text().splitlines())) == [
        {"sample": "A", "fastq_1": str(reads / "a.fq.gz")}
    ]
    run("--repair-bundle", "--verify")
    required = ("manifest.json", "environment.yml", "checksums.sha256")
    repaired = {name: (bundle / name).read_bytes() for name in required}
    assert b'"regenerated_post_hoc": true' in repaired["manifest.json"]
    for line in repaired["checksums.sha256"].decode().splitlines():
        digest, relative = line.split("  ", 1)
        assert (
            hashlib.sha256((bundle.parent / relative).read_bytes()).hexdigest()
            == digest
        )
    assert repaired["checksums.sha256"]
    run("--repair-bundle", "--verify")
    assert {name: (bundle / name).read_bytes() for name in required} == repaired


def test_repair_bundle_recomputes_checksum_after_new_metadata(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    (tmp_path / "commands.sh").write_text("#!/bin/bash\ntrue\n")
    (tmp_path / "checksums.sha256").write_text("old\n")

    def manifest(directory, commands):
        (directory / "manifest.json").write_text('{"regenerated_post_hoc":true}')

    def environment(directory):
        (directory / "environment.yml").write_text("name: recipe\n")

    def checksums(directory):
        assert (directory / "manifest.json").exists()
        assert (directory / "environment.yml").exists()
        (directory / "checksums.sha256").write_text("updated\n")

    assert (
        api.repair_bundle(
            tmp_path,
            required_files=("manifest.json", "environment.yml", "checksums.sha256"),
            commands_finder=lambda **kw: tmp_path / "commands.sh",
            regenerate_manifest=manifest,
            regenerate_environment=environment,
            regenerate_checksums=checksums,
        )
        == 0
    )
    assert (tmp_path / "checksums.sha256").read_text() == "updated\n"
    assert (
        api.repair_bundle(
            tmp_path,
            required_files=("manifest.json", "environment.yml", "checksums.sha256"),
            commands_finder=lambda **kw: tmp_path / "commands.sh",
            regenerate_manifest=manifest,
            regenerate_environment=environment,
            regenerate_checksums=lambda _: pytest.fail(
                "complete bundle must be untouched"
            ),
        )
        == 0
    )


def test_portable_shell_value_preserves_literal_metacharacters():
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    value = '/a b/é/$value/`literal`/"quoted"'
    rendered = api.format_output_value(value, old_value='"old"')
    result = subprocess.run(
        ["bash", "-c", "printf '%s' " + rendered],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == value


# Paths as users name them: spaces, apostrophes, quotes and shell metacharacters.
# Everything a POSIX path may hold on one line: NUL is impossible and \n/\r are
# rejected outright, but tabs, NEL and U+2028 are legal file-name characters.
_PATHS = st.text(
    alphabet=st.characters(blacklist_categories=("Cs",), blacklist_characters="\x00\n\r")
    | st.sampled_from(list(" '\"$`\\;&|()<>*?#~!{}[]:")),
    max_size=24,
)


@given(path=_PATHS, old=st.sampled_from(['"old"', "'old'", "old"]))
@example(path="/Users/Dana's Drive/run 1", old="'old'")
@example(path='/a"b', old="old")
@example(path="/data/run;touch x", old="old")
@example(path="", old="old")
def test_shell_value_round_trips_through_bash_quoting(path, old):
    """What format_output_value writes, the remappers must read back unchanged,
    or a second relocation silently finds no --output and leaves the old path."""
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    word = api.format_output_value(path, old_value=old)
    assert re.fullmatch(api.SHELL_WORD, word)
    assert api.unquote_shell_word(word) == path


@pytest.mark.parametrize(
    "path", ["/Users/Dana's Drive/run 1", '/x/"q" $HOME `id`', "/plain/out", ""]
)
def test_bash_assigns_exactly_the_formatted_value(path):
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    for old in ('"old"', "'old'", "old"):
        word = api.format_output_value(path, old_value=old)
        result = subprocess.run(
            ["bash", "-c", "printf '%s' " + word],
            capture_output=True,
            text=True,
            check=False,
        )
        assert (result.returncode, result.stdout) == (0, path), (old, word, result.stderr)


def test_line_breaks_cannot_be_written_to_a_command_line():
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    with pytest.raises(ValueError, match="line break"):
        api.format_output_value("/a\nb", old_value="old")


@given(path=_PATHS)
@example(path="/data/a #1/ref.fa")
@example(path="/data/it's/ref: fa")
@example(path="true")
@example(path="")
@example(path="/a\x85b\u2028c\tq\ufeff")
def test_yaml_reference_values_round_trip(path):
    """params.yaml is replayed by Nextflow's YAML parser: a remapped reference must
    load as exactly the path written, never truncated at ' #' or retyped."""
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    scalar = api.format_yaml_scalar(path)
    assert yaml.safe_load(f"fasta: {scalar}\n") == {"fasta": path}
    quote = scalar[0] if scalar[:1] in ("'", '"') else ""
    raw = scalar[1:-1] if quote else scalar
    assert api.yaml_scalar_value(quote, raw) == path


_PARAMS_REF = re.compile(
    r"^(?P<prefix>[ \t]*(?:fasta_fai|fasta)[ \t]*:[ \t]+)(?P<q>[\"']?)(?P<value>.*?)(?P=q)(?P<trail>[ \t]*)$",
    re.MULTILINE,
)


@settings(deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(tail=_PATHS, new=st.sampled_from(["/moved", "/moved #2", "/Dana's: refs"]))
def test_yaml_remap_rewrites_only_matching_prefixes_and_stays_loadable(tmp_path, tail, new):
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    params = tmp_path / "params.yaml"
    original = {"fasta": "/old/" + tail, "fasta_fai": "/elsewhere/x.fai", "outdir": "/old/out"}
    serialize = importlib.import_module("clawbio.common.nfcore.parameters").serialize_params_yaml
    params.write_text(serialize(original, pipeline_name="nf-core/x"), encoding="utf-8")
    changes = api.remap_yaml_references(params, "/old", new, dry_run=False, pattern=_PARAMS_REF)
    assert changes == [("fasta", "/old/" + tail, new + "/" + tail)]
    reloaded = yaml.safe_load(params.read_text())
    assert reloaded == {**original, "fasta": new + "/" + tail}
    assert params.read_text(encoding="utf-8").startswith("# 'input' and 'outdir'")


@settings(deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(tail=_PATHS)
@example(tail="a\u2028b\x85c.fq.gz")
def test_csv_remap_keeps_rows_whole(tmp_path, tail):
    """Python's str.splitlines() also breaks on NEL and U+2028, which may appear in
    a file name; the samplesheet must be read as CSV, not as text lines."""
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    path = tmp_path / "samplesheet.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        csv.writer(handle, lineterminator="\n").writerows(
            [["sample", "fastq_1"], ["A", "/old/" + tail], ["B", "/keep/b.fq"]]
        )
    api.remap_csv(path, "/old", "/new", columns=("fastq_1",), dry_run=False)
    with path.open(newline="", encoding="utf-8") as handle:
        assert list(csv.reader(handle)) == [
            ["sample", "fastq_1"], ["A", "/new/" + tail], ["B", "/keep/b.fq"]
        ]


def test_repair_leaves_an_allowlist_checksum_file_alone_when_only_stubs_are_written(tmp_path):
    """scrnaseq's checksums list a fixed allowlist that never includes the stubs, so
    writing a stub must not trigger a recompute that could overwrite a valid file."""
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    (tmp_path / "commands.sh").write_text("#!/bin/bash\ntrue\n")
    (tmp_path / "checksums.sha256").write_text("original\n")
    assert (
        api.repair_bundle(
            tmp_path,
            required_files=("manifest.json", "environment.yml", "checksums.sha256"),
            commands_finder=lambda **kw: tmp_path / "commands.sh",
            regenerate_manifest=lambda d, c: (d / "manifest.json").write_text("{}"),
            regenerate_environment=lambda d: (d / "environment.yml").write_text("x: 1\n"),
            regenerate_checksums=lambda _: pytest.fail("allowlist checksums must stand"),
            checksums_cover_stubs=False,
        )
        == 0
    )
    assert (tmp_path / "checksums.sha256").read_text() == "original\n"


def test_post_hoc_stubs_are_marked_and_name_the_environment(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    api.write_manifest_stub(tmp_path, tmp_path / "commands.sh")
    api.write_environment_stub(tmp_path, name="claw-example")
    manifest = __import__("json").loads((tmp_path / "manifest.json").read_text())
    assert manifest["regenerated_post_hoc"] is True
    assert manifest["commands_sh"] == str(tmp_path / "commands.sh")
    env = yaml.safe_load((tmp_path / "environment.yml").read_text())
    assert env == {"name": "claw-example", "regenerated_post_hoc": True}
