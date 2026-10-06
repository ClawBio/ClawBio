"""End-to-end runs of the entry script, with no network and no real snakemake.

Every test here drives `metagenome_mag_pipeline.py` as a subprocess with the
stub `snakemake` and `conda` on PATH and `--pipeline-dir` pointing at the fake
checkout, so the whole flow — samplesheet, config, launch, parse, report — is
exercised without conda building twenty tool environments.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "metagenome_mag_pipeline.py"
FIXTURES = Path(__file__).resolve().parent / "fixtures"
REPO_ROOT = SKILL_DIR.parent.parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

COMMON = ["--pipeline-dir", "{pipeline}", "--conda-prefix", "{prefix}"]


def run_cli(
    args: list[str],
    *,
    fake_pipeline: Path,
    tmp_path: Path,
    stubs: bool = True,
    extra_env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess:
    """Invoke the script the way AGENTS.md documents, with a hermetic PATH."""
    bindir = tmp_path / "bindir"
    bindir.mkdir(exist_ok=True)
    env = dict(os.environ)
    if stubs:
        for name in ("snakemake", "conda"):
            target = bindir / name
            target.write_text((FIXTURES / f"stub_{name}.py").read_text())
            target.chmod(0o755)
        env["PATH"] = f"{bindir}{os.pathsep}{env['PATH']}"
        env["STUB_SNAKEMAKE_ARGV"] = str(tmp_path / "snakemake_argv.json")
    env["PYTHONPATH"] = str(REPO_ROOT)
    env.pop("CLAWBIO_CACHE_DIR", None)
    if extra_env:
        env.update(extra_env)
    resolved = [
        item.replace("{pipeline}", str(fake_pipeline)).replace("{prefix}", str(tmp_path / "conda"))
        for item in args
    ]
    return subprocess.run(
        [sys.executable, str(SCRIPT), *resolved],
        capture_output=True,
        text=True,
        env=env,
        cwd=str(tmp_path),
        timeout=300,
    )


def demo_args(output: Path, pipeline: Path, tmp_path: Path) -> list[str]:
    return ["--demo", "--output", str(output), *COMMON]


def error_of(result: subprocess.CompletedProcess) -> dict:
    """The structured error the script prints to stderr, which is pretty-printed."""
    decoder = json.JSONDecoder()
    for index, character in enumerate(result.stderr):
        if character != "{":
            continue
        try:
            payload, _ = decoder.raw_decode(result.stderr[index:])
        except json.JSONDecodeError:
            continue
        if isinstance(payload, dict) and "error_code" in payload:
            return payload
    raise AssertionError(f"no JSON error object on stderr:\n{result.stderr}")


class TestDemoEndToEnd:
    def test_produces_the_documented_output_tree(
        self, fake_pipeline: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "demo"
        result = run_cli(demo_args(out, fake_pipeline, tmp_path), fake_pipeline=fake_pipeline, tmp_path=tmp_path)
        assert result.returncode == 0, result.stderr
        for relative in (
            "report.md",
            "result.json",
            "pipeline/samplesheet.tsv",
            "pipeline/config.yaml",
            "pipeline/snakemake.log",
            "reproducibility/commands.sh",
            "reproducibility/environment.yml",
            "reproducibility/checksums.sha256",
            "results/final/qc/read_counts.tsv",
            "results/final/assembly/quast.tsv",
            "results/final/diagnostic/bottleneck.tsv",
            "results/provenance/run_manifest.json",
        ):
            assert (out / relative).exists(), f"missing {relative}"

    def test_result_json_is_ok_and_demo(
        self, fake_pipeline: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "demo"
        run_cli(demo_args(out, fake_pipeline, tmp_path), fake_pipeline=fake_pipeline, tmp_path=tmp_path)
        data = json.loads((out / "result.json").read_text())
        assert data["status"] == "ok"
        assert data["mode"] == "demo"
        assert data["skill"] == "metagenome-mag-pipeline"
        assert data["disclaimer"]
        assert data["pipeline"]["source_kind"] == "local_checkout"

    def test_samplesheet_and_config_are_absolute(
        self, fake_pipeline: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "demo"
        run_cli(demo_args(out, fake_pipeline, tmp_path), fake_pipeline=fake_pipeline, tmp_path=tmp_path)
        sheet = (out / "pipeline" / "samplesheet.tsv").read_text().splitlines()
        for row in sheet[1:]:
            for cell in row.split("\t")[1:3]:
                for path in cell.split(";"):
                    if path:
                        assert path.startswith("/"), path
        import yaml

        config = yaml.safe_load((out / "pipeline" / "config.yaml").read_text())
        for key in ("samplesheet", "results_dir", "database_dir"):
            assert config[key].startswith("/"), key

    def test_recorded_argv_holds_only_absolute_paths(
        self, fake_pipeline: Path, tmp_path: Path
    ) -> None:
        """Upstream resolves relative config paths against the pipeline root."""
        out = tmp_path / "demo"
        run_cli(demo_args(out, fake_pipeline, tmp_path), fake_pipeline=fake_pipeline, tmp_path=tmp_path)
        argv = json.loads((tmp_path / "snakemake_argv.json").read_text())
        for flag in ("-s", "--configfile", "--directory", "--conda-prefix"):
            value = argv[argv.index(flag) + 1]
            assert value.startswith("/"), f"{flag} got a relative path: {value!r}"

    def test_demo_uses_the_upstream_test_configuration(
        self, fake_pipeline: Path, tmp_path: Path
    ) -> None:
        import yaml

        out = tmp_path / "demo"
        run_cli(demo_args(out, fake_pipeline, tmp_path), fake_pipeline=fake_pipeline, tmp_path=tmp_path)
        config = yaml.safe_load((out / "pipeline" / "config.yaml").read_text())
        base_test = yaml.safe_load((FIXTURES / "upstream_config_test.yaml").read_text())
        differing = {
            key for key in set(config) | set(base_test) if config.get(key) != base_test.get(key)
        }
        assert differing == {"samplesheet", "results_dir", "database_dir"}

    def test_checksums_verify(self, fake_pipeline: Path, tmp_path: Path) -> None:
        out = tmp_path / "demo"
        run_cli(demo_args(out, fake_pipeline, tmp_path), fake_pipeline=fake_pipeline, tmp_path=tmp_path)
        result = subprocess.run(
            ["sha256sum", "-c", "reproducibility/checksums.sha256"],
            cwd=out, capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr

    def test_report_carries_the_disclaimer_and_the_demo_caveat(
        self, fake_pipeline: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "demo"
        run_cli(demo_args(out, fake_pipeline, tmp_path), fake_pipeline=fake_pipeline, tmp_path=tmp_path)
        text = (out / "report.md").read_text()
        from magpipe_schemas import DISCLAIMER

        assert DISCLAIMER in " ".join(text.split())
        assert "synthetic" in text.lower()


class TestOutputContract:
    """Every path SKILL.md's Output Structure tree promises must exist after a demo."""

    @staticmethod
    def _parse_output_contract(skill_md: Path) -> list[str]:
        """Extract files promised in SKILL.md '## Output Structure' tree."""
        if not skill_md.exists():
            return []
        text = skill_md.read_text()
        match = re.search(r"##\s*Output Structure\s*\n+```[^\n]*\n(.*?)\n```", text, re.DOTALL)
        if not match:
            return []
        files: list[str] = []
        parents: dict[int, str] = {}
        for raw in match.group(1).splitlines():
            if not raw.strip():
                continue
            parts = re.split(r"\s+#", raw, maxsplit=1)
            entry = parts[0]
            comment = parts[1] if len(parts) > 1 else ""
            found = re.match(r"^([\s│├└─]*)(.*)$", entry)
            prefix, name = found.group(1), found.group(2).strip()
            if not name:
                continue
            depth = len(prefix) // 4
            if depth == 0:
                continue
            if name.endswith("/"):
                parents[depth] = name.rstrip("/")
                for deeper in [d for d in parents if d > depth]:
                    del parents[deeper]
                continue
            if "optional" in comment.lower():
                continue
            prefix_path = "/".join(parents[d] for d in sorted(parents) if d < depth)
            files.append(f"{prefix_path}/{name}" if prefix_path else name)
        return files

    def test_documented_outputs_are_produced(
        self, fake_pipeline: Path, tmp_path: Path
    ) -> None:
        promised = self._parse_output_contract(SKILL_DIR / "SKILL.md")
        assert promised, "SKILL.md's ## Output Structure tree is not parseable"
        out = tmp_path / "demo"
        result = run_cli(
            demo_args(out, fake_pipeline, tmp_path), fake_pipeline=fake_pipeline, tmp_path=tmp_path
        )
        assert result.returncode == 0, f"demo run failed: {result.stderr}"
        missing = [p for p in promised if not (out / p).exists()]
        assert not missing, (
            "SKILL.md promises artifacts the skill did not produce: " + ", ".join(missing)
        )


class TestCheckMode:
    def test_planned_report_and_no_results(
        self, fake_pipeline: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "plan"
        result = run_cli(
            demo_args(out, fake_pipeline, tmp_path) + ["--check"],
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
        )
        assert result.returncode == 0, result.stderr
        assert not (out / "results").exists()
        data = json.loads((out / "result.json").read_text())
        assert data["status"] == "planned"
        assert data["mode"] == "check"
        assert data["findings"] is None
        text = (out / "report.md").read_text()
        assert "Nothing was executed" in text
        assert "snakemake" in text

    def test_dry_run_flag_reaches_snakemake(
        self, fake_pipeline: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "plan"
        run_cli(
            demo_args(out, fake_pipeline, tmp_path) + ["--check"],
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
        )
        argv = json.loads((tmp_path / "snakemake_argv.json").read_text())
        assert "--dry-run" in argv


class TestUserInput:
    def _sheet(self, tmp_path: Path, reads: Path) -> Path:
        from magpipe_samplesheet import write_samplesheet

        samples = []
        from magpipe_samplesheet import Sample

        for name in ("S1", "S2"):
            first = reads / name / f"{name}_R1.fastq.gz"
            second = reads / name / f"{name}_R2.fastq.gz"
            for path in (first, second):
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03")
            samples.append(
                Sample(name, (first,), (second,), name, None, True)
            )
        return write_samplesheet(samples, tmp_path / "sheet.tsv")

    def test_input_samplesheet(self, fake_pipeline: Path, tmp_path: Path) -> None:
        sheet = self._sheet(tmp_path, tmp_path / "reads")
        out = tmp_path / "out"
        result = run_cli(
            ["--input", str(sheet), "--output", str(out), *COMMON],
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
        )
        assert result.returncode == 0, result.stderr
        data = json.loads((out / "result.json").read_text())
        assert data["status"] == "ok"
        assert data["mode"] == "run"
        assert len(data["config"]["samplesheet"]) > 0

    def test_input_reads_directory(self, fake_pipeline: Path, tmp_path: Path) -> None:
        reads = tmp_path / "reads"
        for name in ("S1",):
            for mate in ("R1", "R2"):
                path = reads / name / f"{name}_{mate}.fastq.gz"
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03")
        out = tmp_path / "out"
        result = run_cli(
            ["--input", str(reads), "--output", str(out), *COMMON],
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
        )
        assert result.returncode == 0, result.stderr
        rows = (out / "pipeline" / "samplesheet.tsv").read_text().splitlines()
        assert rows[0].split("\t")[:2] == ["sample", "fastq_1"]
        assert rows[1].split("\t")[0] == "S1"

    def test_full_preset_without_databases_runs_binning_only(
        self, fake_pipeline: Path, tmp_path: Path
    ) -> None:
        reads = tmp_path / "reads"
        for mate in ("R1", "R2"):
            path = reads / "S1" / f"S1_{mate}.fastq.gz"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03")
        out = tmp_path / "out"
        result = run_cli(
            ["--input", str(reads), "--output", str(out), "--preset", "full", *COMMON],
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
        )
        assert result.returncode == 0, result.stderr
        data = json.loads((out / "result.json").read_text())
        assert data["status"] == "ok"
        assert data["findings"]["bins"] is None
        assert data["findings"]["taxonomy"] is None
        assert data["findings"]["annotation"] is None
        skipped = {
            entry["stage"]: entry["reason"]
            for entry in data["stages"]
            if not entry["enabled"]
        }
        assert "--checkm2-db" in skipped["MAG assessment (CheckM2, dRep, GTDB-Tk)"]
        assert "--taxonomy-db" in skipped["contig taxonomy (DIAMOND, LCA)"]
        text = (out / "report.md").read_text()
        assert "--checkm2-db" in text

    def test_single_end_input_rejected_for_assembly(
        self, fake_pipeline: Path, tmp_path: Path
    ) -> None:
        reads = tmp_path / "reads"
        path = reads / "SE1" / "SE1_R1.fastq.gz"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03")
        out = tmp_path / "out"
        result = run_cli(
            ["--input", str(reads), "--output", str(out), "--preset", "assembly", *COMMON],
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
        )
        assert result.returncode == 1
        assert error_of(result)["error_code"] == "ASSEMBLY_REQUIRES_PAIRED"

    def test_qc_preset_accepts_single_end(self, fake_pipeline: Path, tmp_path: Path) -> None:
        reads = tmp_path / "reads"
        path = reads / "SE1" / "SE1_R1.fastq.gz"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03")
        out = tmp_path / "out"
        result = run_cli(
            ["--input", str(reads), "--output", str(out), "--preset", "qc", *COMMON],
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
        )
        assert result.returncode == 0, result.stderr
        assert json.loads((out / "result.json").read_text())["status"] == "ok"


class TestFailures:
    def test_no_input_and_no_demo(self, fake_pipeline: Path, tmp_path: Path) -> None:
        result = run_cli(
            ["--output", str(tmp_path / "out"), *COMMON],
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
        )
        assert result.returncode == 1
        payload = error_of(result)
        assert payload["error_code"] == "MISSING_INPUT"
        assert "--demo" in payload["fix"]
        assert payload["ok"] is False

    def test_non_empty_output_without_force(
        self, fake_pipeline: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "out"
        out.mkdir()
        (out / "stale.txt").write_text("old\n")
        result = run_cli(
            demo_args(out, fake_pipeline, tmp_path), fake_pipeline=fake_pipeline, tmp_path=tmp_path
        )
        assert result.returncode == 1
        payload = error_of(result)
        assert payload["error_code"] == "OUTPUT_DIR_NOT_EMPTY"
        assert "--force" in payload["fix"]

    def test_non_empty_output_with_force(self, fake_pipeline: Path, tmp_path: Path) -> None:
        out = tmp_path / "out"
        out.mkdir()
        (out / "stale.txt").write_text("old\n")
        result = run_cli(
            demo_args(out, fake_pipeline, tmp_path) + ["--force"],
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
        )
        assert result.returncode == 0, result.stderr

    def test_snakemake_missing(self, fake_pipeline: Path, tmp_path: Path) -> None:
        """With no snakemake the run fails, and no findings are invented."""
        out = tmp_path / "out"
        result = run_cli(
            demo_args(out, fake_pipeline, tmp_path),
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
            extra_env={"PATH": "/nonexistent-bin"},
        )
        assert result.returncode == 1
        payload = error_of(result)
        assert payload["error_code"] == "MISSING_SNAKEMAKE"
        assert "snakemake=9.11.2" in payload["fix"]

    def test_conda_missing(self, fake_pipeline: Path, tmp_path: Path) -> None:
        bindir = tmp_path / "only_snakemake"
        bindir.mkdir()
        stub = bindir / "snakemake"
        stub.write_text((FIXTURES / "stub_snakemake.py").read_text())
        stub.chmod(0o755)
        # The stub is `#!/usr/bin/env python3`, so a hermetic PATH still needs one.
        (bindir / "python3").symlink_to(Path(sys.executable))
        out = tmp_path / "out"
        result = run_cli(
            demo_args(out, fake_pipeline, tmp_path),
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
            stubs=False,
            extra_env={"PATH": str(bindir)},
        )
        assert result.returncode == 1
        assert error_of(result)["error_code"] == "MISSING_CONDA"

    def test_stub_failure_exits_nonzero_and_writes_no_findings(
        self, fake_pipeline: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "out"
        result = run_cli(
            demo_args(out, fake_pipeline, tmp_path),
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
            extra_env={"STUB_SNAKEMAKE_EXIT": "1"},
        )
        assert result.returncode == 1
        payload = error_of(result)
        assert payload["error_code"] == "EXECUTION_FAILED"
        assert payload["details"]["log_tail"]
        data = json.loads((out / "result.json").read_text())
        assert data["status"] == "error"
        assert data["findings"] is None

    def test_bad_pipeline_dir(self, tmp_path: Path) -> None:
        empty = tmp_path / "not-a-pipeline"
        empty.mkdir()
        out = tmp_path / "out"
        result = run_cli(
            ["--demo", "--output", str(out), "--pipeline-dir", str(empty),
             "--conda-prefix", str(tmp_path / "conda")],
            fake_pipeline=empty,
            tmp_path=tmp_path,
        )
        assert result.returncode == 1
        assert error_of(result)["error_code"] == "PIPELINE_SOURCE_INVALID"

    def test_remote_sample_needs_the_flag(self, fake_pipeline: Path, tmp_path: Path) -> None:
        sheet = tmp_path / "remote.tsv"
        sheet.write_text(
            "sample\tfastq_1\tfastq_2\tgroup\tsra_run\tlayout\n"
            "SRR1761673\t\t\tSRR1761673\tSRR1761673\tpaired\n"
        )
        out = tmp_path / "out"
        result = run_cli(
            ["--input", str(sheet), "--output", str(out), *COMMON],
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
        )
        assert result.returncode == 1
        payload = error_of(result)
        assert payload["error_code"] == "REMOTE_INPUT_NOT_ALLOWED"
        assert "--allow-remote-inputs" in payload["fix"]

    def test_error_result_json_is_written_when_the_output_dir_exists(
        self, fake_pipeline: Path, tmp_path: Path
    ) -> None:
        out = tmp_path / "out"
        run_cli(
            demo_args(out, fake_pipeline, tmp_path),
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
            extra_env={"STUB_SNAKEMAKE_EXIT": "1"},
        )
        data = json.loads((out / "result.json").read_text())
        assert data["status"] == "error"
        assert data["warnings"]

    def test_missing_input_path(self, fake_pipeline: Path, tmp_path: Path) -> None:
        result = run_cli(
            ["--input", str(tmp_path / "nope"), "--output", str(tmp_path / "out"), *COMMON],
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
        )
        assert result.returncode == 1
        assert error_of(result)["error_code"] == "MISSING_INPUT"


class TestDemoFlagConflicts:
    """--demo is upstream's test configuration, so flags it cannot honour are refused.

    The alternative — accept the flag and quietly not apply it — leaves a user
    believing they ran, say, a binned analysis when no binning rule was enabled.
    A structured error naming every offending flag is the only honest outcome.
    """

    def _run(self, args, fake_pipeline: Path, tmp_path: Path):
        out = tmp_path / "out"
        return run_cli(
            ["--demo", "--output", str(out), *args, *COMMON],
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
        )

    @pytest.mark.parametrize(
        "args",
        [
            ["--preset", "full"],
            ["--preset", "qc"],
            ["--preset", "profile"],
            ["--coassemble"],
            ["--host-fasta", "/db/host.fa"],
            ["--skip-host-removal"],
            ["--checkm2-db", "/db/checkm2.dmnd"],
            ["--gtdbtk-db", "/db/gtdbtk"],
            ["--taxonomy-db", "/db/nr.dmnd"],
            ["--taxdump", "/db/taxdump"],
            ["--kegg-db", "/db/kegg.dmnd"],
            ["--cog-db", "/db/cog.dmnd"],
            ["--pfam-db", "/db/Pfam-A.hmm"],
            ["--binners", "metabat,maxbin"],
            ["--allow-remote-inputs"],
        ],
    )
    def test_refused(self, args, fake_pipeline: Path, tmp_path: Path) -> None:
        result = self._run(args, fake_pipeline, tmp_path)
        assert result.returncode == 1, result.stderr
        payload = error_of(result)
        assert payload["error_code"] == "INVALID_CONFIG"
        assert args[0] in payload["message"] or args[0] in payload["details"]["conflicts"]
        assert "--demo" in payload["fix"]
        assert "config_test.yaml" in payload["message"] or (
            "config_test.yaml" in payload["details"]["base_config"]
        )

    def test_input_and_demo_together_are_refused(self, fake_pipeline: Path, tmp_path: Path) -> None:
        sheet = tmp_path / "s.tsv"
        sheet.write_text("sample\tfastq_1\tfastq_2\tgroup\tsra_run\tlayout\n")
        result = run_cli(
            ["--demo", "--input", str(sheet), "--output", str(tmp_path / "out"), *COMMON],
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
        )
        assert result.returncode == 1
        payload = error_of(result)
        assert payload["error_code"] == "INVALID_CONFIG"
        assert "--input" in payload["message"]

    def test_every_conflicting_flag_is_named_at_once(self, fake_pipeline: Path, tmp_path: Path) -> None:
        result = self._run(
            ["--preset", "full", "--coassemble", "--binners", "metabat,vamb"],
            fake_pipeline,
            tmp_path,
        )
        assert result.returncode == 1
        conflicts = error_of(result)["details"]["conflicts"]
        assert {"--preset", "--coassemble", "--binners"} <= set(conflicts)

    def test_plain_demo_still_runs(self, fake_pipeline: Path, tmp_path: Path) -> None:
        result = self._run([], fake_pipeline, tmp_path)
        assert result.returncode == 0, result.stderr
        data = json.loads((tmp_path / "out" / "result.json").read_text())
        assert data["status"] == "ok"
        assert data["warnings"] == []

    def test_demo_check_still_plans(self, fake_pipeline: Path, tmp_path: Path) -> None:
        result = self._run(["--check"], fake_pipeline, tmp_path)
        assert result.returncode == 0, result.stderr
        assert json.loads((tmp_path / "out" / "result.json").read_text())["status"] == "planned"

    def test_resource_flags_are_honoured_in_demo(self, fake_pipeline: Path, tmp_path: Path) -> None:
        """--cores and --mem-mb reach the scheduler even though demo keeps its own config."""
        result = self._run(["--cores", "2", "--mem-mb", "4000"], fake_pipeline, tmp_path)
        assert result.returncode == 0, result.stderr
        argv = json.loads((tmp_path / "snakemake_argv.json").read_text())
        assert argv[argv.index("--cores") + 1] == "2"
        assert argv[argv.index("--resources") + 1] == "mem_mb=4000"

    def test_refusal_happens_before_any_clone(self, fake_pipeline: Path, tmp_path: Path) -> None:
        """A conflicting flag must not cost the user a clone of the workflow."""
        cache = tmp_path / "cache"
        result = run_cli(
            ["--demo", "--preset", "full", "--output", str(tmp_path / "out"), *COMMON],
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
            extra_env={"CLAWBIO_CACHE_DIR": str(cache)},
        )
        assert result.returncode == 1
        assert not cache.exists()

    def test_base_config_is_recorded_in_result_json(self, fake_pipeline: Path, tmp_path: Path) -> None:
        result = self._run([], fake_pipeline, tmp_path)
        assert result.returncode == 0, result.stderr
        data = json.loads((tmp_path / "out" / "result.json").read_text())
        assert data["config_source"]["base_config"] == "workflow/config/config_test.yaml"
        assert data["config_source"]["demo"] is True
        # The pin is what the wrapper asked for; the commit is what the checkout
        # reports, which is "unknown" for a directory that is not a git repo.
        assert data["config_source"]["pipeline_ref"] == "7351a702d29857801963cae7aeff6600d59ebe77"
        assert data["config_source"]["pipeline_commit"]

    def test_database_flag_for_a_disabled_stage_is_refused(
        self, fake_pipeline: Path, tmp_path: Path
    ) -> None:
        reads = tmp_path / "reads"
        for mate in ("R1", "R2"):
            path = reads / "S1" / f"S1_{mate}.fastq.gz"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03")
        database = tmp_path / "checkm2.dmnd"
        database.write_text("x\n")
        result = run_cli(
            ["--input", str(reads), "--preset", "assembly", "--checkm2-db", str(database),
             "--output", str(tmp_path / "out"), *COMMON],
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
        )
        assert result.returncode == 1
        payload = error_of(result)
        assert payload["error_code"] == "INVALID_CONFIG"
        assert "--checkm2-db" in payload["message"]
        assert "--preset full" in payload["fix"]

    def test_a_real_run_records_the_main_config(self, fake_pipeline: Path, tmp_path: Path) -> None:
        reads = tmp_path / "reads"
        for mate in ("R1", "R2"):
            path = reads / "S1" / f"S1_{mate}.fastq.gz"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03")
        out = tmp_path / "out"
        result = run_cli(
            ["--input", str(reads), "--output", str(out), *COMMON],
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
        )
        assert result.returncode == 0, result.stderr
        data = json.loads((out / "result.json").read_text())
        assert data["config_source"]["base_config"] == "workflow/config/config.yaml"
        assert data["config_source"]["demo"] is False

    def test_report_names_the_base_config(self, fake_pipeline: Path, tmp_path: Path) -> None:
        result = self._run([], fake_pipeline, tmp_path)
        assert result.returncode == 0, result.stderr
        text = (tmp_path / "out" / "report.md").read_text()
        assert "workflow/config/config_test.yaml" in text


class TestResourceAndFlagPlumbing:
    """Flags that change what snakemake is asked to do must reach the config."""

    def _sheet(self, tmp_path: Path) -> Path:
        reads = tmp_path / "reads"
        for mate in ("R1", "R2"):
            path = reads / "S1" / f"S1_{mate}.fastq.gz"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03")
        from magpipe_samplesheet import Sample, write_samplesheet

        return write_samplesheet(
            [
                Sample("S1", (reads / "S1" / "S1_R1.fastq.gz",), (reads / "S1" / "S1_R2.fastq.gz",),
                       "S1", None, True)
            ],
            tmp_path / "sheet.tsv",
        )

    def test_cores_and_memory_reach_the_config(self, fake_pipeline: Path, tmp_path: Path) -> None:
        import yaml

        sheet = self._sheet(tmp_path)
        out = tmp_path / "out"
        result = run_cli(
            ["--input", str(sheet), "--output", str(out), "--cores", "2", "--mem-mb", "4000", *COMMON],
            fake_pipeline=fake_pipeline, tmp_path=tmp_path,
        )
        assert result.returncode == 0, result.stderr
        config = yaml.safe_load((out / "pipeline" / "config.yaml").read_text())
        assert max(config["threads"].values()) <= 2
        assert max(config["mem_mb"].values()) <= 4000
        argv = json.loads((tmp_path / "snakemake_argv.json").read_text())
        assert argv[argv.index("--cores") + 1] == "2"
        assert argv[argv.index("--resources") + 1] == "mem_mb=4000"

    def test_single_binner_is_refused(self, fake_pipeline: Path, tmp_path: Path) -> None:
        sheet = self._sheet(tmp_path)
        result = run_cli(
            ["--input", str(sheet), "--preset", "full", "--binners", "metabat",
             "--output", str(tmp_path / "out"), *COMMON],
            fake_pipeline=fake_pipeline, tmp_path=tmp_path,
        )
        assert result.returncode == 1
        payload = error_of(result)
        assert payload["error_code"] == "INVALID_CONFIG"
        assert "at least two" in payload["message"]

    def test_unknown_binner_is_refused(self, fake_pipeline: Path, tmp_path: Path) -> None:
        sheet = self._sheet(tmp_path)
        result = run_cli(
            ["--input", str(sheet), "--preset", "full", "--binners", "metabat,concoct",
             "--output", str(tmp_path / "out"), *COMMON],
            fake_pipeline=fake_pipeline, tmp_path=tmp_path,
        )
        assert result.returncode == 1
        assert error_of(result)["error_code"] == "INVALID_CONFIG"

    def test_output_path_with_spaces(self, fake_pipeline: Path, tmp_path: Path) -> None:
        out = tmp_path / "my output dir"
        result = run_cli(
            ["--demo", "--output", str(out), *COMMON],
            fake_pipeline=fake_pipeline, tmp_path=tmp_path,
        )
        assert result.returncode == 0, result.stderr
        assert (out / "report.md").is_file()
        checks = subprocess.run(
            ["sha256sum", "-c", "reproducibility/checksums.sha256"],
            cwd=out, capture_output=True, text=True,
        )
        assert checks.returncode == 0, checks.stdout + checks.stderr


class TestCliSurface:
    def test_help_exits_zero(self) -> None:
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--help"],
            capture_output=True, text=True,
            env={**os.environ, "PYTHONPATH": str(REPO_ROOT)},
            timeout=120,
        )
        assert result.returncode == 0
        assert "--preset" in result.stdout

    def test_help_runs_without_pythonpath(self, tmp_path: Path) -> None:
        env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--help"],
            capture_output=True, text=True, env=env, cwd=str(tmp_path), timeout=120,
        )
        assert "ModuleNotFoundError" not in result.stderr, result.stderr
        assert result.returncode == 0

    def test_network_requirements_are_printed_before_launch(
        self, fake_pipeline: Path, tmp_path: Path, capsys
    ) -> None:
        out = tmp_path / "demo"
        result = run_cli(
            demo_args(out, fake_pipeline, tmp_path), fake_pipeline=fake_pipeline, tmp_path=tmp_path
        )
        assert result.returncode == 0
        assert "network" in result.stderr.lower()
        assert "conda" in result.stderr.lower()

    def test_output_is_required(self) -> None:
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--demo"],
            capture_output=True, text=True,
            env={**os.environ, "PYTHONPATH": str(REPO_ROOT)},
            timeout=120,
        )
        assert result.returncode == 2
        assert "--output" in result.stderr

    @pytest.mark.parametrize("preset", ["qc", "profile", "assembly", "full"])
    def test_every_preset_validates_in_check_mode(
        self, preset: str, fake_pipeline: Path, tmp_path: Path
    ) -> None:
        """Every preset plans against real paired-end input, not the demo."""
        reads = tmp_path / f"reads-{preset}"
        for mate in ("R1", "R2"):
            path = reads / "S1" / f"S1_{mate}.fastq.gz"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03")
        out = tmp_path / preset
        result = run_cli(
            ["--input", str(reads), "--output", str(out), "--preset", preset, "--check", *COMMON],
            fake_pipeline=fake_pipeline,
            tmp_path=tmp_path,
        )
        assert result.returncode == 0, result.stderr
        assert json.loads((out / "result.json").read_text())["status"] == "planned"
