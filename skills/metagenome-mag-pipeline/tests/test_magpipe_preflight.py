"""Preflight: pipeline source, snakemake, conda, output directory, network use.

Nothing here launches anything. The point is to fail in a second, naming the
flag that fixes it, rather than after conda has built twenty tool environments.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest
import yaml

SKILL_DIR = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).resolve().parent / "fixtures"
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

import magpipe_preflight as preflight
from magpipe_errors import ErrorCode, SkillError
from magpipe_schemas import (
    MIN_SNAKEMAKE,
    PINNED_REF,
    SNAKEMAKE_PIN,
    UPSTREAM_REPO,
)

BASE = yaml.safe_load((FIXTURES / "upstream_config.yaml").read_text())
BASE_TEST = yaml.safe_load((FIXTURES / "upstream_config_test.yaml").read_text())


def local(fake_pipeline: Path) -> dict:
    return preflight.resolve_pipeline_source(
        fake_pipeline, cache_dir=fake_pipeline.parent / "cache", ref=PINNED_REF, allow_fetch=False
    )


def cache_dir_for(fake_pipeline: Path, tmp_path: Path) -> Path:
    return tmp_path / "cache"


def cached_clone(fake_pipeline: Path, cache: Path) -> Path:
    """Materialise the fake checkout at the cache location the resolver looks for."""
    target = cache / f"metagenomics-workflow-{PINNED_REF[:12]}"
    target.mkdir(parents=True)
    for name in ("workflow/Snakefile", "workflow/config/config.yaml",
                 "workflow/config/config_test.yaml", "test/make_test_dataset.py"):
        destination = target / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text((fake_pipeline / name).read_text())
    return target


class TestLocalCheckout:
    def test_local_checkout_is_accepted(self, fake_pipeline: Path) -> None:
        source = local(fake_pipeline)
        assert source["path"] == fake_pipeline
        assert source["source_kind"] == "local_checkout"
        assert source["ref"] == PINNED_REF
        assert "commit" in source and "dirty" in source

    def test_missing_required_file(self, fake_pipeline: Path) -> None:
        (fake_pipeline / "workflow" / "Snakefile").unlink()
        with pytest.raises(SkillError) as excinfo:
            local(fake_pipeline)
        assert excinfo.value.error_code == ErrorCode.PIPELINE_SOURCE_INVALID
        assert "workflow/Snakefile" in excinfo.value.message
        assert "--pipeline-dir" in excinfo.value.fix

    def test_directory_that_is_not_a_checkout(self, tmp_path: Path) -> None:
        empty = tmp_path / "not-a-pipeline"
        empty.mkdir()
        with pytest.raises(SkillError) as excinfo:
            preflight.resolve_pipeline_source(
                empty, cache_dir=tmp_path / "c", ref=PINNED_REF, allow_fetch=False
            )
        assert excinfo.value.error_code == ErrorCode.PIPELINE_SOURCE_INVALID
        assert len(excinfo.value.details["missing"]) == 4

    def test_missing_directory(self, tmp_path: Path) -> None:
        with pytest.raises(SkillError) as excinfo:
            preflight.resolve_pipeline_source(
                tmp_path / "gone", cache_dir=tmp_path / "c", ref=PINNED_REF, allow_fetch=False
            )
        assert excinfo.value.error_code == ErrorCode.PIPELINE_SOURCE_INVALID


class TestCachedAndFetched:
    def test_cache_hit(self, fake_pipeline: Path, tmp_path: Path) -> None:
        cache = cache_dir_for(fake_pipeline, tmp_path)
        cached_clone(fake_pipeline, cache)
        source = preflight.resolve_pipeline_source(
            None, cache_dir=cache, ref=PINNED_REF, allow_fetch=False
        )
        assert source["source_kind"] == "cached_clone"
        assert source["path"] == cache / f"metagenomics-workflow-{PINNED_REF[:12]}"
        assert source["ref"] == PINNED_REF

    def test_no_cache_and_fetch_disallowed(self, tmp_path: Path) -> None:
        with pytest.raises(SkillError) as excinfo:
            preflight.resolve_pipeline_source(
                None, cache_dir=tmp_path / "c", ref=PINNED_REF, allow_fetch=False
            )
        assert excinfo.value.error_code == ErrorCode.PIPELINE_FETCH_FAILED
        assert "--pipeline-dir" in excinfo.value.fix

    def test_incomplete_cache_is_not_a_hit(self, fake_pipeline: Path, tmp_path: Path) -> None:
        cache = cache_dir_for(fake_pipeline, tmp_path)
        cached_clone(fake_pipeline, cache)
        (cache / f"metagenomics-workflow-{PINNED_REF[:12]}" / "workflow" / "Snakefile").unlink()
        with pytest.raises(SkillError) as excinfo:
            preflight.resolve_pipeline_source(
                None, cache_dir=cache, ref=PINNED_REF, allow_fetch=False
            )
        assert excinfo.value.error_code == ErrorCode.PIPELINE_FETCH_FAILED

    def test_clone_when_fetch_allowed(self, fake_pipeline: Path, tmp_path: Path, monkeypatch) -> None:
        """The clone is faked: no test in this suite touches the network."""
        cache = cache_dir_for(fake_pipeline, tmp_path)
        commands: list[list[str]] = []

        def fake_run(argv, **kwargs):
            commands.append(list(argv))
            if "clone" in argv:
                destination = Path(argv[-1])
                for name in (
                    "workflow/Snakefile",
                    "workflow/config/config.yaml",
                    "workflow/config/config_test.yaml",
                    "test/make_test_dataset.py",
                ):
                    path = destination / name
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_text((fake_pipeline / name).read_text())
            return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

        monkeypatch.setattr(preflight.subprocess, "run", fake_run)
        source = preflight.resolve_pipeline_source(
            None, cache_dir=cache, ref=PINNED_REF, allow_fetch=True
        )
        assert source["source_kind"] == "cached_clone"
        clone = [c for c in commands if "clone" in c]
        assert clone and UPSTREAM_REPO in clone[0]
        checkout = [c for c in commands if "checkout" in c]
        assert checkout and PINNED_REF in checkout[0]

    def test_clone_failure_is_reported(self, tmp_path: Path, monkeypatch) -> None:
        def fake_run(argv, **kwargs):
            raise subprocess.CalledProcessError(128, argv, stderr="network is unreachable")

        monkeypatch.setattr(preflight.subprocess, "run", fake_run)
        with pytest.raises(SkillError) as excinfo:
            preflight.resolve_pipeline_source(
                None, cache_dir=tmp_path / "c", ref=PINNED_REF, allow_fetch=True
            )
        assert excinfo.value.error_code == ErrorCode.PIPELINE_FETCH_FAILED
        assert "network is unreachable" in excinfo.value.message


class TestSnakemake:
    def test_snakemake_missing(self, no_pipeline_tools: None) -> None:
        with pytest.raises(SkillError) as excinfo:
            preflight.check_snakemake()
        assert excinfo.value.error_code == ErrorCode.MISSING_SNAKEMAKE
        assert f"snakemake={SNAKEMAKE_PIN}" in excinfo.value.fix
        assert "conda create" in excinfo.value.fix

    def test_snakemake_ok(self, stub_snakemake: Path) -> None:
        assert preflight.check_snakemake() == SNAKEMAKE_PIN

    def test_snakemake_too_old(self, stub_snakemake: Path, monkeypatch) -> None:
        monkeypatch.setattr(preflight.subprocess, "run", _version_output("8.27.0"))
        with pytest.raises(SkillError) as excinfo:
            preflight.check_snakemake()
        assert excinfo.value.error_code == ErrorCode.SNAKEMAKE_VERSION_TOO_OLD
        assert "8.27.0" in excinfo.value.message
        assert ".".join(str(part) for part in MIN_SNAKEMAKE) in excinfo.value.message

    def test_snakemake_exactly_at_the_floor(self, stub_snakemake: Path, monkeypatch) -> None:
        monkeypatch.setattr(preflight.subprocess, "run", _version_output("9.0.0"))
        assert preflight.check_snakemake() == "9.0.0"

    def test_unparseable_version(self, stub_snakemake: Path, monkeypatch) -> None:
        monkeypatch.setattr(preflight.subprocess, "run", _version_output("unknown"))
        with pytest.raises(SkillError) as excinfo:
            preflight.check_snakemake()
        assert excinfo.value.error_code == ErrorCode.SNAKEMAKE_VERSION_TOO_OLD


def _version_output(version: str):
    def fake_run(argv, **kwargs):
        return subprocess.CompletedProcess(argv, 0, stdout=f"{version}\n", stderr="")

    return fake_run


class TestConda:
    def test_conda_found(self, stub_conda: Path) -> None:
        assert preflight.check_conda() == "conda"

    def test_mamba_accepted(self, stub_conda: Path, monkeypatch) -> None:
        """A machine with only mamba is fine: snakemake drives both frontends."""
        mamba = stub_conda.parent / "mamba_only" / "mamba"
        mamba.parent.mkdir()
        mamba.write_text((FIXTURES / "stub_conda.py").read_text())
        mamba.chmod(0o755)
        monkeypatch.setenv("PATH", str(mamba.parent))
        assert preflight.check_conda() == "mamba"

    def test_conda_missing(self, no_pipeline_tools: None) -> None:
        with pytest.raises(SkillError) as excinfo:
            preflight.check_conda()
        assert excinfo.value.error_code == ErrorCode.MISSING_CONDA
        assert "conda" in excinfo.value.fix


class TestOutputDir:
    def test_empty_directory_is_fine(self, tmp_path: Path) -> None:
        out = tmp_path / "out"
        out.mkdir()
        preflight.check_output_dir(out, force=False)

    def test_absent_directory_is_fine(self, tmp_path: Path) -> None:
        preflight.check_output_dir(tmp_path / "out", force=False)

    def test_non_empty_is_refused(self, tmp_path: Path) -> None:
        out = tmp_path / "out"
        out.mkdir()
        (out / "report.md").write_text("old\n")
        with pytest.raises(SkillError) as excinfo:
            preflight.check_output_dir(out, force=False)
        assert excinfo.value.error_code == ErrorCode.OUTPUT_DIR_NOT_EMPTY
        assert "--force" in excinfo.value.fix

    def test_non_empty_with_force_is_allowed(self, tmp_path: Path) -> None:
        out = tmp_path / "out"
        out.mkdir()
        (out / "report.md").write_text("old\n")
        preflight.check_output_dir(out, force=True)

    def test_only_a_directory_is_accepted(self, tmp_path: Path) -> None:
        target = tmp_path / "out"
        target.write_text("a file, not a directory\n")
        with pytest.raises(SkillError) as excinfo:
            preflight.check_output_dir(target, force=True)
        assert excinfo.value.error_code == ErrorCode.OUTPUT_DIR_NOT_WRITABLE

    def test_unwritable_directory(self, tmp_path: Path) -> None:
        out = tmp_path / "out"
        out.mkdir()
        out.chmod(0o500)
        try:
            with pytest.raises(SkillError) as excinfo:
                preflight.check_output_dir(out, force=True)
            assert excinfo.value.error_code == ErrorCode.OUTPUT_DIR_NOT_WRITABLE
        finally:
            out.chmod(0o700)

    def test_unwritable_parent(self, tmp_path: Path) -> None:
        parent = tmp_path / "locked"
        parent.mkdir()
        parent.chmod(0o500)
        try:
            with pytest.raises(SkillError) as excinfo:
                preflight.check_output_dir(parent / "out", force=False)
            assert excinfo.value.error_code == ErrorCode.OUTPUT_DIR_NOT_WRITABLE
        finally:
            parent.chmod(0o700)


class TestNetworkRequirements:
    def test_conda_environments_always(self) -> None:
        requirements = preflight.network_requirements(BASE, [])
        assert any("conda" in line for line in requirements)

    def test_demo_config_downloads_only_conda_envs(self) -> None:
        config = yaml.safe_load(yaml.safe_dump(BASE_TEST))
        requirements = preflight.network_requirements(config, [])
        joined = " ".join(requirements).lower()
        assert "grch38" not in joined
        assert "metaphlan" not in joined
        assert not any("sra" in line.lower() for line in requirements)

    def test_host_download_listed(self) -> None:
        config = yaml.safe_load(yaml.safe_dump(BASE))
        config["qc"]["decontaminate"]["enabled"] = True
        config["qc"]["decontaminate"]["fasta"] = None
        requirements = preflight.network_requirements(config, [])
        assert any("GRCh38" in line for line in requirements)

    def test_local_host_fasta_removes_the_download(self) -> None:
        config = yaml.safe_load(yaml.safe_dump(BASE))
        config["qc"]["decontaminate"]["enabled"] = True
        config["qc"]["decontaminate"]["fasta"] = "/db/host.fa"
        requirements = preflight.network_requirements(config, [])
        assert not any("GRCh38" in line for line in requirements)

    def test_sra_download_listed_when_enabled(self) -> None:
        config = yaml.safe_load(yaml.safe_dump(BASE))
        config["fetch"]["sra"]["enabled"] = True
        requirements = preflight.network_requirements(config, [])
        assert any("SRA" in line or "run" in line for line in requirements)

    def test_metaphlan_index_listed(self) -> None:
        config = yaml.safe_load(yaml.safe_dump(BASE))
        config["profile"]["profiler"] = "metaphlan4"
        requirements = preflight.network_requirements(config, [])
        assert any("MetaPhlAn" in line for line in requirements)

    def test_no_profiler_no_index(self) -> None:
        config = yaml.safe_load(yaml.safe_dump(BASE))
        config["profile"]["profiler"] = "none"
        requirements = preflight.network_requirements(config, [])
        assert not any("MetaPhlAn" in line for line in requirements)

    def test_minpath_script_listed(self) -> None:
        config = yaml.safe_load(yaml.safe_dump(BASE))
        config["annotation"]["minpath"]["data_dir"] = "/db/minpath"
        requirements = preflight.network_requirements(config, [])
        assert any("MinPath" in line for line in requirements)


class TestRunPreflight:
    def test_aggregates_every_check(
        self, fake_pipeline: Path, tmp_path: Path, stub_snakemake: Path, stub_conda: Path
    ) -> None:
        summary = preflight.run_preflight(
            output_dir=tmp_path / "out",
            config=BASE,
            samples=[],
            pipeline_dir=fake_pipeline,
            cache_dir=tmp_path / "cache",
            force=False,
        )
        assert summary["snakemake"] == SNAKEMAKE_PIN
        assert summary["conda"] == "conda"
        assert summary["pipeline_source"]["source_kind"] == "local_checkout"
        assert any("conda" in line for line in summary["network"])

    def test_output_dir_checked_first(
        self, fake_pipeline: Path, tmp_path: Path, stub_snakemake: Path, stub_conda: Path
    ) -> None:
        out = tmp_path / "out"
        out.mkdir()
        (out / "stale.txt").write_text("x\n")
        with pytest.raises(SkillError) as excinfo:
            preflight.run_preflight(
                output_dir=out,
                config=BASE,
                samples=[],
                pipeline_dir=fake_pipeline,
                cache_dir=tmp_path / "cache",
                force=False,
            )
        assert excinfo.value.error_code == ErrorCode.OUTPUT_DIR_NOT_EMPTY

    def test_tools_can_be_skipped(self, fake_pipeline: Path, tmp_path: Path) -> None:
        summary = preflight.run_preflight(
            output_dir=tmp_path / "out",
            config=BASE,
            samples=[],
            pipeline_dir=fake_pipeline,
            cache_dir=tmp_path / "cache",
            check_tools=False,
        )
        assert summary["snakemake"] is None
        assert summary["conda"] is None
