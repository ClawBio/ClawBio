"""report.md, result.json and the reproducibility bundle.

The report states what ran, what did not and why, and never upgrades a number
into a claim. Three failure modes are tested explicitly: a blank upstream cell
rendered as `not measured` rather than 0, a threshold echoed from the config the
run actually used rather than a value typed into the template, and the demo
caveat appearing only in demo mode.
"""

from __future__ import annotations

import copy
import json
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

SKILL_DIR = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).resolve().parent / "fixtures"
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

from magpipe_config import RunOptions, build_config
from magpipe_outputs import parse_outputs
from magpipe_report import (
    write_report_md,
    write_reproducibility,
    write_result_json,
)
from magpipe_samplesheet import Sample
from magpipe_schemas import DISCLAIMER

BASE = yaml.safe_load((FIXTURES / "upstream_config.yaml").read_text())
BASE_TEST = yaml.safe_load((FIXTURES / "upstream_config_test.yaml").read_text())
COMMAND = "snakemake -s /pipeline/workflow/Snakefile --configfile /out/config.yaml"
SOURCE = {
    "path": "/pipeline",
    "source_kind": "local_checkout",
    "commit": "7351a702d29857801963cae7aeff6600d59ebe77",
    "dirty": False,
    "ref": "7351a702d29857801963cae7aeff6600d59ebe77",
}


def full_config(tmp_path: Path) -> dict:
    db = tmp_path / "dbs"
    (db / "gtdbtk").mkdir(parents=True, exist_ok=True)
    (db / "taxdump").mkdir(parents=True, exist_ok=True)
    for name in ("checkm2.dmnd", "nr.dmnd", "kegg.dmnd", "cog.dmnd", "Pfam-A.hmm"):
        (db / name).write_text("x\n")
    return build_config(
        copy.deepcopy(BASE),
        [Sample("S1", (Path("/r/R1.fq.gz"),), (Path("/r/R2.fq.gz"),), "S1", None, True)],
        RunOptions(
            preset="full",
            checkm2_db=db / "checkm2.dmnd",
            gtdbtk_db=db / "gtdbtk",
            taxonomy_db=db / "nr.dmnd",
            taxdump=db / "taxdump",
            kegg_db=db / "kegg.dmnd",
            cog_db=db / "cog.dmnd",
            pfam_db=db / "Pfam-A.hmm",
        ),
        samplesheet=tmp_path / "samplesheet.tsv",
        results_dir=tmp_path / "results",
        database_dir=db,
    )


def assembly_config(tmp_path: Path) -> dict:
    return build_config(
        copy.deepcopy(BASE),
        [Sample("S1", (Path("/r/R1.fq.gz"),), (Path("/r/R2.fq.gz"),), "S1", None, True)],
        RunOptions(preset="assembly"),
        samplesheet=tmp_path / "samplesheet.tsv",
        results_dir=tmp_path / "results",
        database_dir=tmp_path / "db",
    )


@pytest.fixture
def report(tmp_path: Path, fake_results: Path) -> tuple[str, dict, dict]:
    config = full_config(tmp_path)
    findings = parse_outputs(fake_results, config)
    path = write_report_md(
        tmp_path,
        status="ok",
        mode="run",
        findings=findings,
        config=config,
        pipeline_source=dict(SOURCE),
        warnings=[],
    )
    return path.read_text(), findings, config


class TestReportMd:
    def test_written_to_report_md(self, tmp_path: Path, report) -> None:
        assert (tmp_path / "report.md").is_file()

    def test_contains_the_disclaimer(self, report) -> None:
        text, _, _ = report
        assert DISCLAIMER in " ".join(text.split())

    def test_states_the_pipeline_commit_and_mode(self, report) -> None:
        text, _, _ = report
        assert "7351a702d29857801963cae7aeff6600d59ebe77" in text
        assert "run" in text

    def test_pipeline_line_does_not_repeat_the_commit(self, report) -> None:
        text, _, _ = report
        line = next(line for line in text.splitlines() if line.startswith("**Pipeline**"))
        assert line.count("7351a702d29857801963cae7aeff6600d59ebe77") == 1

    def test_lists_stages_that_ran_and_stages_skipped_with_a_reason(self, report) -> None:
        text, _, config = report
        assert "Read QC" in text and "Assembly" in text
        assert "Stages skipped" in text or "skipped" in text
        config["annotation"]["enabled"] = False
        skipped = write_report_md(
            Path("/tmp/unused"), status="ok", mode="run", findings=None, config=config,
            pipeline_source=dict(SOURCE), warnings=[],
        )
        assert "--kegg-db" in skipped.read_text()

    def test_blank_cells_render_as_not_measured(self, report) -> None:
        text, _, _ = report
        assert "not measured" in text
        assert "S2 | single | 2500 | 750000 | not measured" in text

    def test_no_zero_is_printed_for_an_unmeasured_cell(self, report) -> None:
        text, _, _ = report
        for line in text.splitlines():
            if "retained_pairs" in line or "assembly_capture_pct" in line:
                assert "| 0" not in line, line

    def test_thresholds_come_from_the_config(self, tmp_path: Path, fake_results: Path) -> None:
        config = full_config(tmp_path)
        findings = parse_outputs(fake_results, config)
        config["mag"]["checkm2"]["min_completeness"] = 71.0
        config["mag"]["checkm2"]["max_contamination"] = 3.5
        text = write_report_md(
            tmp_path,
            status="ok",
            mode="run",
            findings=findings,
            config=config,
            pipeline_source=dict(SOURCE),
            warnings=[],
        ).read_text()
        assert "71" in text and "3.5" in text
        assert "50.0" not in text and "10.0" not in text

    def test_relative_abundance_is_labelled(self, report) -> None:
        text, _, _ = report
        assert "%" in text
        assert "relative abundance" in text.lower()

    def test_threshold_passing_is_not_called_a_high_quality_mag(self, report) -> None:
        text, _, _ = report
        lowered = text.lower()
        for claim in ("is a high-quality mag", "qualifies as a high-quality",
                      "hq-mag status"):
            assert claim not in lowered, f"report makes a HQ-MAG claim: {claim!r}"
        assert "does not by itself establish a high-quality mag" in lowered

    def test_links_to_multiqc_and_quast_when_present(self, report) -> None:
        text, _, _ = report
        assert "results/reports/multiqc/multiqc_report.html" in text
        assert "results/reports/assembly/quast.html" in text

    def test_reproducibility_pointer(self, report) -> None:
        text, _, _ = report
        assert "reproducibility/" in text
        assert "commands.sh" in text

    def test_warnings_are_rendered(self, tmp_path: Path, fake_results: Path) -> None:
        config = full_config(tmp_path)
        findings = parse_outputs(fake_results, config)
        text = write_report_md(
            tmp_path,
            status="ok",
            mode="run",
            findings=findings,
            config=config,
            pipeline_source=dict(SOURCE),
            warnings=["Reference database X was not checksummed"],
        ).read_text()
        assert "Reference database X was not checksummed" in text

    def test_demo_caveat_only_in_demo_mode(self, tmp_path: Path, fake_results: Path) -> None:
        config = full_config(tmp_path)
        findings = parse_outputs(fake_results, config)
        run_mode = write_report_md(
            tmp_path / "run", status="ok", mode="run", findings=findings, config=config,
            pipeline_source=dict(SOURCE), warnings=[],
        ).read_text()
        demo_mode = write_report_md(
            tmp_path / "demo", status="ok", mode="demo", findings=findings, config=config,
            pipeline_source=dict(SOURCE), warnings=[],
        ).read_text()
        assert "synthetic" not in run_mode.lower()
        assert "synthetic" in demo_mode.lower()
        assert "no biological truth" in demo_mode.lower()

    def test_check_mode_says_nothing_executed(self, tmp_path: Path) -> None:
        config = assembly_config(tmp_path)
        text = write_report_md(
            tmp_path,
            status="planned",
            mode="check",
            findings=None,
            config=config,
            pipeline_source=dict(SOURCE, command=COMMAND),
            warnings=[],
        ).read_text()
        assert "Nothing was executed" in text
        assert "dry run" in text.lower()
        assert COMMAND in text
        assert "Planned stages" in text

    def test_check_mode_lists_the_network_requirements(self, tmp_path: Path) -> None:
        config = assembly_config(tmp_path)
        text = write_report_md(
            tmp_path,
            status="planned",
            mode="check",
            findings=None,
            config=config,
            pipeline_source=dict(SOURCE, command=COMMAND),
            warnings=[],
        ).read_text()
        assert "Conda environments" in text
        assert "GRCh38" in text

    def test_check_mode_has_no_result_tables(self, tmp_path: Path) -> None:
        config = assembly_config(tmp_path)
        text = write_report_md(
            tmp_path,
            status="planned",
            mode="check",
            findings=None,
            config=config,
            pipeline_source=dict(SOURCE, command=COMMAND),
            warnings=[],
        ).read_text()
        assert "Read QC" not in text
        assert "retained_reads" not in text

    def test_error_status_is_reported_without_findings(self, tmp_path: Path) -> None:
        text = write_report_md(
            tmp_path,
            status="error",
            mode="run",
            findings=None,
            config=assembly_config(tmp_path),
            pipeline_source=dict(SOURCE),
            warnings=["snakemake exited with status 1"],
        ).read_text()
        assert "failed" in text.lower()
        assert "snakemake exited with status 1" in text
        assert DISCLAIMER in " ".join(text.split())

    def test_report_is_utf8_and_ends_with_a_newline(self, tmp_path: Path, report) -> None:
        assert report[0].endswith("\n")


class TestResultJson:
    def test_schema(self, tmp_path: Path, fake_results: Path) -> None:
        config = full_config(tmp_path)
        findings = parse_outputs(fake_results, config)
        path = write_result_json(
            tmp_path,
            status="ok",
            mode="run",
            findings=findings,
            config=config,
            pipeline_source=dict(SOURCE),
            command=COMMAND,
            warnings=[],
        )
        data = json.loads(path.read_text())
        for key in (
            "status",
            "mode",
            "skill",
            "version",
            "pipeline",
            "findings",
            "warnings",
            "disclaimer",
            "config",
            "command",
            "generated_at",
        ):
            assert key in data, key
        assert data["status"] == "ok"
        assert data["mode"] == "run"
        assert data["skill"] == "metagenome-mag-pipeline"
        assert data["disclaimer"] == DISCLAIMER
        assert data["pipeline"]["commit"] == SOURCE["commit"]
        assert data["pipeline"]["ref"] == SOURCE["ref"]
        assert data["findings"]["bins"]["count"] == 3

    def test_config_is_the_effective_config(self, tmp_path: Path) -> None:
        config = assembly_config(tmp_path)
        data = json.loads(
            write_result_json(
                tmp_path, status="planned", mode="check", findings=None, config=config,
                pipeline_source=dict(SOURCE), command=COMMAND, warnings=[],
            ).read_text()
        )
        assert data["config"] == config
        assert data["findings"] is None

    def test_json_is_serialisable_with_paths_as_strings(self, tmp_path: Path, fake_results: Path) -> None:
        config = full_config(tmp_path)
        findings = parse_outputs(fake_results, config)
        data = json.loads(
            write_result_json(
                tmp_path, status="ok", mode="run", findings=findings, config=config,
                pipeline_source=dict(SOURCE, path=Path("/pipeline")), command=COMMAND, warnings=[],
            ).read_text()
        )
        assert data["pipeline"]["path"] == "/pipeline"


class TestReproducibility:
    def test_bundle_contents(self, tmp_path: Path, fake_results: Path) -> None:
        config = full_config(tmp_path)
        findings = parse_outputs(fake_results, config)
        write_result_json(
            tmp_path, status="ok", mode="run", findings=findings, config=config,
            pipeline_source=dict(SOURCE), command=COMMAND, warnings=[],
        )
        write_report_md(
            tmp_path, status="ok", mode="run", findings=findings, config=config,
            pipeline_source=dict(SOURCE), warnings=[],
        )
        files = [
            tmp_path / "report.md",
            tmp_path / "result.json",
            tmp_path / "pipeline" / "config.yaml",
            fake_results / "final" / "qc" / "read_counts.tsv",
        ]
        write_reproducibility(tmp_path, command=COMMAND, files=files)
        assert (tmp_path / "reproducibility" / "commands.sh").is_file()
        assert (tmp_path / "reproducibility" / "environment.yml").is_file()
        assert (tmp_path / "reproducibility" / "checksums.sha256").is_file()

    def test_commands_sh_holds_the_command(self, tmp_path: Path) -> None:
        write_reproducibility(tmp_path, command=COMMAND, files=[tmp_path / "report.md"])
        text = (tmp_path / "reproducibility" / "commands.sh").read_text()
        assert COMMAND in text
        assert text.startswith("#!/usr/bin/env bash")

    def test_environment_pins_snakemake_and_the_channels(self, tmp_path: Path) -> None:
        write_reproducibility(tmp_path, command=COMMAND, files=[])
        env = yaml.safe_load((tmp_path / "reproducibility" / "environment.yml").read_text())
        assert env["name"] == "clawbio-metagenome-mag-pipeline"
        assert "conda-forge" in env["channels"] and "bioconda" in env["channels"]
        assert "snakemake=9.11.2" in env["dependencies"]
        assert "pyyaml" in env["dependencies"]

    def test_checksums_verify_with_sha256sum(self, tmp_path: Path, fake_results: Path) -> None:
        config = full_config(tmp_path)
        findings = parse_outputs(fake_results, config)
        write_result_json(
            tmp_path, status="ok", mode="run", findings=findings, config=config,
            pipeline_source=dict(SOURCE), command=COMMAND, warnings=[],
        )
        write_reproducibility(
            tmp_path,
            command=COMMAND,
            files=[tmp_path / "result.json", fake_results / "final" / "qc" / "read_counts.tsv"],
        )
        result = subprocess.run(
            ["sha256sum", "-c", "reproducibility/checksums.sha256"],
            cwd=tmp_path,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr

    def test_checksums_are_relative_to_the_output_dir(self, tmp_path: Path) -> None:
        nested = tmp_path / "results" / "final" / "qc" / "read_counts.tsv"
        nested.parent.mkdir(parents=True)
        nested.write_text("sample\nS1\n")
        write_reproducibility(tmp_path, command=COMMAND, files=[nested])
        text = (tmp_path / "reproducibility" / "checksums.sha256").read_text()
        assert "results/final/qc/read_counts.tsv" in text
        assert not text.strip().endswith(str(nested))

    def test_missing_files_are_skipped_not_fatal(self, tmp_path: Path) -> None:
        write_reproducibility(
            tmp_path, command=COMMAND, files=[tmp_path / "nope.tsv"]
        )
        assert (tmp_path / "reproducibility" / "checksums.sha256").read_text() == ""
