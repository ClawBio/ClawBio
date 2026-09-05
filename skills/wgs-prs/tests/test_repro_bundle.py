"""Tests for the WGS-PRS shared reproducibility bundle."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

SKILL_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = SKILL_DIR.parents[1]
sys.path.insert(0, str(SKILL_DIR))
sys.path.insert(0, str(PROJECT_ROOT))

import repro_bundle  # noqa: E402
from wgs_prs import BridgeConfig, WgsToPrsBridge  # noqa: E402
from clawbio.common.checksums import sha256_file  # noqa: E402
from clawbio.common.vcf_qc import QcConfig  # noqa: E402


def flat(commands: str) -> str:
    return commands.replace(" \\\n  ", " ")


def make_config(**overrides) -> BridgeConfig:
    cfg = BridgeConfig(
        output_dir="wgs_prs_output",
        sample_id="SAMPLE",
        sex="XX",
        fail_fast=True,
        dry_run=False,
    )
    for key, value in overrides.items():
        setattr(cfg, key, value)
    return cfg


def make_run(
    tmp_path: Path,
    *,
    demo: bool = False,
    input_vcf: Path | None = None,
    fastq_r1: Path | None = None,
    fastq_r2: Path | None = None,
    trait: str | None = None,
    pgs_id: str | None = None,
    panel_id: str = "",
    dry_run: bool = False,
):
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    if input_vcf is None and not demo and fastq_r1 is None:
        input_vcf = tmp_path / "patient-jane-doe.vcf"
        input_vcf.write_text("synthetic vcf content\n", encoding="utf-8")
    if fastq_r1 is not None and not fastq_r1.exists():
        fastq_r1.write_bytes(b"fastq-r1")
    if fastq_r2 is not None and not fastq_r2.exists():
        fastq_r2.write_bytes(b"fastq-r2")

    report_md = output_dir / "bridge_report.md"
    report_json = output_dir / "bridge_report.json"
    report_md.write_text("report\n", encoding="utf-8")
    report_json.write_text("{}\n", encoding="utf-8")

    cfg = make_config(
        demo=demo,
        dry_run=dry_run,
        panel_id=panel_id,
        prs_traits=[trait] if trait else [],
        pgs_ids=[pgs_id] if pgs_id else [],
    )
    paths = repro_bundle.create_reproducibility_bundle(
        output_dir=output_dir,
        config=cfg,
        fastq_r1=fastq_r1,
        fastq_r2=fastq_r2,
        input_vcf=input_vcf,
        output_paths=[report_md, report_json],
    )
    return output_dir, input_vcf, fastq_r1, paths, cfg


def test_bundle_uses_the_shared_reproducibility_layer(tmp_path) -> None:
    output_dir, _vcf, _r1, paths, _cfg = make_run(tmp_path)

    assert paths == {
        "commands": output_dir / "reproducibility" / "commands.sh",
        "environment": output_dir / "reproducibility" / "environment.yml",
        "provenance": output_dir / "reproducibility" / "provenance.json",
        "checksums": output_dir / "reproducibility" / "checksums.sha256",
    }
    assert repro_bundle.write_checksums.__module__ == "clawbio.common.reproducibility"
    assert repro_bundle.write_environment_yml.__module__ == "clawbio.common.reproducibility"
    assert repro_bundle.write_portable_commands_sh.__module__ == (
        "clawbio.common.reproducibility"
    )


def test_vcf_provenance_hashes_without_private_paths(tmp_path) -> None:
    _output_dir, input_vcf, _r1, paths, _cfg = make_run(tmp_path)

    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))
    serialized = json.dumps(provenance, sort_keys=True)

    assert provenance["input"]["entry"] == "vcf"
    assert provenance["input"]["files"] == [
        {"role": "input_vcf", "sha256": sha256_file(input_vcf)},
    ]
    assert str(tmp_path) not in serialized
    assert input_vcf.name not in serialized
    assert "synthetic vcf content" not in serialized
    assert "jane-doe" not in serialized


def test_trait_is_fingerprinted_not_stored(tmp_path) -> None:
    _output_dir, input_vcf, _r1, paths, _cfg = make_run(
        tmp_path, trait="type 2 diabetes",
    )

    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))
    commands = paths["commands"].read_text(encoding="utf-8")

    assert provenance["parameters"]["selection"] == {
        "mode": "trait",
        "query_sha256": "sha256:" + hashlib.sha256(b"type 2 diabetes").hexdigest(),
    }
    assert "type 2 diabetes" not in json.dumps(provenance)
    assert "type 2 diabetes" not in commands
    assert str(tmp_path) not in commands
    assert input_vcf.name not in commands
    assert (
        "${TRAIT_QUERY:?Set TRAIT_QUERY to the trait used for this run}"
        in commands
    )


def test_vcf_commands_require_input_env_and_are_portable(tmp_path) -> None:
    _output_dir, input_vcf, _r1, paths, _cfg = make_run(tmp_path)

    commands = paths["commands"].read_text(encoding="utf-8")
    flat_commands = flat(commands)

    assert "CLAWBIO_ROOT:=/path/to/ClawBio" in commands
    assert '"$CLAWBIO_ROOT/skills/wgs-prs/wgs_prs.py"' in commands
    assert '${INPUT_VCF:?Set INPUT_VCF to the VCF used for this run}' in commands
    assert '--input-vcf "${INPUT_VCF}"' in flat_commands
    assert '--output-dir "$OUTPUT_DIR"' in flat_commands
    assert str(tmp_path) not in commands
    assert input_vcf.name not in commands


def test_fastq_commands_do_not_embed_read_paths(tmp_path) -> None:
    r1 = tmp_path / "private" / "patient_R1.fastq.gz"
    r2 = tmp_path / "private" / "patient_R2.fastq.gz"
    r1.parent.mkdir()
    output_dir, _vcf, _r1, paths, _cfg = make_run(
        tmp_path, input_vcf=None, fastq_r1=r1, fastq_r2=r2,
    )

    commands = paths["commands"].read_text(encoding="utf-8")
    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))
    serialized = json.dumps(provenance)

    assert provenance["input"]["entry"] == "fastq"
    assert {item["role"] for item in provenance["input"]["files"]} == {
        "fastq_r1",
        "fastq_r2",
    }
    assert "patient_R1.fastq.gz" not in commands
    assert "patient_R2.fastq.gz" not in commands
    assert str(tmp_path) not in commands
    assert "patient_R1" not in serialized
    assert '${FASTQ_R1:?Set FASTQ_R1 to the forward reads used for this run}' in commands
    assert '${FASTQ_R2:?Set FASTQ_R2 to the reverse reads used for this run}' in commands
    assert "--input-vcf" not in commands


def test_demo_command_needs_no_private_input_path(tmp_path) -> None:
    _output_dir, _vcf, _r1, paths, _cfg = make_run(tmp_path, demo=True)

    commands = paths["commands"].read_text(encoding="utf-8")
    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))

    assert "--demo" in commands
    assert "INPUT_VCF" not in commands
    assert "FASTQ_R1" not in commands
    assert provenance["input"]["entry"] == "demo"
    assert provenance["parameters"]["selection"]["mode"] == "demo"


def test_pgs_id_is_shell_quoted_in_replay_command(tmp_path) -> None:
    _output_dir, _vcf, _r1, paths, _cfg = make_run(
        tmp_path, pgs_id="PGS000013; echo unsafe",
    )

    commands = paths["commands"].read_text(encoding="utf-8")
    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))

    assert "--pgs-id 'PGS000013; echo unsafe'" in flat(commands)
    assert provenance["parameters"]["selection"] == {
        "mode": "pgs_id",
        "pgs_id": "PGS000013; echo unsafe",
    }


def test_panel_id_run_replays_the_panel(tmp_path) -> None:
    _output_dir, _vcf, _r1, paths, _cfg = make_run(
        tmp_path, panel_id="CLAWBIO-T2D-8",
    )

    commands = paths["commands"].read_text(encoding="utf-8")
    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))

    assert "--panel-id CLAWBIO-T2D-8" in flat(commands)
    assert "--pgs-id" not in flat(commands)
    assert "--trait" not in flat(commands)
    assert provenance["parameters"]["selection"] == {
        "mode": "panel_id",
        "panel_id": "CLAWBIO-T2D-8",
    }


def test_selection_precedence_mirrors_stage3(tmp_path) -> None:
    source = (SKILL_DIR / "wgs_prs.py").read_text(encoding="utf-8")
    stage3 = source[source.index("def _run_stage3("):]
    stage3 = stage3[: stage3.index("def _run_stage4(")]
    order_in_stage3 = [
        name
        for name, _pos in sorted(
            (
                (name, stage3.index(needle))
                for name, needle in (
                    ("pgs_ids", "if self.config.pgs_ids:"),
                    ("prs_traits", "elif self.config.prs_traits:"),
                    ("panel_id", "elif self.config.panel_id:"),
                )
            ),
            key=lambda item: item[1],
        )
    ]
    bundle_source = (SKILL_DIR / "repro_bundle.py").read_text(encoding="utf-8")
    selection_body = bundle_source[bundle_source.index("def _selection("):]
    selection_body = selection_body[: selection_body.index("def _repro_command(")]
    order_in_bundle = [
        name
        for name, _pos in sorted(
            (
                (name, selection_body.index(needle))
                for name, needle in (
                    ("pgs_ids", "pgs_ids = list"),
                    ("prs_traits", "traits = list"),
                    ("panel_id", 'panel_id = str(getattr(config, "panel_id"'),
                )
            ),
            key=lambda item: item[1],
        )
    ]
    assert order_in_stage3 == ["pgs_ids", "prs_traits", "panel_id"]
    assert order_in_bundle == order_in_stage3


def test_environment_declares_eager_runtime_dependencies(tmp_path) -> None:
    _output_dir, _vcf, _r1, paths, _cfg = make_run(tmp_path)
    environment = paths["environment"].read_text(encoding="utf-8")
    for dependency in repro_bundle.REPLAY_PIP_DEPENDENCIES:
        assert dependency in environment


def test_checksum_manifest_covers_outputs_and_bundle_metadata(tmp_path) -> None:
    output_dir, _vcf, _r1, paths, _cfg = make_run(tmp_path)

    entries = {}
    for line in paths["checksums"].read_text(encoding="utf-8").splitlines():
        digest, label = line.split("  ", 1)
        entries[label] = digest

    expected = {
        "bridge_report.md",
        "bridge_report.json",
        "reproducibility/commands.sh",
        "reproducibility/environment.yml",
        "reproducibility/provenance.json",
    }
    assert set(entries) == expected
    for label, digest in entries.items():
        assert digest == sha256_file(output_dir / label)


def test_bundle_text_files_use_lf_line_endings(tmp_path) -> None:
    _output_dir, _vcf, _r1, paths, _cfg = make_run(tmp_path)
    for path in paths.values():
        assert b"\r" not in path.read_bytes(), path


def test_vcf_entry_run_writes_the_bundle(tmp_path) -> None:
    from wgs_prs import write_demo_vcf

    vcf = write_demo_vcf(tmp_path / "test.vcf")
    cfg = BridgeConfig(
        output_dir=str(tmp_path / "out"),
        sample_id="TEST",
        fail_fast=False,
        qc=QcConfig(min_snp_count=50),
    )
    with patch("clawbio.common.vcf_qc.shutil.which", return_value=None):
        bridge = WgsToPrsBridge(cfg)
        report = bridge.run(input_vcf=vcf)

    repro = Path(report.output_dir) / "reproducibility"
    for artefact in ("commands.sh", "environment.yml", "provenance.json", "checksums.sha256"):
        assert (repro / artefact).exists(), artefact
    commands = (repro / "commands.sh").read_text(encoding="utf-8")
    assert "INPUT_VCF" in commands
    assert str(vcf) not in commands
    assert vcf.name not in commands


def test_demo_cli_writes_the_documented_output_contract(tmp_path) -> None:
    output_dir = tmp_path / "demo"
    result = subprocess.run(
        [
            sys.executable,
            str(SKILL_DIR / "wgs_prs.py"),
            "--demo",
            "--output-dir",
            str(output_dir),
        ],
        capture_output=True,
        text=True,
        timeout=120,
        cwd=str(PROJECT_ROOT),
    )
    assert (output_dir / "bridge_report.json").exists(), result.stderr
    assert (output_dir / "bridge_report.md").exists()
    repro = output_dir / "reproducibility"
    for artefact in ("commands.sh", "environment.yml", "provenance.json", "checksums.sha256"):
        assert (repro / artefact).exists(), artefact
    commands = (repro / "commands.sh").read_text(encoding="utf-8")
    assert "--demo" in commands
    assert "INPUT_VCF" not in commands
    provenance = json.loads((repro / "provenance.json").read_text(encoding="utf-8"))
    assert provenance["input"]["entry"] == "demo"
    assert provenance["tool"]["version"] == "0.2.0"
