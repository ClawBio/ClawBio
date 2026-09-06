"""Tests for the GWAS-PRS shared reproducibility bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = SKILL_DIR.parents[1]
sys.path.insert(0, str(SKILL_DIR))
sys.path.insert(0, str(PROJECT_ROOT))

import repro_bundle  # noqa: E402

from clawbio.common.checksums import sha256_file  # noqa: E402


def args_for(
    output_dir: Path,
    input_path: Path,
    *,
    demo: bool = False,
    trait: str | None = None,
    pgs_id: str | None = None,
    panel_id: str | None = None,
):
    """Mirror the argparse.Namespace ``gwas_prs.main()`` builds."""
    if not demo and not any([trait, pgs_id, panel_id]):
        trait = "type 2 diabetes"
    return argparse.Namespace(
        demo=demo,
        input=str(input_path),
        trait=trait,
        pgs_id=pgs_id,
        panel_id=panel_id,
        output=str(output_dir),
        min_overlap=0.5,
        max_variants=50000,
        build="GRCh37",
        no_cache=False,
        cache_dir="/private/cache/path",
    )


def make_run(
    tmp_path: Path,
    *,
    demo: bool = False,
    trait: str | None = None,
    pgs_id: str | None = None,
    panel_id: str | None = None,
    scoring_files: list[dict] | None = None,
):
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    input_path = tmp_path / "patient-jane-doe.txt"
    input_path.write_text("synthetic genotype content\n", encoding="utf-8")
    scoring_path = tmp_path / "CLAWBIO-T2D-8_GRCh37.txt"
    scoring_path.write_text("synthetic scoring content\n", encoding="utf-8")

    output_paths = []
    for name, content in (
        ("prs_report.md", "report\n"),
        ("prs_results.json", "[]\n"),
        ("prs_variants.csv", "header\n"),
        ("result.json", "{}\n"),
    ):
        path = output_dir / name
        path.write_text(content, encoding="utf-8")
        output_paths.append(path)

    args = args_for(
        output_dir, input_path, demo=demo, trait=trait, pgs_id=pgs_id, panel_id=panel_id,
    )
    if scoring_files is None:
        scoring_files = [
            {
                "score_id": "CLAWBIO-T2D-8",
                "pgs_id": None,
                "trait": "Type 2 diabetes",
                "filepath": scoring_path,
                "curated_demo_panel": True,
                "curated_panel_id": "CLAWBIO-T2D-8",
                "legacy_pgs_id": "PGS000013",
                "legacy_pgs_compatibility": False,
                "pgs_catalog_id": None,
            }
        ]
    paths = repro_bundle.create_reproducibility_bundle(
        output_dir=output_dir,
        input_path=input_path,
        input_info={"format": "23andme", "total_snps": 1},
        scoring_files=scoring_files,
        args=args,
        output_paths=output_paths,
    )
    return output_dir, input_path, scoring_path, paths


def flat(commands: str) -> str:
    return commands.replace(" \\\n  ", " ")


def test_bundle_uses_the_shared_reproducibility_layer(tmp_path) -> None:
    output_dir, _input_path, _scoring_path, paths = make_run(tmp_path)

    assert paths == {
        "commands": output_dir / "reproducibility" / "commands.sh",
        "environment": output_dir / "reproducibility" / "environment.yml",
        "provenance": output_dir / "reproducibility" / "provenance.json",
        "checksums": output_dir / "reproducibility" / "checksums.sha256",
    }
    assert repro_bundle.write_checksums.__module__ == ("clawbio.common.reproducibility")
    assert repro_bundle.write_environment_yml.__module__ == (
        "clawbio.common.reproducibility"
    )
    assert repro_bundle.write_portable_commands_sh.__module__ == (
        "clawbio.common.reproducibility"
    )


def test_provenance_contains_hashes_without_private_paths_or_raw_input(
    tmp_path,
) -> None:
    _output_dir, input_path, scoring_path, paths = make_run(tmp_path)

    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))
    serialized = json.dumps(provenance, sort_keys=True)

    assert provenance["input"] == {
        "sha256": sha256_file(input_path),
        "format": "23andme",
        "total_snps": 1,
    }
    assert provenance["scoring_files"] == [
        {
            "score_id": "CLAWBIO-T2D-8",
            "pgs_id": None,
            "filename": scoring_path.name,
            "sha256": sha256_file(scoring_path),
            "curated_demo_panel": True,
            "curated_panel_id": "CLAWBIO-T2D-8",
            "legacy_pgs_id": "PGS000013",
            "legacy_pgs_compatibility": False,
            "pgs_catalog_id": None,
        }
    ]
    assert provenance["parameters"]["selection"] == {
        "mode": "trait",
        "query_sha256": ("sha256:" + hashlib.sha256(b"type 2 diabetes").hexdigest()),
    }
    assert str(tmp_path) not in serialized
    assert input_path.name not in serialized
    assert "synthetic genotype content" not in serialized
    assert "type 2 diabetes" not in serialized
    assert "/private/cache/path" not in serialized


def test_provenance_keeps_a_null_pgs_id_null(tmp_path) -> None:
    """Curated panels carry ``pgs_id: None``; it must not become the string
    "None" (which it did before ``score_id`` was threaded through)."""
    _output_dir, _input_path, _scoring_path, paths = make_run(tmp_path, demo=True)

    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))
    assert provenance["scoring_files"][0]["pgs_id"] is None
    assert '"None"' not in paths["provenance"].read_text(encoding="utf-8")


def test_provenance_sorts_scoring_files_by_score_id(tmp_path) -> None:
    scoring_files = []
    for score_id in ("CLAWBIO-T2D-8", "CLAWBIO-AF-12", "PGS000031"):
        path = tmp_path / f"{score_id}_GRCh37.txt"
        path.write_text(f"{score_id}\n", encoding="utf-8")
        scoring_files.append({
            "score_id": score_id,
            "pgs_id": score_id if score_id.startswith("PGS") else None,
            "filepath": path,
        })
    _output_dir, _input_path, _scoring_path, paths = make_run(
        tmp_path, demo=True, scoring_files=scoring_files,
    )
    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))
    assert [item["score_id"] for item in provenance["scoring_files"]] == [
        "CLAWBIO-AF-12", "CLAWBIO-T2D-8", "PGS000031",
    ]


def test_environment_declares_every_eager_runtime_dependency(tmp_path) -> None:
    _output_dir, _input_path, _scoring_path, paths = make_run(tmp_path)

    environment = paths["environment"].read_text(encoding="utf-8")

    for dependency in repro_bundle.REPLAY_PIP_DEPENDENCIES:
        assert dependency in environment
    assert set(repro_bundle.REPLAY_PIP_DEPENDENCIES) == {
        "requests>=2.31",
        "opentelemetry-sdk>=1.20,<2",
        "numpy>=1.24",
        "pandas>=2.0",
    }


def test_numpy_and_pandas_really_are_eager_imports() -> None:
    """The skill never calls numpy or pandas, but importing ``clawbio.common``
    (for checksums and the reproducibility layer) loads its package
    ``__init__``, which imports ``scrna_io`` and therefore both libraries.
    If this ever stops being true, drop them from REPLAY_PIP_DEPENDENCIES."""
    assert "clawbio.common" in sys.modules
    assert "numpy" in sys.modules
    assert "pandas" in sys.modules
    assert "opentelemetry" in sys.modules


def test_commands_are_portable_and_require_input_file(tmp_path) -> None:
    _output_dir, input_path, _scoring_path, paths = make_run(tmp_path)

    commands = paths["commands"].read_text(encoding="utf-8")
    flat_commands = flat(commands)

    assert "CLAWBIO_ROOT:=/path/to/ClawBio" in commands
    assert '"$CLAWBIO_ROOT/skills/gwas-prs/gwas_prs.py"' in commands
    assert (
        "${INPUT_FILE:?Set INPUT_FILE to the genotype file used for this run}"
        in commands
    )
    assert "${TRAIT_QUERY:?Set TRAIT_QUERY to the trait used for this run}" in commands
    assert '--input "${INPUT_FILE}"' in flat_commands
    assert '--trait "${TRAIT_QUERY}"' in flat_commands
    assert '--output "$OUTPUT_DIR"' in flat_commands
    assert '"${PGS_CACHE_DIR:-$HOME/.clawbio/pgs_cache}"' in flat_commands
    assert "type 2 diabetes" not in commands
    assert str(tmp_path) not in commands
    assert input_path.name not in commands


def test_demo_command_needs_no_private_input_path(tmp_path) -> None:
    _output_dir, _input_path, _scoring_path, paths = make_run(tmp_path, demo=True)

    commands = paths["commands"].read_text(encoding="utf-8")

    assert "--demo" in commands
    assert "INPUT_FILE" not in commands
    assert "--trait" not in commands
    assert "--pgs-id" not in commands
    assert "--panel-id" not in commands


def test_pgs_id_is_shell_quoted_in_replay_command(tmp_path) -> None:
    _output_dir, _input_path, _scoring_path, paths = make_run(
        tmp_path,
        pgs_id="PGS000013; echo unsafe",
    )

    commands = paths["commands"].read_text(encoding="utf-8")

    assert "--pgs-id 'PGS000013; echo unsafe'" in flat(commands)
    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))
    assert provenance["parameters"]["selection"] == {
        "mode": "pgs_id",
        "pgs_id": "PGS000013; echo unsafe",
    }


def test_panel_id_run_replays_the_panel_and_records_it(tmp_path) -> None:
    """``--panel-id`` (added in #380) is the normal way to score a bundled
    panel; the bundle must replay and record it rather than fall through."""
    _output_dir, _input_path, _scoring_path, paths = make_run(
        tmp_path, panel_id="CLAWBIO-T2D-8",
    )

    commands = paths["commands"].read_text(encoding="utf-8")
    flat_commands = flat(commands)
    assert "--panel-id CLAWBIO-T2D-8" in flat_commands
    assert "--pgs-id" not in flat_commands
    assert "--trait" not in flat_commands
    assert "TRAIT_QUERY" not in commands
    assert '--input "${INPUT_FILE}"' in flat_commands

    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))
    assert provenance["parameters"]["selection"] == {
        "mode": "panel_id",
        "panel_id": "CLAWBIO-T2D-8",
    }


@pytest.mark.parametrize(
    "selectors, expected_mode, expected_flag",
    [
        # main() resolves --demo, then --panel-id, then --pgs-id, then --trait.
        # The CLI rejects more than one selector, but the bundle must not
        # depend on that: whichever branch main() would take is the one the
        # replay command and provenance must describe.
        ({"demo": True, "panel_id": "CLAWBIO-T2D-8", "pgs_id": "PGS000031", "trait": "t2d"},
         "demo", "--demo"),
        ({"panel_id": "CLAWBIO-T2D-8", "pgs_id": "PGS000031", "trait": "t2d"},
         "panel_id", "--panel-id CLAWBIO-T2D-8"),
        ({"pgs_id": "PGS000031", "trait": "t2d"}, "pgs_id", "--pgs-id PGS000031"),
        ({"trait": "t2d"}, "trait", '--trait "${TRAIT_QUERY}"'),
    ],
)
def test_selection_precedence_mirrors_main(
    tmp_path, selectors, expected_mode, expected_flag,
) -> None:
    _output_dir, _input_path, _scoring_path, paths = make_run(tmp_path, **selectors)

    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))
    assert provenance["parameters"]["selection"]["mode"] == expected_mode

    flat_commands = flat(paths["commands"].read_text(encoding="utf-8"))
    assert expected_flag in flat_commands
    for other in ("--demo", "--panel-id", "--pgs-id", "--trait"):
        if not expected_flag.startswith(other):
            assert other not in flat_commands, other


def test_selection_precedence_matches_the_order_in_main_source() -> None:
    """Pin the precedence to ``gwas_prs.py`` itself: the ``if/elif`` chain in
    ``main()`` must test the selectors in the same order ``_selection`` does."""
    source = (SKILL_DIR / "gwas_prs.py").read_text(encoding="utf-8")
    main_body = source[source.index("def main():"):]
    order_in_main = [
        name
        for name, _pos in sorted(
            (
                (name, main_body.index(f"{keyword} args.{name}:"))
                for name, keyword in (
                    ("demo", "if"), ("panel_id", "elif"), ("pgs_id", "elif"), ("trait", "elif"),
                )
            ),
            key=lambda item: item[1],
        )
    ]
    bundle_source = (SKILL_DIR / "repro_bundle.py").read_text(encoding="utf-8")
    selection_body = bundle_source[bundle_source.index("def _selection("):]
    order_in_bundle = [
        name
        for name, _pos in sorted(
            ((name, selection_body.index(f"args.{name}")) for name in order_in_main),
            key=lambda item: item[1],
        )
    ]
    assert order_in_main == ["demo", "panel_id", "pgs_id", "trait"]
    assert order_in_bundle == order_in_main


def test_no_selector_is_an_error_not_a_silent_pgs_id_of_none(tmp_path) -> None:
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    args = args_for(output_dir, tmp_path / "in.txt")
    args.trait = None
    with pytest.raises(ValueError, match="no score selector"):
        repro_bundle._safe_parameters(args)


def test_checksum_manifest_covers_outputs_and_bundle_metadata(tmp_path) -> None:
    output_dir, _input_path, _scoring_path, paths = make_run(tmp_path)

    entries = {}
    for line in paths["checksums"].read_text(encoding="utf-8").splitlines():
        digest, label = line.split("  ", 1)
        entries[label] = digest

    expected = {
        "prs_report.md",
        "prs_results.json",
        "prs_variants.csv",
        "result.json",
        "reproducibility/commands.sh",
        "reproducibility/environment.yml",
        "reproducibility/provenance.json",
    }
    assert set(entries) == expected
    for label, digest in entries.items():
        assert digest == sha256_file(output_dir / label)


def test_bundle_text_files_use_lf_line_endings(tmp_path) -> None:
    _output_dir, _input_path, _scoring_path, paths = make_run(tmp_path)

    for path in paths.values():
        assert b"\r" not in path.read_bytes(), path


def test_demo_cli_writes_the_documented_output_contract(tmp_path) -> None:
    output_dir = tmp_path / "demo"

    subprocess.run(
        [
            sys.executable,
            str(SKILL_DIR / "gwas_prs.py"),
            "--demo",
            "--output",
            str(output_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
        timeout=120,
    )

    expected = {
        "prs_report.md",
        "prs_results.json",
        "prs_variants.csv",
        "result.json",
        "reproducibility/commands.sh",
        "reproducibility/environment.yml",
        "reproducibility/provenance.json",
        "reproducibility/checksums.sha256",
    }
    actual = {
        path.relative_to(output_dir).as_posix()
        for path in output_dir.rglob("*")
        if path.is_file()
    }
    assert actual == expected

    results = json.loads((output_dir / "prs_results.json").read_text(encoding="utf-8"))
    assert set(results[0]) == {
        "evidence_assessment",
        "interpretation_scope",
        "score_id",
        "pgs_id",
        "curated_panel_id",
        "legacy_pgs_id",
        "legacy_pgs_compatibility",
        "curated_demo_panel",
        "pgs_catalog_id",
        "trait",
        "raw_score",
        "variants_used",
        "variants_total",
        "overlap_fraction",
        "percentile",
        "risk_category",
        "z_score",
        "method",
        "reference_population",
    }

    skill_text = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
    for path in expected:
        assert path.split("/", 1)[-1] in skill_text
    for field in results[0]:
        assert f"| {field} |" in skill_text, field
    assert "tables/scores.csv" not in skill_text

    provenance = json.loads(
        (output_dir / "reproducibility" / "provenance.json").read_text(encoding="utf-8")
    )
    assert provenance["parameters"]["selection"] == {"mode": "demo"}
    assert {item["score_id"] for item in provenance["scoring_files"]} == {
        result["score_id"] for result in results
    }
    assert all(item["pgs_id"] is None for item in provenance["scoring_files"])


def test_version_and_ci_contracts_stay_in_sync() -> None:
    frontmatter = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8").split("---", 2)[1]
    declared = re.search(r"^\s*version:\s*(\S+)\s*$", frontmatter, re.MULTILINE)
    assert declared is not None
    catalog = json.loads(
        (PROJECT_ROOT / "skills" / "catalog.json").read_text(encoding="utf-8")
    )
    catalog_entry = next(
        item for item in catalog["skills"] if item["name"] == "gwas-prs"
    )
    workflow = (PROJECT_ROOT / ".github" / "workflows" / "ci.yml").read_text(
        encoding="utf-8"
    )

    assert re.fullmatch(r"\d+\.\d+\.\d+", repro_bundle.TOOL_VERSION)
    assert repro_bundle.TOOL_VERSION == declared.group(1)
    assert catalog_entry["version"] == repro_bundle.TOOL_VERSION
    assert "uv run pytest skills/gwas-prs/tests/ -v" in workflow


def test_skill_version_parser_tolerates_reindented_frontmatter(tmp_path, monkeypatch) -> None:
    skill_md = tmp_path / "SKILL.md"
    skill_md.write_text(
        "---\nname: gwas-prs\nmetadata:\n    version: 9.9.9\n---\n# body\nversion: 0.0.0\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(repro_bundle, "__file__", str(tmp_path / "repro_bundle.py"))
    assert repro_bundle._read_skill_version() == "9.9.9"
