from __future__ import annotations

import asyncio
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml

SKILL_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = SKILL_DIR.parents[1]
FIXTURES = Path(__file__).resolve().parent / "fixtures"
DEMO_VCF = SKILL_DIR / "examples" / "demo_patient.vcf"
sys.path.insert(0, str(SKILL_DIR))

import just_prs_mcp_bridge as bridge  # noqa: E402
from mcp_client import McpCallError, to_jsonable, tool_result_data  # noqa: E402


class FakeClient:
    def __init__(self, trait_report: dict[str, Any]) -> None:
        self.trait_report = trait_report
        self.calls: list[tuple[str, dict[str, Any]]] = []

    async def call_tool(self, name: str, arguments: dict[str, Any]) -> Any:
        self.calls.append((name, arguments))
        if name == "search_traits":
            return [
                {
                    "id": "EFO_0001360",
                    "label": "type 2 diabetes mellitus",
                    "n_associated": 20,
                    "n_child_associated": 5,
                }
            ]
        if name == "compute_prs_by_trait":
            return self.trait_report
        if name == "percentile":
            return {
                "pgs_id": arguments["pgs_id"],
                "percentile": 78.0,
                "z_score": 0.772,
                "reliable": True,
                "method": "reference_panel",
                "reference_panel_ancestry": arguments["superpopulation"],
            }
        if name == "absolute_risk" and arguments["pgs_id"] == "PGS000014":
            return {
                "absolute_risk": 0.21,
                "population_prevalence": 0.12,
                "risk_ratio": 1.75,
                "method": "or_per_sd",
                "confidence": "high",
                "prevalence_source": "curated",
                "prevalence_type": "lifetime",
                "caveats": [],
            }
        if name == "absolute_risk":
            raise McpCallError("No absolute-risk estimate available.")
        raise AssertionError(f"Unexpected tool call: {name}")


def load_trait_report() -> dict[str, Any]:
    return json.loads((FIXTURES / "demo_trait_report.json").read_text(encoding="utf-8"))


def test_server_launch_is_pinned_and_local(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("PRS_CACHE_DIR", "/tmp/clawbio-prs-cache")
    monkeypatch.setenv("UNRELATED_SERVICE_API_KEY", "must-not-reach-child")
    command = bridge.server_command()
    env = bridge.server_environment()

    assert command == [
        "uvx",
        "just-prs-mcp@0.3.1",
        "stdio",
        "--mode",
        "essentials",
    ]
    assert env["PRS_MCP_MODE"] == "essentials"
    assert env["PRS_CACHE_DIR"] == "/tmp/clawbio-prs-cache"
    assert "UNRELATED_SERVICE_API_KEY" not in env
    assert all("http://" not in value and "https://" not in value for value in command)


def test_omitted_superpopulation_defaults_to_eur_with_disclosure_flag() -> None:
    resolved, defaulted = bridge.normalize_superpopulation(None)
    assert resolved == "EUR"
    assert defaulted is True

    explicit, explicit_defaulted = bridge.normalize_superpopulation("eur")
    assert explicit == "EUR"
    assert explicit_defaulted is False

    auto, auto_defaulted = bridge.normalize_superpopulation("auto")
    assert auto == "AUTO"
    assert auto_defaulted is False


def test_report_discloses_reference_panel_ancestry_and_eur_default(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    data = bridge.map_trait_report(
        load_trait_report(),
        superpopulation="EUR",
        absolute_risks={},
        superpopulation_defaulted=True,
    )

    result = bridge.write_bundle(
        tmp_path / "ancestry",
        data,
        input_path=DEMO_VCF,
        demo=False,
        superpopulation="EUR",
        superpopulation_defaulted=True,
    )
    report = (tmp_path / "ancestry" / "report.md").read_text(encoding="utf-8")

    assert result["data"]["superpopulation_defaulted"] is True
    assert "Requested ancestry" in report
    assert "EUR (default — not explicitly chosen)" in report
    assert "Reference panel ancestry" in report
    assert "| Reference ancestry |" in report
    assert "| EUR |" in report
    assert "WARNING" in report
    assert "EUR-referenced" in report
    assert "population prevalence baseline" in report.lower()
    assert "**Requested ancestry**: EUR (default — not explicitly chosen)" in report
    assert "**Reference panel ancestry**:" in report

    bridge.emit_superpopulation_default_warning("EUR")
    err = capsys.readouterr().err
    assert "WARNING" in err
    assert "EUR" in err
    assert "explicitly" in err.lower()


def test_explicit_superpopulation_skips_default_warning(tmp_path: Path) -> None:
    data = bridge.map_trait_report(
        load_trait_report(),
        superpopulation="SAS",
        absolute_risks={},
        superpopulation_defaulted=False,
    )
    bridge.write_bundle(
        tmp_path / "explicit",
        data,
        input_path=DEMO_VCF,
        demo=False,
        superpopulation="SAS",
        superpopulation_defaulted=False,
    )
    report = (tmp_path / "explicit" / "report.md").read_text(encoding="utf-8")
    assert "EUR (default — not explicitly chosen)" not in report
    assert "**Requested ancestry**: SAS" in report
    header = report.split("## Interpretation")[0]
    assert "WARNING" not in header


def test_demo_vcf_declares_genotype_format() -> None:
    content = DEMO_VCF.read_text(encoding="utf-8")
    assert "##FORMAT=<ID=GT,Number=1,Type=String" in content
    assert "\tFORMAT\tDEMO\n" in content


def test_fastmcp_root_result_is_unwrapped() -> None:
    class Root:
        def __init__(self) -> None:
            self.root = {"pgs_id": "PGS000014", "score": 0.0}

    assert to_jsonable(Root()) == {"pgs_id": "PGS000014", "score": 0.0}


def test_fastmcp_result_prefers_structured_content_over_dynamic_root() -> None:
    class Result:
        data = object()
        structured_content = {"pgs_id": "PGS000014", "score": 0.0}

    assert tool_result_data(Result()) == {"pgs_id": "PGS000014", "score": 0.0}


def test_fastmcp_root_list_wrapper_is_unwrapped() -> None:
    class Result:
        data = []
        structured_content = {"result": [{"id": "MONDO_0005148"}]}

    assert tool_result_data(Result()) == [{"id": "MONDO_0005148"}]


def test_trait_resolution_rejects_ambiguity() -> None:
    matches = [
        {"id": "EFO_0001360", "label": "type 2 diabetes mellitus"},
        {"id": "EFO_0000400", "label": "diabetes mellitus"},
    ]

    with pytest.raises(bridge.AmbiguousTraitError, match="EFO_0001360"):
        bridge.resolve_trait("diabetes", matches)


def test_trait_report_mapping_preserves_failures_and_honesty_signals() -> None:
    mapped = bridge.map_trait_report(
        load_trait_report(),
        superpopulation="EUR",
        absolute_risks={
            "PGS000014": {
                "status": "available",
                "absolute_risk": 0.21,
                "population_prevalence": 0.12,
                "risk_ratio": 1.75,
            },
            "PGS000013": {
                "status": "unavailable",
                "reason": "No absolute-risk estimate available.",
            },
        },
    )

    assert mapped["trait_id"] == "EFO_0001360"
    assert mapped["filter_summary"].startswith("profile=curated")
    assert mapped["model_agreement"]["reliable_models"] == 1
    assert mapped["model_agreement"]["verdict"] == "insufficient_reliable_models"
    assert mapped["scores"][0]["weight_mass_coverage"] == pytest.approx(0.94)
    assert mapped["scores"][0]["absolute_risk"]["risk_ratio"] == pytest.approx(1.75)
    assert mapped["scores"][1]["absolute_risk"]["status"] == "unavailable"
    assert mapped["scores"][2]["status"] == "failed"
    assert mapped["scores"][2]["error"] == "Scoring file unavailable."


def test_live_workflow_passes_only_local_path_and_public_identifiers() -> None:
    client = FakeClient(load_trait_report())
    result = asyncio.run(
        bridge.analyze_with_client(
            client=client,
            input_path=DEMO_VCF,
            trait="type 2 diabetes",
            trait_id=None,
            pgs_id=None,
            superpopulation="EUR",
            profile="curated",
            top_n=5,
            limit=2,
        )
    )

    by_trait_call = next(call for call in client.calls if call[0] == "compute_prs_by_trait")
    assert by_trait_call[1]["vcf_path"] == str(DEMO_VCF.resolve())
    assert by_trait_call[1]["interpret"] is True
    assert by_trait_call[1]["profile"] == "curated"
    assert by_trait_call[1]["top_n"] == 5
    assert by_trait_call[1]["limit"] == 2
    assert "vcf_content" not in json.dumps(client.calls)
    assert result["scores"][0]["pgs_id"] == "PGS000014"
    assert any(name == "absolute_risk" for name, _ in client.calls)


def test_unreliable_percentile_blocks_absolute_risk() -> None:
    class UnreliablePercentileClient:
        def __init__(self) -> None:
            self.calls: list[str] = []

        async def call_tool(self, name: str, arguments: dict[str, Any]) -> Any:
            self.calls.append(name)
            if name == "percentile":
                return {
                    "pgs_id": arguments["pgs_id"],
                    "percentile": 50.0,
                    "z_score": 0.0,
                    "reliable": False,
                    "caveat": "Weight-mass coverage is below the reliability threshold.",
                }
            if name == "assess_quality":
                return {
                    "quality_label": "Low",
                    "summary": "Percentile is unreliable because coverage is low.",
                }
            if name == "absolute_risk":
                raise AssertionError("Absolute risk must not be computed from an unreliable PRS.")
            raise AssertionError(f"Unexpected tool call: {name}")

    client = UnreliablePercentileClient()
    risks = asyncio.run(
        bridge._absolute_risks_for_rows(
            client,
            [
                {
                    "pgs_id": "PGS000014",
                    "status": "scored",
                    "score": 0.0,
                    "weight_mass_coverage": 0.0,
                }
            ],
            superpopulation="EUR",
        )
    )

    assert risks["PGS000014"]["status"] == "unavailable"
    assert "coverage" in risks["PGS000014"]["reason"].lower()
    assert client.calls == ["percentile", "assess_quality"]


def test_single_score_uses_upstream_reported_trait() -> None:
    class SingleScoreClient:
        async def call_tool(self, name: str, arguments: dict[str, Any]) -> Any:
            if name == "compute_prs":
                return {
                    "pgs_id": arguments["pgs_id"],
                    "score": 0.0,
                    "variants_matched": 0,
                    "variants_total": 100,
                    "match_rate": 0.0,
                    "weight_mass_coverage": 0.0,
                    "detected_genome_build": "GRCh38",
                    "build_mismatch": False,
                    "trait_reported": "Type 2 diabetes",
                    "performance": {
                        "effect_sizes": [
                            {
                                "name_short": "OR",
                                "estimate": 1.55,
                                "ci_lower": 1.4,
                                "ci_upper": 1.7,
                            }
                        ],
                        "class_acc": [{"name_short": "AUROC", "estimate": 0.74}],
                    },
                }
            if name == "percentile":
                return {
                    "pgs_id": arguments["pgs_id"],
                    "percentile": 78.0,
                    "method": "reference_panel",
                    "z_score": 0.772,
                    "reliable": True,
                    "caveat": None,
                    "reference_panel_ancestry": "EUR",
                }
            if name == "assess_quality":
                return {"quality_label": "High", "summary": "Reliable reference percentile."}
            if name == "absolute_risk":
                return {
                    "absolute_risk": 0.21,
                    "population_prevalence": 0.12,
                    "risk_ratio": 1.75,
                }
            raise AssertionError(f"Unexpected tool call: {name}")

    result = asyncio.run(
        bridge.analyze_with_client(
            client=SingleScoreClient(),
            input_path=DEMO_VCF,
            trait=None,
            trait_id=None,
            pgs_id="PGS000014",
            superpopulation="EUR",
            profile="curated",
            top_n=5,
            limit=None,
        )
    )

    assert result["trait"] == "Type 2 diabetes"
    score = result["scores"][0]
    assert score["percentile"] == pytest.approx(78.0)
    assert score["percentile_reliable"] is True
    assert score["reference_panel_ancestry"] == "EUR"
    assert score["effect_size"] == "OR=1.55 [1.40-1.70]"
    assert score["auroc_estimate"] == pytest.approx(0.74)
    assert score["quality_label"] == "High"
    assert score["absolute_risk"]["status"] == "available"


@pytest.mark.parametrize(
    ("superpopulation", "profile"),
    [
        ("AFR", "curated"),
        ("AMR", "curated"),
        ("EAS", "curated"),
        ("EUR", "curated"),
        ("SAS", "curated"),
        ("AFR", "all"),
        ("AMR", "all"),
        ("EAS", "all"),
        ("EUR", "all"),
        ("SAS", "all"),
    ],
)
def test_ten_run_mapping_stress(superpopulation: str, profile: str) -> None:
    client = FakeClient(load_trait_report())
    result = asyncio.run(
        bridge.analyze_with_client(
            client=client,
            input_path=DEMO_VCF,
            trait=None,
            trait_id="EFO_0001360",
            pgs_id=None,
            superpopulation=superpopulation,
            profile=profile,
            top_n=3,
            limit=3,
        )
    )

    request = next(arguments for name, arguments in client.calls if name == "compute_prs_by_trait")
    assert request["superpopulation"] == superpopulation
    assert request["profile"] == profile
    assert all(
        score["requested_superpopulation"] == superpopulation for score in result["scores"]
    )


def test_demo_writes_complete_output_contract(tmp_path: Path) -> None:
    output_dir = tmp_path / "demo"
    result = bridge.run_demo(output_dir)

    expected = {
        output_dir / "report.md",
        output_dir / "result.json",
        output_dir / "tables" / "scores.csv",
        output_dir / "reproducibility" / "commands.sh",
        output_dir / "reproducibility" / "checksums.sha256",
    }
    assert all(path.exists() for path in expected)
    assert result["skill"] == "just-prs-mcp"
    assert result["data"]["scores"]
    report = (output_dir / "report.md").read_text(encoding="utf-8")
    assert "C_wt" in report
    assert "Model agreement" in report
    assert "ClawBio is a research and educational tool" in report
    assert "synthetic" in report.lower()
    commands = (output_dir / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--demo" in commands
    assert "--extra just-prs" in commands
    assert "\n+  " not in commands
    assert "/my/custom/path/" not in commands


def test_replay_script_runs_outside_source_checkout(tmp_path: Path) -> None:
    output_dir = tmp_path / "demo"
    bridge.run_demo(output_dir)
    environment = dict(os.environ)
    environment.pop("VIRTUAL_ENV", None)

    completed = subprocess.run(
        [str(output_dir / "reproducibility" / "commands.sh")],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    commands = (output_dir / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert 'uv run --project "$CLAWBIO_ROOT" --extra just-prs' in commands


def test_all_failed_bundle_reports_failure(tmp_path: Path) -> None:
    report = load_trait_report()
    report["n_scored"] = 0
    report["n_failed"] = len(report["rows"])
    for row in report["rows"]:
        row.update(
            {
                "status": "failed",
                "score": None,
                "percentile": None,
                "percentile_reliable": None,
                "error": "No scoring model completed.",
            }
        )
    data = bridge.map_trait_report(report, superpopulation="EUR", absolute_risks={})

    result = bridge.write_bundle(
        tmp_path / "failed",
        data,
        input_path=DEMO_VCF,
        demo=False,
    )

    assert result["ok"] is False
    assert result["status"] == "failed"


def test_nonempty_output_directory_emits_overwrite_warning(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    output_dir = tmp_path / "existing"
    output_dir.mkdir()
    (output_dir / "keep.txt").write_text("existing", encoding="utf-8")

    bridge.run_demo(output_dir)

    assert "WARNING" in capsys.readouterr().err


def test_live_overwrite_preflight_happens_before_mcp_log(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    events: list[str] = []

    def record_preflight(output_dir: Path) -> None:
        assert output_dir == tmp_path / "live"
        events.append("preflight")

    class StopAfterClientConstruction:
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            assert events == ["preflight"]
            raise RuntimeError("stop after ordering check")

    monkeypatch.setattr(bridge, "_warn_before_overwrite", record_preflight)
    monkeypatch.setattr(bridge, "JustPrsMcpClient", StopAfterClientConstruction)

    with pytest.raises(RuntimeError, match="ordering check"):
        asyncio.run(
            bridge._run_live(
                input_path=DEMO_VCF,
                output_dir=tmp_path / "live",
                trait=None,
                trait_id="EFO_0001360",
                pgs_id=None,
                superpopulation="EUR",
                profile="curated",
                top_n=5,
                limit=1,
                timeout_seconds=300,
            )
        )


def test_skill_metadata_conforms_to_clawbio_template() -> None:
    content = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
    metadata = yaml.safe_load(content.split("---", 2)[1])
    openclaw = metadata["metadata"]["openclaw"]

    assert metadata["name"] == "just-prs-mcp"
    assert metadata["metadata"]["version"] == "0.1.1"
    assert metadata["metadata"]["author"]
    assert metadata["metadata"]["inputs"][0]["required"] is True
    assert metadata["metadata"]["outputs"]
    assert len(openclaw["trigger_keywords"]) >= 3
    for section in (
        "## Trigger",
        "## Scope",
        "## Workflow",
        "## Example Output",
        "## Gotchas",
        "## Safety",
        "## Agent Boundary",
        "## Integration with Bio Orchestrator",
        "## Maintenance",
    ):
        assert section in content
    assert content.count("You will want to") >= 3
    assert len(content.splitlines()) < 500


def test_intents_and_runner_registration() -> None:
    descriptor = json.loads((SKILL_DIR / "INTENTS.json").read_text(encoding="utf-8"))
    from clawbio.cli import SKILLS

    assert descriptor["schema"] == "clawbio.skill_intents.v1"
    assert descriptor["skill"] == "just-prs"
    assert descriptor["routes"][0]["demo_policy"] == "only_when_explicit"
    registered = SKILLS["just-prs"]
    assert registered["script"] == SKILL_DIR / "just_prs_mcp_bridge.py"
    assert registered["default_timeout_seconds"] > 3600
    assert {
        "--trait",
        "--trait-id",
        "--pgs-id",
        "--superpopulation",
        "--profile",
        "--top-n",
        "--limit",
    } <= registered["allowed_extra_flags"]


def test_typer_cli_help_and_validation() -> None:
    script = SKILL_DIR / "just_prs_mcp_bridge.py"
    help_run = subprocess.run(
        [sys.executable, str(script), "--help"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert help_run.returncode == 0
    assert "--trait" in help_run.stdout
    assert "--pgs-id" in help_run.stdout

    invalid = subprocess.run(
        [
            sys.executable,
            str(script),
            "--input",
            str(DEMO_VCF),
            "--output",
            str(SKILL_DIR / "tests" / "_invalid_output"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert invalid.returncode != 0
    assert "trait" in (invalid.stdout + invalid.stderr).lower()
