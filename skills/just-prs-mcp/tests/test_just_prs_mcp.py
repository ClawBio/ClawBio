from __future__ import annotations

import asyncio
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest
from inline_snapshot import snapshot
from clawbio.common import demo_runner
import yaml

pytest.importorskip("typer", reason="needs the just-prs extra (uv sync --extra just-prs)")

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


def test_replay_script_does_not_record_the_checkout_path(tmp_path: Path) -> None:
    """Bundles are shared; the generating machine's checkout path must not be in them."""
    output_dir = tmp_path / "demo"
    bridge.run_demo(output_dir)

    commands = (output_dir / "reproducibility" / "commands.sh").read_text(encoding="utf-8")

    assert str(bridge._PROJECT_ROOT) not in commands


def test_replay_script_runs_outside_source_checkout(tmp_path: Path) -> None:
    output_dir = tmp_path / "demo"
    bridge.run_demo(output_dir)
    environment = dict(os.environ)
    environment.pop("VIRTUAL_ENV", None)
    environment["CLAWBIO_ROOT"] = str(bridge._PROJECT_ROOT)

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


# Rich, which renders typer's help, emits ANSI styling when it sees
# GITHUB_ACTIONS=true or FORCE_COLOR, splitting "--trait" across escape codes.
# Every GitHub runner sets GITHUB_ACTIONS, so this test passed locally and
# failed in CI. Run it under each colour trigger so it fails anywhere.
@pytest.mark.parametrize(
    "colour_env",
    [{}, {"GITHUB_ACTIONS": "true"}, {"FORCE_COLOR": "1"}],
    ids=["plain", "github-actions", "force-color"],
)
def test_typer_cli_help_and_validation(colour_env: dict[str, str]) -> None:
    script = SKILL_DIR / "just_prs_mcp_bridge.py"
    env = {k: v for k, v in os.environ.items() if k not in {"GITHUB_ACTIONS", "FORCE_COLOR", "NO_COLOR"}}
    env.update(colour_env)
    help_run = subprocess.run(
        [sys.executable, str(script), "--help"],
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )
    assert help_run.returncode == 0
    help_text = re.sub(r"\x1b\[[0-9;]*m", "", help_run.stdout)
    assert "--trait" in help_text
    assert "--pgs-id" in help_text

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


# ---------------------------------------------------------------------------
# AUTO ancestry: the silent-EUR path
#
# `--superpopulation AUTO` asks the upstream engine which reference panel a
# model was built against. When it cannot say, the bridge falls back to EUR.
# That fallback must be as loud as the omitted-flag fallback already is:
# an undisclosed EUR reference is the equity failure the skill's own Gotchas
# section forbids.
# ---------------------------------------------------------------------------


def load_ancestry_blind_trait_report() -> dict[str, Any]:
    """Trait report whose models do not declare a reference-panel ancestry."""
    report = load_trait_report()
    for row in report.get("rows", []):
        row["reference_panel_ancestry"] = None
    return report


class AncestryBlindClient(FakeClient):
    """Upstream that never reports which panel a percentile came from."""

    async def call_tool(self, name: str, arguments: dict[str, Any]) -> Any:
        result = await super().call_tool(name, arguments)
        if name == "percentile":
            result = dict(result)
            result["reference_panel_ancestry"] = None
        return result


def _auto_run(client: FakeClient) -> dict[str, Any]:
    return asyncio.run(
        bridge.analyze_with_client(
            client=client,
            input_path=DEMO_VCF,
            trait="type 2 diabetes",
            trait_id=None,
            pgs_id=None,
            superpopulation="AUTO",
            profile="curated",
            top_n=5,
            limit=2,
        )
    )


def test_auto_ancestry_records_the_panel_actually_applied() -> None:
    """When upstream does report a panel, AUTO must record which one it used."""
    result = _auto_run(FakeClient(load_trait_report()))
    scored = [s for s in result["scores"] if s.get("status") != "failed"]
    assert scored, "expected at least one scored model"
    for score in scored:
        assert score["requested_superpopulation"] == "AUTO"
        # AUTO resolved from upstream evidence, so it is not a silent default.
        assert score["superpopulation_applied"] == "EUR"
        assert score["superpopulation_auto_defaulted"] is False


def test_auto_ancestry_fallback_to_eur_is_disclosed() -> None:
    """Upstream cannot say which panel: EUR is applied and must be flagged."""
    result = _auto_run(AncestryBlindClient(load_ancestry_blind_trait_report()))
    scored = [s for s in result["scores"] if s.get("status") != "failed"]
    assert scored
    assert any(s["superpopulation_auto_defaulted"] for s in scored)
    for score in scored:
        if score["superpopulation_auto_defaulted"]:
            assert score["superpopulation_applied"] == "EUR"
    assert result["superpopulation_defaulted"] is True


def test_auto_ancestry_fallback_warns_and_banners(tmp_path: Path) -> None:
    """The undisclosed-EUR fallback reaches stderr and the report body."""
    result = _auto_run(AncestryBlindClient(load_ancestry_blind_trait_report()))
    bundle = tmp_path / "out"
    bridge.write_bundle(
        bundle,
        result,
        input_path=DEMO_VCF,
        demo=False,
        superpopulation="AUTO",
        superpopulation_defaulted=result["superpopulation_defaulted"],
    )
    report = (bundle / "report.md").read_text(encoding="utf-8")
    assert "EUR" in report
    assert "AUTO" in report
    assert "WARNING" in report


def test_auto_ancestry_is_not_presented_as_inferred_from_the_sample() -> None:
    """AUTO reads the model's panel; it never infers the patient's ancestry."""
    text = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
    assert "does not" in text.lower()
    # The CLI help must not claim the sample's ancestry is inferred.
    help_text = bridge.SUPERPOPULATION_HELP
    assert "infers ancestry" not in help_text.lower()


def test_build_mismatch_withholds_absolute_risk() -> None:
    """SKILL.md says a build-mismatched score must not be interpreted.

    A coordinate-build mismatch means the upstream matcher was comparing the
    sample against scoring positions from another assembly, so the raw score,
    its percentile and any absolute risk derived from it are not interpretable.
    """
    report = load_trait_report()
    report["build_mismatch"] = True
    mapped = bridge.map_trait_report(
        report,
        superpopulation="EUR",
        absolute_risks={
            "PGS000014": {
                "status": "available",
                "absolute_risk": 0.21,
                "population_prevalence": 0.12,
                "risk_ratio": 1.75,
            }
        },
    )
    assert mapped["build_mismatch"] is True
    for score in mapped["scores"]:
        assert score["absolute_risk"]["status"] != "available"
        if score["absolute_risk"]["status"] == "withheld":
            assert "build" in score["absolute_risk"]["reason"].lower()


def test_build_mismatch_banner_appears_in_report(tmp_path: Path) -> None:
    report = load_trait_report()
    report["build_mismatch"] = True
    mapped = bridge.map_trait_report(
        report, superpopulation="EUR", absolute_risks={}
    )
    bridge.write_bundle(
        tmp_path / "mismatch",
        mapped,
        input_path=DEMO_VCF,
        demo=False,
        superpopulation="EUR",
        superpopulation_defaulted=False,
    )
    body = (tmp_path / "mismatch" / "report.md").read_text(encoding="utf-8")
    assert "build" in body.lower()
    assert "WARNING" in body


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/just-prs-mcp/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("just-prs-mcp", tmp_path) == snapshot({'report.md': """\
# just-prs Polygenic Risk Report

**Date**: <TIMESTAMP> UTC
**Skill**: just-prs-mcp
**Version**: 0.1.1
**Mode**: Synthetic offline demo
**Trait**: type 2 diabetes mellitus
**Trait ontology ID**: EFO_0001360
**Genome build**: GRCh38
**Requested ancestry**: EUR
**Reference panel ancestry**: EUR
**Input files**:
- `demo_patient.vcf`: `56eba18a865171bf2b87e6ff5f5735616bb40b036990ae9714a1ecdddc3a3835`

---

## Summary

Scored 2/3 PGS ID(s) for type 2 diabetes mellitus on GRCh38; 1 failed.

- Scored models: **2**
- Failed models: **1**
- Curated out: **4**
- Build mismatch: **False**
- Requested ancestry: **EUR**
- Reference panel ancestry: **EUR**

## Model agreement

- Reliable models: **1**
- Verdict: **insufficient_reliable_models**
- Reliable percentile range: **78.00–78.00**
- Descriptive percentile spread: **0.00**

This is a descriptive comparison of reliable models, not a clinical confidence threshold.

## Score details

| PGS ID | Status | Percentile | Reference ancestry | Reliable | C_wt | Match rate | Quality | Absolute risk |
|---|---|---:|---|---|---:|---:|---|---|
| [PGS000014](https://www.pgscatalog.org/score/PGS000014/) | scored | 78.00 | EUR | True | 94.0% | 92.0% | High | unavailable |
| [PGS000013](https://www.pgscatalog.org/score/PGS000013/) | scored | 54.00 | EUR | False | 61.0% | 72.5% | Moderate | unavailable |
| [PGS999999](https://www.pgscatalog.org/score/PGS999999/) | failed | n/a | n/a | None | n/a | n/a | n/a | not_attempted |

## Interpretation boundaries

- A PRS measures genetic predisposition, not the trait itself and not a diagnosis.
- Percentiles are meaningful only when the reference-panel ancestry is appropriate.
- Absolute-risk percentages may use a single-population prevalence baseline upstream, so they can be ancestry-inconsistent with an ancestry-stratified percentile even when both numbers are shown.
- Low C_wt, unreliable percentiles, model disagreement, and build mismatches must not be hidden.
- Lifestyle, environment, age, family history, and genetic factors outside the model also matter.

## Curation provenance

profile=curated dropped 4 scored row(s): 2 no performance evidence, 2 duplicate score family

Upstream engine: `just-prs-mcp==0.3.1` via local stdio.

---

## Disclaimer

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
""", 'result.json': """\
{
  "status": "ok",
  "ok": true,
  "skill": "just-prs-mcp",
  "version": "0.1.1",
  "completed_at": "<TIMESTAMP>",
  "input_checksum": "sha256:56eba18a865171bf2b87e6ff5f5735616bb40b036990ae9714a1ecdddc3a3835",
  "datasets": {
    "just-prs-mcp": "0.3.1",
    "PGS Catalog": "cached synthetic fixture"
  },
  "summary": {
    "trait": "type 2 diabetes mellitus",
    "trait_id": "EFO_0001360",
    "models_scored": 2,
    "models_failed": 1,
    "reliable_models": 1,
    "requested_superpopulation": "EUR",
    "superpopulation_defaulted": false,
    "reference_panel_ancestry": "EUR",
    "demo": true
  },
  "data": {
    "trait_id": "EFO_0001360",
    "trait": "type 2 diabetes mellitus",
    "genome_build": "GRCh38",
    "detected_genome_build": "GRCh38",
    "build_mismatch": false,
    "profile": "curated",
    "requested_superpopulation": "EUR",
    "superpopulation_defaulted": false,
    "superpopulation_auto_defaulted_models": [],
    "reference_panel_ancestry": "EUR",
    "n_requested": 3,
    "n_scored": 2,
    "n_failed": 1,
    "n_skipped": 0,
    "n_filtered": 4,
    "n_omitted": 0,
    "filter_summary": "profile=curated dropped 4 scored row(s): 2 no performance evidence, 2 duplicate score family",
    "upstream_summary": "Scored 2/3 PGS ID(s) for type 2 diabetes mellitus on GRCh38; 1 failed.",
    "model_agreement": {
      "reliable_models": 1,
      "verdict": "insufficient_reliable_models",
      "percentile_min": 78.0,
      "percentile_max": 78.0,
      "percentile_spread": 0.0,
      "percentile_mean": 78.0
    },
    "scores": [
      {
        "pgs_id": "PGS000014",
        "pgs_catalog_url": "https://www.pgscatalog.org/score/PGS000014/",
        "status": "scored",
        "raw_score": 0.42,
        "genome_build": "GRCh38",
        "variants_number": 100,
        "weight_type": "beta",
        "variants_matched": 92,
        "variants_total": 100,
        "match_rate": 0.92,
        "weight_mass_coverage": 0.94,
        "percentile": 78.0,
        "percentile_method": "reference_panel",
        "percentile_reliable": true,
        "percentile_caveat": null,
        "requested_superpopulation": "EUR",
        "superpopulation_applied": "EUR",
        "superpopulation_auto_defaulted": false,
        "reference_panel_ancestry": "EUR",
        "quality_label": "High",
        "quality_summary": "High-quality model with reliable percentile.",
        "effect_size": "OR=1.55 [1.40-1.70]",
        "auroc_estimate": 0.74,
        "absolute_risk": {
          "status": "unavailable",
          "reason": "Offline demo does not call prevalence services."
        },
        "error": null
      },
      {
        "pgs_id": "PGS000013",
        "pgs_catalog_url": "https://www.pgscatalog.org/score/PGS000013/",
        "status": "scored",
        "raw_score": 0.11,
        "genome_build": "GRCh38",
        "variants_number": 80,
        "weight_type": "beta",
        "variants_matched": 58,
        "variants_total": 80,
        "match_rate": 0.725,
        "weight_mass_coverage": 0.61,
        "percentile": 54.0,
        "percentile_method": "theoretical",
        "percentile_reliable": false,
        "percentile_caveat": "Coverage is too low for a reliable percentile.",
        "requested_superpopulation": "EUR",
        "superpopulation_applied": "EUR",
        "superpopulation_auto_defaulted": false,
        "reference_panel_ancestry": "EUR",
        "quality_label": "Moderate",
        "quality_summary": "Moderate model; interpret cautiously.",
        "effect_size": "OR=1.20",
        "auroc_estimate": 0.65,
        "absolute_risk": {
          "status": "unavailable",
          "reason": "Offline demo does not call prevalence services."
        },
        "error": null
      },
      {
        "pgs_id": "PGS999999",
        "pgs_catalog_url": "https://www.pgscatalog.org/score/PGS999999/",
        "status": "failed",
        "raw_score": null,
        "genome_build": null,
        "variants_number": null,
        "weight_type": null,
        "variants_matched": null,
        "variants_total": null,
        "match_rate": null,
        "weight_mass_coverage": null,
        "percentile": null,
        "percentile_method": null,
        "percentile_reliable": null,
        "percentile_caveat": null,
        "requested_superpopulation": "EUR",
        "superpopulation_applied": "EUR",
        "superpopulation_auto_defaulted": false,
        "reference_panel_ancestry": null,
        "quality_label": null,
        "quality_summary": null,
        "effect_size": null,
        "auroc_estimate": null,
        "absolute_risk": {
          "status": "not_attempted",
          "reason": "Score failed, so absolute risk was not attempted."
        },
        "error": "Scoring file unavailable."
      }
    ]
  }
}\
""", 'tables/scores.csv': """\
pgs_id,status,raw_score,variants_matched,variants_total,match_rate,weight_mass_coverage,percentile,percentile_method,percentile_reliable,reference_panel_ancestry,quality_label,auroc_estimate,absolute_risk_status,absolute_risk,risk_ratio,error
PGS000014,scored,0.42,92,100,0.92,0.94,78.0,reference_panel,True,EUR,High,0.74,unavailable,,,
PGS000013,scored,0.11,58,80,0.725,0.61,54.0,theoretical,False,EUR,Moderate,0.65,unavailable,,,
PGS999999,failed,,,,,,,,,,,,not_attempted,,,Scoring file unavailable.
"""})
