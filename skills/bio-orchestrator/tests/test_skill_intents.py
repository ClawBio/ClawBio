import json
from pathlib import Path

from clawbio.skill_intents import (
    SCHEMA,
    load_skill_intent_descriptors,
    plan_skill_intent,
    skill_intent_tool_summary,
    skill_names_for_tool_schema,
)


def _fixture_registry(tmp_path: Path) -> dict:
    skill_dir = tmp_path / "skills" / "fixture-skill"
    examples_dir = skill_dir / "examples"
    examples_dir.mkdir(parents=True)
    script = skill_dir / "fixture_skill.py"
    script.write_text("print('fixture')\n", encoding="utf-8")
    (examples_dir / "status.json").write_text("{}", encoding="utf-8")
    (examples_dir / "prepare.json").write_text("{}", encoding="utf-8")
    (examples_dir / "finish.json").write_text("{}", encoding="utf-8")
    (skill_dir / "INTENTS.json").write_text(
        json.dumps(
            {
                "schema": SCHEMA,
                "skill": "fixture-skill",
                "aliases": ["fixture"],
                "routes": [
                    {
                        "intent_id": "runtime_version",
                        "description": "Check runtime status/version.",
                        "trigger_terms": ["version", "runtime version", "status"],
                        "demo_policy": "never_unless_explicit",
                        "plan": [
                            {
                                "kind": "skill_run",
                                "skill": "fixture-skill",
                                "input": "examples/status.json",
                            }
                        ],
                    },
                    {
                        "intent_id": "demo_report",
                        "description": "Run fixture demo.",
                        "trigger_terms": ["demo", "example"],
                        "demo_policy": "only_when_explicit",
                        "plan": [
                            {"kind": "skill_run", "skill": "fixture-skill", "demo": True}
                        ],
                    },
                    {
                        "intent_id": "two_step",
                        "description": "Run two related fixture steps.",
                        "trigger_terms": ["multi step", "pipeline"],
                        "plan": [
                            {
                                "id": "prepare",
                                "kind": "skill_run",
                                "skill": "fixture-skill",
                                "input": "examples/prepare.json",
                            },
                            {
                                "id": "finish",
                                "kind": "skill_run",
                                "skill": "fixture-skill",
                                "input": "examples/finish.json",
                                "confirmation": {
                                    "required": True,
                                    "reason": "Final step mutates cached fixture state.",
                                },
                            },
                        ],
                    },
                ],
            }
        ),
        encoding="utf-8",
    )
    other_dir = tmp_path / "skills" / "other-skill"
    other_dir.mkdir(parents=True)
    other_script = other_dir / "other_skill.py"
    other_script.write_text("print('other')\n", encoding="utf-8")
    return {
        "fixture-skill": {"script": script, "demo_args": ["--demo"], "allowed_extra_flags": set()},
        "other-skill": {"script": other_script, "demo_args": ["--demo"], "allowed_extra_flags": set()},
    }


def test_descriptor_routes_version_request(tmp_path: Path):
    registry = _fixture_registry(tmp_path)

    plan = plan_skill_intent(
        user_text="What runtime version is installed for fixture?",
        requested_skill="fixture-skill",
        requested_mode=None,
        attachments=None,
        skill_registry=registry,
    )

    assert plan.status == "planned"
    assert plan.intent_id == "runtime_version"
    assert plan.confidence == "high"
    assert len(plan.executions) == 1
    assert "--input" in plan.executions[0].argv
    assert plan.executions[0].argv[-1].endswith("examples/status.json")


def test_demo_mode_requires_explicit_demo_text(tmp_path: Path):
    registry = _fixture_registry(tmp_path)

    weak_demo = plan_skill_intent(
        user_text="What is the fixture status?",
        requested_skill="fixture-skill",
        requested_mode="demo",
        attachments=None,
        skill_registry=registry,
    )

    assert weak_demo.intent_id == "runtime_version"
    assert "--demo" not in weak_demo.executions[0].argv

    explicit_demo = plan_skill_intent(
        user_text="Please run the fixture demo.",
        requested_skill="fixture-skill",
        requested_mode="demo",
        attachments=None,
        skill_registry=registry,
    )

    assert explicit_demo.status == "planned"
    assert explicit_demo.intent_id == "demo_report"
    assert explicit_demo.executions[0].argv[-1] == "--demo"


def test_multistep_route_can_require_confirmation(tmp_path: Path):
    registry = _fixture_registry(tmp_path)

    plan = plan_skill_intent(
        user_text="Run the fixture multi step pipeline.",
        requested_skill="fixture-skill",
        requested_mode=None,
        attachments=None,
        skill_registry=registry,
    )

    assert plan.status == "needs_confirmation"
    assert plan.intent_id == "two_step"
    assert len(plan.executions) == 2
    assert plan.executions[1].requires_confirmation is True
    assert plan.executions[1].route_step_id == "finish"


def test_raw_text_can_override_weak_requested_skill(tmp_path: Path):
    registry = _fixture_registry(tmp_path)

    plan = plan_skill_intent(
        user_text="For fixture, check runtime version.",
        requested_skill="other-skill",
        requested_mode=None,
        attachments=None,
        skill_registry=registry,
    )

    assert plan.skill == "fixture-skill"
    assert plan.intent_id == "runtime_version"
    assert plan.requested_skill == "other-skill"


def test_execution_root_can_differ_from_symlinked_descriptor_root(tmp_path: Path):
    clawbio_root = tmp_path / "ClawBio"
    external_root = tmp_path / "gentle_rs" / "integrations" / "clawbio"
    external_skill = external_root / "skills" / "gentle-cloning"
    (external_skill / "examples").mkdir(parents=True)
    (external_skill / "examples" / "request_runtime_version.json").write_text("{}", encoding="utf-8")
    script = external_skill / "gentle_cloning.py"
    script.write_text("print('gentle')\n", encoding="utf-8")
    (external_skill / "INTENTS.json").write_text(
        json.dumps(
            {
                "schema": SCHEMA,
                "skill": "gentle-cloning",
                "routes": [
                    {
                        "intent_id": "runtime_version",
                        "trigger_terms": ["version"],
                        "plan": [
                            {
                                "kind": "skill_run",
                                "skill": "gentle-cloning",
                                "input": "examples/request_runtime_version.json",
                            }
                        ],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    registry_skill = clawbio_root / "skills" / "gentle-cloning"
    registry_skill.parent.mkdir(parents=True)
    registry_skill.symlink_to(external_skill, target_is_directory=True)
    registry = {"gentle-cloning": {"script": registry_skill / "gentle_cloning.py"}}

    plan = plan_skill_intent(
        user_text="check gentle version",
        requested_skill="auto",
        requested_mode=None,
        attachments=None,
        skill_registry=registry,
        project_root=clawbio_root,
    )

    assert plan.intent_id == "runtime_version"
    assert plan.executions[0].argv[1] == str(clawbio_root / "clawbio.py")
    assert plan.executions[0].input_path.endswith(
        "gentle_rs/integrations/clawbio/skills/gentle-cloning/examples/request_runtime_version.json"
    )


def test_confirmed_demo_tool_call_is_allowed(tmp_path: Path):
    registry = _fixture_registry(tmp_path)

    plan = plan_skill_intent(
        user_text="yes, go ahead",
        requested_skill="fixture-skill",
        requested_mode="demo",
        attachments=None,
        skill_registry=registry,
    )

    assert plan.status == "planned"
    assert plan.executions[0].argv[-1] == "--demo"


def test_drugphoto_keeps_demo_genotype_exception(tmp_path: Path):
    skill_dir = tmp_path / "skills" / "pharmgx-reporter"
    skill_dir.mkdir(parents=True)
    script = skill_dir / "pharmgx_reporter.py"
    script.write_text("print('drugphoto')\n", encoding="utf-8")
    registry = {"drugphoto": {"script": script, "demo_args": ["--demo"], "summary_default": True}}

    plan = plan_skill_intent(
        user_text="Medication photo shows clopidogrel 75mg",
        requested_skill="drugphoto",
        requested_mode="demo",
        attachments=[{"drug_name": "clopidogrel"}, {"visible_dose": "75mg"}],
        skill_registry=registry,
    )

    assert plan.status == "planned"
    assert "--demo" in plan.executions[0].argv
    assert "--drug" in plan.executions[0].argv


def test_unregistered_skill_directory_descriptor_is_discovered_but_not_executable(tmp_path: Path):
    skill_dir = tmp_path / "skills" / "gentle-cloning"
    (skill_dir / "examples").mkdir(parents=True)
    (skill_dir / "examples" / "request_runtime_version.json").write_text("{}", encoding="utf-8")
    (skill_dir / "INTENTS.json").write_text(
        json.dumps(
            {
                "schema": SCHEMA,
                "skill": "gentle-cloning",
                "aliases": ["gentle"],
                "routes": [
                    {
                        "intent_id": "runtime_version",
                        "description": "Check runtime version.",
                        "trigger_terms": ["version", "runtime version"],
                        "plan": [
                            {
                                "kind": "skill_run",
                                "skill": "gentle-cloning",
                                "input": "examples/request_runtime_version.json",
                            }
                        ],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    registry = {}

    descriptors = load_skill_intent_descriptors(registry, tmp_path)
    names = skill_names_for_tool_schema(registry, tmp_path)
    summary = skill_intent_tool_summary(registry, tmp_path)
    plan = plan_skill_intent(
        user_text="gentle runtime version",
        requested_skill="auto",
        requested_mode=None,
        attachments=None,
        skill_registry=registry,
        project_root=tmp_path,
    )

    assert [descriptor["skill"] for descriptor in descriptors] == ["gentle-cloning"]
    assert "gentle-cloning" not in names
    assert "runtime_version" not in summary
    assert plan.status == "needs_registration"
    assert plan.skill == "gentle-cloning"
    assert plan.executions == []
