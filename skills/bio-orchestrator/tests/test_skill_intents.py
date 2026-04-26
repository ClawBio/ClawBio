import json
from pathlib import Path

from clawbio.skill_intents import SCHEMA, plan_skill_intent


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
