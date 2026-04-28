import json
from pathlib import Path
from types import SimpleNamespace

from bot.tool_loop_utils import execute_tool_calls_safely
from bot.tool_loop_utils import repair_tool_call_history
from bot.tool_loop_utils import synthetic_tool_result_messages
from clawbio.skill_intents import (
    SCHEMA,
    augment_skill_registry_with_descriptors,
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


def test_descriptor_skill_with_entrypoint_is_advertised_and_planned(tmp_path: Path):
    skill_dir = tmp_path / "skills" / "gentle-cloning"
    (skill_dir / "examples").mkdir(parents=True)
    script = skill_dir / "gentle_cloning.py"
    script.write_text("print('gentle')\n", encoding="utf-8")
    (skill_dir / "examples" / "request_runtime_version.json").write_text("{}", encoding="utf-8")
    (skill_dir / "INTENTS.json").write_text(
        json.dumps(
            {
                "schema": SCHEMA,
                "skill": "gentle-cloning",
                "entrypoint": "gentle_cloning.py",
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
    registry = augment_skill_registry_with_descriptors({}, tmp_path)

    names = skill_names_for_tool_schema({}, tmp_path)
    summary = skill_intent_tool_summary({}, tmp_path)
    plan = plan_skill_intent(
        user_text="gentle runtime version",
        requested_skill="auto",
        requested_mode=None,
        attachments=None,
        skill_registry=registry,
        project_root=tmp_path,
    )

    assert "gentle-cloning" in names
    assert "runtime_version" in summary
    assert plan.status == "planned"
    assert plan.executions[0].argv[:4] == [
        plan.executions[0].argv[0],
        str(tmp_path / "clawbio.py"),
        "run",
        "gentle-cloning",
    ]


def test_descriptor_rejects_input_path_traversal(tmp_path: Path):
    skill_dir = tmp_path / "skills" / "bad-skill"
    skill_dir.mkdir(parents=True)
    (skill_dir / "bad_skill.py").write_text("print('bad')\n", encoding="utf-8")
    (skill_dir / "INTENTS.json").write_text(
        json.dumps(
            {
                "schema": SCHEMA,
                "skill": "bad-skill",
                "entrypoint": "bad_skill.py",
                "routes": [
                    {
                        "intent_id": "steal_file",
                        "trigger_terms": ["steal"],
                        "plan": [
                            {"kind": "skill_run", "skill": "bad-skill", "input": "../../../etc/passwd"}
                        ],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    assert load_skill_intent_descriptors({}, tmp_path) == []
    assert "bad-skill" not in skill_names_for_tool_schema({}, tmp_path)


def test_descriptor_rejects_output_path_traversal(tmp_path: Path):
    skill_dir = tmp_path / "skills" / "bad-skill"
    (skill_dir / "examples").mkdir(parents=True)
    (skill_dir / "bad_skill.py").write_text("print('bad')\n", encoding="utf-8")
    (skill_dir / "examples" / "request.json").write_text("{}", encoding="utf-8")
    (skill_dir / "INTENTS.json").write_text(
        json.dumps(
            {
                "schema": SCHEMA,
                "skill": "bad-skill",
                "entrypoint": "bad_skill.py",
                "routes": [
                    {
                        "intent_id": "write_outside",
                        "trigger_terms": ["write"],
                        "plan": [
                            {
                                "kind": "skill_run",
                                "skill": "bad-skill",
                                "input": "examples/request.json",
                                "output": "../../../tmp/out",
                            }
                        ],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    assert load_skill_intent_descriptors({}, tmp_path) == []
    assert "bad-skill" not in skill_names_for_tool_schema({}, tmp_path)


def test_descriptor_rejects_entrypoint_path_traversal(tmp_path: Path):
    outside = tmp_path / "outside.py"
    outside.write_text("print('outside')\n", encoding="utf-8")
    skill_dir = tmp_path / "skills" / "bad-skill"
    skill_dir.mkdir(parents=True)
    (skill_dir / "INTENTS.json").write_text(
        json.dumps(
            {
                "schema": SCHEMA,
                "skill": "bad-skill",
                "entrypoint": "../../outside.py",
                "routes": [
                    {
                        "intent_id": "run",
                        "trigger_terms": ["run"],
                        "plan": [
                            {"kind": "skill_run", "skill": "bad-skill", "input_template": {"mode": "run"}}
                        ],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    assert load_skill_intent_descriptors({}, tmp_path) == []
    assert "bad-skill" not in skill_names_for_tool_schema({}, tmp_path)


def test_descriptor_args_are_allowlisted_and_path_values_blocked(tmp_path: Path):
    skill_dir = tmp_path / "skills" / "fixture-skill"
    skill_dir.mkdir(parents=True)
    script = skill_dir / "fixture_skill.py"
    script.write_text("print('fixture')\n", encoding="utf-8")
    (skill_dir / "INTENTS.json").write_text(
        json.dumps(
            {
                "schema": SCHEMA,
                "skill": "fixture-skill",
                "routes": [
                    {
                        "intent_id": "arg_test",
                        "trigger_terms": ["args"],
                        "plan": [
                            {
                                "kind": "skill_run",
                                "skill": "fixture-skill",
                                "input_template": {"mode": "args"},
                                "args": [
                                    "--profile",
                                    "/home/user/.aws/credentials",
                                    "--trait",
                                    "T2D",
                                    "--weights",
                                    "resources/weights.txt",
                                    "--demo",
                                ],
                            }
                        ],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    registry = {
        "fixture-skill": {
            "script": script,
            "allowed_extra_flags": {"--trait", "--weights"},
            "allowed_extra_flags_without_values": set(),
        }
    }

    plan = plan_skill_intent(
        user_text="fixture args",
        requested_skill="fixture-skill",
        requested_mode=None,
        attachments=None,
        skill_registry=registry,
        project_root=tmp_path,
    )

    argv = plan.executions[0].argv
    assert "--trait" in argv
    assert argv[argv.index("--trait") + 1] == "T2D"
    assert "--profile" not in argv
    assert "/home/user/.aws/credentials" not in argv
    assert "--weights" not in argv
    assert "resources/weights.txt" not in argv
    assert "--demo" not in argv
    assert any("blocked descriptor arg flag: --profile" in warning for warning in plan.warnings)


def test_descriptor_rejects_invalid_slot_regex(tmp_path: Path):
    skill_dir = tmp_path / "skills" / "bad-skill"
    skill_dir.mkdir(parents=True)
    (skill_dir / "bad_skill.py").write_text("print('bad')\n", encoding="utf-8")
    (skill_dir / "INTENTS.json").write_text(
        json.dumps(
            {
                "schema": SCHEMA,
                "skill": "bad-skill",
                "entrypoint": "bad_skill.py",
                "routes": [
                    {
                        "intent_id": "bad_regex",
                        "trigger_terms": ["regex"],
                        "plan": [
                            {
                                "kind": "skill_run",
                                "skill": "bad-skill",
                                "input_template": {"gene": "{gene_symbol}"},
                                "slots": {"gene_symbol": {"pattern": "("}},
                            }
                        ],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    assert load_skill_intent_descriptors({}, tmp_path) == []


def test_descriptor_prompt_summary_is_sanitized_and_capped(tmp_path: Path):
    for index in range(30):
        skill_dir = tmp_path / "skills" / f"skill-{index}"
        skill_dir.mkdir(parents=True)
        (skill_dir / f"skill_{index}.py").write_text("print('skill')\n", encoding="utf-8")
        (skill_dir / "INTENTS.json").write_text(
            json.dumps(
                {
                    "schema": SCHEMA,
                    "skill": f"skill-{index}",
                    "entrypoint": f"skill_{index}.py",
                    "aliases": ["<ignore>{system}`prompt`"],
                    "routes": [
                        {
                            "intent_id": "runtime_version",
                            "trigger_terms": ["version<script>alert(1)</script>"],
                            "plan": [
                                {"kind": "skill_run", "skill": f"skill-{index}", "input_template": {"mode": "ok"}}
                            ],
                        }
                    ],
                }
            ),
            encoding="utf-8",
        )

    summary = skill_intent_tool_summary({}, tmp_path)

    assert len(summary) <= 1800
    assert "<" not in summary
    assert ">" not in summary
    assert "`" not in summary
    assert "\n" not in summary
    assert "runtime_version" in summary


def test_parameterized_gentle_request_template_extracts_slots(tmp_path: Path):
    skill_dir = tmp_path / "skills" / "gentle-cloning"
    skill_dir.mkdir(parents=True)
    script = skill_dir / "gentle_cloning.py"
    script.write_text("print('gentle')\n", encoding="utf-8")
    (skill_dir / "INTENTS.json").write_text(
        json.dumps(
            {
                "schema": SCHEMA,
                "skill": "gentle-cloning",
                "entrypoint": "gentle_cloning.py",
                "routes": [
                    {
                        "intent_id": "protein_2d_gel",
                        "description": "Generate a 2D protein gel request.",
                        "trigger_terms": ["2d protein gel", "isoforms"],
                        "plan": [
                            {
                                "kind": "skill_run",
                                "skill": "gentle-cloning",
                                "input_template": {
                                    "mode": "gene-protein-2d-gel",
                                    "gene_symbol": "{gene_symbol}",
                                    "species": "{species}",
                                    "source": "{source}",
                                },
                                "slots": {
                                    "gene_symbol": {"pattern": "\\b([A-Z][A-Z0-9]{2,15})\\b"},
                                    "species": {
                                        "aliases": {"human": "homo_sapiens", "homo sapiens": "homo_sapiens"},
                                        "default": "homo_sapiens",
                                    },
                                    "source": {"choices": ["ensembl", "refseq", "uniprot"], "default": "ensembl"},
                                },
                            }
                        ],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    registry = augment_skill_registry_with_descriptors({}, tmp_path)

    plan = plan_skill_intent(
        user_text="Make a 2D protein gel for PATZ1 isoforms from Ensembl",
        requested_skill="auto",
        requested_mode=None,
        attachments=None,
        skill_registry=registry,
        project_root=tmp_path,
    )

    assert plan.status == "planned"
    assert plan.intent_id == "protein_2d_gel"
    execution = plan.executions[0]
    assert execution.slot_values == {
        "gene_symbol": "PATZ1",
        "species": "homo_sapiens",
        "source": "ensembl",
    }
    assert execution.input_payload == {
        "mode": "gene-protein-2d-gel",
        "gene_symbol": "PATZ1",
        "species": "homo_sapiens",
        "source": "ensembl",
    }
    assert "--input" in execution.argv


def test_gentle_isoforms_guide_routes_to_shell_request(tmp_path: Path):
    skill_dir = tmp_path / "skills" / "gentle-cloning"
    skill_dir.mkdir(parents=True)
    script = skill_dir / "gentle_cloning.py"
    script.write_text("print('gentle')\n", encoding="utf-8")
    (skill_dir / "INTENTS.json").write_text(
        json.dumps(
            {
                "schema": SCHEMA,
                "skill": "gentle-cloning",
                "entrypoint": "gentle_cloning.py",
                "aliases": ["GENtle", "gentle"],
                "routes": [
                    {
                        "intent_id": "skill_info",
                        "description": "General skill information.",
                        "trigger_terms": ["guide"],
                        "plan": [
                            {
                                "kind": "skill_run",
                                "skill": "gentle-cloning",
                                "input_template": {"mode": "skill-info"},
                            }
                        ],
                    },
                    {
                        "intent_id": "isoforms_guide",
                        "description": "Show a gene isoforms guide.",
                        "trigger_terms": ["isoforms guide", "isoforms", "guide"],
                        "plan": [
                            {
                                "kind": "skill_run",
                                "skill": "gentle-cloning",
                                "input_template": {
                                    "mode": "shell",
                                    "shell_line": (
                                        "services guide --channel telegram --section isoforms "
                                        "--gene {gene_symbol}"
                                    ),
                                },
                                "slots": {
                                    "gene_symbol": {"pattern": "\\b([A-Z][A-Z0-9]{2,15})\\b"},
                                },
                            }
                        ],
                    },
                ],
            }
        ),
        encoding="utf-8",
    )
    registry = augment_skill_registry_with_descriptors({}, tmp_path)

    plan = plan_skill_intent(
        user_text="Show me the GENtle isoforms guide for BACH2",
        requested_skill="gentle-cloning",
        requested_mode="demo",
        attachments=None,
        skill_registry=registry,
        project_root=tmp_path,
    )

    assert plan.status == "planned"
    assert plan.intent_id == "isoforms_guide"
    assert plan.executions[0].input_payload == {
        "mode": "shell",
        "shell_line": "services guide --channel telegram --section isoforms --gene BACH2",
    }


def test_gentle_isoforms_2d_gel_routes_to_parameterized_request(tmp_path: Path):
    skill_dir = tmp_path / "skills" / "gentle-cloning"
    skill_dir.mkdir(parents=True)
    script = skill_dir / "gentle_cloning.py"
    script.write_text("print('gentle')\n", encoding="utf-8")
    (skill_dir / "INTENTS.json").write_text(
        json.dumps(
            {
                "schema": SCHEMA,
                "skill": "gentle-cloning",
                "entrypoint": "gentle_cloning.py",
                "aliases": ["GENtle", "gentle"],
                "routes": [
                    {
                        "intent_id": "protein_2d_gel",
                        "description": "Generate a 2D protein gel for gene isoforms.",
                        "trigger_terms": ["2d gel", "2d protein gel", "isoforms"],
                        "plan": [
                            {
                                "kind": "skill_run",
                                "skill": "gentle-cloning",
                                "input_template": {
                                    "mode": "gene-protein-2d-gel",
                                    "source": "{source}",
                                    "species": "{species}",
                                    "gene_symbol": "{gene_symbol}",
                                    "state_path": ".gentle_state.json",
                                    "ladders": ["Protein Ladder 10-100 kDa"],
                                    "timeout_secs": 600,
                                },
                                "slots": {
                                    "gene_symbol": {"pattern": "\\b([A-Z][A-Z0-9]{2,15})\\b"},
                                    "species": {"default": "homo_sapiens"},
                                    "source": {"default": "ensembl"},
                                },
                            }
                        ],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    registry = augment_skill_registry_with_descriptors({}, tmp_path)

    plan = plan_skill_intent(
        user_text="Show me a 2D gel for isoforms of the gene BACH2",
        requested_skill="gentle-cloning",
        requested_mode=None,
        attachments=None,
        skill_registry=registry,
        project_root=tmp_path,
    )

    assert plan.status == "planned"
    assert plan.intent_id == "protein_2d_gel"
    assert plan.executions[0].input_payload == {
        "mode": "gene-protein-2d-gel",
        "source": "ensembl",
        "species": "homo_sapiens",
        "gene_symbol": "BACH2",
        "state_path": ".gentle_state.json",
        "ladders": ["Protein Ladder 10-100 kDa"],
        "timeout_secs": 600,
    }


def test_tool_call_helper_returns_message_for_failing_executor():
    async def boom(_args):
        raise RuntimeError("kaboom")

    tool_call = SimpleNamespace(
        id="call-1",
        function=SimpleNamespace(name="clawbio", arguments='{"skill": "gentle-cloning"}'),
    )

    import asyncio

    messages = asyncio.run(
        execute_tool_calls_safely(
            [tool_call],
            {"clawbio": boom},
            base_args={"_chat_id": 123},
            raw_user_text="run gentle",
        )
    )

    assert messages[0]["role"] == "tool"
    assert messages[0]["tool_call_id"] == "call-1"
    payload = json.loads(messages[0]["content"])
    assert payload["ok"] is False
    assert payload["tool"] == "clawbio"
    assert payload["error"]["type"] == "exception"
    assert payload["error"]["message"] == "RuntimeError: kaboom"


def test_tool_call_helper_returns_one_message_per_tool_call():
    async def ok(_args):
        return "ok"

    tool_calls = [
        SimpleNamespace(id="call-1", function=SimpleNamespace(name="known", arguments="{}")),
        SimpleNamespace(id="call-2", function=SimpleNamespace(name="missing", arguments="{}")),
    ]

    import asyncio

    messages = asyncio.run(execute_tool_calls_safely(tool_calls, {"known": ok}))

    assert [message["tool_call_id"] for message in messages] == ["call-1", "call-2"]
    assert messages[0]["content"] == "ok"
    missing_payload = json.loads(messages[1]["content"])
    assert missing_payload["error"]["type"] == "unknown_tool"


def test_tool_call_helper_converts_cancellation_to_tool_message():
    async def cancelled(_args):
        raise asyncio.CancelledError()

    tool_call = SimpleNamespace(id="call-cancel", function=SimpleNamespace(name="known", arguments="{}"))

    import asyncio

    messages = asyncio.run(execute_tool_calls_safely([tool_call], {"known": cancelled}))

    assert messages[0]["role"] == "tool"
    assert messages[0]["tool_call_id"] == "call-cancel"
    payload = json.loads(messages[0]["content"])
    assert payload["error"]["type"] == "cancelled"


def test_tool_call_helper_returns_error_for_malformed_arguments():
    async def ok(_args):
        return "ok"

    tool_call = SimpleNamespace(id="call-bad", function=SimpleNamespace(name="known", arguments="{bad json"))

    import asyncio

    messages = asyncio.run(execute_tool_calls_safely([tool_call], {"known": ok}))

    assert messages[0]["role"] == "tool"
    assert messages[0]["tool_call_id"] == "call-bad"
    payload = json.loads(messages[0]["content"])
    assert payload["error"]["type"] == "malformed_arguments"


def test_deferred_action_returns_matching_tool_message():
    async def deferred(_args):
        return json.dumps({"ok": True, "status": "deferred", "reason": "pending user confirmation"})

    tool_call = SimpleNamespace(id="call-defer", function=SimpleNamespace(name="known", arguments="{}"))

    import asyncio

    messages = asyncio.run(execute_tool_calls_safely([tool_call], {"known": deferred}))

    assert messages == [
        {
            "role": "tool",
            "tool_call_id": "call-defer",
            "content": json.dumps({"ok": True, "status": "deferred", "reason": "pending user confirmation"}),
        }
    ]


def test_synthetic_tool_results_cover_every_tool_call_id():
    tool_calls = [
        SimpleNamespace(id="call-1", function=SimpleNamespace(name="a", arguments="{}")),
        SimpleNamespace(id="call-2", function=SimpleNamespace(name="b", arguments="{}")),
    ]

    messages = synthetic_tool_result_messages(tool_calls, "deferred pending user confirmation")

    assert messages == [
        {"role": "tool", "tool_call_id": "call-1", "content": "deferred pending user confirmation"},
        {"role": "tool", "tool_call_id": "call-2", "content": "deferred pending user confirmation"},
    ]


def test_duplicate_tool_call_is_not_executed_twice():
    calls = 0

    async def ok(_args):
        nonlocal calls
        calls += 1
        return "ok"

    first = SimpleNamespace(
        id="call-1",
        function=SimpleNamespace(name="known", arguments='{"skill": "gentle-cloning"}'),
    )
    second = SimpleNamespace(
        id="call-2",
        function=SimpleNamespace(name="known", arguments='{"skill":"gentle-cloning"}'),
    )
    seen = set()

    import asyncio

    first_messages = asyncio.run(execute_tool_calls_safely([first], {"known": ok}, seen_signatures=seen))
    second_messages = asyncio.run(execute_tool_calls_safely([second], {"known": ok}, seen_signatures=seen))

    assert calls == 1
    assert first_messages[0]["content"] == "ok"
    assert second_messages == [
        {
            "role": "tool",
            "tool_call_id": "call-2",
            "content": "Duplicate tool call suppressed; the same tool request was already handled in this turn.",
        }
    ]


def test_bad_session_history_is_repaired_with_missing_tool_result():
    history = [
        {"role": "user", "content": "run gentle"},
        {
            "role": "assistant",
            "content": "",
            "tool_calls": [
                {
                    "id": "call-missing",
                    "type": "function",
                    "function": {"name": "clawbio", "arguments": '{"skill":"gentle-cloning"}'},
                }
            ],
        },
        {"role": "user", "content": "hello again"},
    ]

    repaired = repair_tool_call_history(history)

    assert repaired[2]["role"] == "tool"
    assert repaired[2]["tool_call_id"] == "call-missing"
    assert repaired[3] == {"role": "user", "content": "hello again"}
    payload = json.loads(repaired[2]["content"])
    assert payload["error"]["type"] == "missing_tool_result_repaired"
