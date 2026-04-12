import pytest
from skills.clawpathy_autoresearch.dispatcher import CallableDispatcher, DispatchRequest
from skills.clawpathy_autoresearch.sandbox import Sandbox
from skills.clawpathy_autoresearch.executor import execute_skill


def _sandbox():
    return Sandbox(
        allowed_tools=["Read", "Bash"],
        network_allowed=False,
        network_hosts=[],
        data_dir="./data",
        python_packages=["pandas"],
        timeout_seconds=300,
    )


def test_executor_passes_whitelist_and_returns_json():
    captured = {}

    def fake(req: DispatchRequest) -> str:
        captured["req"] = req
        return '```json\n{"result": 42}\n```'

    d = CallableDispatcher(fake)
    out = execute_skill(
        dispatcher=d,
        skill_content="# SKILL\nDo X",
        sandbox=_sandbox(),
        data_dir_abs="/tmp/data",
    )
    assert out == {"result": 42}
    assert captured["req"].role == "executor"
    assert captured["req"].allowed_tools == ["Read", "Bash"]
    assert "# SKILL" in captured["req"].prompt
    assert "allowed_tools" in captured["req"].prompt


def test_executor_raises_on_bad_json():
    d = CallableDispatcher(lambda req: "not json at all")
    with pytest.raises(ValueError, match="JSON"):
        execute_skill(
            dispatcher=d,
            skill_content="# SKILL",
            sandbox=_sandbox(),
            data_dir_abs="/tmp/data",
        )


def test_executor_accepts_raw_json_no_fence():
    d = CallableDispatcher(lambda req: '{"a": 1}')
    out = execute_skill(
        dispatcher=d, skill_content="# SKILL", sandbox=_sandbox(), data_dir_abs="/tmp/data",
    )
    assert out == {"a": 1}
