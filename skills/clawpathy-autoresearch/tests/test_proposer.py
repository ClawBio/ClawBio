from skills.clawpathy_autoresearch.dispatcher import CallableDispatcher, DispatchRequest
from skills.clawpathy_autoresearch.proposer import propose_skill


def test_proposer_returns_replacement_skill():
    captured = {}

    def fake(req: DispatchRequest) -> str:
        captured["prompt"] = req.prompt
        captured["role"] = req.role
        return "```markdown\n# NEW SKILL\nStep 1: do thing\n```"

    d = CallableDispatcher(fake)
    out = propose_skill(
        dispatcher=d,
        task_prompt="make the thing",
        current_skill="# OLD\n",
        last_score=0.5,
        last_breakdown={"err_a": 0.5},
        recent_history=[{"iter": 1, "score": 0.6, "kept": False}],
    )
    assert out == "# NEW SKILL\nStep 1: do thing"
    assert captured["role"] == "proposer"
    assert "make the thing" in captured["prompt"]
    assert "# OLD" in captured["prompt"]
    assert "err_a" in captured["prompt"]


def test_proposer_no_history_first_iteration():
    def fake(req):
        return "# SKILL\ncontent"
    d = CallableDispatcher(fake)
    out = propose_skill(
        dispatcher=d, task_prompt="t", current_skill="",
        last_score=None, last_breakdown=None, recent_history=[],
    )
    assert out == "# SKILL\ncontent"
