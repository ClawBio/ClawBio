import pytest
from skills.clawpathy_autoresearch.dispatcher import (
    Dispatcher,
    CallableDispatcher,
    DispatchRequest,
)


def test_callable_dispatcher_passes_through():
    calls = []

    def fake(req: DispatchRequest) -> str:
        calls.append(req)
        return "ok-response"

    d = CallableDispatcher(fake)
    out = d.dispatch(DispatchRequest(role="proposer", prompt="hi", allowed_tools=["Read"]))
    assert out == "ok-response"
    assert calls[0].role == "proposer"
    assert calls[0].prompt == "hi"
    assert calls[0].allowed_tools == ["Read"]


def test_base_dispatcher_raises_with_guidance():
    d = Dispatcher()
    with pytest.raises(RuntimeError, match="outer agent"):
        d.dispatch(DispatchRequest(role="proposer", prompt="x", allowed_tools=[]))
