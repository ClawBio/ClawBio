"""Pluggable subagent dispatcher.

The loop does not call LLMs directly. It asks a Dispatcher to execute a
role (proposer or executor) with a prompt and a tool whitelist, and
receives the response as text. The outer agent provides a concrete
dispatcher; tests inject a callable.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable


@dataclass
class DispatchRequest:
    role: str
    prompt: str
    allowed_tools: list[str]


class Dispatcher:
    """Base dispatcher. Override `dispatch` in subclasses or use CallableDispatcher."""

    def dispatch(self, request: DispatchRequest) -> str:
        raise RuntimeError(
            "No dispatcher configured. Provide a CallableDispatcher "
            "whose callback routes requests to an outer agent's Agent-tool "
            "invocation, or subclass Dispatcher."
        )


class CallableDispatcher(Dispatcher):
    """Adapter that wraps any callable taking a DispatchRequest."""

    def __init__(self, func: Callable[[DispatchRequest], str]):
        self._func = func

    def dispatch(self, request: DispatchRequest) -> str:
        return self._func(request)
