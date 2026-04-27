"""Shared helpers for robust OpenAI-compatible tool-call loops."""

from __future__ import annotations

import asyncio
import json
from typing import Any, Awaitable, Callable


def tool_call_id(tool_call: Any) -> str:
    """Extract a stable tool_call_id from an OpenAI-compatible tool call."""

    return str(getattr(tool_call, "id", ""))


def synthetic_tool_result_messages(
    tool_calls: list[Any],
    content: str,
) -> list[dict[str, str]]:
    """Build one synthetic tool message for every tool call."""

    return [
        {
            "role": "tool",
            "tool_call_id": tool_call_id(tc),
            "content": content,
        }
        for tc in tool_calls
    ]


async def execute_tool_calls_safely(
    tool_calls: list[Any],
    executors: dict[str, Callable[[dict], Awaitable[str]]],
    *,
    base_args: dict[str, Any] | None = None,
    raw_user_text: str = "",
    audit: Callable[..., None] | None = None,
    audit_context: dict[str, Any] | None = None,
    logger: Any | None = None,
) -> list[dict[str, str]]:
    """Execute tool calls and always return matching tool messages."""

    messages: list[dict[str, str]] = []
    base_args = base_args or {}
    audit_context = audit_context or {}
    expected_ids = [tool_call_id(tc) for tc in tool_calls]
    for tc in tool_calls:
        current_tool_call_id = tool_call_id(tc)
        function = getattr(tc, "function", None)
        func_name = str(getattr(function, "name", ""))
        arguments = getattr(function, "arguments", "{}")
        result = f"Unknown tool: {func_name}" if func_name else "Unknown tool call."
        try:
            executor = executors.get(func_name)
            if executor:
                try:
                    args = json.loads(arguments or "{}")
                except json.JSONDecodeError:
                    args = {}
                args.update(base_args)
                if raw_user_text:
                    args["_raw_user_text"] = raw_user_text
                if logger:
                    logger.info(f"Tool call: {func_name}({json.dumps(args)[:200]})")
                if audit:
                    audit(
                        "tool_call",
                        **audit_context,
                        tool=func_name,
                        args_preview=json.dumps(args, default=str)[:300],
                    )
                result = await executor(args)
        except asyncio.CancelledError:
            if logger:
                logger.warning(f"Tool {func_name} was cancelled")
            if audit:
                audit(
                    "tool_cancelled",
                    **audit_context,
                    tool=func_name,
                )
            result = f"Tool execution cancelled before completion: {func_name}"
        except Exception as tool_err:
            if logger:
                logger.error(f"Tool {func_name} raised: {tool_err}", exc_info=True)
            if audit:
                audit(
                    "tool_error",
                    **audit_context,
                    tool=func_name,
                    error=str(tool_err)[:300],
                )
            result = f"Error executing {func_name}: {type(tool_err).__name__}: {tool_err}"
        messages.append({
            "role": "tool",
            "tool_call_id": current_tool_call_id,
            "content": str(result),
        })
    returned_ids = {message["tool_call_id"] for message in messages}
    for missing_id in expected_ids:
        if missing_id not in returned_ids:
            messages.append({
                "role": "tool",
                "tool_call_id": missing_id,
                "content": "Tool execution did not produce a result.",
            })
    return messages
