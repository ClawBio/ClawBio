"""Shared helpers for robust OpenAI-compatible tool-call loops."""

from __future__ import annotations

import asyncio
import json
from typing import Any, Awaitable, Callable


def tool_call_id(tool_call: Any) -> str:
    """Extract a stable tool_call_id from an OpenAI-compatible tool call."""

    return str(getattr(tool_call, "id", ""))


def tool_call_signature(tool_call: Any) -> str:
    """Return a stable signature for duplicate suppression within one turn."""

    function = getattr(tool_call, "function", None)
    func_name = str(getattr(function, "name", ""))
    arguments = getattr(function, "arguments", "{}") or "{}"
    try:
        normalised_args = json.dumps(json.loads(arguments), sort_keys=True)
    except json.JSONDecodeError:
        normalised_args = str(arguments)
    return f"{func_name}:{normalised_args}"


def _history_tool_call_id(tool_call: Any) -> str:
    if isinstance(tool_call, dict):
        return str(tool_call.get("id", ""))
    return tool_call_id(tool_call)


def tool_error_content(
    tool: str,
    error_type: str,
    message: str,
    *,
    tool_call_id_value: str | None = None,
    retryable: bool = False,
    details: dict[str, Any] | None = None,
) -> str:
    """Return a structured tool-result error payload as JSON text."""

    payload: dict[str, Any] = {
        "ok": False,
        "tool": tool,
        "error": {
            "type": error_type,
            "message": message,
            "retryable": retryable,
        },
    }
    if tool_call_id_value:
        payload["tool_call_id"] = tool_call_id_value
    if details:
        payload["details"] = details
    return json.dumps(payload, sort_keys=True)


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


def repair_tool_call_history(
    history: list[dict[str, Any]],
    *,
    repair_content: str | None = None,
    audit: Callable[..., None] | None = None,
    audit_context: dict[str, Any] | None = None,
    logger: Any | None = None,
) -> list[dict[str, Any]]:
    """Repair OpenAI chat history so every assistant tool call has a tool result.

    OpenAI requires an assistant message with ``tool_calls`` to be followed
    immediately by one ``role='tool'`` message for each emitted
    ``tool_call_id``. This function drops orphan/duplicate tool messages and
    inserts synthetic error tool results for missing ids.
    """

    audit_context = audit_context or {}
    repaired: list[dict[str, Any]] = []
    i = 0
    while i < len(history):
        msg = history[i]
        if msg.get("role") == "tool":
            if logger:
                logger.warning("Dropped orphaned tool message from history")
            if audit:
                audit("history_sanitised", **audit_context, detail="orphaned_tool_message_dropped")
            i += 1
            continue

        repaired.append(msg)
        tool_calls = msg.get("tool_calls") if msg.get("role") == "assistant" else None
        if not tool_calls:
            i += 1
            continue

        expected_ids = [call_id for call_id in (_history_tool_call_id(call) for call in tool_calls) if call_id]
        seen_ids: set[str] = set()
        i += 1
        while i < len(history) and history[i].get("role") == "tool":
            tool_msg = history[i]
            tool_call_id_value = str(tool_msg.get("tool_call_id", ""))
            if tool_call_id_value in expected_ids and tool_call_id_value not in seen_ids:
                repaired.append(tool_msg)
                seen_ids.add(tool_call_id_value)
            else:
                if logger:
                    logger.warning("Dropped unexpected or duplicate tool message from history")
                if audit:
                    audit(
                        "history_sanitised",
                        **audit_context,
                        detail="unexpected_or_duplicate_tool_message_dropped",
                        tool_call_id=tool_call_id_value,
                    )
            i += 1

        for expected_id in expected_ids:
            if expected_id in seen_ids:
                continue
            content = repair_content or tool_error_content(
                "unknown",
                "missing_tool_result_repaired",
                "Recovered a previous assistant tool call that did not have a tool result.",
                tool_call_id_value=expected_id,
                retryable=True,
            )
            repaired.append({
                "role": "tool",
                "tool_call_id": expected_id,
                "content": content,
            })
            if logger:
                logger.warning("Inserted synthetic tool result for unresolved tool_call_id %s", expected_id)
            if audit:
                audit(
                    "history_sanitised",
                    **audit_context,
                    detail="missing_tool_message_inserted",
                    tool_call_id=expected_id,
                )

    history[:] = repaired
    return history


async def execute_tool_calls_safely(
    tool_calls: list[Any],
    executors: dict[str, Callable[[dict], Awaitable[str]]],
    *,
    base_args: dict[str, Any] | None = None,
    raw_user_text: str = "",
    audit: Callable[..., None] | None = None,
    audit_context: dict[str, Any] | None = None,
    logger: Any | None = None,
    seen_signatures: set[str] | None = None,
) -> list[dict[str, str]]:
    """Execute tool calls and always return matching tool messages."""

    messages: list[dict[str, str]] = []
    base_args = base_args or {}
    audit_context = audit_context or {}
    expected_ids = [tool_call_id(tc) for tc in tool_calls]
    for tc in tool_calls:
        current_tool_call_id = tool_call_id(tc)
        signature = tool_call_signature(tc)
        function = getattr(tc, "function", None)
        func_name = str(getattr(function, "name", ""))
        arguments = getattr(function, "arguments", "{}")
        result = tool_error_content(
            func_name or "unknown",
            "unknown_tool",
            f"Unknown tool: {func_name}" if func_name else "Unknown tool call.",
            tool_call_id_value=current_tool_call_id,
        )
        try:
            if seen_signatures is not None:
                if signature in seen_signatures:
                    if logger:
                        logger.warning(f"Duplicate tool call suppressed: {func_name}")
                    if audit:
                        audit(
                            "tool_duplicate_suppressed",
                            **audit_context,
                            tool=func_name,
                        )
                    messages.append({
                        "role": "tool",
                        "tool_call_id": current_tool_call_id,
                        "content": "Duplicate tool call suppressed; the same tool request was already handled in this turn.",
                    })
                    continue
                seen_signatures.add(signature)
            executor = executors.get(func_name)
            if not executor:
                result = tool_error_content(
                    func_name or "unknown",
                    "unknown_tool",
                    f"Unknown tool: {func_name}" if func_name else "Unknown tool call.",
                    tool_call_id_value=current_tool_call_id,
                )
            else:
                try:
                    args = json.loads(arguments or "{}")
                except json.JSONDecodeError as parse_err:
                    messages.append({
                        "role": "tool",
                        "tool_call_id": current_tool_call_id,
                        "content": tool_error_content(
                            func_name,
                            "malformed_arguments",
                            f"Tool arguments were not valid JSON: {parse_err}",
                            tool_call_id_value=current_tool_call_id,
                        ),
                    })
                    continue
                if not isinstance(args, dict):
                    messages.append({
                        "role": "tool",
                        "tool_call_id": current_tool_call_id,
                        "content": tool_error_content(
                            func_name,
                            "malformed_arguments",
                            "Tool arguments must be a JSON object.",
                            tool_call_id_value=current_tool_call_id,
                        ),
                    })
                    continue
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
            result = tool_error_content(
                func_name,
                "cancelled",
                f"Tool execution cancelled before completion: {func_name}",
                tool_call_id_value=current_tool_call_id,
                retryable=True,
            )
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
            result = tool_error_content(
                func_name,
                "exception",
                f"{type(tool_err).__name__}: {tool_err}",
                tool_call_id_value=current_tool_call_id,
            )
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
                "content": tool_error_content(
                    "unknown",
                    "missing_result",
                    "Tool execution did not produce a result.",
                    tool_call_id_value=missing_id,
                    retryable=True,
                ),
            })
    return messages
