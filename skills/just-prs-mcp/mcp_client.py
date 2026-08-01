"""Typed local-stdio client for the pinned just-prs MCP server."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any


class McpCallError(RuntimeError):
    """A just-prs MCP tool failed or returned an unusable payload."""


def to_jsonable(value: Any) -> Any:
    """Convert FastMCP/Pydantic results into JSON-compatible Python values."""
    if hasattr(value, "model_dump"):
        return to_jsonable(value.model_dump(mode="json"))
    if hasattr(value, "root"):
        return to_jsonable(value.root)
    if isinstance(value, Mapping):
        return {str(key): to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(item) for item in value]
    return value


def tool_result_data(result: Any) -> Any:
    """Prefer FastMCP's transport-stable structured content over dynamic Root objects."""
    structured = getattr(result, "structured_content", None)
    if structured is not None:
        data = to_jsonable(structured)
        if isinstance(data, dict) and set(data) == {"result"}:
            return data["result"]
        return data
    return to_jsonable(getattr(result, "data", None))


class JustPrsMcpClient:
    """Async context manager for a local ``uvx`` just-prs-mcp process."""

    def __init__(
        self,
        command: list[str],
        environment: dict[str, str],
        *,
        timeout_seconds: int = 3600,
        log_file: Path | None = None,
    ) -> None:
        if len(command) < 2:
            raise ValueError("The MCP server command must include an executable and arguments.")
        self._command = command
        self._environment = environment
        self._timeout_seconds = timeout_seconds
        self._log_file = log_file
        self._client: Any = None

    async def __aenter__(self) -> JustPrsMcpClient:
        try:
            from fastmcp import Client
            from fastmcp.client.transports import StdioTransport
        except ImportError as exc:
            raise McpCallError(
                "The just-prs bridge requires the optional dependencies. "
                "Run `uv sync --extra just-prs` in the ClawBio checkout."
            ) from exc

        transport = StdioTransport(
            command=self._command[0],
            args=self._command[1:],
            env=self._environment,
            log_file=self._log_file,
        )
        self._client = Client(transport, timeout=self._timeout_seconds)
        try:
            await self._client.__aenter__()
        except Exception as exc:
            self._client = None
            raise McpCallError(f"Could not start local just-prs-mcp: {exc}") from exc
        return self

    async def __aexit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: Any,
    ) -> None:
        if self._client is not None:
            await self._client.__aexit__(exc_type, exc, traceback)
            self._client = None

    async def call_tool(self, name: str, arguments: dict[str, Any]) -> Any:
        """Call one tool and return a JSON-compatible structured payload."""
        if self._client is None:
            raise McpCallError("The MCP client must be entered before calling tools.")
        try:
            result = await self._client.call_tool(name, arguments)
        except Exception as exc:
            raise McpCallError(f"{name} failed: {exc}") from exc
        data = tool_result_data(result)
        if data is None:
            raise McpCallError(f"{name} returned no structured data.")
        return data
