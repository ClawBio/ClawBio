"""Explicit OpenAI and Ollama providers for synchronous text generation.

Task prompts and fallback behaviour belong to the calling skill. These clients
only handle configuration, requests, and conversion to a common result.
"""

from __future__ import annotations

import json
import math
import os
from collections.abc import Mapping
from dataclasses import dataclass
from urllib.parse import urlsplit


class ProviderError(RuntimeError):
    """Generation failed or returned no complete text; safe to show to users."""


@dataclass(frozen=True)
class GenerationResult:
    """Generated text and provenance, without credentials or the source prompt."""

    text: str
    provider: str
    model: str
    finish_reason: str
    prompt_tokens: int | None = None
    completion_tokens: int | None = None


def _setting(explicit: str | None, *env_names: str, default: str = "") -> str:
    if explicit is not None:
        return explicit.strip()
    return next((os.environ[name].strip() for name in env_names if os.environ.get(name, "").strip()), default)


_RESERVED_PARAMS = frozenset({
    "model", "messages", "stream", "stream_options", "n", "tools", "tool_choice",
    "functions", "function_call", "modalities", "audio", "api_key", "base_url",
    "timeout", "max_retries", "extra_body", "extra_headers", "extra_query",
    "max_tokens", "max_completion_tokens", "max_output_tokens",
})


def validate_model_params(model_params: Mapping[str, object] | None) -> dict:
    """Copy JSON generation settings without accepting transport/contract overrides."""
    if model_params is None:
        return {}
    if not isinstance(model_params, Mapping) or any(not isinstance(key, str) for key in model_params):
        raise ValueError("model_params must be a JSON object with string keys.")
    reserved = _RESERVED_PARAMS.intersection(model_params)
    if reserved:
        raise ValueError(f"model_params contains reserved fields: {', '.join(sorted(reserved))}. Use the dedicated arguments where available.")
    try:
        # Copy nested values too; never modify the caller's configuration.
        return json.loads(json.dumps(dict(model_params), allow_nan=False))
    except (TypeError, ValueError):
        raise ValueError("model_params must contain JSON-serializable values and finite numbers.") from None


class _ChatProvider:
    """Shared transport for the subset of Chat Completions used by both services."""

    name: str
    default_base_url: str
    token_limit_field: str

    def __init__(
        self,
        model: str | None = None,
        *,
        api_key: str | None = None,
        base_url: str | None = None,
        timeout: float = 120.0,
    ) -> None:
        prefix = self.name.upper()
        self.model = _setting(model, f"{prefix}_MODEL", "CLAWBIO_MODEL")
        if not self.model:
            raise ValueError(f"A model is required: pass model or set {prefix}_MODEL / CLAWBIO_MODEL.")
        self.base_url = _setting(base_url, f"{prefix}_BASE_URL", default=self.default_base_url)
        parsed = urlsplit(self.base_url)
        if (parsed.scheme not in {"http", "https"} or not parsed.hostname
                or parsed.username is not None or parsed.password is not None
                or parsed.query or parsed.fragment):
            raise ValueError("base_url must be an HTTP(S) API base URL without credentials, query, or fragment.")
        if not math.isfinite(timeout) or timeout <= 0:
            raise ValueError("timeout must be a finite positive number of seconds.")

        if self.name == "openai":
            key = _setting(api_key, "OPENAI_API_KEY", "LLM_API_KEY")
            if not key:
                raise ValueError("OpenAI requires an API key: set OPENAI_API_KEY or pass api_key.")
        else:
            # The SDK requires a nonempty key; local Ollama ignores this dummy
            # value. Never reuse an OpenAI or bot credential for a local server.
            key = _setting(api_key, default="ollama") or "ollama"

        try:
            from openai import OpenAI
        except ImportError:
            raise ImportError("Text generation requires the openai package; install the ClawBio dependencies.") from None

        options = {}
        if self.name == "ollama":
            # Prevent SDK defaults from forwarding cloud account identifiers.
            options.update(organization="", project="")
        self._client = OpenAI(
            api_key=key,
            base_url=self.base_url,
            timeout=timeout,
            max_retries=0,
            **options,
        )

    def generate(
        self,
        prompt: str,
        *,
        system: str | None = None,
        max_output_tokens: int | None = None,
        model_params: Mapping[str, object] | None = None,
    ) -> GenerationResult:
        """Generate complete text, or raise ProviderError without switching services.

        model_params forwards explicit generation settings (e.g. temperature,
        top_p or reasoning_effort); omitted settings use the server defaults.
        The provider/model validates which settings it supports. Connection,
        message, streaming and token-limit overrides are not accepted here.
        Output limits are optional; OpenAI reasoning models count reasoning
        tokens toward their completion limit as well as visible output.
        """
        from openai import APIConnectionError, APIError, APIStatusError, APITimeoutError

        if not isinstance(prompt, str) or not prompt.strip():
            raise ValueError("prompt must be nonempty text.")
        if max_output_tokens is not None and (
            type(max_output_tokens) is not int or max_output_tokens <= 0
        ):
            raise ValueError("max_output_tokens must be a positive integer.")

        messages = []
        if system:
            messages.append({"role": "system", "content": system})
        messages.append({"role": "user", "content": prompt})
        body = validate_model_params(model_params)
        if max_output_tokens is not None:
            # extra_body works with older SDKs too, which do not expose the
            # newer max_completion_tokens keyword in their Python signature.
            body[self.token_limit_field] = max_output_tokens
        options = {"extra_body": body} if body else {}
        try:
            response = self._client.chat.completions.create(
                model=self.model, messages=messages, **options,
            )
        except APITimeoutError:
            raise ProviderError(f"{self.name} request timed out; check the server or increase timeout.") from None
        except APIConnectionError:
            raise ProviderError(f"Could not connect to {self.name}; check the configured endpoint and server.") from None
        except APIStatusError as exc:
            hints = {
                400: "Check model_params and the settings supported by the selected model.",
                401: "Check credentials.",
                403: "Check credentials and model access.",
                404: "Check the endpoint and model name; for Ollama, ensure the model is installed.",
                429: "Check quota or retry later.",
            }
            hint = hints.get(exc.status_code, "Check the provider server and request settings.")
            # Provider response bodies may echo credentials or source text.
            raise ProviderError(f"{self.name} request failed (HTTP {exc.status_code}). {hint}") from None
        except APIError:
            raise ProviderError(f"{self.name} returned an invalid API response.") from None

        choices = getattr(response, "choices", None)
        if not choices:
            raise ProviderError(f"{self.name} returned no text choices.")
        choice = choices[0]
        finish_reason = getattr(choice, "finish_reason", None)
        if finish_reason != "stop":
            # A partial sentence must not silently become a saved summary.
            reason = finish_reason if finish_reason in {"length", "content_filter", "tool_calls", "function_call"} else "unknown"
            raise ProviderError(f"{self.name} did not finish a text response (finish_reason={reason}).")
        message = getattr(choice, "message", None)
        content = getattr(message, "content", None)
        if getattr(message, "refusal", None) or not isinstance(content, str) or not content.strip():
            raise ProviderError(f"{self.name} returned no usable text (empty output or refusal).")
        usage = getattr(response, "usage", None)
        return GenerationResult(
            text=content.strip(),
            provider=self.name,
            model=getattr(response, "model", None) or self.model,
            finish_reason=finish_reason,
            prompt_tokens=getattr(usage, "prompt_tokens", None),
            completion_tokens=getattr(usage, "completion_tokens", None),
        )

    def close(self) -> None:
        self._client.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()


class OpenAIProvider(_ChatProvider):
    """OpenAI-hosted chat models; requires explicit or environment credentials."""

    name = "openai"
    default_base_url = "https://api.openai.com/v1"
    token_limit_field = "max_completion_tokens"


class OllamaProvider(_ChatProvider):
    """An Ollama server, defaulting to localhost with no user API key required."""

    name = "ollama"
    default_base_url = "http://localhost:11434/v1"
    token_limit_field = "max_tokens"


def create_provider(
    provider: str,
    model: str | None = None,
    *,
    api_key: str | None = None,
    base_url: str | None = None,
    timeout: float = 120.0,
) -> OpenAIProvider | OllamaProvider:
    """Select a provider explicitly; never infer or fall back to another service."""
    providers = {"openai": OpenAIProvider, "ollama": OllamaProvider}
    if provider not in providers:
        raise ValueError("Unsupported provider; choose 'openai' or 'ollama'.")
    return providers[provider](model, api_key=api_key, base_url=base_url, timeout=timeout)
