"""Provider contract tests using the real SDK with intercepted HTTP requests."""

from __future__ import annotations

import json

import pytest
from openai import _base_client

# SDK 2.x (the repository lock) uses httpx; SDK 3.x uses httpx2.
# Intercept the transport actually used by the installed SDK in either case.
httpx = getattr(_base_client, "httpx", None) or _base_client.httpx2

from clawbio.providers import (
    OllamaProvider,
    OpenAIProvider,
    ProviderError,
    create_provider,
)


@pytest.fixture(autouse=True)
def clean_environment(monkeypatch):
    for name in (
        "OPENAI_API_KEY", "OPENAI_BASE_URL", "OPENAI_MODEL", "OPENAI_ORG_ID",
        "OPENAI_PROJECT_ID", "LLM_API_KEY", "LLM_BASE_URL", "CLAWBIO_MODEL",
        "OLLAMA_BASE_URL", "OLLAMA_MODEL",
    ):
        monkeypatch.delenv(name, raising=False)


@pytest.fixture
def endpoint(monkeypatch):
    state = {
        "requests": [],
        "status": 200,
        "body": {
            "id": "chatcmpl-test",
            "object": "chat.completion",
            "created": 0,
            "model": "resolved-model",
            "choices": [{
                "index": 0,
                "message": {"role": "assistant", "content": " A short answer. "},
                "finish_reason": "stop",
            }],
            "usage": {"prompt_tokens": 20, "completion_tokens": 5, "total_tokens": 25},
        },
    }

    def send(client, request, **kwargs):
        state["requests"].append(request)
        if state.get("error"):
            raise state["error"]("transport failed", request=request)
        return httpx.Response(state["status"], json=state["body"], request=request)

    # Intercept every SDK request: these tests cannot call a real model.
    monkeypatch.setattr(httpx.Client, "send", send)
    return state


def test_openai_sends_key_and_messages_to_openai(endpoint, monkeypatch):
    monkeypatch.setenv("OPENAI_API_KEY", "test-openai-key")
    with create_provider("openai", model="test-model") as provider:
        result = provider.generate("Source text", system="Summarise faithfully.", max_output_tokens=400)
    request = endpoint["requests"][0]
    assert str(request.url) == "https://api.openai.com/v1/chat/completions"
    assert request.headers["authorization"] == "Bearer test-openai-key"
    body = json.loads(request.content)
    assert body["model"] == "test-model"
    assert body["messages"] == [
        {"role": "system", "content": "Summarise faithfully."},
        {"role": "user", "content": "Source text"},
    ]
    assert body["max_completion_tokens"] == 400
    assert "max_tokens" not in body
    assert "temperature" not in body  # not supported uniformly across models
    assert result.text == "A short answer."
    assert result.provider == "openai"
    assert result.model == "resolved-model"
    assert result.prompt_tokens == 20
    assert result.completion_tokens == 5


def test_ollama_needs_no_key_and_ignores_cloud_configuration(endpoint, monkeypatch):
    monkeypatch.setenv("OPENAI_API_KEY", "never-send-this")
    monkeypatch.setenv("LLM_API_KEY", "never-send-this-either")
    monkeypatch.setenv("OPENAI_BASE_URL", "https://cloud.example/v1")
    monkeypatch.setenv("LLM_BASE_URL", "https://another-cloud.example/v1")
    monkeypatch.setenv("OPENAI_ORG_ID", "private-org")
    monkeypatch.setenv("OPENAI_PROJECT_ID", "private-project")
    with create_provider("ollama", model="local-model") as provider:
        result = provider.generate("Hello", max_output_tokens=100)
    request = endpoint["requests"][0]
    assert str(request.url) == "http://localhost:11434/v1/chat/completions"
    assert request.headers["authorization"] == "Bearer ollama"
    assert request.headers.get("openai-organization", "") == ""
    assert request.headers.get("openai-project", "") == ""
    body = json.loads(request.content)
    assert body["max_tokens"] == 100
    assert "max_completion_tokens" not in body
    assert body["messages"] == [{"role": "user", "content": "Hello"}]
    assert result.provider == "ollama"


@pytest.mark.parametrize("name,cls", [("openai", OpenAIProvider), ("ollama", OllamaProvider)])
def test_provider_specific_environment(name, cls, monkeypatch, endpoint):
    monkeypatch.setenv(f"{name.upper()}_MODEL", "configured-model")
    monkeypatch.setenv(f"{name.upper()}_BASE_URL", "http://localhost:8000/v1/")
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")
    with create_provider(name) as provider:
        assert isinstance(provider, cls)
        provider.generate("Hello")
    assert str(endpoint["requests"][0].url) == "http://localhost:8000/v1/chat/completions"
    assert json.loads(endpoint["requests"][0].content)["model"] == "configured-model"


def test_explicit_settings_override_environment(monkeypatch, endpoint):
    monkeypatch.setenv("OPENAI_MODEL", "wrong-model")
    monkeypatch.setenv("OPENAI_API_KEY", "wrong-key")
    monkeypatch.setenv("OPENAI_BASE_URL", "https://wrong.example/v1")
    with OpenAIProvider(model="chosen-model", api_key="chosen-key", base_url="http://localhost:9000/v1") as provider:
        provider.generate("Hello")
    request = endpoint["requests"][0]
    assert request.headers["authorization"] == "Bearer chosen-key"
    assert request.url.host == "localhost"
    assert json.loads(request.content)["model"] == "chosen-model"


def test_existing_bot_key_and_model_conventions(monkeypatch, endpoint):
    monkeypatch.setenv("LLM_API_KEY", "bot-key")
    monkeypatch.setenv("CLAWBIO_MODEL", "bot-model")
    with create_provider("openai") as provider:
        provider.generate("Hello")
    assert endpoint["requests"][0].headers["authorization"] == "Bearer bot-key"
    assert json.loads(endpoint["requests"][0].content)["model"] == "bot-model"


def test_missing_openai_key_fails_before_request(endpoint):
    with pytest.raises(ValueError, match="OPENAI_API_KEY"):
        create_provider("openai", model="test-model")
    assert endpoint["requests"] == []


@pytest.mark.parametrize("name", ["openai", "ollama"])
def test_model_is_required(name, endpoint):
    with pytest.raises(ValueError, match="model"):
        create_provider(name, api_key="test-key")
    assert endpoint["requests"] == []


def test_unknown_provider_is_not_silently_substituted(endpoint):
    with pytest.raises(ValueError, match="Unsupported provider"):
        create_provider("unknown", model="test-model")
    assert endpoint["requests"] == []


@pytest.mark.parametrize("timeout", [0, -1, float("nan"), float("inf")])
def test_invalid_timeout_is_rejected(timeout, endpoint):
    with pytest.raises(ValueError, match="timeout"):
        OllamaProvider(model="local-model", timeout=timeout)
    assert endpoint["requests"] == []


@pytest.mark.parametrize("url", ["localhost:11434", "file:///tmp/model", "https://host/v1?key=secret"])
def test_invalid_base_url_is_rejected(url, endpoint):
    with pytest.raises(ValueError, match="base_url"):
        OllamaProvider(model="local-model", base_url=url)
    assert endpoint["requests"] == []


def test_timeout_is_forwarded_and_client_closed(endpoint):
    with OllamaProvider(model="local-model", timeout=90) as provider:
        provider.generate("Hello")
        client = provider._client
    assert endpoint["requests"][0].extensions["timeout"]["read"] == 90
    assert client.is_closed()


@pytest.mark.parametrize("limit", [0, -5, 1.5, True])
def test_invalid_output_limit_does_not_send_request(limit, endpoint):
    with OllamaProvider(model="local-model") as provider:
        with pytest.raises(ValueError, match="max_output_tokens"):
            provider.generate("Hello", max_output_tokens=limit)
    assert endpoint["requests"] == []


def test_blank_prompt_does_not_send_request(endpoint):
    with OllamaProvider(model="local-model") as provider:
        with pytest.raises(ValueError, match="prompt"):
            provider.generate("  ")
    assert endpoint["requests"] == []


@pytest.mark.parametrize("status", [401, 404, 429, 500])
def test_api_errors_are_safe_and_never_switch_provider(status, endpoint):
    endpoint["status"] = status
    endpoint["body"] = {"error": {"message": "secret-key source-abstract", "type": "test_error"}}
    with OllamaProvider(model="local-model") as provider:
        with pytest.raises(ProviderError) as error:
            provider.generate("source-abstract")
    assert str(status) in str(error.value)
    assert "secret-key" not in str(error.value)
    assert "source-abstract" not in str(error.value)
    assert len(endpoint["requests"]) == 1


@pytest.mark.parametrize("error_type,match", [(httpx.ReadTimeout, "timed out"), (httpx.ConnectError, "connect")])
def test_transport_failure_is_actionable(error_type, match, endpoint):
    endpoint["error"] = error_type
    with OllamaProvider(model="local-model") as provider:
        with pytest.raises(ProviderError, match=match):
            provider.generate("Hello")
    assert len(endpoint["requests"]) == 1


@pytest.mark.parametrize("content", [None, "", "   "])
def test_empty_model_output_is_not_a_summary(content, endpoint):
    endpoint["body"]["choices"][0]["message"]["content"] = content
    with OllamaProvider(model="local-model") as provider:
        with pytest.raises(ProviderError, match="text"):
            provider.generate("Hello")


@pytest.mark.parametrize("reason", ["length", "content_filter", "tool_calls"])
def test_incomplete_or_non_text_completion_is_rejected(reason, endpoint):
    endpoint["body"]["choices"][0]["finish_reason"] = reason
    with OllamaProvider(model="local-model") as provider:
        with pytest.raises(ProviderError, match=reason):
            provider.generate("Hello")


def test_missing_choices_is_a_provider_error(endpoint):
    endpoint["body"]["choices"] = []
    with OllamaProvider(model="local-model") as provider:
        with pytest.raises(ProviderError):
            provider.generate("Hello")


def test_usage_is_optional(endpoint):
    endpoint["body"].pop("usage")
    with OllamaProvider(model="local-model") as provider:
        result = provider.generate("Hello")
    assert result.prompt_tokens is None
    assert result.completion_tokens is None


def test_cli_prints_text(endpoint, capsys):
    from clawbio.providers.__main__ import main

    assert main(["--provider", "ollama", "--model", "local-model", "--prompt", "Hello"]) == 0
    assert capsys.readouterr().out == "A short answer.\n"


def test_cli_json_records_provider_and_model(endpoint, capsys):
    from clawbio.providers.__main__ import main

    assert main(["--provider", "ollama", "--model", "local-model", "--prompt", "Hello", "--json"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["provider"] == "ollama"
    assert result["model"] == "resolved-model"
    assert result["text"] == "A short answer."


def test_cli_configuration_failure_is_concise(endpoint, capsys):
    from clawbio.providers.__main__ import main

    assert main(["--provider", "openai", "--model", "test-model", "--prompt", "Hello"]) == 1
    captured = capsys.readouterr()
    assert "OPENAI_API_KEY" in captured.err
    assert "Traceback" not in captured.err
    assert captured.out == ""
    assert endpoint["requests"] == []


@pytest.mark.parametrize("name", ["openai", "ollama"])
def test_model_params_are_forwarded_without_mutating_input(name, endpoint):
    params = {"temperature": 0.2, "top_p": 0.9, "seed": 42, "stop": ["END"]}
    with create_provider(name, model="test-model", api_key="test-key") as provider:
        provider.generate("Hello", max_output_tokens=300, model_params=params)
    body = json.loads(endpoint["requests"][0].content)
    for key, value in params.items():
        assert body[key] == value
    token_field = "max_completion_tokens" if name == "openai" else "max_tokens"
    assert body[token_field] == 300
    assert params == {"temperature": 0.2, "top_p": 0.9, "seed": 42, "stop": ["END"]}


def test_provider_specific_generation_fields_can_be_forwarded(endpoint):
    with OllamaProvider(model="local-model") as provider:
        provider.generate("Hello", model_params={"reasoning_effort": "low"})
    assert json.loads(endpoint["requests"][0].content)["reasoning_effort"] == "low"


@pytest.mark.parametrize("key", [
    "model", "messages", "stream", "stream_options", "n", "tools", "tool_choice",
    "functions", "function_call", "modalities", "audio", "api_key", "base_url",
    "timeout", "max_retries", "extra_body", "extra_headers", "extra_query",
    "max_tokens", "max_completion_tokens", "max_output_tokens",
])
def test_model_params_cannot_override_connection_or_text_contract(key, endpoint):
    with OllamaProvider(model="local-model") as provider:
        with pytest.raises(ValueError, match=key):
            provider.generate("Hello", model_params={key: "override"})
    assert endpoint["requests"] == []


@pytest.mark.parametrize("params", [
    ["temperature"], {1: "value"}, {"temperature": float("nan")},
    {"temperature": float("inf")}, {"stop": {"not", "json"}},
])
def test_model_params_must_be_a_json_object(params, endpoint):
    with OllamaProvider(model="local-model") as provider:
        with pytest.raises(ValueError, match="model_params"):
            provider.generate("Hello", model_params=params)
    assert endpoint["requests"] == []


def test_invalid_model_setting_reports_parameter_hint(endpoint):
    endpoint["status"] = 400
    endpoint["body"] = {"error": {"message": "Do not echo source-text or secret-key"}}
    with OllamaProvider(model="local-model") as provider:
        with pytest.raises(ProviderError, match="model_params"):
            provider.generate("source-text", model_params={"temperature": 0.2})
    assert len(endpoint["requests"]) == 1


def test_cli_accepts_model_params_json(endpoint, capsys):
    from clawbio.providers.__main__ import main

    assert main([
        "--provider", "ollama", "--model", "local-model", "--prompt", "Hello",
        "--model-params", '{"temperature": 0.2}',
    ]) == 0
    assert json.loads(endpoint["requests"][0].content)["temperature"] == 0.2
    assert capsys.readouterr().out == "A short answer.\n"


@pytest.mark.parametrize("params", ['{"temperature":', '[]', '{"stream": true}'])
def test_cli_rejects_invalid_model_params_without_network(params, endpoint, capsys):
    from clawbio.providers.__main__ import main

    assert main([
        "--provider", "ollama", "--model", "local-model", "--prompt", "Hello",
        "--model-params", params,
    ]) == 1
    assert "Error:" in capsys.readouterr().err
    assert endpoint["requests"] == []
