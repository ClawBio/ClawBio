# Text-generation providers

OpenAI and Ollama expose the same Python interface for generating text. The
selected provider determines the destination, authentication and supported
request parameters. Skill-specific prompts, summary length rules and fallback
behaviour belong to the calling skill.

The existing `FlockRouter` in `flock.py` remains a separate skill-routing client.
`chat.py` contains the shared chat transport, OpenAI/Ollama configuration,
generation result and provider factory.
The PubMed summariser uses these providers with `--summary-method llm` and
`--provider openai` or `--provider ollama`. Its default `first-sentence` mode
copies the abstract opening without calling a model. See the
[PubMed examples](../../skills/pubmed-summariser/examples/README.md) for setup.

## Try a provider

Run commands from the repository root using a Python environment with ClawBio's
dependencies installed. The `openai` SDK is already a project dependency. To
test only this provider interface in a small environment, install `openai`.

### OpenAI

Set `OPENAI_API_KEY` in your environment. Replace `YOUR_OPENAI_MODEL` with a
chat-capable model available to your account:

```sh
python -m clawbio.providers --provider openai --model YOUR_OPENAI_MODEL --prompt "Explain an abstract in one sentence."
```

This makes a real request to OpenAI. The command does not read `.env` files;
export environment variables in your shell or load them in your application.

### Ollama

Start your Ollama server and choose an installed chat model (`ollama list`).
Replace `YOUR_LOCAL_MODEL` with its exact name, including the tag if present:

```sh
python -m clawbio.providers --provider ollama --model YOUR_LOCAL_MODEL --prompt "Explain an abstract in one sentence."
```

No API key is required for the default local server. The provider uses Ollama's
OpenAI-compatible API at `http://localhost:11434/v1`. It supplies the dummy key
`ollama` required by the client library, not a cloud credential. It does not
start the server or download a model automatically.

To use a different server address:

```sh
python -m clawbio.providers --provider ollama --model YOUR_LOCAL_MODEL --base-url http://localhost:11435/v1 --prompt "Say hello." --json
```

`--json` prints the generated text, provider, returned model identifier, finish
reason, and token counts when supplied by the server. Nothing is saved to disk
automatically. Use `--timeout 240` if loading a local model needs more time.

## Python interface

```python
from dataclasses import asdict
from clawbio.providers import ProviderError, create_provider

with create_provider("ollama", model="YOUR_LOCAL_MODEL") as provider:
    try:
        result = provider.generate(
            "Source text to summarise.",
            system="Summarise only the supplied text in two sentences.",
            model_params={"temperature": 0.2},  # only if your model supports it
        )
    except ProviderError as exc:
        # The caller decides whether to stop or record an excerpt fallback.
        print(exc)
    else:
        print(result.text)
        record = asdict(result)
```

Use `create_provider("openai", model="YOUR_OPENAI_MODEL")` for OpenAI, or
instantiate `OpenAIProvider` / `OllamaProvider` directly. `api_key`, `base_url`
and `timeout` are optional keyword arguments. The context manager closes the
HTTP client; a long-lived application can reuse a provider and call `close()`.

## Configuration

Explicit Python arguments or CLI options take precedence over environment
variables. A provider must always be chosen explicitly.

| Setting | OpenAI | Ollama |
|---|---|---|
| Model | `OPENAI_MODEL`, then `CLAWBIO_MODEL` | `OLLAMA_MODEL`, then `CLAWBIO_MODEL` |
| API key | `OPENAI_API_KEY`, then `LLM_API_KEY` | No environment key; dummy `ollama` unless explicitly supplied in Python |
| API base URL | `OPENAI_BASE_URL`, otherwise `https://api.openai.com/v1` | `OLLAMA_BASE_URL`, otherwise `http://localhost:11434/v1` |
| Timeout | 120 seconds by default | 120 seconds by default |

There is no default model. The existing bot's `CLAWBIO_MODEL` and `LLM_API_KEY`
conventions are supported, but its generic `LLM_BASE_URL` is deliberately not
read: a shared cloud endpoint must not redirect a request selected for Ollama.
Use provider-specific endpoint variables or `--base-url`, including the `/v1`
path for standard servers. OpenAI credentials and cloud account identifiers
are not inherited by Ollama.

The selected endpoint receives the supplied text. For literature skills, send
public titles and abstracts; ClawBio's prohibition on uploading patient/genetic
data still applies. A custom Ollama endpoint can be remote, so choosing Ollama
alone is not a guarantee of local processing.

## Behaviour and limits

- This interface supports synchronous, non-streaming text chat only. It does
  not provide embeddings, tools, images, model downloads or structured-output
  guarantees.
- `max_output_tokens` / `--max-output-tokens` is optional. OpenAI receives it as
  `max_completion_tokens`; Ollama receives `max_tokens`. With OpenAI reasoning
  models, the budget also includes reasoning tokens, not just visible text.
- `model_params` forwards explicit JSON generation settings such as
  `temperature`, `top_p`, `seed`, `stop` and `reasoning_effort`. Omitted settings
  use server defaults. The selected model validates support; unsupported
  parameters are not silently dropped or retried with different values.
  The equivalent CLI option is `--model-params '{"temperature": 0.2}'` (shell
  quoting varies). Keep credentials/endpoints in provider configuration and
  output limits in `max_output_tokens`. Messages, streaming, multiple choices,
  tools and other non-text modes cannot be overridden through this dictionary.
- No automatic retries or switching providers are performed. A failed local
  request never triggers a cloud request.
- Timeouts, HTTP errors, refusals, empty responses and incomplete completions
  raise `ProviderError`. Truncated text (`finish_reason=length`) is not returned
  as a successful result; increase the budget if needed. The CLI exits with
  status 1 and a concise message on failure.
- The utility does not decide whether to generate a 300-character excerpt.
  The PubMed summariser handles that fallback: if generation raises
  `ProviderError`, it saves and displays the abstract opening with the failure
  reason. The complete source abstract is preserved in either mode.

## Tests

```sh
python -m pytest clawbio/tests/test_providers.py -q
```

Tests intercept the SDK's HTTP transport. They verify request destinations,
authentication, provider-specific parameters, metadata, failure handling and
the manual-test command without calling a live model or requiring credentials.

API references: [OpenAI Chat Completions](https://developers.openai.com/api/reference/resources/chat/subresources/completions/methods/create),
[Ollama OpenAI compatibility](https://docs.ollama.com/api/openai-compatibility).
