"""Small manual test command: python -m clawbio.providers --help."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict

from . import ProviderError, create_provider


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Test a configured OpenAI or Ollama text provider.")
    parser.add_argument("--provider", required=True, choices=("openai", "ollama"))
    parser.add_argument("--model", help="Model name; otherwise use provider-specific MODEL or CLAWBIO_MODEL.")
    parser.add_argument("--prompt", required=True, help="Text to send to the selected model.")
    parser.add_argument("--system", help="Optional instructions for the model.")
    parser.add_argument("--base-url", help="Full API base URL, including /v1 for a standard endpoint.")
    parser.add_argument("--timeout", type=float, default=120.0, help="Request timeout in seconds (default 120).")
    parser.add_argument("--max-output-tokens", type=int, help="Optional completion token limit.")
    parser.add_argument("--model-params", help='Generation settings as JSON, e.g. {"temperature": 0.2}.')
    parser.add_argument("--json", action="store_true", help="Print text plus provider/model/usage metadata as JSON.")
    args = parser.parse_args(argv)
    try:
        model_params = json.loads(args.model_params) if args.model_params is not None else None
        with create_provider(
            args.provider, model=args.model, base_url=args.base_url, timeout=args.timeout,
        ) as provider:
            result = provider.generate(
                args.prompt, system=args.system, max_output_tokens=args.max_output_tokens,
                model_params=model_params,
            )
    except (ValueError, ImportError, ProviderError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(asdict(result), ensure_ascii=False, indent=2) if args.json else result.text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
