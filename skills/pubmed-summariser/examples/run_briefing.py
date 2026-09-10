"""Edit the settings below, then run this file with your ClawBio interpreter.

Fetches public PubMed records and runs the actual summarizer. See README.md
for environment setup and first-sentence/Ollama/OpenAI settings.
"""
from __future__ import annotations

import importlib.util
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

# Settings to experiment with.
QUERY = "PARP inhibitor resistance in BRCA1-mutated ovarian cancer"
MAX_RESULTS = 2
SUMMARY_METHOD = "first-sentence"  # copy abstract opening (no AI), or "llm" for an AI summary
PROVIDER = "ollama"         # "ollama" or "openai"; only used in LLM mode
MODEL = None               # exact model name, or provider's MODEL environment setting
MODEL_PARAMS = {}          # for the tested Ollama Qwen: {"reasoning_effort": "none"}
LLM_TIMEOUT = 120

REPO_ROOT = Path(__file__).resolve().parents[3]
SKILL_SCRIPT = Path(__file__).resolve().parents[1] / "pubmed_summariser.py"


def main() -> None:
    # Load the same entry point used by the CLI; no separate search or summary logic.
    spec = importlib.util.spec_from_file_location("pubmed_example", SKILL_SCRIPT)
    summarizer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(summarizer)

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S")
    output_dir = REPO_ROOT / "output" / "pubmed-examples" / f"{stamp}-{uuid4().hex[:8]}"
    arguments = ["--query", QUERY, "--max-results", str(MAX_RESULTS),
                 "--summary-method", SUMMARY_METHOD, "--output", str(output_dir)]
    if SUMMARY_METHOD == "llm":
        arguments += ["--provider", PROVIDER, "--llm-timeout", str(LLM_TIMEOUT),
                      "--model-params", json.dumps(MODEL_PARAMS)]
        if MODEL:
            arguments += ["--model", MODEL]
        if PROVIDER == "ollama":
            # Preserve proxy exclusions while keeping localhost traffic local.
            hosts = [host.strip() for key in ("NO_PROXY", "no_proxy")
                     for host in os.environ.get(key, "").split(",") if host.strip()]
            os.environ["NO_PROXY"] = os.environ["no_proxy"] = ",".join(
                dict.fromkeys([*hosts, "localhost", "127.0.0.1", "::1"])
            )

    summarizer.main(arguments)

    # Inspect the saved contract, including the actual method after any fallback.
    result = json.loads((output_dir / "result.json").read_text(encoding="utf-8"))
    for paper in result["data"]["papers"]:
        print(f"PMID {paper['pmid']}: {len(paper['abstract'])} source characters, "
              f"{len(paper['abstract_sections'])} sections, "
              f"summary method={paper['summary']['method']}")
        if paper["summary"]["fallback_reason"]:
            print("Fallback:", paper["summary"]["fallback_reason"])
    print("Inspect full abstracts and metadata:", output_dir / "result.json")
    print("Open the readable report:", output_dir / "report.html")


if __name__ == "__main__":
    main()
