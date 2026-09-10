# Run a PubMed briefing

[`run_briefing.py`](run_briefing.py) demonstrates calling the actual summarizer
from Python and reading its saved JSON. It searches for **PARP inhibitor resistance
in BRCA1-mutated ovarian cancer**, retrieves up to two papers, preserves their
complete abstracts, and copies each abstract's first sentence (up to 300 characters).
Change `QUERY` to your own research topic.

## Use your existing environment

From the repository root, activate the environment containing your ClawBio
dependencies. Adjust the activation path if your environment is outside the repo.
The example does not automatically read `.env` files or install a dotenv package.
Copy `.env.example` to `.env`, replace the placeholders, and load the variables
into the current shell before running an LLM example. Never commit `.env` or put
real keys in `.env.example`.

PowerShell:

```powershell
& .venv\Scripts\Activate.ps1
# First time only: Copy-Item skills/pubmed-summariser/examples/.env.example .env
# Load simple KEY=value lines into this PowerShell process:
Get-Content .env | ForEach-Object { if ($_ -match '^\s*([^#=]+)\s*=\s*(.*)\s*$') { Set-Item -Path "Env:$($matches[1].Trim())" -Value $matches[2].Trim() } }
python -c "import sys; print(sys.executable)"
python skills/pubmed-summariser/examples/run_briefing.py
```

Bash:

```bash
source .venv/bin/activate
# First time only: cp skills/pubmed-summariser/examples/.env.example .env
set -a; source .env; set +a
python -c "import sys; print(sys.executable)"
python skills/pubmed-summariser/examples/run_briefing.py
```

You can also run the script directly in PyCharm with your existing interpreter;
shell activation is unnecessary in that case. The example installs nothing and
requires live PubMed access. Keep patient information out of search queries.

## Try a different summary method

Edit the settings at the top of `run_briefing.py` and rerun it.

**Abstract opening, no AI:**

```python
SUMMARY_METHOD = "first-sentence"
```

**Ollama:** start the server and install your selected model first. For the
locally tested Qwen model:

```python
SUMMARY_METHOD = "llm"
PROVIDER = "ollama"
MODEL = "qwen3.5:4b"
MODEL_PARAMS = {"reasoning_effort": "none"}
```

Use parameters supported by your model. An empty dictionary uses server defaults.
The provider defaults to `http://localhost:11434/v1`; `OLLAMA_BASE_URL` overrides it.

**OpenAI:** load `OPENAI_API_KEY` as above or configure it in your PyCharm run
configuration. The provider reads `OPENAI_API_KEY` (and accepts `LLM_API_KEY` as
a fallback); `OPENAI_KEY` is not recognized. Do not put keys in this example file.
Select a model available to your account; requests send public titles/abstracts
and can incur API charges.

```python
SUMMARY_METHOD = "llm"
PROVIDER = "openai"
MODEL = "YOUR_MODEL"
MODEL_PARAMS = {}
```

## Inspect the result

Each run creates a unique directory under ignored `output/pubmed-examples/`:

- `result.json`: full abstracts, structured sections, summaries and provenance.
- `report.html`: short summaries with expandable full abstracts.
- `report.md`: Markdown briefing.
- `reproducibility/`: command, suggested environment and output checksums.

The script prints the actual summary method for each paper. A failed LLM request
shows the abstract opening labeled **Abstract opening — no AI**; missing abstracts have method `unavailable`. Verify
generated summaries against the preserved source abstracts.

For all CLI options and output field definitions, see [SKILL.md](../SKILL.md).
