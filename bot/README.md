# RoboTerri ClawBio Bot

A Telegram bot that runs ClawBio bioinformatics skills using any LLM as the reasoning engine. Send genetic data, medication photos, or natural language questions -- get personalised genomic reports back.

## Features

- Works with **any LLM provider**: OpenAI, Anthropic (via OpenRouter), Google, Mistral, Groq, Together, Ollama, LM Studio, etc.
- Runs all ClawBio skills: pharmgx, equity, nutrigx, metagenomics, compare, drugphoto
- Handles text messages, genetic file uploads (.txt, .csv, .vcf, .fastq), and medication photos
- All genetic data stays local -- nothing leaves your machine
- Reports and figures sent directly in Telegram

## Prerequisites

- Python 3.11+
- A Telegram account
- An API key from any OpenAI-compatible LLM provider
- ClawBio cloned and working (`python3 clawbio.py run pharmgx --demo`)

## Setup

### 1. Install dependencies

```bash
pip3 install -r bot/requirements.txt
```

### 2. Create a Telegram bot

1. Open Telegram and search for **@BotFather**
2. Send `/newbot`, choose a name and username
3. Save the **bot token** BotFather gives you

### 3. Get your chat ID

1. Start a conversation with your new bot
2. Send any message
3. Visit `https://api.telegram.org/bot<YOUR_TOKEN>/getUpdates`
4. Find `"chat":{"id":123456789}` -- that number is your chat ID

### 4. Configure environment

Create a `.env` file in the ClawBio root directory:

```
TELEGRAM_BOT_TOKEN=your-bot-token-here
TELEGRAM_CHAT_ID=your-chat-id-here
LLM_API_KEY=your-api-key-here
CLAWBIO_MODEL=gpt-4o
```

### Provider examples

Any provider that speaks the OpenAI chat completions API works. Set `LLM_BASE_URL` to point to your provider:

```bash
# OpenAI (default -- no LLM_BASE_URL needed)
LLM_API_KEY=sk-...
CLAWBIO_MODEL=gpt-4o

# Anthropic via OpenRouter
LLM_BASE_URL=https://openrouter.ai/api/v1
LLM_API_KEY=sk-or-...
CLAWBIO_MODEL=anthropic/claude-sonnet-4-5-20250929

# Google Gemini via OpenRouter
LLM_BASE_URL=https://openrouter.ai/api/v1
LLM_API_KEY=sk-or-...
CLAWBIO_MODEL=google/gemini-2.5-pro

# Groq
LLM_BASE_URL=https://api.groq.com/openai/v1
LLM_API_KEY=gsk_...
CLAWBIO_MODEL=llama-3.3-70b-versatile

# Together AI
LLM_BASE_URL=https://api.together.xyz/v1
LLM_API_KEY=...
CLAWBIO_MODEL=meta-llama/Meta-Llama-3.1-70B-Instruct-Turbo

# Ollama (local, free)
LLM_BASE_URL=http://localhost:11434/v1
LLM_API_KEY=ollama
CLAWBIO_MODEL=llama3.1

# LM Studio (local, free)
LLM_BASE_URL=http://localhost:1234/v1
LLM_API_KEY=lm-studio
CLAWBIO_MODEL=local-model
```

> **Note**: Photo/drug detection requires a model with vision capabilities (e.g. gpt-4o, claude-sonnet, gemini-pro). Tool calling requires a model that supports function calling. Most major providers support both.

### 5. Run

```bash
python3 bot/roboterri.py
```

## Commands

| Command | Description |
|---|---|
| `/start` | Show welcome message and available commands |
| `/skills` | List all available ClawBio skills |
| `/demo <skill>` | Run a skill with demo data (e.g. `/demo pharmgx`) |

## Usage

- **Text**: Ask any bioinformatics question -- the LLM routes to the right skill
- **File upload**: Send a 23andMe .txt, AncestryDNA .csv, or VCF file for analysis
- **Photo**: Send a photo of medication packaging for personalised drug guidance
- **Demo**: Type `/demo pharmgx` to see a pharmacogenomics report with synthetic data

## Security

- The bot only responds to your configured `TELEGRAM_CHAT_ID`
- All genetic data is processed locally
- Never commit your `.env` file (already in `.gitignore`)

---

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
