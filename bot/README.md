# RoboTerri ClawBio Bot

A Telegram and Discord bot that runs ClawBio bioinformatics skills using any LLM as the reasoning engine. Send genetic data, medication photos, or natural language questions -- get personalised genomic reports back.

## Features

- Works with **any LLM provider**: OpenAI, Anthropic (via OpenRouter), Google, Mistral, Groq, Together, Ollama, LM Studio, etc.
- Runs all ClawBio skills: pharmgx, equity, nutrigx, metagenomics, compare, drugphoto
- Handles text messages, genetic file uploads (.txt, .csv, .vcf, .fastq), and medication photos
- All genetic data stays local -- nothing leaves your machine
- Reports and figures sent directly in Telegram

## Prerequisites

- Python 3.11+
- A Telegram or Discord account
- An API key from any OpenAI-compatible LLM provider
- ClawBio cloned and working (`python3 clawbio.py run pharmgx --demo`)

## Setup

### 1. Install dependencies

```bash
pip3 install -r bot/requirements.txt
```

### 2a. Create a Telegram bot

1. Open Telegram and search for **@BotFather**
2. Send `/newbot`, choose a name and username
3. Save the **bot token** BotFather gives you

### 3a. Get your Telegram chat ID (optional)

The bot is open to all users by default. To identify yourself as admin (bypasses rate limits):

1. Start a conversation with your new bot
2. Send any message
3. Visit `https://api.telegram.org/bot<YOUR_TOKEN>/getUpdates`
4. Find `"chat":{"id":123456789}` -- that number is your chat ID

### 2b. Create a Discord bot

1. Go to **discord.com/developers/applications** and click **New Application**
2. Go to the **Bot** tab, click **Reset Token**, and save the **bot token**
3. Under **Privileged Gateway Intents**, enable **Message Content Intent**
4. Go to **OAuth2 > URL Generator**, select the `bot` scope
5. Under **Bot Permissions**, select: Send Messages, Attach Files, Read Message History
6. Copy the generated URL and open it to invite the bot to your server

### 3b. Get your Discord channel ID

1. In Discord, go to **User Settings > Advanced** and enable **Developer Mode**
2. Right-click the channel you want the bot to respond in
3. Click **Copy Channel ID**

### 4. Configure environment

Create a `.env` file in the ClawBio root directory:

```
# --- For Telegram ---
TELEGRAM_BOT_TOKEN=your-bot-token-here

# --- For Discord ---
DISCORD_BOT_TOKEN=your-discord-bot-token-here
DISCORD_CHANNEL_ID=your-channel-id-here

# --- LLM (shared by both) ---
LLM_API_KEY=your-api-key-here
CLAWBIO_MODEL=gemini-2.5-flash

# Optional:
TELEGRAM_CHAT_ID=your-chat-id-here    # Admin chat ID (bypasses rate limits)
RATE_LIMIT_PER_HOUR=10                 # Max messages per user per hour (default: 10)
```

### Provider examples

Any provider that speaks the OpenAI chat completions API works. Set `LLM_BASE_URL` to point to your provider:

```bash
# Google Gemini (default -- free tier)
LLM_BASE_URL=https://generativelanguage.googleapis.com/v1beta/openai/
LLM_API_KEY=your-google-api-key
CLAWBIO_MODEL=gemini-2.5-flash

# OpenAI
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
# Telegram bot
python3 bot/roboterri.py

# Discord bot
python3 bot/roboterri_discord.py

# Both at the same time (separate terminals)
python3 bot/roboterri.py &
python3 bot/roboterri_discord.py &
```

## Commands

### Telegram

| Command | Description |
|---|---|
| `/start` | Show welcome message and available commands |
| `/skills` | List all available ClawBio skills |
| `/demo <skill>` | Run a skill with demo data (e.g. `/demo pharmgx`) |

### Discord

| Command | Description |
|---|---|
| `!start` | Show welcome message and available commands |
| `!skills` | List all available ClawBio skills |
| `!demo <skill>` | Run a skill with demo data (e.g. `!demo pharmgx`) |

## Usage

- **Text**: Ask any bioinformatics question -- the LLM routes to the right skill
- **File upload**: Send a 23andMe .txt, AncestryDNA .csv, or VCF file for analysis
- **Photo**: Send a photo of medication packaging for personalised drug guidance
- **Demo**: Type `/demo pharmgx` (Telegram) or `!demo pharmgx` (Discord) to see a pharmacogenomics report with synthetic data

### Zero-cost setup (recommended for public bots)

```bash
# Google Gemini free tier (1,500 requests/day)
LLM_BASE_URL=https://generativelanguage.googleapis.com/v1beta/openai/
LLM_API_KEY=your-google-api-key
CLAWBIO_MODEL=gemini-2.5-flash
```

Get a free Google API key at https://aistudio.google.com/apikey

## Security

- All genetic data is processed locally -- nothing leaves your machine
- Per-user rate limiting prevents abuse (default: 10 messages/hour)
- Admin chat ID is optional -- set it to bypass rate limits for yourself
- Never commit your `.env` file (already in `.gitignore`)

---

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
