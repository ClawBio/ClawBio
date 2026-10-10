#!/usr/bin/env python3
"""
roboterri_whatsapp.py — RoboTerri ClawBio WhatsApp Bot
=======================================================
A WhatsApp bot using Meta's Cloud API that runs ClawBio bioinformatics
skills using any LLM as the reasoning engine. Handles text messages,
genetic file uploads, and medication photos.

Works with any OpenAI-compatible provider: OpenAI, Anthropic (via proxy),
Google, Mistral, Groq, Together, OpenRouter, Ollama, LM Studio, etc.

Prerequisites:
    uv sync --group bot

Usage:
    # Set environment variables in .env (see bot/README.md)
    uv run python bot/roboterri_whatsapp.py

Setup (Meta WhatsApp Cloud API):
    1. Go to https://developers.facebook.com and create an app (type: Business)
    2. Add the WhatsApp product to your app
    3. In WhatsApp > API Setup, get your:
       - Phone Number ID
       - Permanent access token (System User token recommended)
    4. Configure a webhook:
       - URL: https://your-domain/webhook  (use ngrok for local dev)
       - Verify token: set WHATSAPP_VERIFY_TOKEN in .env
       - Subscribe to: messages
    5. Set env vars in .env:
       WHATSAPP_TOKEN=your_permanent_access_token
       WHATSAPP_PHONE_NUMBER_ID=your_phone_number_id
       WHATSAPP_VERIFY_TOKEN=your_chosen_verify_token
       WHATSAPP_ADMIN_PHONE=your_phone_number  (optional, e.g. 447700900000)
"""

import asyncio
import base64
import json
import logging
import os
import re
import shutil
import sys
import tempfile
import threading
import time
from datetime import datetime
from pathlib import Path

import requests
from dotenv import load_dotenv
from flask import Flask, request, jsonify
from openai import AsyncOpenAI, APIError

# Shared bot security helpers (identity isolation + webhook signature check).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from security import (
    UnsafePath,
    is_allowed_extension,
    safe_write_path,
    sanitize_filename,
    scoped_get,
    is_sender_allowed,
    verify_whatsapp_signature,
    upload_tmp_path,
)

# --------------------------------------------------------------------------- #
# Config
# --------------------------------------------------------------------------- #

load_dotenv()

WHATSAPP_TOKEN = os.environ.get("WHATSAPP_TOKEN", "")
WHATSAPP_PHONE_NUMBER_ID = os.environ.get("WHATSAPP_PHONE_NUMBER_ID", "")
WHATSAPP_VERIFY_TOKEN = os.environ.get("WHATSAPP_VERIFY_TOKEN", "")
WHATSAPP_APP_SECRET = os.environ.get("WHATSAPP_APP_SECRET", "")
WHATSAPP_ADMIN_PHONE = os.environ.get("WHATSAPP_ADMIN_PHONE", "")
# Sender allow-list: only these numbers (plus the admin) may use the bot.
# Public access is an explicit opt-in via WHATSAPP_PUBLIC=1, never the default.
WHATSAPP_ALLOWED_NUMBERS = {
    n.strip() for n in os.environ.get("WHATSAPP_ALLOWED_NUMBERS", "").split(",") if n.strip()
}
WHATSAPP_PUBLIC = os.environ.get("WHATSAPP_PUBLIC", "").lower() in ("1", "true", "yes")
WHATSAPP_PORT = int(os.environ.get("WHATSAPP_PORT", "5001"))

LLM_API_KEY = os.environ.get("LLM_API_KEY", "")
LLM_BASE_URL = os.environ.get("LLM_BASE_URL", "")
CLAWBIO_MODEL = os.environ.get("CLAWBIO_MODEL", "gpt-4.1-mini")

# Rate limiting: messages per user per hour (0 = unlimited)
RATE_LIMIT_PER_HOUR = int(os.environ.get("RATE_LIMIT_PER_HOUR", "10"))

if not WHATSAPP_TOKEN:
    print("Error: WHATSAPP_TOKEN not set. See bot/README.md for setup.")
    sys.exit(1)
if not WHATSAPP_PHONE_NUMBER_ID:
    print("Error: WHATSAPP_PHONE_NUMBER_ID not set. See bot/README.md for setup.")
    sys.exit(1)
if not LLM_API_KEY:
    print("Error: LLM_API_KEY not set. See bot/README.md for setup.")
    sys.exit(1)
if not WHATSAPP_APP_SECRET:
    print("WARNING: WHATSAPP_APP_SECRET not set -- all webhook POSTs are rejected "
          "(fail closed). Set it to your Meta app secret to enable signature checks.")

CLAWBIO_DIR = Path(__file__).resolve().parent.parent
CLAWBIO_PY = CLAWBIO_DIR / "clawbio.py"
SOUL_MD = CLAWBIO_DIR / "SOUL.md"
OUTPUT_DIR = CLAWBIO_DIR / "output"
DATA_DIR = CLAWBIO_DIR / "data"

# Owner's genome — used as default when admin asks about their own PGx/nutrition/risk
OWNER_GENOME = CLAWBIO_DIR / "skills" / "genome-compare" / "data" / "manuel_corpas_23andme.txt.gz"

# Security limits
MAX_UPLOAD_BYTES = 50 * 1024 * 1024  # 50 MB
MAX_PHOTO_BYTES = 20 * 1024 * 1024   # 20 MB

# WhatsApp API base URL
WA_API_BASE = f"https://graph.facebook.com/v21.0/{WHATSAPP_PHONE_NUMBER_ID}"
WA_HEADERS = {
    "Authorization": f"Bearer {WHATSAPP_TOKEN}",
    "Content-Type": "application/json",
}

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
)
logger = logging.getLogger("roboterri-whatsapp")


# ---------------------------------------------------------------------------
# Redact tokens from log output
# ---------------------------------------------------------------------------
class _TokenRedactFilter(logging.Filter):
    def __init__(self, token: str):
        super().__init__()
        self._token = token

    def filter(self, record: logging.LogRecord) -> bool:
        if self._token and self._token in record.getMessage():
            record.msg = str(record.msg).replace(self._token, "[REDACTED]")
            if isinstance(record.args, tuple):
                record.args = tuple(
                    str(a).replace(self._token, "[REDACTED]")
                    for a in record.args
                )
            elif isinstance(record.args, dict):
                record.args = {
                    k: str(v).replace(self._token, "[REDACTED]")
                    for k, v in record.args.items()
                }
        return True


if WHATSAPP_TOKEN:
    _redact = _TokenRedactFilter(WHATSAPP_TOKEN)
    logging.getLogger("urllib3").addFilter(_redact)
    logging.getLogger("requests").addFilter(_redact)


# ---------------------------------------------------------------------------
# Structured audit log (JSONL)
# ---------------------------------------------------------------------------
_AUDIT_LOG_DIR = CLAWBIO_DIR / "bot" / "logs"
_AUDIT_LOG_DIR.mkdir(parents=True, exist_ok=True)
_AUDIT_LOG_PATH = _AUDIT_LOG_DIR / "audit_whatsapp.jsonl"


def _audit(event: str, **kwargs):
    """Append a structured JSON event to the audit log."""
    from datetime import timezone as _tz
    entry = {"ts": datetime.now(_tz.utc).isoformat(), "event": event, **kwargs}
    try:
        with open(_AUDIT_LOG_PATH, "a", encoding="utf-8") as f:
            f.write(json.dumps(entry, default=str) + "\n")
    except OSError:
        pass


def is_admin(phone: str) -> bool:
    """Check if the phone number is the admin."""
    return bool(WHATSAPP_ADMIN_PHONE) and phone == WHATSAPP_ADMIN_PHONE


# --------------------------------------------------------------------------- #
# System prompt
# --------------------------------------------------------------------------- #

if SOUL_MD.exists():
    _soul = SOUL_MD.read_text(encoding="utf-8")
    logger.info(f"Loaded SOUL.md ({len(_soul)} chars)")
else:
    _soul = (
        "You are RoboTerri, an AI agent inspired by Professor Teresa K. Attwood. "
        "Respond in Terri's warm, direct style with characteristic dashes and emoticons."
    )
    logger.warning("SOUL.md not found, using fallback prompt")

ROLE_GUARDRAILS = """
Operational constraints:
1. You are a bioinformatics assistant powered by ClawBio skills.
2. Keep outputs concise, evidence-led, and explicit about confidence and gaps.
3. When the user sends a genetic data file (23andMe .txt, AncestryDNA .csv, VCF, FASTQ) or asks about pharmacogenomics, nutrigenomics, equity scoring, metagenomics, or genome comparison, use the clawbio tool. When the user asks about disease risk, polygenic risk scores, or "what am I at risk for", use skill='prs'. For a unified profile report use skill='profile'. For gene-drug database lookups use skill='clinpgx'. For variant lookups (rsID, "look up rs...") use skill='gwas'. For quick demos say "run pharmgx demo", "run prs demo", "run profile demo" etc. Reports and figures are sent automatically after your summary.
4. TOOL OUTPUT RELAY (STRICT): When the clawbio tool returns results, relay the output VERBATIM. Do not paraphrase, summarise, or rewrite tool results. Tool outputs contain precise data (IBS scores, percentages, gene-drug interactions) that must not be altered. You may add a brief intro line before the verbatim output but never replace or condense it.
5. OWNER GENOME: The bot owner (admin) has their genome pre-loaded. When the admin asks about "my pharmacogenomics", "my risk", "my nutrition", "my genome", or similar personal queries WITHOUT uploading a file, use mode='file' — the system will automatically use the owner's genome. Do NOT ask the admin to upload a file.
6. DEMO FALLBACK: When a non-admin user asks about pharmacogenomics, nutrigenomics, risk scores, or any skill that needs genetic data but has NOT uploaded a file, do NOT just ask for a file and stop. Instead, offer to run the demo with built-in synthetic data (mode='demo') so they can see the skill in action. Example: "I can run a demo with synthetic data so you can see what the report looks like — shall I go ahead?" If they agree (or if the request is clearly exploratory), run it immediately.
"""

SYSTEM_PROMPT = f"{_soul}\n\n{ROLE_GUARDRAILS}"

# --------------------------------------------------------------------------- #
# State
# --------------------------------------------------------------------------- #

_client_kwargs = {"api_key": LLM_API_KEY}
if LLM_BASE_URL:
    _client_kwargs["base_url"] = LLM_BASE_URL
llm = AsyncOpenAI(**_client_kwargs)

conversations: dict[str, list] = {}  # keyed by phone number
MAX_HISTORY = 20

# Per-user received file storage
_received_files: dict[str, dict] = {}

# Pending media queue: phone -> list of {"type": "document"|"photo", "path": str}
_pending_media: dict[str, list[dict]] = {}

# Pending text queue: bypass LLM paraphrasing for compare/drugphoto
_pending_text: dict[str, list[str]] = {}

# Per-user voice reply toggle
_voice_enabled: dict[str, bool] = {}

# Dedup: track recently processed message IDs
_processed_messages: set[str] = set()
_processed_messages_max = 1000

BOT_START_TIME = time.time()

# Async event loop for running coroutines from Flask threads
_loop = asyncio.new_event_loop()
_loop_thread = threading.Thread(target=_loop.run_forever, daemon=True)
_loop_thread.start()


def run_async(coro):
    """Run an async coroutine from a sync Flask context."""
    future = asyncio.run_coroutine_threadsafe(coro, _loop)
    return future.result(timeout=180)


# --------------------------------------------------------------------------- #
# WhatsApp API helpers
# --------------------------------------------------------------------------- #


def wa_send_text(to: str, text: str):
    """Send a text message via WhatsApp Cloud API."""
    # WhatsApp has a 4096 char limit per message
    MAX_LEN = 4096
    text = strip_markup(text)
    if not text:
        return

    chunks = []
    if len(text) <= MAX_LEN:
        chunks = [text]
    else:
        remaining = text
        while remaining:
            if len(remaining) <= MAX_LEN:
                chunks.append(remaining)
                break
            split_at = remaining.rfind("\n\n", 0, MAX_LEN)
            if split_at == -1:
                split_at = remaining.rfind("\n", 0, MAX_LEN)
            if split_at == -1:
                split_at = MAX_LEN
            chunks.append(remaining[:split_at])
            remaining = remaining[split_at:].lstrip("\n")

    for chunk in chunks:
        if not chunk.strip():
            continue
        payload = {
            "messaging_product": "whatsapp",
            "to": to,
            "type": "text",
            "text": {"body": chunk},
        }
        try:
            resp = requests.post(
                f"{WA_API_BASE}/messages",
                headers=WA_HEADERS,
                json=payload,
                timeout=30,
            )
            if resp.status_code != 200:
                logger.error(f"WhatsApp send failed ({resp.status_code}): {resp.text[:300]}")
        except Exception as e:
            logger.error(f"WhatsApp send error: {e}")


def wa_send_document(to: str, filepath: str, caption: str = ""):
    """Upload and send a document via WhatsApp Cloud API."""
    path = Path(filepath)
    if not path.exists():
        return

    # Step 1: Upload media
    mime_map = {
        ".md": "text/markdown",
        ".txt": "text/plain",
        ".csv": "text/csv",
        ".html": "text/html",
        ".json": "application/json",
        ".png": "image/png",
        ".jpg": "image/jpeg",
        ".jpeg": "image/jpeg",
        ".pdf": "application/pdf",
        ".mp3": "audio/mpeg",
    }
    mime = mime_map.get(path.suffix.lower(), "application/octet-stream")

    try:
        upload_resp = requests.post(
            f"{WA_API_BASE}/media",
            headers={"Authorization": f"Bearer {WHATSAPP_TOKEN}"},
            files={"file": (path.name, open(path, "rb"), mime)},
            data={"messaging_product": "whatsapp"},
            timeout=60,
        )
        if upload_resp.status_code != 200:
            logger.error(f"Media upload failed: {upload_resp.text[:300]}")
            return
        media_id = upload_resp.json().get("id")
        if not media_id:
            logger.error("No media ID returned from upload")
            return
    except Exception as e:
        logger.error(f"Media upload error: {e}")
        return

    # Step 2: Send media message
    if path.suffix.lower() in (".png", ".jpg", ".jpeg"):
        msg_type = "image"
        media_obj = {"id": media_id}
        if caption:
            media_obj["caption"] = caption
    else:
        msg_type = "document"
        media_obj = {"id": media_id, "filename": path.name}
        if caption:
            media_obj["caption"] = caption

    payload = {
        "messaging_product": "whatsapp",
        "to": to,
        "type": msg_type,
        msg_type: media_obj,
    }
    try:
        resp = requests.post(
            f"{WA_API_BASE}/messages",
            headers=WA_HEADERS,
            json=payload,
            timeout=30,
        )
        if resp.status_code != 200:
            logger.error(f"Media send failed: {resp.text[:300]}")
    except Exception as e:
        logger.error(f"Media send error: {e}")


def wa_download_media(media_id: str) -> bytes | None:
    """Download media from WhatsApp Cloud API by media ID."""
    try:
        # Get media URL
        resp = requests.get(
            f"https://graph.facebook.com/v21.0/{media_id}",
            headers={"Authorization": f"Bearer {WHATSAPP_TOKEN}"},
            timeout=30,
        )
        if resp.status_code != 200:
            logger.error(f"Media URL fetch failed: {resp.text[:300]}")
            return None
        media_url = resp.json().get("url")
        if not media_url:
            return None

        # Download media content
        dl_resp = requests.get(
            media_url,
            headers={"Authorization": f"Bearer {WHATSAPP_TOKEN}"},
            timeout=60,
        )
        if dl_resp.status_code != 200:
            logger.error(f"Media download failed: {dl_resp.status_code}")
            return None
        return dl_resp.content
    except Exception as e:
        logger.error(f"Media download error: {e}")
        return None


def wa_mark_read(message_id: str):
    """Mark a message as read."""
    try:
        requests.post(
            f"{WA_API_BASE}/messages",
            headers=WA_HEADERS,
            json={
                "messaging_product": "whatsapp",
                "status": "read",
                "message_id": message_id,
            },
            timeout=10,
        )
    except Exception:
        pass


# --------------------------------------------------------------------------- #
# Tool definitions (OpenAI function-calling format)
# --------------------------------------------------------------------------- #

TOOLS = [
    {
        "type": "function",
        "function": {
            "name": "clawbio",
            "description": (
                "Run a ClawBio bioinformatics skill. Available skills: "
                "pharmgx (pharmacogenomics report from 23andMe/AncestryDNA data), "
                "equity (HEIM equity score from VCF or ancestry CSV), "
                "nutrigx (nutrigenomics dietary advice from genetic data), "
                "metagenomics (metagenomic profiling from FASTQ), "
                "compare (genome comparison: IBS vs George Church + ancestry estimation), "
                "drugphoto (identify a drug from a photo and get personalised dosage guidance "
                "using demo genotype data -- always use mode='demo'), "
                "prs (polygenic risk scores from GWAS -- disease risk: T2D, atrial fibrillation, CAD, etc.), "
                "clinpgx (gene-drug interaction database lookup via PharmGKB/CPIC), "
                "gwas (federated variant lookup across 9 genomic databases by rsID), "
                "profile (unified genomic profile report combining all skill results). "
                "Use mode='demo' to run with built-in demo data. "
                "Use mode='file' when the user has sent a genetic data file. "
                "Use skill='auto' to let the orchestrator detect the right skill. "
                "IMPORTANT: When this tool returns results, relay the output VERBATIM. "
                "Do not paraphrase, summarise, or rewrite. The output contains exact numerical "
                "results (IBS scores, percentages, gene-drug interactions) that must be shown unchanged."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "skill": {
                        "type": "string",
                        "enum": ["pharmgx", "equity", "nutrigx", "metagenomics",
                                 "compare", "drugphoto", "prs", "clinpgx",
                                 "gwas", "profile", "auto"],
                        "description": (
                            "Which bioinformatics skill to run. Use 'auto' to let "
                            "the orchestrator detect from the file type or query."
                        ),
                    },
                    "mode": {
                        "type": "string",
                        "enum": ["file", "demo"],
                        "description": (
                            "file: use a file the user sent via WhatsApp. "
                            "demo: run with built-in demo/synthetic data."
                        ),
                    },
                    "query": {
                        "type": "string",
                        "description": (
                            "Natural language query for auto-routing via the "
                            "orchestrator (only used when skill='auto' and no file)."
                        ),
                    },
                    "extra_args": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": (
                            "Additional CLI arguments for power users "
                            "(e.g. ['--weights', '0.4,0.3,0.15,0.15'])."
                        ),
                    },
                    "drug_name": {
                        "type": "string",
                        "description": (
                            "Drug name identified from a photo (brand or generic, "
                            "e.g. 'Plavix' or 'clopidogrel'). Required when skill='drugphoto'."
                        ),
                    },
                    "visible_dose": {
                        "type": "string",
                        "description": (
                            "Dosage visible on the packaging (e.g. '50mg', '75mg'). "
                            "Optional -- enriches the recommendation."
                        ),
                    },
                    "trait": {
                        "type": "string",
                        "description": (
                            "Disease/trait to assess risk for (e.g. 'type 2 diabetes', "
                            "'atrial fibrillation'). Used with prs skill."
                        ),
                    },
                    "gene": {
                        "type": "string",
                        "description": (
                            "Gene symbol for clinpgx lookup (e.g. 'CYP2D6'). "
                            "Used with clinpgx skill."
                        ),
                    },
                    "rsid": {
                        "type": "string",
                        "description": (
                            "rsID for GWAS variant lookup (e.g. 'rs3798220'). "
                            "Used with gwas skill."
                        ),
                    },
                },
                "required": ["skill", "mode"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "save_file",
            "description": (
                "Save a file that was sent via WhatsApp to a specific folder. "
                "The file is temporarily stored after download; use this tool to "
                "move it to the requested destination. Only works for the most "
                "recently received file. Default: saves to ClawBio data/ directory."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "destination_folder": {
                        "type": "string",
                        "description": (
                            "The folder path to save the file in (absolute path). "
                            "Default: ClawBio data/ directory."
                        ),
                    },
                    "filename": {
                        "type": "string",
                        "description": (
                            "Optional filename to save as. If not provided, uses "
                            "the original filename from WhatsApp."
                        ),
                    },
                },
                "required": ["destination_folder"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "write_file",
            "description": (
                "Create or overwrite a file on the filesystem with the given content. "
                "Use this to write reports, markdown documents, text files, etc. "
                "Default destination: ClawBio data/ directory."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "content": {
                        "type": "string",
                        "description": "The full text content to write to the file.",
                    },
                    "filename": {
                        "type": "string",
                        "description": "Filename including extension (e.g. 'report.md', 'notes.txt').",
                    },
                    "destination_folder": {
                        "type": "string",
                        "description": (
                            "Folder path (absolute). Default: ClawBio data/ directory."
                        ),
                    },
                },
                "required": ["content", "filename"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "generate_audio",
            "description": (
                "Generate an MP3 audio file from text using edge-tts (Microsoft Edge "
                "text-to-speech). Good for converting reports into accessible audio. "
                "Available voices: en-GB-RyanNeural (British male, default), "
                "en-GB-SoniaNeural (British female), en-US-GuyNeural (American male). "
                "Typical speed: ~150 words/minute."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "text": {
                        "type": "string",
                        "description": "The text to convert to speech.",
                    },
                    "filename": {
                        "type": "string",
                        "description": "Output MP3 filename (e.g. 'report-audio.mp3').",
                    },
                    "voice": {
                        "type": "string",
                        "description": "TTS voice. Default: 'en-GB-RyanNeural'.",
                    },
                    "rate": {
                        "type": "string",
                        "description": "Speech rate adjustment (e.g. '-5%', '+10%'). Default: '-5%'.",
                    },
                    "destination_folder": {
                        "type": "string",
                        "description": "Folder to save the MP3 (absolute path). Default: ClawBio data/.",
                    },
                },
                "required": ["text", "filename"],
            },
        },
    },
]


# --------------------------------------------------------------------------- #
# Security helpers
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
# execute_clawbio
# --------------------------------------------------------------------------- #


async def execute_clawbio(args: dict, phone: str) -> str:
    """Execute a ClawBio bioinformatics skill via subprocess."""
    skill_key = args.get("skill", "auto")
    mode = args.get("mode", "demo")
    query = args.get("query", "")

    # Auto-routing via orchestrator
    if skill_key == "auto":
        orch_script = CLAWBIO_DIR / "skills" / "bio-orchestrator" / "orchestrator.py"
        if not orch_script.exists():
            return "Error: bio-orchestrator not found."

        orch_input = query
        if mode == "file":
            info = scoped_get(_received_files, phone)
            if info:
                orch_input = info["path"]
        if not orch_input:
            return "Error: skill='auto' requires either a file or a query to route."

        try:
            proc = await asyncio.create_subprocess_exec(
                sys.executable, str(orch_script),
                "--input", orch_input,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=str(orch_script.parent),
            )
            stdout, stderr = await asyncio.wait_for(proc.communicate(), timeout=30)
            if proc.returncode != 0:
                return f"Orchestrator error: {stderr.decode()[-500:]}"
            routing = json.loads(stdout.decode())
            detected = routing.get("detected_skill", "")
            orch_to_key = {
                "pharmgx-reporter": "pharmgx",
                "equity-scorer": "equity",
                "nutrigx": "nutrigx",
                "claw-metagenomics": "metagenomics",
                "genome-compare": "compare",
                "gwas-prs": "prs",
                "clinpgx": "clinpgx",
                "gwas-lookup": "gwas",
                "profile-report": "profile",
            }
            skill_key = orch_to_key.get(detected, "")
            if not skill_key:
                avail = list(orch_to_key.values())
                return (
                    f"Orchestrator detected skill '{detected}' which is not "
                    f"available via WhatsApp. Available: {avail}"
                )
            logger.info(f"Auto-routed to: {skill_key} (via {routing.get('detection_method', '?')})")
        except asyncio.TimeoutError:
            return "Error: orchestrator timed out."
        except json.JSONDecodeError:
            return "Error: could not parse orchestrator output."
        except Exception as e:
            return f"Error running orchestrator: {e}"

    # Resolve input and profile for file mode
    input_path = None
    profile_path = None
    info = scoped_get(_received_files, phone)
    if info:
        input_path = info.get("path")
        profile_path = info.get("profile_path")

    if mode == "file" and not input_path and not profile_path:
        if OWNER_GENOME.exists():
            input_path = str(OWNER_GENOME)
            logger.info(f"No file uploaded — using owner genome: {OWNER_GENOME.name}")
        else:
            return "Error: no file received. Send a genetic data file first, then run the skill."

    # Build output directory
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = OUTPUT_DIR / f"{skill_key}_{ts}"

    # Build command
    cmd = [sys.executable, str(CLAWBIO_PY), "run", skill_key]

    if skill_key == "profile":
        if mode == "demo":
            cmd.append("--demo")
        elif profile_path:
            cmd.extend(["--profile", profile_path])
        else:
            return "Error: no profile available. Send a genetic data file first to create a profile."
    elif skill_key == "prs":
        if mode == "demo":
            cmd.append("--demo")
        elif profile_path:
            cmd.extend(["--profile", profile_path])
        elif input_path:
            cmd.extend(["--input", str(input_path)])
        trait = args.get("trait", "")
        if trait:
            cmd.extend(["--trait", trait])
    elif skill_key == "clinpgx":
        if mode == "demo":
            cmd.append("--demo")
        else:
            gene = args.get("gene", "")
            if gene:
                cmd.extend(["--gene", gene])
            else:
                cmd.append("--demo")
    elif skill_key == "gwas":
        if mode == "demo":
            cmd.append("--demo")
        else:
            rsid = args.get("rsid", "")
            if rsid:
                cmd.extend(["--rsid", rsid])
            else:
                cmd.append("--demo")
    elif mode == "demo":
        cmd.append("--demo")
    elif input_path:
        cmd.extend(["--input", str(input_path)])

    if skill_key not in ("compare", "drugphoto"):
        cmd.extend(["--output", str(out_dir)])

    if skill_key == "drugphoto":
        drug_name = args.get("drug_name", "")
        visible_dose = args.get("visible_dose", "")
        if drug_name:
            cmd.extend(["--drug", drug_name])
        if visible_dose:
            cmd.extend(["--dose", visible_dose])

    try:
        proc = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        stdout_bytes, stderr_bytes = await asyncio.wait_for(
            proc.communicate(), timeout=120,
        )
        stdout_str = stdout_bytes.decode(errors="replace")
        stderr_str = stderr_bytes.decode(errors="replace")
    except asyncio.TimeoutError:
        return f"{skill_key} timed out after 120 seconds."
    except Exception as e:
        import traceback as _tb
        return f"{skill_key} crashed:\n{_tb.format_exc()[-1500:]}"

    if proc.returncode != 0:
        err = stderr_str[-1500:] if stderr_str else stdout_str[-1500:] if stdout_str else "unknown error"
        return f"{skill_key} failed (exit {proc.returncode}):\n{err}"

    if skill_key in ("compare", "drugphoto", "profile"):
        raw_output = stdout_str.strip()
        if raw_output:
            _pending_text.setdefault(phone, []).append(raw_output)
        return "Result sent directly to chat. Do not repeat or paraphrase it."

    if out_dir.exists():
        media_items = []
        for f in sorted(out_dir.rglob("*")):
            if not f.is_file():
                continue
            if f.suffix in (".md", ".html"):
                media_items.append({"type": "document", "path": str(f)})
            elif f.suffix == ".png":
                media_items.append({"type": "photo", "path": str(f)})
        if media_items:
            _pending_media["__current__"] = _pending_media.get("__current__", []) + media_items

    report_text = ""
    if out_dir.exists():
        for pattern in ["report.md", "*_report.md", "*.md"]:
            for md_file in sorted(out_dir.glob(pattern)):
                if md_file.name.startswith("."):
                    continue
                report_text = md_file.read_text(encoding="utf-8")
                break
            if report_text:
                break

    if not report_text:
        return stdout_str if stdout_str else f"{skill_key} completed. Output: {out_dir}"

    keep_lines = []
    skip = False
    for line in report_text.split("\n"):
        if line.startswith("## Chromosome Breakdown"):
            skip = True
        elif line.startswith("## Ancestry Composition"):
            skip = False
        elif line.startswith("## Methods"):
            skip = True
        elif line.startswith("## About"):
            skip = False
        elif line.startswith("## Reproducibility"):
            skip = True
        elif line.startswith("## Disclaimer"):
            skip = False
        if line.startswith("!["):
            continue
        if not skip:
            keep_lines.append(line)

    return "\n".join(keep_lines).strip()


# --------------------------------------------------------------------------- #
# execute_save_file
# --------------------------------------------------------------------------- #


async def execute_save_file(args: dict, phone: str) -> str:
    """Save the most recently received file to the requested destination."""
    file_info = scoped_get(_received_files, phone)

    if not file_info:
        return "No recently received file to save. Send a file first."

    src_path = Path(file_info["path"])
    if not src_path.exists():
        return "The temporary file has expired. Please send it again."

    try:
        final_path = safe_write_path(
            args.get("destination_folder"),
            args.get("filename") or file_info["filename"],
            root=DATA_DIR, require_allowed_extension=True,
        )
    except UnsafePath as exc:
        return f"Error: {exc}"

    shutil.copy2(str(src_path), str(final_path))
    logger.info(f"Saved file: {final_path}")

    try:
        src_path.unlink()
    except OSError:
        pass

    return f"File saved to {final_path}"


# --------------------------------------------------------------------------- #
# execute_write_file
# --------------------------------------------------------------------------- #


async def execute_write_file(args: dict, phone: str) -> str:
    """Create or overwrite a file with the given content."""
    content = args.get("content")
    filename = args.get("filename")
    if not content:
        return "Error: 'content' is required. Provide the full text to write."
    if not filename:
        return "Error: 'filename' is required (e.g. 'report.md')."

    try:
        filepath = safe_write_path(
            args.get("destination_folder"), filename, root=DATA_DIR,
        )
    except UnsafePath as exc:
        return f"Error: {exc}"

    filepath.write_text(content, encoding="utf-8")
    logger.info(f"Wrote file: {filepath} ({len(content)} chars)")
    return f"File written to {filepath} ({len(content)} chars)"


# --------------------------------------------------------------------------- #
# execute_generate_audio
# --------------------------------------------------------------------------- #


async def execute_generate_audio(args: dict, phone: str) -> str:
    """Generate MP3 audio from text using edge-tts."""
    text = args.get("text")
    filename = args.get("filename")
    if not text:
        return "Error: 'text' is required. Provide the text to convert to speech."
    if not filename:
        return "Error: 'filename' is required (e.g. 'report.mp3')."
    if not filename.endswith(".mp3"):
        filename += ".mp3"

    voice = args.get("voice", "en-GB-RyanNeural")
    rate = args.get("rate", "-5%")
    try:
        filepath = safe_write_path(
            args.get("destination_folder"), filename, root=DATA_DIR,
        )
    except UnsafePath as exc:
        return f"Error: {exc}"
    dest = filepath.parent

    text_path = dest / f".tmp_{filename}.txt"
    text_path.write_text(text, encoding="utf-8")

    edge_tts_bin = Path.home() / "Library" / "Python" / "3.9" / "bin" / "edge-tts"
    if not edge_tts_bin.exists():
        edge_tts_bin = "edge-tts"

    try:
        proc = await asyncio.create_subprocess_exec(
            str(edge_tts_bin),
            "--voice", voice,
            f"--rate={rate}",
            "--file", str(text_path),
            "--write-media", str(filepath),
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        stdout, stderr = await asyncio.wait_for(proc.communicate(), timeout=300)

        try:
            text_path.unlink()
        except OSError:
            pass

        if proc.returncode != 0:
            err = stderr.decode()[-300:] if stderr else "unknown error"
            return f"Audio generation failed (exit {proc.returncode}): {err}"

        size_mb = filepath.stat().st_size / (1024 * 1024)
        word_count = len(text.split())
        est_minutes = word_count / 150

        logger.info(f"Generated audio: {filepath} ({size_mb:.1f} MB, ~{est_minutes:.0f} min)")
        return (
            f"Audio saved to {filepath} ({size_mb:.1f} MB, "
            f"~{word_count} words, ~{est_minutes:.0f} min estimated)"
        )

    except asyncio.TimeoutError:
        try:
            text_path.unlink()
        except OSError:
            pass
        return "Audio generation timed out after 5 minutes."
    except FileNotFoundError:
        try:
            text_path.unlink()
        except OSError:
            pass
        return "edge-tts not found. Install with: pip3 install edge-tts"


# --------------------------------------------------------------------------- #
# LLM tool loop
# --------------------------------------------------------------------------- #

TOOL_EXECUTORS = {
    "clawbio": execute_clawbio,
    "save_file": execute_save_file,
    "write_file": execute_write_file,
    "generate_audio": execute_generate_audio,
}

MAX_TOOL_ITERATIONS = 10


async def llm_tool_loop(phone: str, user_content: str | list) -> str:
    """Run the LLM tool-use loop."""
    history = conversations.setdefault(phone, [])

    if isinstance(user_content, str):
        history.append({"role": "user", "content": user_content})
    else:
        oai_parts = []
        for block in user_content:
            if block.get("type") == "text":
                oai_parts.append({"type": "text", "text": block["text"]})
            elif block.get("type") == "image":
                src = block.get("source", {})
                data_uri = f"data:{src['media_type']};base64,{src['data']}"
                oai_parts.append({
                    "type": "image_url",
                    "image_url": {"url": data_uri},
                })
        history.append({"role": "user", "content": oai_parts})

    if len(history) > MAX_HISTORY:
        history[:] = history[-MAX_HISTORY:]

    # Sanitise orphaned tool messages
    sanitised: list[dict] = []
    for msg in history:
        if msg.get("role") == "tool":
            if sanitised and sanitised[-1].get("role") == "assistant":
                if sanitised[-1].get("tool_calls"):
                    sanitised.append(msg)
                    continue
            logger.warning("Dropped orphaned tool message from history")
            _audit("history_sanitised", phone=phone,
                   detail="orphaned_tool_message_dropped")
            continue
        sanitised.append(msg)
    history[:] = sanitised

    last_message = None
    for _iteration in range(MAX_TOOL_ITERATIONS):
        try:
            response = await llm.chat.completions.create(
                model=CLAWBIO_MODEL,
                max_tokens=8192,
                messages=[{"role": "system", "content": SYSTEM_PROMPT}] + history,
                tools=TOOLS,
            )
        except APIError as e:
            logger.error(f"LLM API error: {e}")
            return f"Sorry, I'm having trouble thinking right now -- API error: {e}"

        choice = response.choices[0]
        last_message = choice.message

        assistant_msg = {"role": "assistant", "content": last_message.content or ""}
        if last_message.tool_calls:
            assistant_msg["tool_calls"] = [
                {
                    "id": tc.id,
                    "type": "function",
                    "function": {
                        "name": tc.function.name,
                        "arguments": tc.function.arguments,
                    },
                }
                for tc in last_message.tool_calls
            ]
        history.append(assistant_msg)

        if not last_message.tool_calls:
            return last_message.content or "(no response)"

        for tc in last_message.tool_calls:
            func_name = tc.function.name
            executor = TOOL_EXECUTORS.get(func_name)
            if executor:
                try:
                    args = json.loads(tc.function.arguments)
                except json.JSONDecodeError:
                    args = {}
                logger.info(f"Tool call: {func_name}({json.dumps(args)[:200]})")
                _audit("tool_call", phone=phone, tool=func_name,
                       args_preview=json.dumps(args, default=str)[:300])
                try:
                    result = await executor(args, phone)
                except Exception as tool_err:
                    logger.error(f"Tool {func_name} raised: {tool_err}", exc_info=True)
                    _audit("tool_error", phone=phone, tool=func_name,
                           error=str(tool_err)[:300])
                    result = f"Error executing {func_name}: {type(tool_err).__name__}: {tool_err}"
            else:
                result = f"Unknown tool: {func_name}"

            history.append({
                "role": "tool",
                "tool_call_id": tc.id,
                "content": result,
            })

    return last_message.content if last_message and last_message.content else "(max tool iterations reached)"


# --------------------------------------------------------------------------- #
# Rate limiting
# --------------------------------------------------------------------------- #

_rate_buckets: dict[str, list[float]] = {}


def _check_rate_limit(phone: str) -> bool:
    """Return True if the user is within rate limits (or is admin)."""
    if RATE_LIMIT_PER_HOUR <= 0 or is_admin(phone):
        return True
    now = time.time()
    window = 3600
    bucket = _rate_buckets.setdefault(phone, [])
    bucket[:] = [t for t in bucket if now - t < window]
    if len(bucket) >= RATE_LIMIT_PER_HOUR:
        return False
    bucket.append(now)
    return True


# --------------------------------------------------------------------------- #
# Text helpers
# --------------------------------------------------------------------------- #

IMAGE_EXTENSIONS = {".png", ".jpg", ".jpeg", ".gif", ".webp", ".bmp"}
GENETIC_EXTENSIONS = {".txt", ".csv", ".vcf", ".fastq", ".fq", ".gz"}


def strip_markup(text: str) -> str:
    """Remove markdown/emoji formatting -- SOUL.md mandates plain text only."""
    text = re.sub(r"```.*?```", "", text, flags=re.DOTALL)
    text = re.sub(r"\*\*(.+?)\*\*", r"\1", text)
    text = re.sub(r"__(.+?)__", r"\1", text)
    text = re.sub(r"(?<!\w)\*(.+?)\*(?!\w)", r"\1", text)
    text = re.sub(r"(?<!\w)_(.+?)_(?!\w)", r"\1", text)
    text = re.sub(r"`(.+?)`", r"\1", text)
    text = re.sub(r"\[(.+?)\]\(.+?\)", r"\1", text)
    text = re.sub(r"^#{1,6}\s+", "", text, flags=re.MULTILINE)
    text = re.sub(r"^>\s?", "", text, flags=re.MULTILINE)
    text = re.sub(r"^[\s]*[-*]\s+", "", text, flags=re.MULTILINE)
    text = re.sub(
        r"[\U0001F300-\U0001F9FF\U00002702-\U000027B0\U0000FE00-\U0000FE0F"
        r"\U0001FA00-\U0001FA6F\U0001FA70-\U0001FAFF\U00002600-\U000026FF"
        r"\U0000200D\U00002B50\U00002B55\U000023CF\U000023E9-\U000023F3"
        r"\U000023F8-\U000023FA\U0000231A\U0000231B\U00003030\U000000A9"
        r"\U000000AE\U00002122\U00002139\U00002194-\U00002199"
        r"\U000021A9-\U000021AA\U0000FE0F]+",
        "",
        text,
    )
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip()


# --------------------------------------------------------------------------- #
# Voice reply helper
# --------------------------------------------------------------------------- #


async def _send_voice_reply(phone: str, text: str) -> bool:
    """Generate voice reply and send via WhatsApp as audio."""
    with tempfile.TemporaryDirectory() as tmpdir:
        text_file = os.path.join(tmpdir, "reply.txt")
        aiff_file = os.path.join(tmpdir, "reply.aiff")
        mp3_file = os.path.join(tmpdir, "reply.mp3")

        with open(text_file, "w", encoding="utf-8") as f:
            f.write(text)

        proc = await asyncio.create_subprocess_exec(
            "say", "-v", "Daniel", "-r", "170",
            "-f", text_file, "-o", aiff_file,
            stdout=asyncio.subprocess.DEVNULL,
            stderr=asyncio.subprocess.DEVNULL,
        )
        await proc.wait()
        if proc.returncode != 0:
            logger.warning("say command failed for voice reply")
            return False

        proc = await asyncio.create_subprocess_exec(
            "ffmpeg", "-y", "-i", aiff_file,
            "-codec:a", "libmp3lame", "-b:a", "128k",
            mp3_file,
            stdout=asyncio.subprocess.DEVNULL,
            stderr=asyncio.subprocess.DEVNULL,
        )
        await proc.wait()
        if proc.returncode != 0:
            logger.warning("ffmpeg MP3 conversion failed for voice reply")
            return False

        wa_send_document(phone, mp3_file, caption="Voice reply")

    return True


# --------------------------------------------------------------------------- #
# Drain pending media
# --------------------------------------------------------------------------- #


def drain_pending_media(phone: str):
    """Send any queued ClawBio media after the text reply."""
    items = _pending_media.pop("__current__", [])
    for item in items:
        try:
            path = Path(item["path"])
            if not path.exists():
                continue
            caption = path.stem.replace("_", " ").title() if item["type"] == "photo" else ""
            wa_send_document(phone, str(path), caption=caption)
        except Exception as e:
            logger.warning(f"Failed to send media {item['path']}: {e}")


# --------------------------------------------------------------------------- #
# Message handler (called from Flask webhook)
# --------------------------------------------------------------------------- #


def handle_whatsapp_message(phone: str, msg: dict):
    """Process an incoming WhatsApp message."""
    if not is_sender_allowed(phone, WHATSAPP_ALLOWED_NUMBERS, WHATSAPP_ADMIN_PHONE, WHATSAPP_PUBLIC):
        _audit("sender_rejected", phone=phone)
        return
    msg_id = msg.get("id", "")
    msg_type = msg.get("type", "")

    # Dedup
    if msg_id in _processed_messages:
        return
    _processed_messages.add(msg_id)
    if len(_processed_messages) > _processed_messages_max:
        # Trim oldest (sets are unordered, but this prevents unbounded growth)
        while len(_processed_messages) > _processed_messages_max // 2:
            _processed_messages.pop()

    wa_mark_read(msg_id)

    # Rate limit
    if not _check_rate_limit(phone):
        _audit("rate_limited", phone=phone)
        wa_send_text(phone,
                     f"You've reached the limit of {RATE_LIMIT_PER_HOUR} messages per hour. "
                     "Please try again later.")
        return

    _audit("message", phone=phone, msg_type=msg_type)

    # ----- Text messages ----- #
    if msg_type == "text":
        text = msg.get("text", {}).get("body", "").strip()
        if not text:
            return

        logger.info(f"Message from {phone}: {text[:100]}")

        # Commands
        text_lower = text.lower().strip()

        if text_lower in ("!start", "start", "hi", "hello", "help"):
            wa_send_text(phone,
                "Welcome to ClawBio -- open-source bioinformatics at your fingertips!\n\n"
                "I can analyse genetic data, check drug interactions, assess nutritional "
                "genomics, estimate polygenic risk scores, and more.\n\n"
                "Commands:\n"
                "  !skills  -- list available bioinformatics skills\n"
                "  !demo <skill>  -- run a demo (pharmgx, equity, nutrigx, compare, prs, profile)\n"
                "  !voice  -- toggle voice replies on/off\n"
                "  !status  -- bot info\n"
                "  !health  -- system health check\n\n"
                "Or just chat -- ask any bioinformatics question.\n"
                "Send a genetic data file (.txt, .csv, .vcf) to analyse it.\n"
                "Send a photo of a medication for personalised drug guidance.\n\n"
                "ClawBio is a research tool, not a medical device. "
                "Consult a healthcare professional before making medical decisions."
            )
            return

        if text_lower == "!skills":
            try:
                import subprocess
                proc = subprocess.run(
                    [sys.executable, str(CLAWBIO_PY), "list"],
                    capture_output=True, text=True, timeout=15,
                    cwd=str(CLAWBIO_DIR),
                )
                wa_send_text(phone, proc.stdout.strip() or "No skills found.")
            except Exception as e:
                wa_send_text(phone, f"Error listing skills: {e}")
            return

        if text_lower == "!voice":
            current = _voice_enabled.get(phone, False)
            _voice_enabled[phone] = not current
            state = "ON" if not current else "OFF"
            wa_send_text(phone,
                f"Voice replies toggled {state}.\n"
                f"{'I will now send voice memos alongside text replies.' if not current else 'Back to text-only replies.'}"
            )
            return

        if text_lower == "!status":
            uptime_secs = int(time.time() - BOT_START_TIME)
            hours, remainder = divmod(uptime_secs, 3600)
            minutes, secs = divmod(remainder, 60)
            uptime_str = f"{hours}h {minutes}m {secs}s"

            skills_dir = CLAWBIO_DIR / "skills"
            skill_count = sum(
                1 for d in skills_dir.iterdir()
                if d.is_dir() and (d / "SKILL.md").exists()
            ) if skills_dir.exists() else 0

            wa_send_text(phone,
                f"RoboTerri ClawBio Status\n"
                f"========================\n"
                f"Bot uptime: {uptime_str}\n"
                f"LLM model: {CLAWBIO_MODEL}\n"
                f"Skills available: {skill_count}\n"
                f"Platform: WhatsApp\n"
            )
            return

        if text_lower == "!health":
            checks = []
            if CLAWBIO_PY.exists():
                checks.append("ClawBio CLI: OK")
            else:
                checks.append("ClawBio CLI: MISSING")
            if SOUL_MD.exists():
                checks.append(f"SOUL.md: OK ({len(_soul)} chars)")
            else:
                checks.append("SOUL.md: MISSING (using fallback)")
            skills_dir = CLAWBIO_DIR / "skills"
            if skills_dir.exists():
                implemented = sum(1 for d in skills_dir.iterdir()
                                  if d.is_dir() and (d / "SKILL.md").exists() and any(d.glob("*.py")))
                stub_only = sum(1 for d in skills_dir.iterdir()
                                if d.is_dir() and (d / "SKILL.md").exists() and not any(d.glob("*.py")))
                checks.append(f"Skills (implemented): {implemented}")
                checks.append(f"Skills (stub/planned): {stub_only}")
            if OUTPUT_DIR.exists():
                checks.append(f"Output runs: {sum(1 for _ in OUTPUT_DIR.iterdir())}")
            wa_send_text(phone,
                "ClawBio Health Check\n"
                "====================\n" + "\n".join(checks)
            )
            return

        if text_lower.startswith("!demo"):
            parts = text.split(maxsplit=1)
            skill = parts[1].strip() if len(parts) > 1 else "pharmgx"
            wa_send_text(phone, f"Running {skill} demo -- this may take a moment...")
            try:
                reply = run_async(llm_tool_loop(phone,
                    f"Run the {skill} demo using the clawbio tool with mode='demo'."))
                _pt = _pending_text.pop(phone, None)
                if _pt:
                    reply = "\n\n".join(_pt)
                wa_send_text(phone, reply)
                drain_pending_media(phone)
                if _voice_enabled.get(phone):
                    try:
                        run_async(_send_voice_reply(phone, reply))
                    except Exception as ve:
                        logger.warning(f"Voice reply failed: {ve}")
            except Exception as e:
                logger.error(f"Demo error: {e}", exc_info=True)
                wa_send_text(phone, f"Demo failed: {e}")
            return

        # Regular text message -> LLM
        try:
            reply = run_async(llm_tool_loop(phone, text))
            _pt = _pending_text.pop(phone, None)
            if _pt:
                reply = "\n\n".join(_pt)
            wa_send_text(phone, reply)
            drain_pending_media(phone)
            if _voice_enabled.get(phone):
                try:
                    run_async(_send_voice_reply(phone, reply))
                except Exception as ve:
                    logger.warning(f"Voice reply failed: {ve}")
        except Exception as e:
            logger.error(f"Error handling message: {e}", exc_info=True)
            wa_send_text(phone, f"Sorry, something went wrong -- {type(e).__name__}: {e}")

    # ----- Image messages ----- #
    elif msg_type == "image":
        image_info = msg.get("image", {})
        media_id = image_info.get("id", "")
        caption = image_info.get("caption", "")
        mime_type = image_info.get("mime_type", "image/jpeg")

        logger.info(f"Image from {phone}: media_id={media_id}, mime={mime_type}")

        img_bytes = wa_download_media(media_id)
        if not img_bytes:
            wa_send_text(phone, "Sorry, I couldn't download that image. Please try again.")
            return

        if len(img_bytes) > MAX_PHOTO_BYTES:
            wa_send_text(phone,
                f"Photo too large ({len(img_bytes) / (1024*1024):.1f} MB). "
                f"Maximum: {MAX_PHOTO_BYTES / (1024*1024):.0f} MB.")
            return

        _audit("photo", phone=phone, size_bytes=len(img_bytes), media_type=mime_type)

        img_b64 = base64.standard_b64encode(img_bytes).decode("ascii")

        filename = sanitize_filename(f"whatsapp_image_{datetime.now().strftime('%H%M%S')}.jpg")
        tmp_path = upload_tmp_path(phone, filename)
        tmp_path.write_bytes(img_bytes)
        _received_files[phone] = {"path": str(tmp_path), "filename": filename}

        content_blocks = [
            {
                "type": "image",
                "source": {
                    "type": "base64",
                    "media_type": mime_type,
                    "data": img_b64,
                },
            },
        ]
        if caption:
            content_blocks.append({"type": "text", "text": caption})
        else:
            content_blocks.append({
                "type": "text",
                "text": (
                    "[Image sent without caption. Look at this image. "
                    "If it shows a medication, drug packaging, pill bottle, blister pack, or "
                    "any pharmaceutical product: immediately identify the drug name and any "
                    "visible dosage, then call the clawbio tool with skill='drugphoto', "
                    "mode='demo', drug_name=<identified drug>, and visible_dose=<dose if readable>. "
                    "Do NOT ask what is needed -- just run the lookup automatically. "
                    "If the image is not a medication, describe what you see and ask if "
                    "anything specific is needed.]"
                ),
            })

        try:
            reply = run_async(llm_tool_loop(phone, content_blocks))
            _pt = _pending_text.pop(phone, None)
            if _pt:
                reply = "\n\n".join(_pt)
            wa_send_text(phone, reply)
            if _voice_enabled.get(phone):
                try:
                    run_async(_send_voice_reply(phone, reply))
                except Exception as ve:
                    logger.warning(f"Voice reply failed: {ve}")
        except Exception as e:
            logger.error(f"Photo handling error: {e}", exc_info=True)
            wa_send_text(phone, f"Sorry, I couldn't process that image -- {type(e).__name__}: {e}")

    # ----- Document messages ----- #
    elif msg_type == "document":
        doc_info = msg.get("document", {})
        media_id = doc_info.get("id", "")
        filename = sanitize_filename(doc_info.get("filename", "document"))
        mime_type = doc_info.get("mime_type", "application/octet-stream")
        caption = doc_info.get("caption", "")

        logger.info(f"Document from {phone}: {filename}, mime={mime_type}")

        # Check if it's an image sent as document
        ext = Path(filename).suffix.lower()
        if mime_type.startswith("image/") or ext in IMAGE_EXTENSIONS:
            # Treat as image
            img_bytes = wa_download_media(media_id)
            if not img_bytes:
                wa_send_text(phone, "Sorry, I couldn't download that file. Please try again.")
                return
            if len(img_bytes) > MAX_PHOTO_BYTES:
                wa_send_text(phone,
                    f"Photo too large ({len(img_bytes) / (1024*1024):.1f} MB). "
                    f"Maximum: {MAX_PHOTO_BYTES / (1024*1024):.0f} MB.")
                return

            img_b64 = base64.standard_b64encode(img_bytes).decode("ascii")
            tmp_path = upload_tmp_path(phone, filename)
            tmp_path.write_bytes(img_bytes)
            _received_files[phone] = {"path": str(tmp_path), "filename": filename}

            content_blocks = [
                {"type": "image", "source": {"type": "base64", "media_type": mime_type, "data": img_b64}},
                {"type": "text", "text": caption or "[Image sent as document. Identify any medication and run drugphoto if applicable.]"},
            ]
            try:
                reply = run_async(llm_tool_loop(phone, content_blocks))
                _pt = _pending_text.pop(phone, None)
                if _pt:
                    reply = "\n\n".join(_pt)
                wa_send_text(phone, reply)
            except Exception as e:
                wa_send_text(phone, f"Sorry, couldn't process that image -- {e}")
            return

        # Genetic data file
        if ext not in GENETIC_EXTENSIONS:
            wa_send_text(phone,
                f"I received '{filename}' but I can only process genetic data files "
                "(.txt, .csv, .vcf, .fastq, .fq, .gz). Please send the correct file type.")
            return

        file_bytes = wa_download_media(media_id)
        if not file_bytes:
            wa_send_text(phone, "Sorry, I couldn't download that file. Please try again.")
            return

        if len(file_bytes) > MAX_UPLOAD_BYTES:
            wa_send_text(phone,
                f"File too large ({len(file_bytes) / (1024*1024):.1f} MB). "
                f"Maximum: {MAX_UPLOAD_BYTES / (1024*1024):.0f} MB.")
            return

        _audit("document", phone=phone, filename=filename, size_bytes=len(file_bytes))

        tmp_path = upload_tmp_path(phone, filename)
        tmp_path.write_bytes(file_bytes)
        _received_files[phone] = {"path": str(tmp_path), "filename": filename}

        # Auto-create patient profile
        profile_path = None
        try:
            import subprocess
            upload_proc = subprocess.run(
                [sys.executable, str(CLAWBIO_PY), "upload", "--input", str(tmp_path)],
                capture_output=True, text=True, timeout=30,
            )
            for line in upload_proc.stdout.splitlines():
                if "profile" in line.lower() and ("/" in line or "\\" in line):
                    for token in line.split():
                        if token.endswith(".json"):
                            profile_path = token
                            break
            if profile_path:
                _received_files[phone]["profile_path"] = profile_path
                logger.info(f"Auto-created profile: {profile_path}")
        except Exception as prof_err:
            logger.warning(f"Auto-profile creation failed (non-fatal): {prof_err}")

        parts = [f"[Document received: {filename} ({len(file_bytes)} bytes)]"]
        if profile_path:
            parts.append(f"[Patient profile auto-created: {profile_path}]")
        if caption:
            parts.append(caption)
        else:
            profile_note = (
                " A patient profile has been created -- the user can now ask "
                "follow-up questions like 'what am I at risk for?' (prs) or "
                "'show my full profile' (profile) without re-uploading."
            ) if profile_path else ""
            parts.append(
                "The user sent this genetic data file. Detect the file type and "
                "run the appropriate ClawBio skill using mode='file'. For .txt "
                "files (23andMe format) use pharmgx. For .csv (AncestryDNA) use "
                "pharmgx. For .vcf use equity. For .fastq use metagenomics. "
                "If unsure, use skill='auto'." + profile_note
            )

        try:
            reply = run_async(llm_tool_loop(phone, "\n\n".join(parts)))
            _pt = _pending_text.pop(phone, None)
            if _pt:
                reply = "\n\n".join(_pt)
            wa_send_text(phone, reply)
            drain_pending_media(phone)
            if _voice_enabled.get(phone):
                try:
                    run_async(_send_voice_reply(phone, reply))
                except Exception as ve:
                    logger.warning(f"Voice reply failed: {ve}")
        except Exception as e:
            logger.error(f"Document handling error: {e}", exc_info=True)
            wa_send_text(phone, f"Sorry, couldn't process that document -- {e}")

    else:
        logger.info(f"Unsupported message type from {phone}: {msg_type}")


# --------------------------------------------------------------------------- #
# Flask webhook
# --------------------------------------------------------------------------- #

app = Flask(__name__)


@app.route("/webhook", methods=["GET"])
def webhook_verify():
    """Handle WhatsApp webhook verification (GET)."""
    mode = request.args.get("hub.mode", "")
    token = request.args.get("hub.verify_token", "")
    challenge = request.args.get("hub.challenge", "")

    if mode == "subscribe" and token == WHATSAPP_VERIFY_TOKEN:
        logger.info("Webhook verified successfully")
        return challenge, 200
    else:
        logger.warning(f"Webhook verification failed: mode={mode}, token={token[:10]}...")
        return "Forbidden", 403


@app.route("/webhook", methods=["POST"])
def webhook_receive():
    """Handle incoming WhatsApp messages (POST)."""
    raw_body = request.get_data()
    signature = request.headers.get("X-Hub-Signature-256")
    if not verify_whatsapp_signature(raw_body, signature, WHATSAPP_APP_SECRET):
        _audit("webhook_signature_rejected", sig_present=bool(signature))
        return "forbidden", 403
    data = request.get_json(silent=True)
    if not data:
        return "OK", 200

    try:
        for entry in data.get("entry", []):
            for change in entry.get("changes", []):
                value = change.get("value", {})
                messages = value.get("messages", [])
                for msg in messages:
                    phone = msg.get("from", "")
                    if phone:
                        # Process in a thread to avoid blocking the webhook
                        threading.Thread(
                            target=handle_whatsapp_message,
                            args=(phone, msg),
                            daemon=True,
                        ).start()
    except Exception as e:
        logger.error(f"Webhook processing error: {e}", exc_info=True)

    # Always return 200 quickly to avoid WhatsApp retries
    return "OK", 200


@app.route("/health", methods=["GET"])
def health_check():
    """Health check endpoint."""
    uptime_secs = int(time.time() - BOT_START_TIME)
    return jsonify({
        "status": "ok",
        "uptime_seconds": uptime_secs,
        "model": CLAWBIO_MODEL,
        "platform": "whatsapp",
    })


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #


def main():
    """Start the bot."""
    logger.info(f"Starting RoboTerri WhatsApp bot (model: {CLAWBIO_MODEL})")
    logger.info(f"ClawBio directory: {CLAWBIO_DIR}")
    logger.info(f"Webhook port: {WHATSAPP_PORT}")
    if LLM_BASE_URL:
        logger.info(f"LLM base URL: {LLM_BASE_URL}")
    logger.info(f"Admin phone: {WHATSAPP_ADMIN_PHONE or 'not set (public mode)'}")
    logger.info(f"Rate limit: {RATE_LIMIT_PER_HOUR} msgs/hour per user (0=unlimited)")
    _audit("bot_start", model=CLAWBIO_MODEL,
           admin_phone=WHATSAPP_ADMIN_PHONE or None, rate_limit=RATE_LIMIT_PER_HOUR)

    print(f"RoboTerri WhatsApp bot is running on port {WHATSAPP_PORT}. Press Ctrl+C to stop.")
    print(f"Webhook URL: http://localhost:{WHATSAPP_PORT}/webhook")
    print("Use ngrok or a reverse proxy to expose this endpoint to the internet.")

    app.run(host="0.0.0.0", port=WHATSAPP_PORT, debug=False)


if __name__ == "__main__":
    main()
