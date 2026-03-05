#!/usr/bin/env python3
"""
roboterri.py — RoboTerri ClawBio Telegram Bot
==============================================
A Telegram bot that runs ClawBio bioinformatics skills using any LLM
as the reasoning engine. Handles text messages, genetic file uploads,
and medication photos.

Works with any OpenAI-compatible provider: OpenAI, Anthropic (via proxy),
Google, Mistral, Groq, Together, OpenRouter, Ollama, LM Studio, etc.

Prerequisites:
    pip3 install python-telegram-bot[job-queue] openai python-dotenv

Usage:
    # Set environment variables in .env (see bot/README.md)
    python3 bot/roboterri.py
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
import time
from datetime import datetime
from pathlib import Path

from dotenv import load_dotenv
from openai import AsyncOpenAI, APIError
from telegram import Update
from telegram.ext import (
    Application,
    CommandHandler,
    ContextTypes,
    MessageHandler,
    filters,
)

# --------------------------------------------------------------------------- #
# Config
# --------------------------------------------------------------------------- #

load_dotenv()

TELEGRAM_BOT_TOKEN = os.environ.get("TELEGRAM_BOT_TOKEN", "")
TELEGRAM_CHAT_ID = int(os.environ.get("TELEGRAM_CHAT_ID", os.environ.get("AUTHORISED_CHAT_ID", "0")))
LLM_API_KEY = os.environ.get("LLM_API_KEY", os.environ.get("OPENAI_API_KEY", ""))
LLM_BASE_URL = os.environ.get("LLM_BASE_URL", "")
CLAWBIO_MODEL = os.environ.get("CLAWBIO_MODEL", "gpt-4o")

if not TELEGRAM_BOT_TOKEN:
    print("Error: TELEGRAM_BOT_TOKEN not set. See bot/README.md for setup.")
    sys.exit(1)
if not LLM_API_KEY:
    print("Error: LLM_API_KEY not set. See bot/README.md for setup.")
    sys.exit(1)
if not TELEGRAM_CHAT_ID:
    print("Error: TELEGRAM_CHAT_ID not set. See bot/README.md for setup.")
    sys.exit(1)

CLAWBIO_DIR = Path(__file__).resolve().parent.parent
CLAWBIO_PY = CLAWBIO_DIR / "clawbio.py"
SOUL_MD = CLAWBIO_DIR / "SOUL.md"
OUTPUT_DIR = CLAWBIO_DIR / "output"
DATA_DIR = CLAWBIO_DIR / "data"

# Security limits (TG-004)
MAX_UPLOAD_BYTES = 50 * 1024 * 1024  # 50 MB (Telegram document limit)
MAX_PHOTO_BYTES = 20 * 1024 * 1024   # 20 MB (Telegram photo limit)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
)
logger = logging.getLogger("roboterri")

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
"""

SYSTEM_PROMPT = f"{_soul}\n\n{ROLE_GUARDRAILS}"

# --------------------------------------------------------------------------- #
# State
# --------------------------------------------------------------------------- #

_client_kwargs = {"api_key": LLM_API_KEY}
if LLM_BASE_URL:
    _client_kwargs["base_url"] = LLM_BASE_URL
llm = AsyncOpenAI(**_client_kwargs)

conversations: dict[int, list] = {}
MAX_HISTORY = 20

# Per-chat received file storage
_received_files: dict[int, dict] = {}

# Pending media queue: chat_id -> list of {"type": "document"|"photo", "path": str}
_pending_media: dict[int, list[dict]] = {}

# Pending text queue: bypass LLM paraphrasing for compare/drugphoto
_pending_text: list[str] = []

BOT_START_TIME = time.time()

# --------------------------------------------------------------------------- #
# Tool definition (OpenAI function-calling format)
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
                            "file: use a file the user sent via Telegram. "
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
                "Save a file that was sent via Telegram to a specific folder. "
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
                            "the original filename from Telegram."
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
# execute_clawbio
# --------------------------------------------------------------------------- #


async def execute_clawbio(args: dict) -> str:
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
            for _cid, info in _received_files.items():
                orch_input = info["path"]
                break
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
                "nutrigx_advisor": "nutrigx",
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
                    f"available via Telegram. Available: {avail}"
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
    for _cid, info in _received_files.items():
        input_path = info.get("path")
        profile_path = info.get("profile_path")
        break

    if mode == "file" and not input_path and not profile_path:
        return "Error: no file received. Send a genetic data file first, then run the skill."

    # Build output directory
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = OUTPUT_DIR / f"{skill_key}_{ts}"

    # Build command
    cmd = [sys.executable, str(CLAWBIO_PY), "run", skill_key]

    # Profile-based skills: prefer --profile over --input
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

    # Skills with summary_default (compare, drugphoto) skip --output
    if skill_key not in ("compare", "drugphoto"):
        cmd.extend(["--output", str(out_dir)])

    # Pass drug_name and visible_dose for drugphoto
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

    # For compare / drugphoto / profile: send stdout directly (bypass LLM paraphrasing)
    if skill_key in ("compare", "drugphoto", "profile"):
        raw_output = stdout_str.strip()
        if raw_output:
            _pending_text.append(raw_output)
        return "Result sent directly to chat. Do not repeat or paraphrase it."

    # For other skills: collect report + figures from output directory
    output_files = sorted([f.name for f in out_dir.rglob("*") if f.is_file()]) if out_dir.exists() else []

    # Queue figures and reports for Telegram delivery
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
            _pending_media[0] = _pending_media.get(0, []) + media_items

    # Read report for chat display
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

    # Extract key sections (drop chromosome table, methods, reproducibility, disclaimer)
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
        elif line.startswith("## Disclaimer"):
            skip = True
        elif line.startswith("## Reproducibility"):
            skip = True
        if line.startswith("!["):
            continue
        if not skip:
            keep_lines.append(line)

    return "\n".join(keep_lines).strip()


# --------------------------------------------------------------------------- #
# Security helpers (TG-002)
# --------------------------------------------------------------------------- #


def _sanitize_filename(filename: str) -> str:
    """Strip path traversal components and dangerous characters from a filename."""
    # Take only the basename (no directory components)
    filename = Path(filename).name
    # Remove null bytes and control characters
    filename = re.sub(r"[\x00-\x1f]", "", filename)
    # Collapse path traversal attempts
    filename = filename.replace("..", "").replace("/", "").replace("\\", "")
    if not filename:
        filename = "unnamed_file"
    return filename


def _resolve_dest(folder: str | None) -> Path:
    """Resolve a destination folder, defaulting to ClawBio data/."""
    dest = Path(folder) if folder else DATA_DIR
    if not dest.is_absolute():
        dest = CLAWBIO_DIR / dest
    dest.mkdir(parents=True, exist_ok=True)
    return dest


def _validate_path(filepath: Path, allowed_root: Path) -> bool:
    """Ensure filepath is under allowed_root (path traversal defense)."""
    try:
        filepath.resolve().relative_to(allowed_root.resolve())
        return True
    except ValueError:
        return False


# --------------------------------------------------------------------------- #
# execute_save_file
# --------------------------------------------------------------------------- #


async def execute_save_file(args: dict) -> str:
    """Save the most recently received file to the requested destination."""
    file_info = None
    for _cid, info in _received_files.items():
        file_info = info
        break

    if not file_info:
        return "No recently received file to save. Send a file first."

    src_path = Path(file_info["path"])
    if not src_path.exists():
        return "The temporary file has expired. Please send it again."

    dest_path = _resolve_dest(args.get("destination_folder"))
    filename = _sanitize_filename(args.get("filename") or file_info["filename"])
    final_path = dest_path / filename

    if not _validate_path(final_path, dest_path):
        return f"Error: filename '{filename}' would escape the destination directory."

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


async def execute_write_file(args: dict) -> str:
    """Create or overwrite a file with the given content."""
    content = args.get("content")
    filename = args.get("filename")
    if not content:
        return "Error: 'content' is required. Provide the full text to write."
    if not filename:
        return "Error: 'filename' is required (e.g. 'report.md')."

    dest = _resolve_dest(args.get("destination_folder"))
    filename = _sanitize_filename(filename)
    filepath = dest / filename

    if not _validate_path(filepath, dest):
        return f"Error: filename '{filename}' would escape the destination directory."

    filepath.write_text(content, encoding="utf-8")
    logger.info(f"Wrote file: {filepath} ({len(content)} chars)")
    return f"File written to {filepath} ({len(content)} chars)"


# --------------------------------------------------------------------------- #
# execute_generate_audio
# --------------------------------------------------------------------------- #


async def execute_generate_audio(args: dict) -> str:
    """Generate MP3 audio from text using edge-tts."""
    text = args.get("text")
    filename = args.get("filename")
    if not text:
        return "Error: 'text' is required. Provide the text to convert to speech."
    if not filename:
        return "Error: 'filename' is required (e.g. 'report.mp3')."
    if not filename.endswith(".mp3"):
        filename += ".mp3"

    filename = _sanitize_filename(filename)
    voice = args.get("voice", "en-GB-RyanNeural")
    rate = args.get("rate", "-5%")
    dest = _resolve_dest(args.get("destination_folder"))
    filepath = dest / filename

    if not _validate_path(filepath, dest):
        return f"Error: filename '{filename}' would escape the destination directory."

    # Write text to temp file for edge-tts
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
# Drain pending media
# --------------------------------------------------------------------------- #


async def _drain_pending_media(update: Update, context) -> None:
    """Send any queued ClawBio media (documents + figures) after the text reply."""
    items = _pending_media.pop(0, [])
    if not items:
        return
    chat_id = update.effective_chat.id
    for item in items:
        try:
            path = Path(item["path"])
            if not path.exists():
                continue
            if item["type"] == "document":
                await context.bot.send_document(
                    chat_id=chat_id,
                    document=open(path, "rb"),
                    filename=path.name,
                )
            elif item["type"] == "photo":
                await context.bot.send_photo(
                    chat_id=chat_id,
                    photo=open(path, "rb"),
                    caption=path.stem.replace("_", " ").title(),
                )
        except Exception as e:
            logger.warning(f"Failed to send media {item['path']}: {e}")


# --------------------------------------------------------------------------- #
# LLM tool loop (OpenAI-compatible chat completions + function calling)
# --------------------------------------------------------------------------- #

TOOL_EXECUTORS = {
    "clawbio": execute_clawbio,
    "save_file": execute_save_file,
    "write_file": execute_write_file,
    "generate_audio": execute_generate_audio,
}

MAX_TOOL_ITERATIONS = 10


async def llm_tool_loop(chat_id: int, user_content: str | list) -> str:
    """
    Run the LLM tool-use loop (OpenAI chat completions format):
    1. Append user message to history
    2. Call LLM with system prompt + history + tools
    3. If tool_calls -> execute -> append results -> call again
    4. Return final text
    """
    history = conversations.setdefault(chat_id, [])

    # Build user message in OpenAI format
    if isinstance(user_content, str):
        history.append({"role": "user", "content": user_content})
    else:
        # Multimodal content blocks — convert to OpenAI format
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

        # Append assistant message to history
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

        # No tool calls — return text
        if not last_message.tool_calls:
            return last_message.content or "(no response)"

        # Execute tool calls and append results
        for tc in last_message.tool_calls:
            func_name = tc.function.name
            executor = TOOL_EXECUTORS.get(func_name)
            if executor:
                try:
                    args = json.loads(tc.function.arguments)
                except json.JSONDecodeError:
                    args = {}
                logger.info(f"Tool call: {func_name}({json.dumps(args)[:200]})")
                try:
                    result = await executor(args)
                except Exception as tool_err:
                    logger.error(f"Tool {func_name} raised: {tool_err}", exc_info=True)
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
# Telegram helpers
# --------------------------------------------------------------------------- #


def is_authorised(update: Update) -> bool:
    """Check if the message is from the authorised chat."""
    return update.effective_chat.id == TELEGRAM_CHAT_ID


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


async def send_long_message(update: Update, text: str):
    """Send a message, splitting at 4096 chars if needed. Strips markup."""
    text = strip_markup(text)
    MAX_LEN = 4096
    if len(text) <= MAX_LEN:
        await update.message.reply_text(text)
        return
    chunks = []
    while text:
        if len(text) <= MAX_LEN:
            chunks.append(text)
            break
        split_at = text.rfind("\n\n", 0, MAX_LEN)
        if split_at == -1:
            split_at = text.rfind("\n", 0, MAX_LEN)
        if split_at == -1:
            split_at = MAX_LEN
        chunks.append(text[:split_at])
        text = text[split_at:].lstrip("\n")
    for chunk in chunks:
        if chunk.strip():
            await update.message.reply_text(chunk)


# --------------------------------------------------------------------------- #
# Command handlers
# --------------------------------------------------------------------------- #


async def cmd_start(update: Update, context: ContextTypes.DEFAULT_TYPE):
    """Handle /start command."""
    if not is_authorised(update):
        return
    await update.message.reply_text(
        "Hi there! RoboTerri here -- your ClawBio bioinformatics assistant ;-)\n\n"
        "Commands:\n"
        "  /skills  -- list available ClawBio skills\n"
        "  /demo <skill>  -- run a demo (pharmgx, equity, nutrigx, compare, prs, profile)\n"
        "  /voice  -- toggle voice replies on/off\n"
        "  /status  -- bot uptime and model info\n"
        "  /health  -- system health check\n\n"
        "Or just chat -- I can answer bioinformatics questions.\n"
        "Send a genetic data file to get your pharmacogenomics report.\n"
        "Then ask: \"what diseases am I at risk for?\" or \"show my full profile\"\n"
        "Send a photo of a medication for personalised drug guidance."
    )


async def cmd_skills(update: Update, context: ContextTypes.DEFAULT_TYPE):
    """Handle /skills command -- list available ClawBio skills."""
    if not is_authorised(update):
        return
    try:
        proc = await asyncio.create_subprocess_exec(
            sys.executable, str(CLAWBIO_PY), "list",
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            cwd=str(CLAWBIO_DIR),
        )
        stdout, _ = await asyncio.wait_for(proc.communicate(), timeout=15)
        output = stdout.decode(errors="replace").strip()
        await send_long_message(update, output or "No skills found.")
    except Exception as e:
        await update.message.reply_text(f"Error listing skills: {e}")


async def cmd_demo(update: Update, context: ContextTypes.DEFAULT_TYPE):
    """Handle /demo <skill> command -- run a skill with demo data."""
    if not is_authorised(update):
        return
    skill = context.args[0] if context.args else "pharmgx"
    await update.message.reply_text(f"Running {skill} demo -- this may take a moment...")
    await context.bot.send_chat_action(
        chat_id=update.effective_chat.id, action="typing"
    )
    try:
        reply = await llm_tool_loop(
            update.effective_chat.id,
            f"Run the {skill} demo using the clawbio tool with mode='demo'."
        )
        if _pending_text:
            reply = "\n\n".join(_pending_text)
            _pending_text.clear()
        await send_long_message(update, reply)
        await _drain_pending_media(update, context)
    except Exception as e:
        logger.error(f"Demo error: {e}", exc_info=True)
        await update.message.reply_text(f"Demo failed: {e}")


async def cmd_voice(update: Update, context: ContextTypes.DEFAULT_TYPE):
    """Handle /voice command -- toggle voice replies on/off."""
    if not is_authorised(update):
        return
    current = context.user_data.get("voice_replies", False)
    context.user_data["voice_replies"] = not current
    state = "ON" if not current else "OFF"
    await update.message.reply_text(
        f"Voice replies toggled {state}.\n"
        f"{'I will now send voice memos alongside text replies.' if not current else 'Back to text-only replies.'}"
    )


async def cmd_status(update: Update, context: ContextTypes.DEFAULT_TYPE):
    """Handle /status command -- report uptime and model info."""
    if not is_authorised(update):
        return

    uptime_secs = int(time.time() - BOT_START_TIME)
    hours, remainder = divmod(uptime_secs, 3600)
    minutes, secs = divmod(remainder, 60)
    uptime_str = f"{hours}h {minutes}m {secs}s"

    # Count available skills
    skills_dir = CLAWBIO_DIR / "skills"
    skill_count = sum(
        1 for d in skills_dir.iterdir()
        if d.is_dir() and (d / "SKILL.md").exists()
    ) if skills_dir.exists() else 0

    status_msg = (
        f"RoboTerri ClawBio Status\n"
        f"========================\n"
        f"Bot uptime: {uptime_str}\n"
        f"LLM model: {CLAWBIO_MODEL}\n"
        f"Skills available: {skill_count}\n"
        f"ClawBio dir: {CLAWBIO_DIR}\n"
    )
    if LLM_BASE_URL:
        status_msg += f"LLM endpoint: {LLM_BASE_URL}\n"

    await update.message.reply_text(status_msg)


async def cmd_health(update: Update, context: ContextTypes.DEFAULT_TYPE):
    """Handle /health command -- system health check."""
    if not is_authorised(update):
        return

    checks = []

    # ClawBio CLI
    if CLAWBIO_PY.exists():
        checks.append("ClawBio CLI: OK")
    else:
        checks.append("ClawBio CLI: MISSING")

    # SOUL.md
    if SOUL_MD.exists():
        checks.append(f"SOUL.md: OK ({len(_soul)} chars)")
    else:
        checks.append("SOUL.md: MISSING (using fallback)")

    # Skills
    skills_dir = CLAWBIO_DIR / "skills"
    if skills_dir.exists():
        implemented = []
        stub_only = []
        for d in sorted(skills_dir.iterdir()):
            if not d.is_dir() or not (d / "SKILL.md").exists():
                continue
            has_py = any(d.glob("*.py"))
            if has_py:
                implemented.append(d.name)
            else:
                stub_only.append(d.name)
        checks.append(f"Skills (implemented): {len(implemented)}")
        checks.append(f"Skills (stub/planned): {len(stub_only)}")
    else:
        checks.append("Skills directory: MISSING")

    # Output directory
    if OUTPUT_DIR.exists():
        output_count = sum(1 for _ in OUTPUT_DIR.iterdir()) if OUTPUT_DIR.exists() else 0
        checks.append(f"Output runs: {output_count}")
    else:
        checks.append("Output directory: not yet created")

    # edge-tts availability
    edge_tts_bin = Path.home() / "Library" / "Python" / "3.9" / "bin" / "edge-tts"
    if edge_tts_bin.exists():
        checks.append("edge-tts: OK")
    else:
        checks.append("edge-tts: not found (audio generation unavailable)")

    await update.message.reply_text(
        "ClawBio Health Check\n"
        "====================\n" + "\n".join(checks)
    )


# --------------------------------------------------------------------------- #
# Voice reply helper
# --------------------------------------------------------------------------- #


async def _send_voice_reply(bot, chat_id: int, text: str) -> bool:
    """Convert text to OGG/Opus voice message and send via Telegram.

    Uses macOS say (Daniel voice) -> ffmpeg -> OGG/Opus.
    Returns True if voice was sent, False on failure.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        text_file = os.path.join(tmpdir, "reply.txt")
        aiff_file = os.path.join(tmpdir, "reply.aiff")
        ogg_file = os.path.join(tmpdir, "reply.ogg")

        with open(text_file, "w", encoding="utf-8") as f:
            f.write(text)

        # Generate speech with macOS say (Daniel = British male)
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

        # Convert to OGG/Opus (Telegram voice format)
        proc = await asyncio.create_subprocess_exec(
            "ffmpeg", "-y", "-i", aiff_file,
            "-codec:a", "libopus", "-b:a", "48k",
            ogg_file,
            stdout=asyncio.subprocess.DEVNULL,
            stderr=asyncio.subprocess.DEVNULL,
        )
        await proc.wait()
        if proc.returncode != 0:
            logger.warning("ffmpeg OGG conversion failed for voice reply")
            return False

        with open(ogg_file, "rb") as audio:
            await bot.send_voice(chat_id=chat_id, voice=audio)

    return True


# --------------------------------------------------------------------------- #
# Message handlers
# --------------------------------------------------------------------------- #


async def handle_message(update: Update, context: ContextTypes.DEFAULT_TYPE):
    """Handle incoming text messages via the LLM tool loop."""
    if not is_authorised(update):
        return
    if not update.message or not update.message.text:
        return

    user_text = update.message.text
    logger.info(f"Message from {update.effective_user.first_name}: {user_text[:100]}")

    try:
        await context.bot.send_chat_action(
            chat_id=update.effective_chat.id, action="typing"
        )
        reply = await llm_tool_loop(update.effective_chat.id, user_text)
        if _pending_text:
            reply = "\n\n".join(_pending_text)
            _pending_text.clear()
        await send_long_message(update, reply)
        await _drain_pending_media(update, context)

        # Voice reply if toggled on
        if context.user_data.get("voice_replies"):
            try:
                await _send_voice_reply(
                    context.bot, update.effective_chat.id, reply
                )
            except Exception as ve:
                logger.warning(f"Voice reply failed: {ve}")
    except Exception as e:
        logger.error(f"Error handling message: {e}", exc_info=True)
        await update.message.reply_text(
            f"Sorry, something went wrong -- {type(e).__name__}: {e}"
        )


async def handle_photo(update: Update, context: ContextTypes.DEFAULT_TYPE):
    """Handle photos: download -> base64 -> LLM vision (drug detection)."""
    if not is_authorised(update):
        return
    if not update.message:
        return

    await context.bot.send_chat_action(
        chat_id=update.effective_chat.id, action="typing"
    )

    try:
        photo = update.message.photo[-1] if update.message.photo else None
        doc = update.message.document if not photo else None

        if not photo and not doc:
            return

        if doc:
            mime = doc.mime_type or ""
            if not mime.startswith("image/"):
                return
            file = await doc.get_file()
            media_type = mime
            filename = doc.file_name or "image.jpg"
        else:
            file = await photo.get_file()
            media_type = "image/jpeg"
            filename = "photo.jpg"

        img_bytes = await file.download_as_bytearray()

        # File size check (TG-004)
        if len(img_bytes) > MAX_PHOTO_BYTES:
            await update.message.reply_text(
                f"Photo too large ({len(img_bytes) / (1024*1024):.1f} MB). "
                f"Maximum: {MAX_PHOTO_BYTES / (1024*1024):.0f} MB."
            )
            return

        img_b64 = base64.standard_b64encode(bytes(img_bytes)).decode("ascii")
        logger.info(f"Photo received: {len(img_bytes)} bytes, type={media_type}")

        # Sanitize filename (TG-002)
        filename = _sanitize_filename(filename)

        # Store for potential file-based skill use
        tmp_path = Path(tempfile.gettempdir()) / f"roboterri_{filename}"
        tmp_path.write_bytes(bytes(img_bytes))
        _received_files[update.effective_chat.id] = {
            "path": str(tmp_path), "filename": filename,
        }

        caption = update.message.caption or ""
        # Use internal format; llm_tool_loop converts to OpenAI image_url format
        content_blocks = [
            {
                "type": "image",
                "source": {
                    "type": "base64",
                    "media_type": media_type,
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

        reply = await llm_tool_loop(update.effective_chat.id, content_blocks)
        if _pending_text:
            reply = "\n\n".join(_pending_text)
            _pending_text.clear()
        await send_long_message(update, reply)

    except Exception as e:
        logger.error(f"Photo handling error: {e}", exc_info=True)
        await update.message.reply_text(
            f"Sorry, I couldn't process that image -- {type(e).__name__}: {e}"
        )


async def handle_document(update: Update, context: ContextTypes.DEFAULT_TYPE):
    """Handle documents: download -> detect genetic file -> route to skill."""
    if not is_authorised(update):
        return
    if not update.message or not update.message.document:
        return

    doc = update.message.document
    mime = doc.mime_type or ""

    # Images handled by handle_photo
    if mime.startswith("image/"):
        return

    await context.bot.send_chat_action(
        chat_id=update.effective_chat.id, action="typing"
    )

    try:
        file = await doc.get_file()
        filename = _sanitize_filename(doc.file_name or "document")
        file_size = doc.file_size or 0

        # File size check (TG-004)
        if file_size > MAX_UPLOAD_BYTES:
            await update.message.reply_text(
                f"File too large ({file_size / (1024*1024):.1f} MB). "
                f"Maximum: {MAX_UPLOAD_BYTES / (1024*1024):.0f} MB."
            )
            return

        tmp_path = Path(tempfile.gettempdir()) / f"roboterri_{filename}"
        await file.download_to_drive(str(tmp_path))
        logger.info(f"Document received: {filename} ({file_size} bytes, {mime})")

        _received_files[update.effective_chat.id] = {
            "path": str(tmp_path), "filename": filename,
        }

        # Auto-create a patient profile for follow-up skill calls
        profile_path = None
        try:
            upload_proc = await asyncio.create_subprocess_exec(
                sys.executable, str(CLAWBIO_PY), "upload",
                "--input", str(tmp_path),
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            up_stdout, up_stderr = await asyncio.wait_for(
                upload_proc.communicate(), timeout=30,
            )
            up_out = up_stdout.decode(errors="replace")
            # Parse profile path from upload output
            for line in up_out.splitlines():
                if "profile" in line.lower() and ("/" in line or "\\" in line):
                    # Extract path-like token from the line
                    for token in line.split():
                        if token.endswith(".json"):
                            profile_path = token
                            break
            if profile_path:
                _received_files[update.effective_chat.id]["profile_path"] = profile_path
                logger.info(f"Auto-created profile: {profile_path}")
            else:
                logger.info(f"Profile upload output (no path parsed): {up_out[:200]}")
        except Exception as prof_err:
            logger.warning(f"Auto-profile creation failed (non-fatal): {prof_err}")

        caption = update.message.caption or ""
        parts = [f"[Document received: {filename} ({mime}, {file_size} bytes)]"]
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

        reply = await llm_tool_loop(
            update.effective_chat.id, "\n\n".join(parts)
        )
        if _pending_text:
            reply = "\n\n".join(_pending_text)
            _pending_text.clear()
        await send_long_message(update, reply)
        await _drain_pending_media(update, context)

    except Exception as e:
        logger.error(f"Document handling error: {e}", exc_info=True)
        await update.message.reply_text(
            f"Sorry, I couldn't process that document -- {type(e).__name__}: {e}"
        )


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #


def main():
    """Start the bot."""
    logger.info(f"Starting RoboTerri ClawBio bot (model: {CLAWBIO_MODEL})")
    logger.info(f"ClawBio directory: {CLAWBIO_DIR}")
    if LLM_BASE_URL:
        logger.info(f"LLM base URL: {LLM_BASE_URL}")
    logger.info(f"Authorised chat ID: {TELEGRAM_CHAT_ID}")

    app = Application.builder().token(TELEGRAM_BOT_TOKEN).build()

    # Commands
    app.add_handler(CommandHandler("start", cmd_start))
    app.add_handler(CommandHandler("skills", cmd_skills))
    app.add_handler(CommandHandler("demo", cmd_demo))
    app.add_handler(CommandHandler("voice", cmd_voice))
    app.add_handler(CommandHandler("status", cmd_status))
    app.add_handler(CommandHandler("health", cmd_health))

    # Message handlers
    app.add_handler(MessageHandler(
        filters.PHOTO | (filters.Document.IMAGE & ~filters.COMMAND),
        handle_photo,
    ))
    app.add_handler(MessageHandler(
        filters.Document.ALL & ~filters.Document.IMAGE & ~filters.COMMAND,
        handle_document,
    ))
    app.add_handler(MessageHandler(
        filters.TEXT & ~filters.COMMAND,
        handle_message,
    ))

    print("RoboTerri ClawBio bot is running. Press Ctrl+C to stop.")
    app.run_polling()


if __name__ == "__main__":
    main()
