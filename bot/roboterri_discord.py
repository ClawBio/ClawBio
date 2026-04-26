#!/usr/bin/env python3
"""
roboterri_discord.py — RoboTerri ClawBio Discord Bot
=====================================================
A Discord bot that runs ClawBio bioinformatics skills using any LLM
as the reasoning engine. Handles text messages, genetic file uploads,
and medication photos.

Works with any OpenAI-compatible provider: OpenAI, Anthropic (via proxy),
Google, Mistral, Groq, Together, OpenRouter, Ollama, LM Studio, etc.

Prerequisites:
    pip3 install discord.py openai python-dotenv

Usage:
    # Set environment variables in .env (see bot/README.md)
    python3 bot/roboterri_discord.py
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

import discord
from dotenv import load_dotenv
from openai import AsyncOpenAI, APIError

_PROJECT_ROOT_FOR_IMPORT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT_FOR_IMPORT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT_FOR_IMPORT))

from clawbio.skill_intents import (
    load_default_skill_registry,
    plan_skill_intent,
    skill_intent_tool_summary,
    skill_names_for_tool_schema,
)
from bot.tool_loop_utils import execute_tool_calls_safely

# --------------------------------------------------------------------------- #
# Config
# --------------------------------------------------------------------------- #

load_dotenv()

DISCORD_BOT_TOKEN = os.environ.get("DISCORD_BOT_TOKEN", "")
LLM_API_KEY = os.environ.get("LLM_API_KEY", "")
LLM_BASE_URL = os.environ.get("LLM_BASE_URL", "")
CLAWBIO_MODEL = os.environ.get("CLAWBIO_MODEL", "gpt-4.1-mini")
ADMIN_USER_ID = int(os.environ.get("DISCORD_ADMIN_USER_ID", "0") or "0")

# Rate limiting: messages per user per hour (0 = unlimited)
RATE_LIMIT_PER_HOUR = int(os.environ.get("RATE_LIMIT_PER_HOUR", "10"))

CHANNELS_FILE = Path(__file__).resolve().parent / ".channels.json"


def load_channels() -> list[dict]:
    """Load authorised channels from .channels.json."""
    if not CHANNELS_FILE.exists():
        # Fall back to env var for backwards compatibility
        env_id = os.environ.get("DISCORD_CHANNEL_ID", "0")
        if env_id and env_id != "0":
            return [{"id": int(env_id), "name": "default", "skills": "all"}]
        return []
    with open(CHANNELS_FILE, encoding="utf-8") as f:
        return json.load(f)


CHANNELS = load_channels()
AUTHORISED_CHANNEL_IDS = {ch["id"] for ch in CHANNELS}


def reload_channels():
    """Reload channels from .channels.json in place."""
    CHANNELS.clear()
    CHANNELS.extend(load_channels())
    AUTHORISED_CHANNEL_IDS.clear()
    AUTHORISED_CHANNEL_IDS.update(ch["id"] for ch in CHANNELS)


def get_channel_config(channel_id: int) -> dict | None:
    """Return config for a channel, or None if not authorised."""
    for ch in CHANNELS:
        if ch["id"] == channel_id:
            return ch
    return None


if not DISCORD_BOT_TOKEN:
    print("Error: DISCORD_BOT_TOKEN not set. See bot/README.md for setup.")
    sys.exit(1)
if not LLM_API_KEY:
    print("Error: LLM_API_KEY not set. See bot/README.md for setup.")
    sys.exit(1)
if not AUTHORISED_CHANNEL_IDS:
    print("Error: No channels configured. Add channels to bot/.channels.json or set DISCORD_CHANNEL_ID in .env.")
    sys.exit(1)

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

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
)
logger = logging.getLogger("roboterri-discord")


# ---------------------------------------------------------------------------
# Redact bot token from log output
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


if DISCORD_BOT_TOKEN:
    _redact = _TokenRedactFilter(DISCORD_BOT_TOKEN)
    for _ln in ("discord", "discord.http", "discord.gateway"):
        logging.getLogger(_ln).addFilter(_redact)


# ---------------------------------------------------------------------------
# Structured audit log (JSONL)
# ---------------------------------------------------------------------------
_AUDIT_LOG_DIR = CLAWBIO_DIR / "bot" / "logs"
_AUDIT_LOG_DIR.mkdir(parents=True, exist_ok=True)
_AUDIT_LOG_PATH = _AUDIT_LOG_DIR / "audit_discord.jsonl"


def _audit(event: str, **kwargs):
    """Append a structured JSON event to the audit log."""
    from datetime import timezone as _tz
    entry = {"ts": datetime.now(_tz.utc).isoformat(), "event": event, **kwargs}
    try:
        with open(_AUDIT_LOG_PATH, "a", encoding="utf-8") as f:
            f.write(json.dumps(entry, default=str) + "\n")
    except OSError:
        pass


def _user_ctx(message: discord.Message) -> dict:
    """Extract user identity for audit logging."""
    return {
        "user_id": message.author.id,
        "username": str(message.author),
        "display_name": message.author.display_name,
        "channel_id": message.channel.id,
        "channel_name": getattr(message.channel, "name", "DM"),
        "is_admin": is_admin(message),
    }


def is_admin(message: discord.Message) -> bool:
    """Check if the message is from the admin user."""
    return bool(ADMIN_USER_ID) and message.author.id == ADMIN_USER_ID


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

conversations: dict[int, list] = {}
MAX_HISTORY = 20

# Per-channel received file storage
_received_files: dict[int, dict] = {}

# Pending media queue: channel_id -> list of {"type": "document"|"photo", "path": str}
_pending_media: dict[int, list[dict]] = {}

# Pending text queue: bypass LLM paraphrasing for compare/drugphoto
_pending_text: list[str] = []

# Per-user voice reply toggle: user_id -> bool
_voice_enabled: dict[int, bool] = {}

BOT_START_TIME = time.time()

_SKILL_REGISTRY = load_default_skill_registry(CLAWBIO_DIR)
_SKILL_TOOL_ENUM = skill_names_for_tool_schema(_SKILL_REGISTRY, CLAWBIO_DIR)
_DESCRIPTOR_TOOL_SUMMARY = skill_intent_tool_summary(_SKILL_REGISTRY, CLAWBIO_DIR)

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
                + (f"Descriptor-provided skill intents: {_DESCRIPTOR_TOOL_SUMMARY}. " if _DESCRIPTOR_TOOL_SUMMARY else "")
                + (
                "Use mode='demo' to run with built-in demo data. "
                "Use mode='file' when the user has sent a genetic data file. "
                "Use skill='auto' to let the orchestrator detect the right skill. "
                "IMPORTANT: When this tool returns results, relay the output VERBATIM. "
                "Do not paraphrase, summarise, or rewrite. The output contains exact numerical "
                "results (IBS scores, percentages, gene-drug interactions) that must be shown unchanged."
                )
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "skill": {
                        "type": "string",
                        "enum": _SKILL_TOOL_ENUM,
                        "description": (
                            "Which bioinformatics skill to run. Use 'auto' to let "
                            "the orchestrator detect from the file type or query."
                        ),
                    },
                    "mode": {
                        "type": "string",
                        "enum": ["file", "demo"],
                        "description": (
                            "file: use a file the user sent via Discord. "
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
                "Save a file that was sent via Discord to a specific folder. "
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
                            "the original filename from Discord."
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
                "Generate an MP3 audio file from text using OpenAI TTS. "
                "Produces natural, human-sounding speech. Good for converting reports "
                "into accessible audio. "
                "Available voices: nova (warm female, default), shimmer (smooth female), "
                "alloy (neutral), echo (male), fable (British), onyx (deep male). "
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
                        "description": "TTS voice. Default: 'nova'.",
                        "enum": ["nova", "shimmer", "alloy", "echo", "fable", "onyx"],
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


def _sanitize_filename(filename: str) -> str:
    """Strip path traversal components and dangerous characters from a filename."""
    filename = Path(filename).name
    filename = re.sub(r"[\x00-\x1f]", "", filename)
    filename = filename.replace("..", "").replace("/", "").replace("\\", "")
    if not filename:
        filename = "unnamed_file"
    return filename


def _resolve_dest(folder: str | None) -> Path:
    """Resolve a destination folder, restricted to CLAWBIO_DIR."""
    dest = Path(folder) if folder else DATA_DIR
    if not dest.is_absolute():
        dest = CLAWBIO_DIR / dest
    try:
        dest.resolve().relative_to(CLAWBIO_DIR.resolve())
    except ValueError:
        logger.warning(f"Path escape blocked: {dest}")
        _audit("security", severity="HIGH", detail="path_escape_blocked",
               attempted_path=str(dest), function="_resolve_dest")
        dest = DATA_DIR
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
# execute_clawbio
# --------------------------------------------------------------------------- #


async def execute_clawbio(args: dict) -> str:
    """Execute a ClawBio bioinformatics skill via subprocess."""
    skill_key = args.get("skill", "auto")
    mode = args.get("mode", "demo")
    query = args.get("query", "")
    raw_user_text = args.get("_raw_user_text") or query or ""
    skill_registry = _SKILL_REGISTRY
    preplanned_plan = None

    # Auto-routing via orchestrator
    if skill_key == "auto":
        descriptor_plan = plan_skill_intent(
            user_text=raw_user_text,
            requested_skill=skill_key,
            requested_mode=mode,
            attachments=[],
            skill_registry=skill_registry,
            project_root=CLAWBIO_DIR,
        )
        if descriptor_plan.intent_id != "legacy_fallback":
            preplanned_plan = descriptor_plan
        else:
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
                        f"available via Discord. Available: {avail}"
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
        # Fall back to owner's genome for admin users
        if OWNER_GENOME.exists():
            input_path = str(OWNER_GENOME)
            logger.info(f"No file uploaded — using owner genome: {OWNER_GENOME.name}")
        else:
            return "Error: no file received. Send a genetic data file first, then run the skill."

    attachments = []
    if input_path or profile_path:
        attachments.append({"path": str(input_path) if input_path else None, "profile_path": profile_path})
    for key in ("trait", "gene", "rsid", "drug_name", "visible_dose"):
        if args.get(key):
            attachments.append({key: args[key]})

    plan = preplanned_plan or plan_skill_intent(
        user_text=raw_user_text,
        requested_skill=skill_key,
        requested_mode=mode,
        attachments=attachments,
        skill_registry=skill_registry,
        project_root=CLAWBIO_DIR,
    )
    _audit(
        "skill_intent_plan",
        channel_id=args.get("_channel_id"),
        raw_user_text_sha256=plan.raw_user_text_sha256,
        raw_user_text_preview=plan.raw_user_text[:200],
        selected_skill=plan.skill,
        selected_intent=plan.intent_id,
        matched_route=plan.matched_route,
        commands=[item.argv for item in plan.executions],
    )
    logger.info(
        "Skill intent plan: skill=%s intent=%s status=%s reason=%s",
        plan.skill, plan.intent_id, plan.status, plan.reason,
    )
    if plan.status == "needs_confirmation":
        return f"Confirmation required before running {plan.skill}: {plan.reason}"
    if plan.status == "needs_input" or not plan.executions:
        return plan.reason or "I need an input file or a clearer skill request before running ClawBio."

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    stdout_parts = []
    output_dirs: list[Path] = []
    executed_skills: list[str] = []

    for index, execution in enumerate(plan.executions):
        cmd = list(execution.argv)
        run_skill = execution.skill
        executed_skills.append(run_skill)
        if "--output" in cmd:
            out_dir = Path(cmd[cmd.index("--output") + 1])
        elif run_skill not in ("compare", "drugphoto"):
            suffix = f"_{index + 1}" if len(plan.executions) > 1 else ""
            out_dir = OUTPUT_DIR / f"{run_skill}_{ts}{suffix}"
            cmd.extend(["--output", str(out_dir)])
        else:
            out_dir = None
        if out_dir:
            output_dirs.append(out_dir)
        _audit(
            "skill_execution_command",
            channel_id=args.get("_channel_id"),
            selected_skill=run_skill,
            selected_intent=plan.intent_id,
            command=cmd,
            output_bundle_path=str(out_dir) if out_dir else None,
        )

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
            return f"{run_skill} timed out after 120 seconds."
        except Exception:
            import traceback as _tb
            return f"{run_skill} crashed:\n{_tb.format_exc()[-1500:]}"

        stdout_parts.append(stdout_str)
        if proc.returncode != 0:
            err = stderr_str[-1500:] if stderr_str else stdout_str[-1500:] if stdout_str else "unknown error"
            return f"{run_skill} failed (exit {proc.returncode}):\n{err}"

    skill_key = executed_skills[-1] if executed_skills else skill_key
    stdout_str = "\n".join(part for part in stdout_parts if part)
    out_dir = output_dirs[-1] if output_dirs else OUTPUT_DIR / f"{skill_key}_{ts}"

    # For compare / drugphoto / profile: send stdout directly (bypass LLM paraphrasing)
    if any(item in ("compare", "drugphoto", "profile") for item in executed_skills):
        raw_output = stdout_str.strip()
        if raw_output:
            _pending_text.append(raw_output)
        return "Result sent directly to chat. Do not repeat or paraphrase it."

    # For other skills: collect report + figures from output directory
    media_items = []
    for bundle_dir in output_dirs:
        if not bundle_dir.exists():
            continue
        for f in sorted(bundle_dir.rglob("*")):
            if not f.is_file():
                continue
            if f.suffix == ".md":
                media_items.append({"type": "document", "path": str(f)})
            elif f.suffix == ".png":
                media_items.append({"type": "photo", "path": str(f)})
    if media_items:
        _pending_media[0] = _pending_media.get(0, []) + media_items

    # Read report for chat display
    report_text = ""
    for bundle_dir in output_dirs:
        if not bundle_dir.exists():
            continue
        for pattern in ["report.md", "*_report.md", "*.md"]:
            for md_file in sorted(bundle_dir.glob(pattern)):
                if md_file.name.startswith("."):
                    continue
                report_text = md_file.read_text(encoding="utf-8")
                break
            if report_text:
                break
        if report_text:
            break

    if not report_text:
        return stdout_str if stdout_str else f"{skill_key} completed. Output: {out_dir}"

    # Trim verbose sections for readability but ALWAYS keep disclaimer.
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
            skip = False  # always show disclaimer
        if line.startswith("!["):
            continue
        if not skip:
            keep_lines.append(line)

    return "\n".join(keep_lines).strip()


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
    """Generate MP3 audio from text using OpenAI TTS API."""
    text = args.get("text")
    filename = args.get("filename")
    if not text:
        return "Error: 'text' is required. Provide the text to convert to speech."
    if not filename:
        return "Error: 'filename' is required (e.g. 'report.mp3')."
    if not filename.endswith(".mp3"):
        filename += ".mp3"

    filename = _sanitize_filename(filename)
    voice = args.get("voice", "nova")
    dest = _resolve_dest(args.get("destination_folder"))
    filepath = dest / filename

    if not _validate_path(filepath, dest):
        return f"Error: filename '{filename}' would escape the destination directory."

    # OpenAI TTS has a 4096-char input limit — split if needed
    MAX_CHUNK = 4096
    chunks = [text[i:i + MAX_CHUNK] for i in range(0, len(text), MAX_CHUNK)]

    try:
        # Use a direct OpenAI client for TTS (not the LLM proxy)
        tts_client = AsyncOpenAI(api_key=LLM_API_KEY)

        if len(chunks) == 1:
            response = await asyncio.wait_for(
                tts_client.audio.speech.create(
                    model="tts-1",
                    voice=voice,
                    input=chunks[0],
                ),
                timeout=300,
            )
            response.stream_to_file(str(filepath))
        else:
            # Multiple chunks: generate and concatenate
            part_files = []
            for i, chunk in enumerate(chunks):
                part_path = dest / f".tmp_{filename}_part{i}.mp3"
                response = await asyncio.wait_for(
                    tts_client.audio.speech.create(
                        model="tts-1",
                        voice=voice,
                        input=chunk,
                    ),
                    timeout=300,
                )
                response.stream_to_file(str(part_path))
                part_files.append(part_path)

            # Concatenate with ffmpeg
            list_file = dest / f".tmp_{filename}_list.txt"
            list_file.write_text(
                "\n".join(f"file '{p}'" for p in part_files),
                encoding="utf-8",
            )
            proc = await asyncio.create_subprocess_exec(
                "ffmpeg", "-y", "-f", "concat", "-safe", "0",
                "-i", str(list_file), "-c", "copy", str(filepath),
                stdout=asyncio.subprocess.DEVNULL,
                stderr=asyncio.subprocess.DEVNULL,
            )
            await proc.wait()

            # Cleanup temp files
            for p in part_files:
                try:
                    p.unlink()
                except OSError:
                    pass
            try:
                list_file.unlink()
            except OSError:
                pass

        size_mb = filepath.stat().st_size / (1024 * 1024)
        word_count = len(text.split())
        est_minutes = word_count / 150

        logger.info(f"Generated audio: {filepath} ({size_mb:.1f} MB, ~{est_minutes:.0f} min)")
        return (
            f"Audio saved to {filepath} ({size_mb:.1f} MB, "
            f"~{word_count} words, ~{est_minutes:.0f} min estimated)"
        )

    except asyncio.TimeoutError:
        return "Audio generation timed out after 5 minutes."
    except APIError as e:
        return f"OpenAI TTS API error: {e}"


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


async def llm_tool_loop(channel_id: int, user_content: str | list) -> str:
    """
    Run the LLM tool-use loop (OpenAI chat completions format):
    1. Append user message to history
    2. Call LLM with system prompt + history + tools
    3. If tool_calls -> execute -> append results -> call again
    4. Return final text
    """
    history = conversations.setdefault(channel_id, [])

    # Build user message in OpenAI format
    raw_user_text = user_content if isinstance(user_content, str) else ""
    if isinstance(user_content, str):
        history.append({"role": "user", "content": user_content})
    else:
        # Multimodal content blocks -- convert to OpenAI format
        oai_parts = []
        for block in user_content:
            if block.get("type") == "text":
                raw_user_text = f"{raw_user_text}\n{block['text']}".strip()
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

    # Sanitise: strip orphaned tool messages that lack a preceding
    # assistant message with tool_calls (prevents API 400 errors).
    sanitised: list[dict] = []
    for msg in history:
        if msg.get("role") == "tool":
            # Only keep if previous message is assistant with tool_calls
            if sanitised and sanitised[-1].get("role") == "assistant":
                if sanitised[-1].get("tool_calls"):
                    sanitised.append(msg)
                    continue
            logger.warning("Dropped orphaned tool message from history")
            _audit("history_sanitised", channel_id=channel_id,
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

        # No tool calls -- return text
        if not last_message.tool_calls:
            return last_message.content or "(no response)"

        tool_messages = await execute_tool_calls_safely(
            last_message.tool_calls,
            TOOL_EXECUTORS,
            base_args={"_channel_id": channel_id},
            raw_user_text=raw_user_text,
            audit=_audit,
            audit_context={"channel_id": channel_id},
            logger=logger,
        )
        history.extend(tool_messages)

    return last_message.content if last_message and last_message.content else "(max tool iterations reached)"


# --------------------------------------------------------------------------- #
# Rate limiting
# --------------------------------------------------------------------------- #

_rate_buckets: dict[int, list[float]] = {}


def _check_rate_limit(message: discord.Message) -> bool:
    """Return True if the user is within rate limits (or is admin)."""
    if RATE_LIMIT_PER_HOUR <= 0 or is_admin(message):
        return True
    uid = message.author.id
    now = time.time()
    window = 3600  # 1 hour
    bucket = _rate_buckets.setdefault(uid, [])
    bucket[:] = [t for t in bucket if now - t < window]
    if len(bucket) >= RATE_LIMIT_PER_HOUR:
        return False
    bucket.append(now)
    return True


# --------------------------------------------------------------------------- #
# Discord helpers
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


async def send_long_message(channel: discord.abc.Messageable, text: str):
    """Send a message, splitting at 2000 chars (Discord limit). Strips markup."""
    text = strip_markup(text)
    if not text:
        return
    MAX_LEN = 2000
    if len(text) <= MAX_LEN:
        await channel.send(text)
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
            await channel.send(chunk)


async def drain_pending_media(channel: discord.abc.Messageable) -> None:
    """Send any queued ClawBio media (documents + figures) after the text reply."""
    items = _pending_media.pop(0, [])
    if not items:
        return
    for item in items:
        try:
            path = Path(item["path"])
            if not path.exists():
                continue
            caption = path.stem.replace("_", " ").title() if item["type"] == "photo" else ""
            await channel.send(
                content=caption or None,
                file=discord.File(str(path), filename=path.name),
            )
        except Exception as e:
            logger.warning(f"Failed to send media {item['path']}: {e}")


# --------------------------------------------------------------------------- #
# Voice reply helper
# --------------------------------------------------------------------------- #


async def _send_voice_reply(channel: discord.abc.Messageable, text: str) -> bool:
    """Convert text to MP3 voice message and send via Discord.

    Uses OpenAI TTS API (nova voice) for natural-sounding speech.
    Returns True if voice was sent, False on failure.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        mp3_file = os.path.join(tmpdir, "voice_reply.mp3")

        try:
            tts_client = AsyncOpenAI(api_key=LLM_API_KEY)
            # Truncate to OpenAI's 4096-char limit for voice replies
            tts_text = text[:4096]
            response = await asyncio.wait_for(
                tts_client.audio.speech.create(
                    model="tts-1",
                    voice="nova",
                    input=tts_text,
                ),
                timeout=120,
            )
            response.stream_to_file(mp3_file)
        except Exception as e:
            logger.warning(f"OpenAI TTS failed for voice reply: {e}")
            return False

        await channel.send(
            file=discord.File(mp3_file, filename="voice_reply.mp3"),
        )

    return True


# --------------------------------------------------------------------------- #
# Discord client
# --------------------------------------------------------------------------- #

intents = discord.Intents.default()
intents.message_content = True
client = discord.Client(intents=intents)


@client.event
async def on_ready():
    logger.info(f"Logged in as {client.user} (id: {client.user.id})")
    logger.info(f"Authorised channels: {[ch['name'] for ch in CHANNELS]} ({len(CHANNELS)})")
    logger.info(f"LLM model: {CLAWBIO_MODEL}")
    logger.info(f"Admin user ID: {ADMIN_USER_ID or 'not set (public mode)'}")
    logger.info(f"Rate limit: {RATE_LIMIT_PER_HOUR} msgs/hour per user (0=unlimited)")
    if LLM_BASE_URL:
        logger.info(f"LLM base URL: {LLM_BASE_URL}")
    _audit("bot_start", model=CLAWBIO_MODEL,
           admin_user=ADMIN_USER_ID, rate_limit=RATE_LIMIT_PER_HOUR)
    print(f"RoboTerri Discord bot is running as {client.user}. Press Ctrl+C to stop.")


@client.event
async def on_message(message: discord.Message):
    # Ignore own messages
    if message.author == client.user:
        return

    # Only respond in authorised channels
    if message.channel.id not in AUTHORISED_CHANNEL_IDS:
        return

    # ----- Commands ----- #

    content = message.content.strip()

    if content == "!reload":
        reload_channels()
        await message.channel.send(
            f"Reloaded channel config -- {len(CHANNELS)} channel(s) authorised."
        )
        logger.info(f"Reloaded .channels.json: {AUTHORISED_CHANNEL_IDS}")
        return

    if content == "!start":
        await message.channel.send(
            "Welcome to ClawBio -- open-source bioinformatics at your fingertips!\n\n"
            "I can analyse genetic data, check drug interactions, assess nutritional "
            "genomics, estimate polygenic risk scores, and more.\n\n"
            "Commands:\n"
            "  `!skills`  -- list available bioinformatics skills\n"
            "  `!demo <skill>`  -- run a demo (pharmgx, equity, nutrigx, compare, prs, profile)\n"
            "  `!voice`  -- toggle voice replies on/off\n"
            "  `!status`  -- bot info\n"
            "  `!health`  -- system health check\n\n"
            "Or just chat -- ask any bioinformatics question.\n"
            "Attach a genetic data file (.txt, .csv, .vcf) to analyse it.\n"
            "Attach a photo of a medication for personalised drug guidance.\n\n"
            "ClawBio is a research tool, not a medical device. "
            "Consult a healthcare professional before making medical decisions."
        )
        return

    if content == "!skills":
        try:
            proc = await asyncio.create_subprocess_exec(
                sys.executable, str(CLAWBIO_PY), "list",
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=str(CLAWBIO_DIR),
            )
            stdout, _ = await asyncio.wait_for(proc.communicate(), timeout=15)
            output = stdout.decode(errors="replace").strip()
            await send_long_message(message.channel, output or "No skills found.")
        except Exception as e:
            await message.channel.send(f"Error listing skills: {e}")
        return

    if content == "!voice":
        uid = message.author.id
        current = _voice_enabled.get(uid, False)
        _voice_enabled[uid] = not current
        state = "ON" if not current else "OFF"
        await message.channel.send(
            f"Voice replies toggled {state}.\n"
            f"{'I will now send voice memos alongside text replies.' if not current else 'Back to text-only replies.'}"
        )
        return

    if content == "!status":
        uptime_secs = int(time.time() - BOT_START_TIME)
        hours, remainder = divmod(uptime_secs, 3600)
        minutes, secs = divmod(remainder, 60)
        uptime_str = f"{hours}h {minutes}m {secs}s"

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

        await message.channel.send(status_msg)
        return

    if content == "!health":
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
            output_count = sum(1 for _ in OUTPUT_DIR.iterdir())
            checks.append(f"Output runs: {output_count}")
        else:
            checks.append("Output directory: not yet created")

        # TTS availability
        if LLM_API_KEY:
            checks.append("TTS: OpenAI TTS (nova voice)")
        else:
            checks.append("TTS: unavailable (no LLM_API_KEY)")

        await message.channel.send(
            "ClawBio Health Check\n"
            "====================\n" + "\n".join(checks)
        )
        return

    if content.startswith("!demo"):
        if not _check_rate_limit(message):
            _audit("rate_limited", **_user_ctx(message))
            await message.channel.send(
                f"You've reached the limit of {RATE_LIMIT_PER_HOUR} messages per hour. "
                "Please try again later."
            )
            return
        parts = content.split(maxsplit=1)
        skill = parts[1].strip() if len(parts) > 1 else "pharmgx"
        await message.channel.send(f"Running {skill} demo -- this may take a moment...")
        async with message.channel.typing():
            try:
                reply = await llm_tool_loop(
                    message.channel.id,
                    f"Run the {skill} demo using the clawbio tool with mode='demo'.",
                )
                if _pending_text:
                    reply = "\n\n".join(_pending_text)
                    _pending_text.clear()
                await send_long_message(message.channel, reply)
                await drain_pending_media(message.channel)
                # Voice reply if toggled on
                if _voice_enabled.get(message.author.id):
                    try:
                        await _send_voice_reply(message.channel, reply)
                    except Exception as ve:
                        logger.warning(f"Voice reply failed: {ve}")
            except Exception as e:
                logger.error(f"Demo error: {e}", exc_info=True)
                await message.channel.send(f"Demo failed: {e}")
        return

    # ----- Attachments: images and genetic data files ----- #

    has_image = False
    has_genetic_file = False

    for attachment in message.attachments:
        ext = Path(attachment.filename).suffix.lower()
        content_type = attachment.content_type or ""

        if content_type.startswith("image/") or ext in IMAGE_EXTENSIONS:
            if not _check_rate_limit(message):
                _audit("rate_limited", **_user_ctx(message))
                await message.channel.send(
                    f"You've reached the limit of {RATE_LIMIT_PER_HOUR} messages per hour. "
                    "Please try again later."
                )
                return

            has_image = True
            # Download image and encode to base64
            img_bytes = await attachment.read()

            # File size check
            if len(img_bytes) > MAX_PHOTO_BYTES:
                await message.channel.send(
                    f"Photo too large ({len(img_bytes) / (1024*1024):.1f} MB). "
                    f"Maximum: {MAX_PHOTO_BYTES / (1024*1024):.0f} MB."
                )
                return

            img_b64 = base64.standard_b64encode(img_bytes).decode("ascii")

            media_type = content_type if content_type.startswith("image/") else "image/jpeg"
            filename = _sanitize_filename(attachment.filename)
            logger.info(f"Image received: {filename} ({len(img_bytes)} bytes, {media_type})")
            _audit("photo", **_user_ctx(message), size_bytes=len(img_bytes),
                   media_type=media_type)

            # Store for potential file-based skill use
            tmp_path = Path(tempfile.gettempdir()) / f"roboterri_{filename}"
            tmp_path.write_bytes(img_bytes)
            _received_files[message.channel.id] = {
                "path": str(tmp_path), "filename": filename,
            }

            caption = message.content.strip() if message.content else ""
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

            async with message.channel.typing():
                try:
                    reply = await llm_tool_loop(message.channel.id, content_blocks)
                    if _pending_text:
                        reply = "\n\n".join(_pending_text)
                        _pending_text.clear()
                    await send_long_message(message.channel, reply)
                    # Voice reply if toggled on
                    if _voice_enabled.get(message.author.id):
                        try:
                            await _send_voice_reply(message.channel, reply)
                        except Exception as ve:
                            logger.warning(f"Voice reply failed: {ve}")
                except Exception as e:
                    logger.error(f"Photo handling error: {e}", exc_info=True)
                    await message.channel.send(
                        f"Sorry, I couldn't process that image -- {type(e).__name__}: {e}"
                    )

        elif ext in GENETIC_EXTENSIONS:
            if not _check_rate_limit(message):
                _audit("rate_limited", **_user_ctx(message))
                await message.channel.send(
                    f"You've reached the limit of {RATE_LIMIT_PER_HOUR} messages per hour. "
                    "Please try again later."
                )
                return

            has_genetic_file = True
            file_bytes = await attachment.read()

            # File size check
            if len(file_bytes) > MAX_UPLOAD_BYTES:
                await message.channel.send(
                    f"File too large ({len(file_bytes) / (1024*1024):.1f} MB). "
                    f"Maximum: {MAX_UPLOAD_BYTES / (1024*1024):.0f} MB."
                )
                return

            filename = _sanitize_filename(attachment.filename)
            tmp_path = Path(tempfile.gettempdir()) / f"roboterri_{filename}"
            tmp_path.write_bytes(file_bytes)
            logger.info(f"Document received: {filename} ({len(file_bytes)} bytes)")
            _audit("document", **_user_ctx(message), filename=filename,
                   size_bytes=len(file_bytes))

            _received_files[message.channel.id] = {
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
                        for token in line.split():
                            if token.endswith(".json"):
                                profile_path = token
                                break
                if profile_path:
                    _received_files[message.channel.id]["profile_path"] = profile_path
                    logger.info(f"Auto-created profile: {profile_path}")
                else:
                    logger.info(f"Profile upload output (no path parsed): {up_out[:200]}")
            except Exception as prof_err:
                logger.warning(f"Auto-profile creation failed (non-fatal): {prof_err}")

            caption = message.content.strip() if message.content else ""
            parts_list = [f"[Document received: {filename} ({len(file_bytes)} bytes)]"]
            if profile_path:
                parts_list.append(f"[Patient profile auto-created: {profile_path}]")
            if caption:
                parts_list.append(caption)
            else:
                profile_note = (
                    " A patient profile has been created -- the user can now ask "
                    "follow-up questions like 'what am I at risk for?' (prs) or "
                    "'show my full profile' (profile) without re-uploading."
                ) if profile_path else ""
                parts_list.append(
                    "The user sent this genetic data file. Detect the file type and "
                    "run the appropriate ClawBio skill using mode='file'. For .txt "
                    "files (23andMe format) use pharmgx. For .csv (AncestryDNA) use "
                    "pharmgx. For .vcf use equity. For .fastq use metagenomics. "
                    "If unsure, use skill='auto'." + profile_note
                )

            async with message.channel.typing():
                try:
                    reply = await llm_tool_loop(
                        message.channel.id, "\n\n".join(parts_list)
                    )
                    if _pending_text:
                        reply = "\n\n".join(_pending_text)
                        _pending_text.clear()
                    await send_long_message(message.channel, reply)
                    await drain_pending_media(message.channel)
                    # Voice reply if toggled on
                    if _voice_enabled.get(message.author.id):
                        try:
                            await _send_voice_reply(message.channel, reply)
                        except Exception as ve:
                            logger.warning(f"Voice reply failed: {ve}")
                except Exception as e:
                    logger.error(f"Document handling error: {e}", exc_info=True)
                    await message.channel.send(
                        f"Sorry, I couldn't process that document -- {type(e).__name__}: {e}"
                    )

    # If there were only attachments (no text beyond them), we're done
    if has_image or has_genetic_file:
        return

    # ----- Plain text messages ----- #

    if not content:
        return

    # Ignore bot commands already handled above
    if content.startswith("!"):
        return

    if not _check_rate_limit(message):
        _audit("rate_limited", **_user_ctx(message))
        await message.channel.send(
            f"You've reached the limit of {RATE_LIMIT_PER_HOUR} messages per hour. "
            "Please try again later."
        )
        return

    user_text = content
    logger.info(f"Message from {message.author.display_name}: {user_text[:100]}")
    _audit("message", **_user_ctx(message), text_preview=user_text[:200],
           text_len=len(user_text))

    async with message.channel.typing():
        try:
            reply = await llm_tool_loop(message.channel.id, user_text)
            if _pending_text:
                reply = "\n\n".join(_pending_text)
                _pending_text.clear()
            await send_long_message(message.channel, reply)
            await drain_pending_media(message.channel)
            # Voice reply if toggled on
            if _voice_enabled.get(message.author.id):
                try:
                    await _send_voice_reply(message.channel, reply)
                except Exception as ve:
                    logger.warning(f"Voice reply failed: {ve}")
        except Exception as e:
            logger.error(f"Error handling message: {e}", exc_info=True)
            await message.channel.send(
                f"Sorry, something went wrong -- {type(e).__name__}: {e}"
            )


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #


def main():
    """Start the bot."""
    logger.info(f"Starting RoboTerri Discord bot (model: {CLAWBIO_MODEL})")
    logger.info(f"ClawBio directory: {CLAWBIO_DIR}")
    if LLM_BASE_URL:
        logger.info(f"LLM base URL: {LLM_BASE_URL}")
    logger.info(f"Admin user ID: {ADMIN_USER_ID or 'not set (public mode)'}")
    logger.info(f"Rate limit: {RATE_LIMIT_PER_HOUR} msgs/hour per user (0=unlimited)")
    logger.info(f"Authorised channels: {[ch['name'] for ch in CHANNELS]} ({len(CHANNELS)})")

    client.run(DISCORD_BOT_TOKEN, log_handler=None)


if __name__ == "__main__":
    main()
