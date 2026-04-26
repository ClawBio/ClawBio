"""Shared deterministic skill-intent planner for ClawBio chat adapters.

The planner reads optional ``INTENTS.json`` descriptors from skill directories
and turns user text plus a weak requested skill/mode hint into one or more
``clawbio.py run ...`` executions. It never imports or executes skill-local
code.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path
from typing import Any


SCHEMA = "clawbio.skill_intents.v1"
DESCRIPTOR_FILENAMES = ("INTENTS.json", "skill_intents.json")
CONFIDENCE_HIGH = "high"
CONFIDENCE_MEDIUM = "medium"
CONFIDENCE_LOW = "low"

_DEMO_TERMS = (
    "demo",
    "demonstration",
    "synthetic",
    "example data",
    "sample data",
    "test data",
)
_CONFIRM_TERMS = ("yes", "confirm", "confirmed", "go ahead", "proceed", "run it")


@dataclass
class SkillIntentExecution:
    """One planned command-line execution."""

    kind: str
    skill: str
    argv: list[str]
    output_dir: str | None = None
    input_path: str | None = None
    requires_confirmation: bool = False
    confirmation_reason: str | None = None
    route_step_id: str | None = None


@dataclass
class SkillExecutionPlan:
    """Structured result returned to chat adapters and other callers."""

    status: str
    raw_user_text: str
    raw_user_text_sha256: str
    skill: str | None = None
    intent_id: str | None = None
    confidence: str = CONFIDENCE_LOW
    reason: str = ""
    matched_route: dict[str, Any] | None = None
    executions: list[SkillIntentExecution] = field(default_factory=list)
    requested_skill: str | None = None
    requested_mode: str | None = None
    descriptor_path: str | None = None
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def load_default_skill_registry(project_root: str | Path | None = None) -> dict[str, dict[str, Any]]:
    """Load ``SKILLS`` from the repository's top-level ``clawbio.py`` script."""

    root = Path(project_root) if project_root else Path(__file__).resolve().parent.parent
    script = root / "clawbio.py"
    spec = importlib.util.spec_from_file_location("_clawbio_cli_registry", script)
    if spec is None or spec.loader is None:
        return {}
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return getattr(module, "SKILLS", {})


def plan_skill_intent(
    user_text: str,
    requested_skill: str | None,
    requested_mode: str | None,
    attachments: list | None,
    skill_registry: dict,
    project_root: str | Path | None = None,
) -> SkillExecutionPlan:
    """Plan one or more ClawBio skill executions from platform-neutral inputs.

    ``requested_skill`` and ``requested_mode`` are treated as hints, typically
    from an LLM tool call. The raw text wins when it strongly matches an
    intent descriptor, so adapters can recover from weak tool-call choices.
    """

    text = user_text or ""
    explicit_demo = _demo_allowed(text, requested_mode)
    effective_mode = _normalise_mode(requested_mode, explicit_demo)
    execution_root = Path(project_root) if project_root else _project_root_from_registry(skill_registry)
    descriptors = load_skill_intent_descriptors(skill_registry, execution_root)
    requested_skill = _resolve_skill_alias(requested_skill, skill_registry, descriptors)

    matches = _score_descriptor_routes(text, requested_skill, descriptors, explicit_demo)
    if matches:
        best = matches[0]
        descriptor, route, score, matched_terms = best
        return _plan_descriptor_route(
            text=text,
            requested_skill=requested_skill,
            requested_mode=requested_mode,
            explicit_demo=explicit_demo,
            descriptor=descriptor,
            route=route,
            score=score,
            matched_terms=matched_terms,
            project_root=execution_root,
            skill_registry=skill_registry,
        )

    return _plan_legacy_fallback(
        text=text,
        requested_skill=requested_skill,
        requested_mode=requested_mode,
        effective_mode=effective_mode,
        explicit_demo=explicit_demo,
        attachments=attachments or [],
        skill_registry=skill_registry,
        project_root=execution_root,
    )


def load_skill_intent_descriptors(
    skill_registry: dict,
    project_root: str | Path | None = None,
) -> list[dict[str, Any]]:
    """Return validated descriptors from registered scripts and ``skills/*`` dirs."""

    descriptors: list[dict[str, Any]] = []
    seen_paths: set[str] = set()
    for alias, info in skill_registry.items():
        skill_dir = _skill_dir(info)
        if not skill_dir:
            continue
        data = _read_descriptor(skill_dir, alias)
        if data:
            seen_paths.add(str(data["_descriptor_path"]))
            descriptors.append(data)

    root = Path(project_root) if project_root else _project_root_from_registry(skill_registry)
    skills_root = root / "skills"
    if skills_root.exists():
        for skill_dir in sorted(p for p in skills_root.iterdir() if p.is_dir()):
            data = _read_descriptor(skill_dir.resolve(), skill_dir.name)
            if not data or str(data["_descriptor_path"]) in seen_paths:
                continue
            seen_paths.add(str(data["_descriptor_path"]))
            descriptors.append(data)
    return descriptors


def skill_names_for_tool_schema(
    skill_registry: dict,
    project_root: str | Path | None = None,
) -> list[str]:
    """Return registry and descriptor skill names suitable for chat tool enums."""

    names = set(skill_registry.keys())
    for descriptor in load_skill_intent_descriptors(skill_registry, project_root):
        if descriptor.get("skill") and _descriptor_has_executable_route(descriptor, skill_registry):
            names.add(str(descriptor["skill"]))
    names.add("auto")
    return sorted(names)


def skill_intent_tool_summary(
    skill_registry: dict,
    project_root: str | Path | None = None,
) -> str:
    """Compact human-readable descriptor route summary for LLM tool descriptions."""

    summaries = []
    for descriptor in load_skill_intent_descriptors(skill_registry, project_root):
        if not _descriptor_has_executable_route(descriptor, skill_registry):
            continue
        skill = descriptor.get("skill") or descriptor.get("_registry_alias")
        intents = [
            str(route.get("intent_id"))
            for route in descriptor.get("routes", [])
            if isinstance(route, dict) and route.get("intent_id")
        ]
        if intents:
            summaries.append(f"{skill} intents: {', '.join(intents)}")
    return "; ".join(summaries)


def _read_descriptor(skill_dir: Path, alias: str) -> dict[str, Any] | None:
    for filename in DESCRIPTOR_FILENAMES:
        path = skill_dir / filename
        if not path.exists():
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if data.get("schema") != SCHEMA:
            continue
        routes = data.get("routes")
        if not isinstance(routes, list):
            continue
        skill_name = str(data.get("skill") or skill_dir.name)
        data["_descriptor_path"] = str(path)
        data["_skill_dir"] = str(skill_dir)
        data["_registry_alias"] = alias
        data["_skill_name"] = skill_name
        return data
    return None


def _descriptor_has_executable_route(descriptor: dict[str, Any], skill_registry: dict) -> bool:
    descriptor_skill = str(descriptor.get("skill") or descriptor.get("_registry_alias"))
    for route in descriptor.get("routes", []):
        if not isinstance(route, dict):
            continue
        for step in route.get("plan") or []:
            if not isinstance(step, dict) or step.get("kind", "skill_run") != "skill_run":
                continue
            step_skill = str(step.get("skill") or descriptor_skill)
            if step_skill in skill_registry:
                return True
    return False


def _plan_descriptor_route(
    text: str,
    requested_skill: str | None,
    requested_mode: str | None,
    explicit_demo: bool,
    descriptor: dict[str, Any],
    route: dict[str, Any],
    score: int,
    matched_terms: list[str],
    project_root: Path,
    skill_registry: dict,
) -> SkillExecutionPlan:
    descriptor_skill = str(descriptor.get("skill") or descriptor.get("_registry_alias"))
    intent_id = str(route.get("intent_id") or "default")
    reason = (
        f"Matched descriptor route '{intent_id}' for skill '{descriptor_skill}' "
        f"using terms: {', '.join(matched_terms) or 'skill hint'}."
    )
    confidence = CONFIDENCE_HIGH if score >= 7 else CONFIDENCE_MEDIUM
    skill_dir = Path(str(descriptor["_skill_dir"]))
    executions: list[SkillIntentExecution] = []
    warnings: list[str] = []
    missing_skills: set[str] = set()

    for index, step in enumerate(route.get("plan") or []):
        if not isinstance(step, dict):
            warnings.append(f"Skipped non-object plan step at index {index}.")
            continue
        if step.get("kind", "skill_run") != "skill_run":
            warnings.append(f"Skipped unsupported plan step kind at index {index}.")
            continue
        step_demo = bool(step.get("demo", False))
        if step_demo and not explicit_demo:
            warnings.append(f"Skipped demo step {index}; user did not request a demo.")
            continue
        step_skill = str(step.get("skill") or descriptor_skill)
        if step_skill not in skill_registry:
            missing_skills.add(step_skill)
            continue
        argv = [sys.executable, str(project_root / "clawbio.py"), "run", step_skill]
        input_path = _resolve_descriptor_input(step.get("input"), skill_dir)
        if step_demo:
            argv.append("--demo")
        elif input_path:
            argv.extend(["--input", str(input_path)])
        for arg in _safe_argv_list(step.get("args")):
            argv.append(arg)
        output_dir = None
        if step.get("output"):
            output_dir = str(_resolve_descriptor_input(step.get("output"), skill_dir))
            argv.extend(["--output", output_dir])
        confirmation = step.get("confirmation") or {}
        requires_confirmation = bool(
            step.get("requires_confirmation")
            or route.get("requires_confirmation")
            or (isinstance(confirmation, dict) and confirmation.get("required"))
        )
        confirmation_reason = None
        if isinstance(confirmation, dict):
            confirmation_reason = confirmation.get("reason")
        executions.append(
            SkillIntentExecution(
                kind="skill_run",
                skill=step_skill,
                argv=argv,
                output_dir=output_dir,
                input_path=str(input_path) if input_path else None,
                requires_confirmation=requires_confirmation,
                confirmation_reason=confirmation_reason,
                route_step_id=str(step.get("id") or index),
            )
        )

    if missing_skills:
        missing = ", ".join(sorted(missing_skills))
        return SkillExecutionPlan(
            status="needs_registration",
            raw_user_text=text,
            raw_user_text_sha256=_sha(text),
            skill=descriptor_skill,
            intent_id=intent_id,
            confidence=confidence,
            reason=(
                f"Matched descriptor route '{intent_id}', but skill(s) {missing} "
                "are not registered in clawbio.py SKILLS yet."
            ),
            matched_route={
                "intent_id": intent_id,
                "description": route.get("description", ""),
                "matched_terms": matched_terms,
                "score": score,
                "demo_policy": route.get("demo_policy", "never_unless_explicit"),
            },
            executions=[],
            requested_skill=requested_skill,
            requested_mode=requested_mode,
            descriptor_path=descriptor.get("_descriptor_path"),
            warnings=[*warnings, f"Register {missing} before exposing it for execution."],
        )

    status = "planned"
    if any(item.requires_confirmation for item in executions) and not _contains_any(text, _CONFIRM_TERMS):
        status = "needs_confirmation"

    matched_route = {
        "intent_id": intent_id,
        "description": route.get("description", ""),
        "matched_terms": matched_terms,
        "score": score,
        "demo_policy": route.get("demo_policy", "never_unless_explicit"),
    }
    return SkillExecutionPlan(
        status=status,
        raw_user_text=text,
        raw_user_text_sha256=_sha(text),
        skill=descriptor_skill,
        intent_id=intent_id,
        confidence=confidence,
        reason=reason,
        matched_route=matched_route,
        executions=executions,
        requested_skill=requested_skill,
        requested_mode=requested_mode,
        descriptor_path=descriptor.get("_descriptor_path"),
        warnings=warnings,
    )


def _plan_legacy_fallback(
    text: str,
    requested_skill: str | None,
    requested_mode: str | None,
    effective_mode: str | None,
    explicit_demo: bool,
    attachments: list,
    skill_registry: dict,
    project_root: Path,
) -> SkillExecutionPlan:
    skill = requested_skill or _infer_legacy_skill(text) or "auto"
    input_path, profile_path = _attachment_paths(attachments)
    extra_args: list[str] = []
    if skill == "prs":
        trait = _extract_attachment_value(attachments, "trait")
        if trait:
            extra_args.extend(["--trait", trait])
    elif skill == "clinpgx":
        gene = _extract_attachment_value(attachments, "gene")
        if gene:
            extra_args.extend(["--gene", gene])
    elif skill == "gwas":
        rsid = _extract_attachment_value(attachments, "rsid")
        if rsid:
            extra_args.extend(["--rsid", rsid])
    elif skill == "drugphoto":
        drug = _extract_attachment_value(attachments, "drug_name")
        dose = _extract_attachment_value(attachments, "visible_dose")
        if drug:
            extra_args.extend(["--drug", drug])
        if dose:
            extra_args.extend(["--dose", dose])

    argv = [sys.executable, str(project_root / "clawbio.py"), "run", skill]
    if (explicit_demo and effective_mode == "demo") or (skill == "drugphoto" and requested_mode == "demo"):
        argv.append("--demo")
    elif skill in ("profile", "prs") and profile_path:
        argv.extend(["--profile", profile_path])
    elif input_path:
        argv.extend(["--input", input_path])
    elif skill in ("clinpgx", "gwas") and extra_args:
        pass
    elif skill != "auto":
        # Deterministic fallback: do not silently switch to demo for weak tool calls.
        return SkillExecutionPlan(
            status="needs_input",
            raw_user_text=text,
            raw_user_text_sha256=_sha(text),
            skill=skill,
            intent_id="legacy_fallback",
            confidence=CONFIDENCE_LOW,
            reason="No intent descriptor matched and no input/profile was available.",
            requested_skill=requested_skill,
            requested_mode=requested_mode,
            warnings=["Demo mode is only planned when the user explicitly asks for a demo."],
        )
    argv.extend(extra_args)
    return SkillExecutionPlan(
        status="planned",
        raw_user_text=text,
        raw_user_text_sha256=_sha(text),
        skill=skill,
        intent_id="legacy_fallback",
        confidence=CONFIDENCE_MEDIUM if requested_skill else CONFIDENCE_LOW,
        reason="No skill intent descriptor matched; using the legacy requested skill and mode.",
        matched_route={"intent_id": "legacy_fallback", "matched_terms": [], "score": 0},
        executions=[
            SkillIntentExecution(
                kind="skill_run",
                skill=skill,
                argv=argv,
                output_dir=None,
                input_path=input_path,
            )
        ],
        requested_skill=requested_skill,
        requested_mode=requested_mode,
        warnings=[] if explicit_demo or requested_mode != "demo" else [
            "Ignored weak demo mode because the user text did not explicitly request a demo."
        ],
    )


def _score_descriptor_routes(
    text: str,
    requested_skill: str | None,
    descriptors: list[dict[str, Any]],
    explicit_demo: bool,
) -> list[tuple[dict[str, Any], dict[str, Any], int, list[str]]]:
    norm = _normalise_text(text)
    matches: list[tuple[dict[str, Any], dict[str, Any], int, list[str]]] = []
    for descriptor in descriptors:
        descriptor_skill = str(descriptor.get("skill") or descriptor.get("_registry_alias"))
        skill_aliases = [descriptor_skill, descriptor.get("_registry_alias"), *descriptor.get("aliases", [])]
        skill_hint = requested_skill in skill_aliases if requested_skill else False
        for route in descriptor.get("routes", []):
            if not isinstance(route, dict):
                continue
            demo_policy = route.get("demo_policy", "never_unless_explicit")
            if demo_policy == "only_when_explicit" and not explicit_demo:
                continue
            terms = [
                *route.get("trigger_terms", []),
                *route.get("aliases", []),
            ]
            matched_terms = [term for term in terms if _term_matches(norm, str(term))]
            score = len(matched_terms) * 4
            if skill_hint:
                score += 3
            if _contains_any(norm, [str(term).lower() for term in skill_aliases if term]):
                score += 2
            if route.get("intent_id") and _term_matches(norm, str(route["intent_id"]).replace("_", " ")):
                score += 2
            if score >= 4:
                matches.append((descriptor, route, score, matched_terms))
    matches.sort(key=lambda item: item[2], reverse=True)
    return matches


def _normalise_mode(requested_mode: str | None, explicit_demo: bool) -> str | None:
    if requested_mode == "demo" and explicit_demo:
        return "demo"
    if requested_mode == "file":
        return "file"
    return None


def _demo_allowed(text: str, requested_mode: str | None) -> bool:
    if _contains_any(text, _DEMO_TERMS):
        return True
    return requested_mode == "demo" and _contains_any(text, _CONFIRM_TERMS)


def _infer_legacy_skill(text: str) -> str | None:
    norm = _normalise_text(text)
    legacy_terms = [
        ("prs", ("polygenic", "risk score", "disease risk", "at risk")),
        ("profile", ("profile report", "full profile", "unified profile")),
        ("clinpgx", ("gene drug", "cpic", "pharmgkb")),
        ("gwas", ("rsid", "rs number", "variant lookup", "look up rs")),
        ("pharmgx", ("pharmacogen", "drug response", "pgx")),
        ("nutrigx", ("nutrition", "diet", "caffeine", "lactose")),
        ("compare", ("compare genome", "genome comparison", "ibs")),
        ("metagenomics", ("metagenomic", "fastq", "microbiome")),
    ]
    for skill, terms in legacy_terms:
        if _contains_any(norm, terms):
            return skill
    return None


def _resolve_skill_alias(
    requested_skill: str | None,
    skill_registry: dict,
    descriptors: list[dict[str, Any]],
) -> str | None:
    if not requested_skill:
        return None
    if requested_skill in skill_registry:
        return requested_skill
    for descriptor in descriptors:
        aliases = [descriptor.get("skill"), descriptor.get("_registry_alias"), *descriptor.get("aliases", [])]
        if requested_skill in aliases:
            return str(descriptor.get("skill") or descriptor.get("_registry_alias"))
    return requested_skill


def _project_root_from_registry(skill_registry: dict) -> Path:
    for info in skill_registry.values():
        skill_dir = _registry_skill_dir(info)
        if skill_dir:
            return skill_dir.parent.parent
    return Path(__file__).resolve().parent.parent


def _registry_skill_dir(info: dict[str, Any]) -> Path | None:
    script = info.get("script") if isinstance(info, dict) else None
    if not script:
        return None
    return Path(script).parent


def _skill_dir(info: dict[str, Any]) -> Path | None:
    script = info.get("script") if isinstance(info, dict) else None
    if not script:
        return None
    return Path(script).resolve().parent


def _resolve_descriptor_input(value: Any, skill_dir: Path) -> Path | None:
    if not value:
        return None
    raw = Path(str(value))
    if raw.is_absolute():
        return raw
    return (skill_dir / raw).resolve()


def _safe_argv_list(value: Any) -> list[str]:
    if not isinstance(value, list):
        return []
    safe = []
    for item in value:
        text = str(item)
        if "\x00" in text or "\n" in text or "\r" in text:
            continue
        safe.append(text)
    return safe


def _attachment_paths(attachments: list) -> tuple[str | None, str | None]:
    input_path = None
    profile_path = None
    for item in attachments:
        if not isinstance(item, dict):
            continue
        input_path = input_path or item.get("path") or item.get("input_path")
        profile_path = profile_path or item.get("profile_path")
    return input_path, profile_path


def _extract_attachment_value(attachments: list, key: str) -> str | None:
    for item in attachments:
        if isinstance(item, dict) and item.get(key):
            return str(item[key])
    return None


def _normalise_text(text: str) -> str:
    return re.sub(r"\s+", " ", text.lower()).strip()


def _term_matches(normalised_text: str, term: str) -> bool:
    term_norm = _normalise_text(term)
    if not term_norm:
        return False
    if " " in term_norm:
        return term_norm in normalised_text
    return re.search(rf"\b{re.escape(term_norm)}\b", normalised_text) is not None


def _contains_any(text: str, terms: tuple[str, ...] | list[str]) -> bool:
    norm = _normalise_text(text)
    return any(_term_matches(norm, term) for term in terms)


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()
