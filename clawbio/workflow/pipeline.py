"""pipeline.yaml loader, validator, and topological sort."""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


PLACEHOLDER_RE = re.compile(r"\{([a-zA-Z_][a-zA-Z0-9_]*)\}")


class PipelineError(ValueError):
    pass


@dataclass
class InputSpec:
    name: str
    required: bool = True


@dataclass
class Step:
    id: str
    summary: str
    command_template: str
    requires_inputs: list[str] = field(default_factory=list)
    prereqs: list[str] = field(default_factory=list)
    requires: dict[str, Any] = field(default_factory=dict)
    emits: list[str] = field(default_factory=list)
    optional: bool = False


@dataclass
class Pipeline:
    skill: str
    version: int
    inputs: list[InputSpec]
    steps: list[Step]
    skill_dir: Path

    @property
    def step_ids(self) -> list[str]:
        return [s.id for s in self.steps]

    def step(self, step_id: str) -> Step:
        for s in self.steps:
            if s.id == step_id:
                return s
        raise PipelineError(f"unknown step id: {step_id}")


def load(skill_dir: Path) -> Pipeline:
    """Load <skill_dir>/pipeline.yaml, validate, and return a Pipeline."""
    skill_dir = Path(skill_dir).resolve()
    yaml_path = skill_dir / "pipeline.yaml"
    if not yaml_path.is_file():
        raise PipelineError(f"no pipeline.yaml at {yaml_path}")

    with yaml_path.open() as f:
        data = yaml.safe_load(f)

    if not isinstance(data, dict):
        raise PipelineError("pipeline.yaml must be a mapping")

    skill = data.get("skill")
    version = data.get("version")
    if not isinstance(skill, str) or not skill:
        raise PipelineError("`skill` must be a non-empty string")
    if version != 1:
        raise PipelineError(f"unsupported pipeline schema version: {version!r} (expected 1)")

    inputs = _parse_inputs(data.get("inputs", []))
    steps = _parse_steps(data.get("steps", []))

    _validate(steps, inputs)

    return Pipeline(
        skill=skill,
        version=version,
        inputs=inputs,
        steps=steps,
        skill_dir=skill_dir,
    )


def _parse_inputs(raw: Any) -> list[InputSpec]:
    if not isinstance(raw, list):
        raise PipelineError("`inputs` must be a list")
    out: list[InputSpec] = []
    for item in raw:
        if not isinstance(item, dict) or "name" not in item:
            raise PipelineError(f"bad input entry: {item!r}")
        out.append(InputSpec(name=item["name"], required=bool(item.get("required", True))))
    return out


def _parse_steps(raw: Any) -> list[Step]:
    if not isinstance(raw, list) or not raw:
        raise PipelineError("`steps` must be a non-empty list")
    out: list[Step] = []
    for item in raw:
        if not isinstance(item, dict):
            raise PipelineError(f"bad step entry: {item!r}")
        for required_key in ("id", "summary", "command_template"):
            if required_key not in item:
                raise PipelineError(f"step missing {required_key!r}: {item!r}")
        out.append(Step(
            id=item["id"],
            summary=item["summary"],
            command_template=item["command_template"],
            requires_inputs=list(item.get("requires_inputs", []) or []),
            prereqs=list(item.get("prereqs", []) or []),
            requires=dict(item.get("requires", {}) or {}),
            emits=list(item.get("emits", []) or []),
            optional=bool(item.get("optional", False)),
        ))
    return out


def _validate(steps: list[Step], inputs: list[InputSpec]) -> None:
    ids = [s.id for s in steps]
    if len(set(ids)) != len(ids):
        raise PipelineError(f"duplicate step ids: {ids}")

    input_names = {i.name for i in inputs}
    emits_so_far: set[str] = set()
    id_to_index = {sid: i for i, sid in enumerate(ids)}

    for step in steps:
        for inp in step.requires_inputs:
            if inp not in input_names:
                raise PipelineError(f"step {step.id!r} requires undeclared input {inp!r}")

        for prereq in step.prereqs:
            if prereq not in id_to_index:
                raise PipelineError(f"step {step.id!r} has unknown prereq {prereq!r}")
            if id_to_index[prereq] >= id_to_index[step.id]:
                raise PipelineError(
                    f"step {step.id!r} prereq {prereq!r} is not declared earlier — cycles forbidden"
                )

        for key in step.requires:
            if key not in emits_so_far:
                raise PipelineError(
                    f"step {step.id!r} requires {key!r} but no prior step emits it"
                )

        for placeholder in PLACEHOLDER_RE.findall(step.command_template):
            if placeholder in ("out", "skill_dir"):
                continue
            if placeholder in input_names:
                continue
            if placeholder in emits_so_far or placeholder in step.emits:
                continue
            raise PipelineError(
                f"step {step.id!r} command_template references unknown placeholder {{{placeholder}}}"
            )

        emits_so_far.update(step.emits)


def render_command(
    step: Step,
    inputs: dict[str, str],
    emitted: dict[str, str],
    out_dir: str,
    skill_dir: str,
) -> str:
    """Substitute {var} placeholders. Raises if any are unresolved."""
    ctx = {"out": out_dir, "skill_dir": skill_dir, **inputs, **emitted}

    def _sub(match: re.Match[str]) -> str:
        key = match.group(1)
        if key not in ctx:
            raise PipelineError(f"unresolved placeholder {{{key}}} for step {step.id!r}")
        return str(ctx[key])

    return PLACEHOLDER_RE.sub(_sub, step.command_template)
