"""Structured failures shared by nf-core wrappers."""

from __future__ import annotations

import math

from dataclasses import dataclass, field
from typing import Any


@dataclass
class SkillError(Exception):
    stage: str
    error_code: str
    message: str
    fix: str
    details: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        super().__init__(self.message)

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": "error",
            "ok": False,
            "stage": self.stage,
            "error_code": self.error_code,
            "message": self.message,
            "fix": self.fix,
            "details": json_safe_details(self.details),
        }


def json_safe_details(value):
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {key: json_safe_details(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe_details(item) for item in value]
    return value
