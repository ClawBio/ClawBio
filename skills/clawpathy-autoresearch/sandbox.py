"""Sandbox spec loader for clawpathy-autoresearch.

Tasks declare executor constraints in sandbox.yaml: which tools the executor
may call, whether it has network, which hosts it may reach, what packages
are guaranteed importable, and a timeout.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml


@dataclass
class Sandbox:
    allowed_tools: list[str]
    network_allowed: bool = False
    network_hosts: list[str] = field(default_factory=list)
    data_dir: str = "./data"
    python_packages: list[str] = field(default_factory=list)
    timeout_seconds: int = 600

    def as_contract_text(self) -> str:
        """Render as a YAML-like block suitable for inclusion in a prompt."""
        return (
            f"allowed_tools: {self.allowed_tools}\n"
            f"network:\n"
            f"  allowed: {self.network_allowed}\n"
            f"  hosts: {self.network_hosts}\n"
            f"data_dir: {self.data_dir}\n"
            f"python_packages: {self.python_packages}\n"
            f"timeout_seconds: {self.timeout_seconds}\n"
        )


def load_sandbox(path: Path) -> Sandbox:
    raw = yaml.safe_load(Path(path).read_text()) or {}
    if "allowed_tools" not in raw or not raw["allowed_tools"]:
        raise ValueError(f"sandbox.yaml missing 'allowed_tools': {path}")
    net = raw.get("network") or {}
    allowed = net.get("allowed", False)
    if isinstance(allowed, list):
        network_allowed = True
        network_hosts = list(allowed)
    else:
        network_allowed = bool(allowed)
        network_hosts = []
    return Sandbox(
        allowed_tools=list(raw["allowed_tools"]),
        network_allowed=network_allowed,
        network_hosts=network_hosts,
        data_dir=raw.get("data_dir", "./data"),
        python_packages=list(raw.get("python_packages") or []),
        timeout_seconds=int(raw.get("timeout_seconds", 600)),
    )
