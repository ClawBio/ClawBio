"""CLI entry for clawpathy-autoresearch.

Usage:
    python -m skills.clawpathy_autoresearch.autoresearch \\
        --workspace path/to/workspace [--dry-run]

The loop requires a subagent dispatcher. The CLI does not provide a real
one — the outer agent (or a future SDK integration) must wire it. With
--dry-run, a mock dispatcher is used that produces stub responses; useful
for validating workspace layout.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .dispatcher import CallableDispatcher, DispatchRequest
from .loop import run_loop


def _dry_run_dispatcher() -> CallableDispatcher:
    def fake(req: DispatchRequest) -> str:
        if req.role == "proposer":
            return "```markdown\n# stub skill\n```"
        return '```json\n{}\n```'
    return CallableDispatcher(fake)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    if args.dry_run:
        dispatcher = _dry_run_dispatcher()
    else:
        print(
            "error: no dispatcher configured. Use --dry-run, or invoke "
            "run_loop() directly with a CallableDispatcher provided by the "
            "outer agent.",
            file=sys.stderr,
        )
        return 2

    result = run_loop(workspace_dir=args.workspace, dispatcher=dispatcher)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
