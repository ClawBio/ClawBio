"""CLI: python -m clawpathy_autoresearch <workspace> [--resume]"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .dispatcher import ClaudeCLIDispatcher
from .loop import run_loop
from .loop_parallel import run_parallel_loop


def main() -> None:
    p = argparse.ArgumentParser(prog="clawpathy-autoresearch")
    p.add_argument("workspace", type=Path)
    p.add_argument("--fresh", dest="resume", action="store_false", help="Wipe history.jsonl and start over (default: resume)")
    p.set_defaults(resume=True)
    p.add_argument("--proposer-model", default="sonnet")
    p.add_argument("--executor-model", default="sonnet")
    p.add_argument("--judge-model", default="opus")
    p.add_argument("--parallel", type=int, default=0, metavar="K",
                   help="Run K proposals in parallel per iteration (0 = serial loop)")
    args = p.parse_args()

    dispatcher = ClaudeCLIDispatcher()
    common = dict(
        workspace_root=args.workspace,
        dispatcher=dispatcher,
        proposer_model=args.proposer_model,
        executor_model=args.executor_model,
        judge_model=args.judge_model,
        resume=args.resume,
    )
    if args.parallel and args.parallel > 1:
        result = run_parallel_loop(K=args.parallel, **common)
    else:
        result = run_loop(**common)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
