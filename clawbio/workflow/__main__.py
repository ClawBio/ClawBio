"""argparse dispatch for the clawbio-flow CLI."""
from __future__ import annotations

import argparse
import sys

from clawbio.workflow import bundle, commands
from clawbio.workflow.commands import CommandError
from clawbio.workflow.pipeline import PipelineError


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="clawbio-flow", description="Stateful workflow CLI for ClawBio skills.")
    sub = p.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("start", help="start a new run from a skill directory")
    sp.add_argument("skill_path", help="path to a directory containing pipeline.yaml")
    sp.add_argument("--input", action="append", default=[], metavar="KEY=VALUE",
                    help="input value (repeatable)")

    ss = sub.add_parser("status", help="show run status")
    ss.add_argument("--check-bundle", action="store_true",
                    help="exit 2 if a finished run hasn't been bundled (for Stop hooks)")

    sub.add_parser("next", help="print the next ready step id")

    so = sub.add_parser("open", help="open a step (returns command, marks in_progress)")
    so.add_argument("step_id")

    sd = sub.add_parser("done", help="mark a step done with its outputs")
    sd.add_argument("step_id")
    sd.add_argument("--output", action="append", default=[], metavar="KEY=PATH",
                    help="output key=path (repeatable, must match step `emits`)")
    sd.add_argument("--exit-code", type=int, default=0)

    sk = sub.add_parser("skip", help="skip an optional step (only valid when pending)")
    sk.add_argument("step_id")
    sk.add_argument("--reason", required=True, help="why this optional step was skipped")

    sub.add_parser("bundle", help="emit reproducibility bundle (commands.sh, checksums, env)")

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        if args.cmd == "start":
            commands.cmd_start(args.skill_path, args.input)
        elif args.cmd == "status":
            commands.cmd_status(check_bundle=args.check_bundle)
        elif args.cmd == "next":
            commands.cmd_next()
        elif args.cmd == "open":
            commands.cmd_open(args.step_id)
        elif args.cmd == "done":
            commands.cmd_done(args.step_id, args.output, exit_code=args.exit_code)
        elif args.cmd == "skip":
            commands.cmd_skip(args.step_id, args.reason)
        elif args.cmd == "bundle":
            bundle.cmd_bundle()
        else:
            parser.error(f"unknown command {args.cmd!r}")
    except (CommandError, PipelineError, FileNotFoundError, FileExistsError, RuntimeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
