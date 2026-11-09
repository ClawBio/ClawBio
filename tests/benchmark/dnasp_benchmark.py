"""Time the demo analysis and Ka/Ks calls, without timing imports or parsing.

Run from the repository root:
    uv run python tests/benchmark/dnasp_benchmark.py

To compare commits, run this same script with --repo pointing to each checkout.
Use the same Python environment and repeat count for both. Each case warms up
with timeit.autorange(), then reports the median time per call. No speed
threshold is enforced because timings depend on the machine.
"""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import json
import platform
import statistics
import sys
import tempfile
import timeit
from functools import partial
from pathlib import Path


def load_dnasp(repo: Path):
    spec = importlib.util.spec_from_file_location(
        "_dnasp_benchmark_skill", repo / "skills" / "dnasp" / "dnasp.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def make_cases(dn):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "demo.fas"
        path.write_text(dn.DEMO_FASTA, encoding="utf-8")
        aln = dn.parse_fasta(path)
    outgroup_index = aln.names.index("outgroup")
    ingroup = [seq for i, seq in enumerate(aln.seqs) if i != outgroup_index]
    ingroup_aln = dn.Alignment(
        names=[name for i, name in enumerate(aln.names) if i != outgroup_index],
        seqs=ingroup,
    )
    return {
        "demo_analysis": partial(
            dn.run_analysis, aln, analyses=dn.VALID_ANALYSES - {"divergence"}
        ),
        "demo_analysis_outgroup": partial(
            dn.run_analysis,
            ingroup_aln,
            analyses=dn.VALID_ANALYSES - {"divergence"},
            outgroup=aln.seqs[outgroup_index],
        ),
        "demo_kaks_all_pairs": partial(dn.compute_ka_ks, aln.seqs),
        "demo_kaks_outgroup": partial(
            dn.compute_ka_ks, ingroup, outgroup=aln.seqs[outgroup_index]
        ),
        "single_sequence_early_return": partial(dn.compute_ka_ks, [aln.seqs[0]]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="checkout containing skills/dnasp/dnasp.py",
    )
    parser.add_argument("--repeat", type=int, default=7)
    args = parser.parse_args()
    if args.repeat < 1:
        parser.error("--repeat must be positive")
    dn = load_dnasp(args.repo.resolve())
    results = {}
    for name, run in make_cases(dn).items():
        # run_analysis emits diagnostics; give each call a fresh buffer.
        def quiet_run(run=run):
            with (
                contextlib.redirect_stdout(io.StringIO()),
                contextlib.redirect_stderr(io.StringIO()),
            ):
                return run()

        timer = timeit.Timer(quiet_run if name.startswith("demo_analysis") else run)
        number, _ = timer.autorange()
        samples = timer.repeat(repeat=args.repeat, number=number)
        results[name] = {
            "median_ms": statistics.median(samples) * 1000 / number,
            "calls_per_sample": number,
        }
    print(
        json.dumps(
            {
                "python": platform.python_version(),
                "repeat": args.repeat,
                "results": results,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
