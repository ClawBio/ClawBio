#!/usr/bin/env python3
"""A stand-in for the `conda` executable, for tests only.

`snakemake --use-conda` would call this to build tool environments; the tests
only need it on PATH, plus a `--version` that check_conda accepts.
"""

import sys

if "--version" in sys.argv[1:]:
    print("conda 25.3.0")
    sys.exit(0)

print("stub conda: nothing to do", file=sys.stderr)
sys.exit(0)
