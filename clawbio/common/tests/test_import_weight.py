"""A reproducibility recipe must be able to replay the run it sits next to.

Skills that import only the standard library declare `pip_deps=[]`. If importing
clawbio.common.reproducibility drags numpy, pandas or opentelemetry in through
the package __init__, `conda env create` followed by `bash commands.sh` stops at
the import and the bundle is decorative.
"""
from __future__ import annotations

import subprocess
import sys

HEAVY = {"numpy", "pandas", "opentelemetry", "scipy", "matplotlib", "anndata", "scanpy"}


def test_reproducibility_import_stays_stdlib_only():
    code = """import sys, json
import clawbio.common.reproducibility
print(json.dumps(sorted({m.split('.')[0] for m in sys.modules})))"""
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    import json

    loaded = set(json.loads(out.stdout))
    assert not (loaded & HEAVY), f"heavy imports pulled in by the helper: {sorted(loaded & HEAVY)}"
