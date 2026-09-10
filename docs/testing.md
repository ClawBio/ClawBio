# Test Verification

This documents a full test run of the repository after a clean `uv sync`
install, and the methodology needed to run the suite correctly.

## Known issue: don't run a single bare `pytest` / `uv run pytest`

Running the entire suite in one process (`uv run pytest`, no path argument)
fails at collection time with:

```
ValueError: Plugin already registered under a different name:
.../skills/gwas-catalog-region-fetch/tests/conftest.py=<module 'tests.conftest'
from '.../skills/eqtl-catalogue-region-fetch/tests/conftest.py'>
```

Cause: `pytest.ini` sets `--import-mode=importlib`. Each skill's `tests/`
directory has its own `tests/__init__.py`, but the skill directory itself
(e.g. `skills/gwas-catalog-region-fetch/`) does not. Without a package
boundary above `tests/`, pytest's import-mode resolves every skill's test
package to the same module name, `tests.conftest`, and any two skills whose
test folders both load in the same process collide.

This is why `.github/workflows/ci.yml` never invokes pytest once across the
whole repo — it runs pytest **separately per test path**. Do the same
locally:

```bash
# Good: one path per invocation (this is what CI does)
uv run pytest skills/pharmgx-reporter/tests/ -v
uv run pytest bot/tests/test_security.py -v

# Or loop over every path declared in pytest.ini's testpaths:
awk '/^testpaths/{f=1;next}/^[a-z_]+ =/{f=0}f' pytest.ini | \
  while read -r p; do uv run pytest "$p" -q; done
```

A proper fix would add an `__init__.py` to each skill directory (or an
equivalent `rootdir`-scoped package marker) so importlib mode can generate
unique module names; that has not been done as of this writing.

## Verification run (2026-09-04, after `uv sync`)

Ran all 79 `testpaths` entries individually, per the workaround above.

**Result: 4,290 passed · 55 failed · 19 errors · 98 skipped.** 67 of 79
paths were fully clean. 0 paths failed to collect.

### Failures caused by missing optional dependencies (9 fully, 1 partially, of 13 failing paths)

These are not installation defects. `uv sync` only installs the core
dependency set declared in `pyproject.toml`; each skill with heavier or
niche dependencies declares them separately in its own `SKILL.md` `install:`
block, and they are opt-in. Install a skill's declared dependencies only if
you intend to use that skill:

| Suite | Missing module(s) | Install |
|---|---|---|
| `robotary/tests` | `fastapi` | **Undocumented** — not declared in any SKILL.md, requirements file, or `pyproject.toml`. `robotary/server.py` imports it directly. Until this is fixed upstream, install manually: `uv add fastapi` (or `pip install fastapi` in your environment). |
| `skills/affinity-proteomics/tests` | `statsmodels` | Declared in `skills/affinity-proteomics/SKILL.md`: `pip install somadata scipy statsmodels seaborn scikit-learn` |
| `skills/cell-detection/tests` | `skimage` (scikit-image), `tifffile`, `cellpose` | Declared in `skills/cell-detection/SKILL.md`: `pip install "cellpose>=4.0" tifffile "czifile>=2019.7.2.2" "nd2>=0.11.1" Pillow scikit-image` |
| `skills/data-extractor/tests` | `cv2` | Declared in `skills/data-extractor/SKILL.md` as `opencv-python-headless`: `pip install anthropic opencv-python-headless` |
| `skills/eqtl-catalogue-region-fetch/tests` | `pysam` | Declared in `skills/eqtl-catalogue-region-fetch/SKILL.md`: `pip install pysam pandas requests` |
| `skills/gwas-catalog-region-fetch/tests` | `pysam` | Declared in `skills/gwas-catalog-region-fetch/SKILL.md`: `pip install pysam pandas requests` |
| `skills/variant-annotation/tests` | `pysam` | Declared in `skills/variant-annotation/SKILL.md`: `uv add pysam requests` |
| `skills/celltype-specificity-profiler/tests` | `scanpy` | Declared in `skills/celltype-specificity-profiler/SKILL.md`: `uv add scanpy anndata numpy scipy pandas` |
| `skills/drug-repurposing-screen/tests` | `pyarrow` (parquet engine) | Declared in `skills/drug-repurposing-screen/SKILL.md`: `pip install numpy pandas scipy pyyaml pyarrow` |
| `skills/proteomics-clock/tests` (partial — see below) | `seaborn` | Declared in `skills/proteomics-clock/SKILL.md`: `pip install pandas numpy matplotlib seaborn requests` |

### Failures unrelated to missing dependencies (3 fully, 1 partially, of 13 failing paths)

These reproduced with all relevant imports present and look like genuine
test/code issues, not environment gaps. Not yet root-caused further:

- **`clawbio/common/tests`** — `test_generate_report_header_requires_skill_version`
  failed with `Failed: DID NOT RAISE TypeError`. The code path no longer
  raises where the test expects it to.
- **`skills/bgpt-mcp/tests`** — several `KeyError`/`AssertionError` failures
  (`'endpoints'`, `'version'`, `'inputs'` missing) indicate the skill's
  metadata file has drifted from the schema its own tests expect.
- **`skills/bigquery-public/tests`** — `test_runner_security_filter_preserves_allowed_value_starting_with_dash`
  failed with `AttributeError: module 'clawbio_runner' has no attribute 'subprocess'`.
- **`skills/proteomics-clock/tests`** — independent of the `seaborn` import
  failure in the same file, one numeric assertion
  (`assert 5 < np.float64(3.8125...)`) fails on its own. This suite is
  counted in both tables above and here.

### Clean paths

The remaining 67 of 79 test paths, including `pharmgx-reporter`, passed
with no failures or errors using only the core `uv sync` install.
