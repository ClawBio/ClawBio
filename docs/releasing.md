# Packaging & distribution

How ClawBio ships as `pip install clawbio` and `conda install -c bioconda clawbio`.

## What changed to make the package installable

Before this work, the wheel shipped only the thin `clawbio/` bridge package: no
engine, no skills, no console command, so `pip install clawbio` produced an
import error and ran nothing. The fix:

- The CLI engine moved from the root `clawbio.py` into the package at
  `clawbio/cli.py`, so it ships inside the wheel. Root `clawbio.py` is now a thin
  shim that keeps `python clawbio.py ...` and the `spec_from_file_location` test
  hooks working.
- A `clawbio` console command is wired via `[project.scripts]`
  (`clawbio = "clawbio.cli:main"`).
- `clawbio/cli.py` resolves its content root for both layouts: the repo root in a
  source checkout, the package directory in an installed wheel. Output and patient
  profiles are routed to the current working directory once installed, never into
  read-only package data.
- The build backend is `hatchling`. `hatch_build.py` force-includes a curated
  subset of `skills/` and `examples/` into the wheel under `clawbio/` **without
  moving anything on disk**, so the Claude Code plugin marketplace and every
  public `skills/...` URL stay intact.

### Wheel payload policy (`hatch_build.py`)

- Every skill's logic (`.md` / `.py` / `.yaml` / `.sh`) ships, so agents can read
  and call every skill (live count: `skill_count` in `skills/catalog.json`).
- Small data files (< 256 KB) ship; heavy demo data (gzipped genomes,
  expected-output images, big catalogs) does not.
- The headline skills in `HEADLINE_SKILLS` ship their full demo payload so the
  quick-start works offline immediately: `pharmgx-reporter`, `drug-photo`,
  `gwas-lookup`, `just-prs-mcp`.

Result: a wheel of a few MB (4.0 MB at 0.7.1, vs the broken 100 KB one), every
skill present, `clawbio run pharmgx --demo` working straight after install.

## Build locally

```bash
uv build                      # builds dist/clawbio-<v>.tar.gz and the wheel
shasum -a 256 dist/clawbio-*.tar.gz   # sha256 the bioconda recipe pins
```

Verify in a clean environment:

```bash
uv venv /tmp/cb && uv pip install --python /tmp/cb/bin/python dist/clawbio-*.whl
/tmp/cb/bin/clawbio run pharmgx --demo
```

## Publish to PyPI (automated: one tag push)

Releases are published automatically by `.github/workflows/publish.yml` when a
`v*` tag is pushed. It authenticates with PyPI via Trusted Publishing (OIDC), so
there is no API token or secret stored in the repo. PyPI release names and
versions are permanent and cannot be reused, so release deliberately.

To cut a release:

1. Bump the version in `clawbio/__init__.py` (single source of truth; the wheel
   reads it via `[tool.hatch.version]`). Keep `CITATION.cff` / `CHANGELOG.md` in
   step.
2. Commit and push to `main`.
3. Tag and push:
   ```bash
   git tag v0.5.3 && git push origin v0.5.3
   ```
   The workflow then verifies the tag matches `clawbio/__init__.py`, builds the
   sdist + wheel, runs `twine check`, and publishes to PyPI. Watch it under the
   repo's Actions tab.

### Trusted Publisher registration

The publisher is registered on PyPI under the `clawbio` project -> **Manage** ->
**Publishing**, as a GitHub publisher with exactly these values. If any of them
change (repo rename, workflow rename, environment rename), update them on PyPI
too, or publishing fails with `invalid-publisher`:

| Field | Value |
|-------|-------|
| Owner | `ClawBio` |
| Repository name | `ClawBio` |
| Workflow name | `publish.yml` |
| Environment name | `pypi` |

The `pypi` GitHub environment already exists: it accepts `v*` tags only and
needs a reviewer's approval before the publish job runs.

If a release fails with `invalid-publisher`, compare the table above with the
claims printed in the job log.

Once the first release has gone out over OIDC, delete the old `PYPI_API_TOKEN`
secret from the `pypi` environment and revoke the token on PyPI.

Trusted Publishing also turns on attestations: `pypa/gh-action-pypi-publish`
signs every file with Sigstore and uploads a PEP 740 attestation. These show as
"Provenance" on the PyPI file page. Token uploads get none.

### Manual fallback (if you ever need to publish without the Action)

```bash
rm -rf dist && uv build
uvx twine check dist/*
uvx twine upload dist/*          # needs a PyPI token in ~/.pypirc
git tag v0.5.3 && git push origin v0.5.3
```

## Bioconda

The recipe lives in bioconda-recipes, not here:
[`recipes/clawbio/meta.yaml`](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/clawbio/meta.yaml).
Maintainer: @camlloyd.

## Notes / open items

- The sdist is ~12 MB because it carries the full `skills/` tree (the wheel build
  hook reads it to assemble the curated wheel). If a leaner sdist is wanted later,
  add an sdist-target filter mirroring the wheel hook.
