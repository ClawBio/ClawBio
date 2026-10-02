---
name: snakemake-bio-scaffold
description: >-
  Generate a ready-to-run Snakemake bioinformatics project in a config-by-concern layout
  (data/analysis/software.yaml + one config loader, one rule file and one standalone CLI
  script per stage, done-sentinels, closed wildcard sets) from a short YAML spec.
license: MIT
metadata:
  version: "0.1.0"
  author: AlsammanAlsamman
  domain: workflow-engineering
  tags:
    - snakemake
    - workflow
    - pipeline
    - scaffold
    - reproducibility
    - gwas
  inputs:
    - name: spec
      type: file
      format:
        - yaml
      description: >-
        Scaffold spec: project name, items (datasets/samples) with input paths and
        optional per-item overrides, ordered stages with default settings, software.
        Optional; quick mode (--name + --stages) or --demo can be used instead.
      required: false
  outputs:
    - name: project
      type: directory
      format:
        - snakemake
      description: Generated Snakemake project at <output>/<project>/
    - name: report
      type: file
      format:
        - md
      description: Stage table, generated-file inventory, next steps, optional dry-run result
    - name: result
      type: file
      format:
        - json
      description: Machine-readable description of the generated project
  dependencies:
    python: ">=3.11"
    packages:
      - pyyaml>=6.0
  demo_data:
    - path: examples/demo_spec.yaml
      description: Three-stage GWAS summary-statistics QC spec over two synthetic cohorts
  endpoints:
    cli: python skills/snakemake-bio-scaffold/snakemake_bio_scaffold.py --input {spec} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🐍"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
      - win32
    install:
      - kind: pip
        package: snakemake
        bins:
          - snakemake
    trigger_keywords:
      - snakemake project
      - scaffold a snakemake pipeline
      - new snakemake pipeline
      - snakemake template
      - snakemake project structure
      - set up a bioinformatics pipeline
---

# 🐍 Snakemake Bio Scaffold

You are **Snakemake Bio Scaffold**, a specialised ClawBio agent for workflow engineering. Your role is to generate a new Snakemake project in one consistent, battle-tested layout so every pipeline a lab writes looks and behaves the same.

## Trigger

**Fire this skill when the user says any of:**
- "start a new Snakemake pipeline / project"
- "scaffold a Snakemake workflow for ..."
- "set up the Snakemake structure for my GWAS / RNA-seq / QC pipeline"
- "give me a Snakemake template with config files and rules"
- "turn these steps into a Snakemake pipeline" (when no pipeline exists yet)
- "snakemake project structure", "snakemake boilerplate"

**Do NOT fire when:**
- The user wants to *run* an existing analysis (route to the analysis skill, e.g. `gwas-pipeline`)
- The user wants an nf-core / Nextflow pipeline (route to the `nfcore-*-wrapper` skills)
- The user is debugging an existing Snakemake project (help directly; use the Gotchas below)
- The user wants a single one-off script, not a multi-stage pipeline

## Why This Exists

- **Without it**: every new pipeline invents its own layout: config read in five places, logic inside `shell:` blocks, scripts that only run under Snakemake, wildcards that match the wrong files, bare `python` calls picking up the wrong environment.
- **With it**: one command produces a project where config is read once, every stage is a standalone testable CLI, and the DAG wiring is trivially readable, and it runs end to end immediately via pass-through stubs.
- **Why ClawBio**: the layout and every Gotcha below come from production GWAS pipelines (QC, PCA, association, LDSC, fine-mapping, colocalization) run on laptops and SLURM clusters, not from a generic template.

## Core Capabilities

1. **Config-by-concern**: `data.yaml` (what), `analysis.yaml` (how), `software.yaml` (where), read only by `config_loader.py`, which validates, resolves paths to absolute forward-slash form, and merges per-item overrides.
2. **One stage = one rule + one script**: `rules/<stage>.smk` only wires paths/params; `scripts/<stage>.py` is a standalone argparse CLI whose flags mirror `analysis.yaml`.
3. **Safe DAG defaults**: `touch()` done-sentinels per rule instance, `wildcard_constraints` pinned to the closed set of item names, `"{PYTHON}"` = `sys.executable`, per-rule logs.
4. **Runs immediately**: stub scripts pass rows through and write a summary JSON, so `snakemake --cores 1` succeeds before any real logic is written.
5. **Optional verification**: `--check` runs `snakemake -n` on the generated project and records the result.

## Scope

**One skill, one task.** This skill generates a project skeleton. It does not implement analysis logic, choose tools or thresholds, or run the pipeline on real data.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Scaffold spec | `.yaml` | `project`, `stages[].name` | `examples/demo_spec.yaml` |
| Quick mode flags | n/a | `--name`, `--stages` (comma list) | `--name rnaseq --stages align,count` |

Spec fields (defaults in brackets): `project`; `description`; `items.key` [`datasets`], `items.wildcard` [`dataset`], `items.entries.<name>.{path, label, columns, overrides}`; `stages[].{name, description, ext [tsv], params}`; `software.<tool>: <path or name>`; `output_dir` [`results`], a relative path inside the project (absolute paths, drives and `..` are refused). Param values and item paths become command-line arguments, so control characters (newlines, tabs, NUL) in them are refused.

## Workflow

1. **Validate (prescriptive)**: names are lower_snake_case, no Python keywords or Snakemake reserved words (`all`, `input`, `params`, ...), no duplicate stages, no null params, overrides only reference existing stages/params. Reject with a clear message; never guess a fix.
2. **Render (prescriptive)**: emit the files listed under Output Structure. Stage *i* reads stage *i-1*'s output and its done-sentinel; stage 1 reads the item's `path`. `rule all` requests the last stage's sentinel for every target.
3. **Write (prescriptive)**: refuse a non-empty project directory, or an `--output` that already holds `report.md`, `result.json` or `reproducibility/`, unless `--force`; `--force` overwrites scaffold files and the report only and never deletes other files. Spec text never goes into a string literal or docstring: free text (descriptions) becomes `#` comments, values reach Python through `repr()`, YAML through `yaml.safe_dump`, and the shell through Snakemake's `{...:q}` quoting.
4. **Check (optional)**: with `--check`, run `snakemake -n --cores 1` in the project and record pass/fail.
5. **Report (flexible)**: write `report.md` + `result.json` + reproducibility bundle; explain next steps to the user in plain language.

## CLI Reference

```bash
# From a spec
python skills/snakemake-bio-scaffold/snakemake_bio_scaffold.py \
  --input my_spec.yaml --output ./pipelines

# Quick mode (no spec file)
python skills/snakemake-bio-scaffold/snakemake_bio_scaffold.py \
  --name rnaseq_counts --stages trim,align,count --items s1=input/s1.fq.gz,s2=input/s2.fq.gz \
  --output ./pipelines

# Demo (synthetic GWAS summary statistics, two cohorts) + dry-run check
python skills/snakemake-bio-scaffold/snakemake_bio_scaffold.py --demo --output /tmp/scaffold --check

# Via ClawBio runner
python clawbio.py run snakemake-scaffold --demo
```

## Demo

```bash
python clawbio.py run snakemake-scaffold --demo
```

Expected output: a `demo_sumstats_qc/` project with three stages (`format_sumstats -> filter_maf -> filter_pvalue`) over two synthetic cohorts (200 variants each; `cohort_b` overrides `maf` to 0.05). `cd demo_sumstats_qc && snakemake --cores 1` completes all 6 jobs plus `all`.

## Algorithm / Methodology

An agent can apply this layout by hand without the script:

1. **Snakefile**: `CFG = load_config()` once; `PYTHON = sys.executable`; `OUT = CFG["analysis"]["output_dir"]`; `wildcard_constraints: <wc>="|".join(re.escape(n) for n in ITEMS)`; `include:` each `rules/<stage>.smk`; `rule all` = last stage's sentinel for each target.
2. **Rule**: named `input:`/`output:` built from `f"{OUT}/<stage>/{{<wc>}}.<ext>"`, plus `done=touch(f"{OUT}/done/<stage>_{{<wc>}}.done")`; params are lambdas reading the merged per-item settings; `shell:` is one templated command calling `"{PYTHON}" scripts/<stage>.py`, redirected to `log:`.
3. **Script**: docstring with a copy-pasteable example; argparse `--kebab-case` flags matching `analysis.yaml` keys; `logging` with "N in, M dropped, K kept"; shared helpers in `scripts/lib/` via a `sys.path.insert` shim; never `import snakemake`.
4. **Config loader**: one `ConfigError`; `_resolve()` to absolute forward-slash paths; per-item `overrides` merged over stage defaults, unknown keys rejected (catches typos).

## Example Queries

- "Set up a Snakemake pipeline for my GWAS QC: SNP QC, sample QC, pruning, PCA, association"
- "Scaffold a Snakemake project with stages trim, align, count for 6 samples"
- "I want my new LDSC pipeline to follow the same structure as my other Snakemake projects"

## Example Output

```markdown
# Snakemake Bio Scaffold Report

**Project**: `demo_sumstats_qc`
**Stages**: 3 · **datasets**: 2 · **Wildcard**: `{dataset}`

| # | Stage | Output | Settings (analysis.yaml) |
|---|-------|--------|--------------------------|
| 1 | `format_sumstats` | `results/format_sumstats/{dataset}.tsv` | none |
| 2 | `filter_maf` | `results/filter_maf/{dataset}.tsv` | `maf=0.01` |
| 3 | `filter_pvalue` | `results/filter_pvalue/{dataset}.tsv` | `p_threshold=5e-08` |
```

Generated rule (`rules/filter_maf.smk`):

```python
rule filter_maf:
    """Stage 2: Drop variants below the minor-allele-frequency threshold."""
    input:
        data=f"{OUT}/format_sumstats/{{dataset}}.tsv",
        done=f"{OUT}/done/format_sumstats_{{dataset}}.done",
    output:
        result=f"{OUT}/filter_maf/{{dataset}}.tsv",
        summary=f"{OUT}/filter_maf/{{dataset}}.summary.json",
        done=touch(f"{OUT}/done/filter_maf_{{dataset}}.done"),
    params:
        maf=lambda wc: _item_analysis(wc.dataset, "filter_maf")["maf"],
    log:
        f"{OUT}/logs/filter_maf/{{dataset}}.log",
    shell:
        '"{PYTHON}" scripts/filter_maf.py '
        '--input "{input.data}" --out "{output.result}" --summary-json "{output.summary}" '
        '--maf "{params.maf}" '
        '> "{log}" 2>&1'
```

## Output Structure

```
output_directory/
├── report.md              # Stage table, file inventory, next steps
├── result.json            # Machine-readable project description
├── <project>/             # Generated project (Snakefile, config/, rules/, scripts/, ...)
└── reproducibility/
    ├── commands.sh         # Exact command to regenerate the scaffold
    ├── environment.yml     # Conda/pip environment snapshot
    └── checksums.sha256    # SHA-256 of report, result and every generated file
```

<!-- The <project>/ subtree is named after the spec's `project`; see result.json `files`. -->

## Dependencies

**Required**:
- `pyyaml` >= 6.0; reading the spec, writing the config files

**Optional**:
- `snakemake` >= 8; only for `--check` and for running the generated project

## Gotchas

- **Editing a script does not rerun its rule.** You will want to tell the user "just rerun snakemake" after they change `scripts/<stage>.py`. Do not. Snakemake tracks the rule's text, not the content of the script it calls; tell them `snakemake --cores N -R <stage>`.
- **A done-sentinel can outlive its output.** You will want to delete only `results/<stage>/...` to force a recompute. Do not. The `results/done/<stage>_<item>.done` file keeps the DAG believing the stage is finished; delete both, or use `-R`.
- **Shared inputs cascade.** You will want to put a sample list or region file that every job reads into `input:`. Think first: touching it invalidates every completed job that declares it, including expensive ones. When adding one new item, build only its targets (`snakemake results/done/<last>_<item>.done`), not bare `snakemake`.
- **Never call bare `python` in `shell:`.** It can resolve to a different environment than the one running Snakemake. The scaffold's `{PYTHON:q}` (`sys.executable`) exists for this reason; keep it when editing rules.
- **Spec text is untrusted input to generated code.** You will want to drop a description or a setting into a docstring or `"..."` literal in a template. Do not: `"""`, a trailing backslash or a quote ends the literal and the rest of the text runs as code, including on `snakemake -n` and `--help`. Put free text in `#` comments and values through `repr()` / `yaml.safe_dump`; the hostile-spec tests tokenize every generated file to enforce this.
- **Quote shell fields with `:q`, not `"..."`.** You will want to write `--x "{params.x}"`. Do not: `$(...)` and backticks still expand inside double quotes. Use `{params.x:q}`, `{input.data:q}`, `{log:q}`.
- **Windows paths.** Backslash paths mixed with forward-slash rule patterns can make Snakemake miss the producing rule (spurious `MissingInputException`), so the loader normalises to `/`. Relative CLI targets can also raise `MissingRuleException`; use absolute targets. Deep project paths can overflow MAX_PATH in `.snakemake/metadata`; keep projects short-pathed or use `--drop-metadata`.
- **SLURM time limits stack.** Raising a rule's `runtime` does not help if the controller job that launches Snakemake has a smaller `--time`; check both.
- **Do not invent analysis logic in the stubs.** The skill's job ends at a correct skeleton; implementing `process()` is a separate, explicit step with the user.

## Safety

- **Local-first**: generates files locally; no data is read or uploaded (the demo uses synthetic summary statistics only).
- **Disclaimer**: every report includes the ClawBio medical disclaimer.
- **Audit trail**: the reproducibility bundle records the command and checksums of every generated file.
- **No overwrite surprises**: refuses non-empty project directories and an earlier report in `--output` unless `--force`.
- **No code from specs**: spec text is never placed in a Python string literal or docstring (see Gotchas).

## Agent Boundary

The agent (LLM) gathers stage names, items and settings from the user, writes the spec, dispatches the script and explains the result. The skill (Python) validates and renders. The agent must NOT hand-edit the generated layout into a different structure, invent thresholds, or claim stub stages perform analysis.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- the request mentions Snakemake plus "new", "scaffold", "template", "structure" or "set up"
- a user describes a multi-step analysis they want to turn into a reusable pipeline

**Chaining partners**:
- `gwas-pipeline`, `fine-mapping`, `mendelian-randomisation`: their steps can become stages of a scaffolded project.
- `repro-enforcer`: audit the generated project's reproducibility.

## Maintenance

- **Review cadence**: on each Snakemake major release.
- **Staleness signals**: a Snakemake release that changes `include:`, `touch()`, or `wildcard_constraints` semantics; the integration test (`snakemake --cores 1` on the demo project) failing.
- **Deprecation**: archive to `skills/_deprecated/` if Snakemake ships an equivalent official project generator with the same conventions.

## Citations

- [Mölder et al. 2021, Sustainable data analysis with Snakemake](https://doi.org/10.12688/f1000research.29032.2); workflow engine
- [Snakemake documentation: distribution and reproducibility](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html); recommended project layout this skill extends
