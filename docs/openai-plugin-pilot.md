# ClawBio Research plugin pilot

This is a distribution experiment for three existing workflows, packaged for
ChatGPT and Codex using the portable Agent Plugins layout. It does not change
the existing Claude plugin or restore the deprecated MCP server.

The supported entry point runs bundled demonstration data only. It does not
accept private genomic inputs. Python execution and the dependencies listed in
the package are required. Installation in ChatGPT does not itself grant a Python
environment, install dependencies or prove that a workflow can execute there.

## Build and verify

From the repository root, in an environment with the analysis dependencies:

```bash
python -m pytest tests/test_openai_plugin_package.py -q
python scripts/build_openai_plugin.py --output dist/clawbio-research-0.1.0.zip
python -m zipfile -e dist/clawbio-research-0.1.0.zip dist/clawbio-research
python dist/clawbio-research/scripts/run_research_demo.py pharmgx --output output/plugin-pharmgx
python dist/clawbio-research/scripts/run_research_demo.py acmg --output output/plugin-acmg
python dist/clawbio-research/scripts/run_research_demo.py equity --output output/plugin-equity
```

Choose fresh output names; the builder and launcher refuse overwrites. To replay
a completed run into a new directory:

```bash
bash output/plugin-pharmgx/reproducibility/commands.sh output/plugin-pharmgx-replay
cd output/plugin-pharmgx
sha256sum -c reproducibility/checksums.sha256
```

The archive has a root `plugin.json`, three generated `skills/` folders, a
bundled source runtime and a fixed demo launcher. The generated skills preserve
the upstream methodology under `references/` but limit their execution procedure
to the demo launcher. Upstream algorithms are copied without alteration.
Package provenance maps copied bytes to their source paths and hashes and records
the base Git commit. Modified local source remains distinguishable by its hash;
the base commit is not a claim that every included byte was already committed.

Only explicit source/runtime and reference-data paths are included. Tests, eval
inputs, personal genomes, profiles, existing reports and credentials are excluded.
The pharmacogenomic input is a constructed six-locus fixture with deliberately
incomplete coverage, replacing the upstream Corpasome demo in this archive only.
The variant-classification fixture is the upstream public GIAB-derived panel
with cached evidence, not fresh clinical evidence. The equity fixture is the
upstream synthetic VCF and population map.

## Supported behavior and evidence

| Workflow | Demonstration | Boundary |
|---|---|---|
| PharmGx | Research drug/gene report and missing-genotype handling | No private input, free-text diplotype parsing or medication decision |
| ACMG | Cached research classification and existing self-audits | No live patient VEP calls or claim of clinical validation |
| Equity | Representation, heterozygosity, FST, PCA and exploratory HEIM composite | No institutional equity certification or measured clinical fairness |

Each run writes the upstream report and structured results, a package receipt,
the canonical disclaimer, a replay command requiring a fresh output directory,
an environment recipe with observed direct dependency versions, and checksums.
Upstream replay/environment/checksum files are retained under `upstream_*` names
when replaced. The environment recipe is a direct-dependency snapshot, not a
complete transitive lock; use the repository's `uv.lock` for its locked environment.

The tests build and extract the archive into a path containing spaces, execute
all workflows from outside the repository, deny socket connections, inspect
missing CYP2D6 abstention, verify output checksums, and reject private-input
flags, runtime tampering and existing output directories. A successful run is
execution evidence, not independent scientific validation.

The integrity check compares bytes with this package's manifest. It is not a
signature, a security sandbox or proof of publisher authenticity. The underlying
source scripts retain upstream interfaces; unsupported direct invocation is not
prevented by the launcher. The supported plugin instructions never direct an
agent to those interfaces for private data.

## Test installation before submission

Use the [official packaging instructions](https://developers.openai.com/plugins/build/plugins)
to expose the extracted folder through a local marketplace in ChatGPT desktop.
A repo-scoped `.agents/plugins/marketplace.json` can point to
`./dist/clawbio-research` with a local source, `AVAILABLE` installation policy,
`ON_INSTALL` authentication policy and `Productivity` category. This is an
authoring source, not a public listing. Do not commit an entry pointing to an
unbuilt or ignored `dist/` folder.

Validate in a new ChatGPT/Codex session: direct and indirect demo requests,
missing Python/dependencies, requests to upload private data, and requests for
diagnosis. Preserve a transcript, execution receipt and usability problems. The
local package tests do not prove skill activation, host installation or directory
acceptance. Obtain one external research user and a repeat use before adding
more workflows; measure successful completions, friction and human review time.

Public submission has a separate
[verified identity and review process](https://developers.openai.com/plugins/deploy/submission).
Confirm publisher metadata and all required listing assets/URLs in the portal.
This prototype contains no invented privacy/terms URLs, app references, hooks or
remote MCP endpoint. It has not been submitted, reviewed or published. A remote
MCP cannot currently be added later to an existing skills-only listing; decide
on that boundary before first public submission. Private institutional operation
will need a separately justified design rather than a cloud upload workaround.

## Subsequent experiments

1. **Agent infrastructure:** use an existing execution host to invoke the same
   packaged workflows and compare receipts. Add an Agents API adapter only after
   the integration identifies a capability gap; a new generic agent service is
   not required for this package.
2. **Evidence gates:** keep deterministic validation in code. A finite-output
   model classifier can route requests but does not validate a genotype, resolve
   copy number or establish clinical evidence. Existing missing-gene handling
   and self-audits are exercised here. A broader safety claim requires a separate
   adversarial suite covering partial calls, conflicting inputs, unresolved CNV,
   unsupported alleles and free-text negation if that interface is introduced.
3. **Institutional pilot:** first define a named institution, workflow, owner,
   intended population and acceptance contract. Measure raw provenance coverage,
   replay success, demographic coverage and stratified task errors with explicit
   denominators and uncertainty. Do not combine these into a 0-100 trust/equity
   score without a specified, justified and independently evaluated framework.
4. **Personal orchestration:** evaluate Dot against the existing executive loop
   using the same work cases and attention measurements. The SOUL doctrine and
   canonical memory stay authoritative. GitHub delivery is not proof of local
   GBrain ingestion, active orchestration or reduced attention cost.
