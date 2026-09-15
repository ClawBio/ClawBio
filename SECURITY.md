# Security Policy

ClawBio runs bioinformatics skills on a user's own machine, often against that
user's own genome. A vulnerability here can expose genetic data, so reports are
handled privately until a fix or mitigation is available.

## Reporting a vulnerability

Do not open a public issue for a security problem.

1. Preferred: GitHub private vulnerability reporting. Open
   <https://github.com/ClawBio/ClawBio/security/advisories/new> and describe
   the problem. The report is visible only to the maintainers until an
   advisory is published.
2. Fallback: email the maintainer listed in `pyproject.toml`
   (`mc.admin@manuelcorpas.com`) with `ClawBio security` in the subject line.

Include what you can of: the skill or module affected, the version or commit,
a minimal reproduction, and what an attacker gains. A proof of concept against
the bundled demo data is welcome. Never send real patient or personal genomic
data in a report.

## What to expect

ClawBio has one maintainer at present, so these are commitments made in good
faith, not contractual service levels.

- Acknowledgement within 5 working days.
- An initial assessment (accepted, needs more information, or out of scope)
  within 14 days.
- A fix or a documented mitigation targeted within 90 days of acceptance,
  and sooner for anything that exposes genetic data or executes untrusted
  input.
- Coordinated disclosure: a public advisory and a CHANGELOG entry are
  published with the fix. Reporters are credited unless they ask not to be.

## Scope

In scope:

- The `clawbio` Python package, the CLI, and every skill under `skills/`.
- Skill documentation, demo data or fixtures that could steer an AI agent
  into unsafe behaviour (prompt injection through `SKILL.md` or bundled files).
- Secrets, tokens or credentials committed to the repository, or leaked
  through a skill's output or reproducibility bundle.
- Supply chain: dependency pins, `uv.lock`, GitHub Actions workflows, and the
  PyPI release path.
- Any skill that sends data off the machine without saying so in its
  `SKILL.md` and in [docs/data-handling.md](docs/data-handling.md).

Out of scope, please use a normal issue instead:

- Scientific accuracy: wrong annotations, misleading percentiles, incorrect
  guideline mappings. These matter and are handled in the open through issues
  and ClawBench, not through this policy.
- Vulnerabilities in third-party services that a skill calls (Ensembl, NCBI,
  PGS Catalog, Open Targets and others). Report those upstream, and tell us if
  ClawBio should stop calling the service.
- Rate-limit or denial-of-service effects on public APIs caused by running a
  skill in a loop.

## Supported versions

Security fixes go into the latest minor release only.

| Version | Supported |
|---------|-----------|
| 0.7.x   | Yes       |
| < 0.7   | No, upgrade to the current release |

## Good-faith research

Research that follows this policy, avoids accessing or destroying other
people's data, and gives the maintainers reasonable time to respond will be
treated as authorised and welcome.
