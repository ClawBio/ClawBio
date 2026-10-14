# Governance

This document says how ClawBio is run today. It describes the project at its
actual size: one lead maintainer, a small number of people with repository
rights, and community contributors. There is no steering committee, no
foundation and no voting body. If the project grows to need one, this file
changes first.

The people holding each right are listed in [MAINTAINERS.md](MAINTAINERS.md),
which is checked against the GitHub API rather than written from memory.

## Roles

**Lead maintainer.** Manuel Corpas. Holds final say on scope, releases,
security response and this document. Owns the PyPI project, the clawbio.ai
domain and the GitHub organisation.

**Maintainer.** A person with `maintain` or `admin` permission on the
repository. Reviews and merges community pull requests within the limits
below, triages issues, and may cut a release when asked to by the lead
maintainer.

**Contributor.** Anyone who opens a pull request. Contributors keep copyright
in their work and license it under MIT by contributing; see
[CONTRIBUTING.md](CONTRIBUTING.md).

**The daily PR agent.** An unattended job, run from the lead maintainer's own
machine, that audits community pull requests once a day. Its policy is fixed
in code and tested; it is summarised under "Automated merges" below.

## How decisions are made

Most decisions are made in the pull request or issue where they arise, by the
maintainer handling it. Anything that changes what ClawBio promises to users
is decided by the lead maintainer and recorded where the change lands:

- the public interface of the `clawbio` package and CLI, and any deprecation;
- the SKILL.md contract, the catalogue schema and the maturity tiers;
- what CI requires of a contribution;
- release timing and version numbers;
- security response, including whether to disclose and when;
- this document, MAINTAINERS.md and CODEOWNERS.

A decision of that kind is recorded in the CHANGELOG under `Unreleased` when
it changes behaviour, or in the pull request that makes it otherwise. There is
no separate RFC process. If you want to argue for a change, open an issue and
make the case; the lead maintainer answers there.

## Who merges what

| Change | Who may merge | Conditions |
|--------|---------------|------------|
| A new or changed skill under `skills/` | Any maintainer, or the daily PR agent | CI green, audit verdict SAFE, a domain read of the science, no conflict with `main` |
| `clawbio/`, `scripts/`, `.github/`, `pyproject.toml`, `uv.lock` | Lead maintainer | Review by the code owner named in `.github/CODEOWNERS` |
| Documentation and site pages | Any maintainer | CI green |
| A release tag `v*` | Lead maintainer | Tag ruleset: only repository admins can create or move one |
| GOVERNANCE.md, MAINTAINERS.md, CODEOWNERS, SECURITY.md | Lead maintainer | Code-owner review |

Nobody merges their own pull request into `clawbio/` or `.github/` without a
second person's review, and the daily PR agent never acts on a pull request
authored by the lead maintainer.

### Automated merges

The daily PR agent comments on, approves and merges a community pull request
only when all of the following hold: the static audit and a domain read of the
diff both return SAFE, CI is green (a run that never started is not green),
GitHub reports the branch mergeable, and the diff contains no prompt-injection
markers. Anything else is a comment or a request for changes, and anything the
scanner flags but the domain read does not confirm is held for a person. It
merges at most five pull requests per run. Its model cannot reach GitHub; the
model returns a recommendation and a separate, tool-free step acts on it using
state read from the API. A file named `DATA/DISABLED` on the machine that runs
it stops all merging while auditing continues. Every rule has a test; the
policy changes only by changing the code and its tests together.

## Becoming a maintainer

There is no application form. Someone becomes a maintainer when:

1. they have had several pull requests merged over a few months, including at
   least one review of somebody else's contribution that caught something real;
2. the lead maintainer asks them, or they ask and the lead maintainer agrees;
3. they accept in writing, in an issue or a pull request against
   MAINTAINERS.md, which states what they are taking on (the expected load is
   about two hours a week);
4. the lead maintainer grants the repository permission and merges the
   MAINTAINERS.md change.

Maintainer status lapses after six months without activity on the repository.
The person is told first; if they want to stay they say so and remain. It is
withdrawn immediately for a breach of the [Code of Conduct](CONTRIBUTING.md)
or of the security policy in [SECURITY.md](SECURITY.md).

ClawBio is looking for a second maintainer. If the description above fits you,
open an issue.

## Resolving disputes

1. Talk it through in the pull request or issue. Most disagreements are about
   a missing fact, and the fact settles it.
2. If two maintainers disagree, the lead maintainer decides and writes the
   reason in the thread.
3. If the disagreement is with the lead maintainer, say so plainly in the
   thread; the decision stands but the objection is recorded with it, and the
   question is reopened when new evidence arrives.
4. Conduct complaints go to the address in SECURITY.md and are handled by the
   lead maintainer, or by another maintainer if the complaint concerns the lead
   maintainer.

## If the lead maintainer is unavailable

If the lead maintainer has not responded on the repository for thirty days,
any repository admin listed in MAINTAINERS.md may: merge pull requests that
meet the conditions above, cut a patch release for a security fix, and respond
to security reports. They may not change this document, the licence, the
package name or the organisation's membership. A second owner account for the
GitHub organisation, controlled by the lead maintainer and protected by a
hardware key, is planned so that loss of one account does not lock the
project.

## Repository settings

The settings this document relies on, and their state on 11 September 2026:

- `main` is not branch-protected. This document recommends protecting it with
  required CI, required code-owner review for the paths in CODEOWNERS, and no
  force pushes. The `scientific-audit` CI job is excluded from the required
  set, because it runs an external audit suite pinned to a commit and is
  informative rather than blocking.
- Release tags `v*` are protected by a repository ruleset.
- Dependency updates are proposed weekly by Dependabot; CodeQL runs on every
  pull request and weekly on `main`; `uv audit` and `pip-audit` report on the
  lockfile.

## Changing this document

By pull request, reviewed by the lead maintainer, with the reason for the
change in the pull request body.
