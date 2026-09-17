# Maintainers

Who holds which right on ClawBio, as of 11 September 2026. Every list below
was read from the GitHub API on that date (`gh api
repos/ClawBio/ClawBio/collaborators`, `gh api orgs/ClawBio/members`) rather
than written from memory. When you change a permission, change this file in
the same pull request. The roles are defined in [GOVERNANCE.md](GOVERNANCE.md).

## Lead maintainer

| Person | GitHub | ORCID |
|--------|--------|-------|
| Manuel Corpas | [@manuelcorpas](https://github.com/manuelcorpas) | [0000-0002-5765-9827](https://orcid.org/0000-0002-5765-9827) |

## Second maintainer

| Person | GitHub | ORCID | Since |
|--------|--------|-------|-------|
| Cameron Lloyd | [@camlloyd](https://github.com/camlloyd) | [0000-0003-1055-2282](https://orcid.org/0000-0003-1055-2282) | 2026-09-17 |

Listed by name only at his request. Holds repository admin and is named in
`.github/CODEOWNERS` beside the lead maintainer.

## Repository permissions (ClawBio/ClawBio)

| GitHub | Permission | Notes |
|--------|------------|-------|
| [@manuelcorpas](https://github.com/manuelcorpas) | admin | Lead maintainer |
| [@jaymoore-research](https://github.com/jaymoore-research) | admin | Repository admin; backup for the thirty-day rule in GOVERNANCE.md |
| [@camlloyd](https://github.com/camlloyd) | admin | Second maintainer since 17 September 2026; code owner with the lead maintainer; backup for the thirty-day rule |
| [@afonsoguerra](https://github.com/afonsoguerra) | read | Organisation member |

## Organisation (github.com/ClawBio)

| GitHub | Role |
|--------|------|
| [@manuelcorpas](https://github.com/manuelcorpas) | owner |
| [@afonsoguerra](https://github.com/afonsoguerra) | member |

A second owner account for the organisation, controlled by the lead maintainer
and protected by a hardware key, is planned and not yet created.

## Other controls

| Control | Held by |
|---------|---------|
| PyPI project `clawbio` | Manuel Corpas |
| Domain `clawbio.ai` (DNS) and GitHub Pages | Manuel Corpas |
| Release tags `v*` | Repository admins only, by ruleset |
| The daily PR agent | Runs from Manuel Corpas's machine; see GOVERNANCE.md |

## Becoming a maintainer

The expectations and the path are in GOVERNANCE.md, "Becoming a maintainer".
The second-maintainer role was filled on 17 September 2026. Further
maintainers are welcome on the same path; open an issue.
