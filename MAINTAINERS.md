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

## Repository permissions (ClawBio/ClawBio)

| GitHub | Permission | Notes |
|--------|------------|-------|
| [@manuelcorpas](https://github.com/manuelcorpas) | admin | Lead maintainer |
| [@jaymoore-research](https://github.com/jaymoore-research) | admin | Repository admin; backup for the thirty-day rule in GOVERNANCE.md |
| [@camlloyd](https://github.com/camlloyd) | maintain (includes push) | Reviews and merges skill contributions within the limits in GOVERNANCE.md |
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

## Second maintainer wanted

ClawBio has one lead maintainer and is seeking a second maintainer who can
share review of community skill contributions and act as backup for security
response. The expectations and the path are in GOVERNANCE.md, "Becoming a
maintainer". If that is you, open an issue.
