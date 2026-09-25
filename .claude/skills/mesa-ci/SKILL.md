---
name: mesa-ci
description: mesa's CI workflows, devcontainer, and the toolchain version doctrine — DESCRIPTION is the single source of truth and R/Bioconductor versions are derived, never hardcoded. Use when editing .github/workflows, .devcontainer, or pinning/bumping any tool version.
---

# CI, devcontainer and toolchain versions

## The doctrine: DESCRIPTION is the single source of truth

**Never hardcode an R or Bioconductor version anywhere.** `.devcontainer/resolve_versions.sh`
is the one place that maps R → Bioconductor for the whole repo. It parses
`R (>= X.Y.Z)` from `DESCRIPTION`, looks up the latest *released* Bioconductor for that R
from `bioconductor.org/config.yaml`, and prints `KEY=VALUE` lines that CI appends to
`$GITHUB_OUTPUT` / `$GITHUB_ENV`:

```bash
bash .devcontainer/resolve_versions.sh
# R_VERSION=4.6
# R_VERSION_FULL=4.6.0
# BIOC_VERSION=3.23
# BIOC_RELEASE=RELEASE_3_23
# ROXYGEN_VERSION=8.1.0
```

It is deliberately dependency-free — curl, grep, sed and POSIX awk only, no gawk
extensions — so it runs on a bare runner before any container is pulled, and behaves the
same on macOS and Linux. **Keep it that way** when editing it.

Upgrading R is therefore a one-line edit to `DESCRIPTION`. CI and the image build pick it
up automatically. Optional Bioc pinning goes in `.devcontainer/versions.env`, which is
otherwise not edited.

### Version inventory

Before pinning or bumping any tool, check where its version is supposed to live. The
authoritative table is in [`.devcontainer/UPGRADE_GUIDE.md`](../../../.devcontainer/UPGRADE_GUIDE.md).
**If a value cannot be derived, it must be listed there with every place it is duplicated** —
that table is the defence against a version being updated in one file and missed in another.

The one known non-derivable value is the **RSPM Ubuntu codename** (`noble`), hardcoded in
*both* `check-bioc.yml` and `.devcontainer/install.R`. It tracks the Bioconductor base
image's Ubuntu release (focal = 20.04, jammy = 22.04, noble = 24.04), not the R version, so
it only changes when the base image does — and when it does, **both** files need it.

## Workflows

| Workflow | Purpose |
|---|---|
| `check-bioc.yml` | The authoritative check. biocthis-generated (`biocthis::use_bioc_github_action()`), runs `R CMD check` + `BiocCheck` across platforms, plus covr and pkgdown on `main`. Triggered by changes to `R/`, `tests/`, `vignettes/`, `inst/`, `DESCRIPTION`, `NAMESPACE`, and by any PR. |
| `build-image.yml` | Builds and pushes the `slim` and `full` devcontainer images to ghcr.io. Triggers on `.devcontainer/**` or `DESCRIPTION` changes. ~20–40 min. |

Because `check-bioc.yml` is biocthis-generated, prefer regenerating or making surgical
edits over rewriting it — gratuitous divergence from upstream makes future biocthis updates
painful.

Put `/nocache` in a commit message to bypass the CI package cache; bump `cache-version` in
the workflow env to invalidate it for everyone.

## Devcontainer

`.devcontainer/install.R` installs mesa's declared dependencies **read from DESCRIPTION** —
the package list is never hand-maintained. The `slim`/`full` split is the single explicit
`full_only` list in that file; `slim` omits heavy genome/annotation packages, so
data-dependent tests skip there (see `mesa-tests`).

Genuine extras are installed separately. GitHub-only packages (`ggtree`, `immunedeconv`)
are **pinned to explicit SHAs** for reproducibility — keep that pattern for anything
installed from GitHub.

Note that the extras loop installs only when a package is **absent**
(`if (!requireNamespace(pkg))`), so it will not correct an image that already carries a
wrong version. Anything that needs a specific version must compare `packageVersion()` and
reinstall on mismatch, not merely check presence — that is what the roxygen2 block just
below it does.

**roxygen2 is pinned.** `resolve_versions.sh` emits `ROXYGEN_VERSION` from `DESCRIPTION`'s
`Config/roxygen2/version` (falling back to the pre-8.0 `RoxygenNote`), `install.R` installs
exactly that version, and the `roxygen-drift` job in `check-bioc.yml` regenerates the docs
with it and fails on any diff. The job deliberately does **not** use the latest roxygen —
otherwise every upstream release would turn it red and it would become the churn it exists
to prevent. Upgrading roxygen2 is a one-line `DESCRIPTION` edit, in its own PR.

## Upgrading R or Bioconductor

Follow [`.devcontainer/UPGRADE_GUIDE.md`](../../../.devcontainer/UPGRADE_GUIDE.md) —
it covers compatibility checking (`check_bioc_compat.R`), the one-line DESCRIPTION edit,
HPC library updates, Codespaces verification and a PR checklist.
