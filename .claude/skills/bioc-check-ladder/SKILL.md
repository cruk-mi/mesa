---
name: bioc-check-ladder
description: How to verify mesa — devtools::test, R CMD check, BiocCheck and coverage — and crucially WHERE each can actually run, since a typical laptop cannot run the full ladder. Use before declaring work done, when a check fails, or when interpreting BiocCheck output.
---

# The mesa check ladder

Run in this order. Each rung is cheaper than the one after it, so fail fast.

```r
devtools::test()                 # 1. unit tests
devtools::check()                # 2. R CMD check
BiocCheck::BiocCheck(".")        # 3. Bioconductor-specific checks
covr::package_coverage()         # 4. coverage, target >= 80%
```

`BiocCheck` is **complementary to, not a replacement for**, `R CMD check`. Always run it
*after*, never instead.

## Where each rung can actually run — check before promising

**Do not claim a check passed without running it, and do not tell the human to run a check
their machine cannot execute.** Verify the environment first:

```r
R.version.string                              # vs DESCRIPTION's R (>= X.Y.Z)
requireNamespace("BiocCheck", quietly = TRUE)
requireNamespace("covr", quietly = TRUE)
```

Common blockers on a local machine:

| Blocker | Effect |
|---|---|
| Installed R older than `DESCRIPTION`'s `R (>= X.Y.Z)` | `devtools::check()` cannot install the package; the whole ladder above rung 1 is unavailable |
| `BiocCheck` / `covr` not installed | rungs 3 and 4 unavailable (they are not default installs) |
| Docker not running | the devcontainer fallback is unavailable too |

`devtools::test()` and `lintr` generally work anywhere.

### The real venues

- **Devcontainer / Codespaces** — `ghcr.io/cruk-mi/mesa/devcontainer:codespaces-slim`.
  Carries the right R and the full Bioc stack. The `slim` variant omits heavy genome data
  packages, so data-dependent tests skip; `full` has them.
- **CI** — `.github/workflows/check-bioc.yml` runs the full ladder across platforms. This
  is the authoritative result.

When a check cannot run locally, say so plainly and push the branch so CI runs it. That is
a legitimate outcome, not a failure to report.

## Reading BiocCheck output

Bioconductor requires maintainers to address **all errors and warnings, and most notes**.

| Level | Treatment |
|---|---|
| `ERROR` | Blocking. Submission is rejected. Fix before the PR. |
| `WARNING` | Blocking in practice. Reviewers will require a fix. |
| `NOTE` | Justify or fix. Some are unavoidable; say which and why in the PR. |

Recurring ones in this package:

- **Tarball size** — BiocCheck's limit is 5,000,000 bytes. Stray files swept into the build
  are the usual cause; check `.Rbuildignore` before blaming the data objects. See PR #101.
- **Non-standard top-level files** — anything at the package root that is not a standard R
  package file gets flagged. Agent config (`CLAUDE.md`, `AGENTS.md`, `.claude/`) must stay
  in `.Rbuildignore`.
- **Coverage below 80 %** — see `mesa-tests`.

Verify the tarball directly rather than assuming:

```bash
R CMD build . && tar tzf mesa_*.tar.gz | grep -vE '^mesa/(R|man|tests|data|vignettes|inst)/'
```

## Precedence

The `r-lib` plugin's `cran-extrachecks` skill covers CRAN's ad-hoc requirements. Useful
background, but **CRAN is not Bioconductor** — its advice does not substitute for
`BiocCheck`, and its release conventions do not apply here.
