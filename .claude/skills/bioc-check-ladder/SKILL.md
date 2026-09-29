---
name: bioc-check-ladder
description: How to verify mesa — devtools::test, R CMD check, BiocCheck and coverage — where each can run, and how to set a macOS machine up to run the whole ladder locally. Use before declaring work done, when a check fails, or when interpreting BiocCheck output.
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

- **A local macOS machine** — set up once, this runs the whole ladder and reproduces CI's
  BiocCheck result exactly. See below.
- **Devcontainer / Codespaces** — `ghcr.io/cruk-mi/mesa/devcontainer:codespaces-slim`.
  Carries the right R and the full Bioc stack. The `slim` variant omits heavy genome data
  packages, so data-dependent tests skip; `full` has them.
- **CI** — `.github/workflows/check-bioc.yml` runs the full ladder across platforms. This
  is the authoritative result.

When a check cannot run locally, say so plainly and push the branch so CI runs it. That is
a legitimate outcome, not a failure to report.

## Setting a macOS machine up to run the whole ladder

Done once on an arm64 Mac, verified reproducing CI's result for `97dd9f2`. Versions come
from `.devcontainer/resolve_versions.sh` — never pick them by hand. Ask it which R you
need, since that is the one step needing sudo:

```bash
bash .devcontainer/resolve_versions.sh
# R_VERSION=4.6  R_VERSION_FULL=4.6.0  BIOC_VERSION=3.23  BIOC_RELEASE=RELEASE_3_23
# (a ROXYGEN_VERSION=... line is added once the resolver learns it, in #105)
```

```bash
brew install r-rig                       # a formula; `--cask r-rig` does not exist (that is r-rig-app, the GUI)
rig add 4.6.0 && rig default 4.6.0       # use the resolved R_VERSION_FULL; needs sudo, so an agent cannot do this step
Rscript .claude/scripts/setup-r-toolchain.R
```

`setup-r-toolchain.R` calls the resolver itself, so nothing is passed to it and nothing is
hardcoded in it. It pins Bioconductor to the resolved release, installs every declared
dependency, installs the check tooling (`rcmdcheck`, `BiocCheck`, `devtools`, `covr`), and
installs the **pinned** roxygen2 — the version `DESCRIPTION` records as having generated
the committed `man/`, not the current release. It refuses to run against a different R
than the resolver names, and its closing line restates the three versions so a mismatch is
visible before you regenerate anything.

On a branch predating `ROXYGEN_VERSION` in the resolver, it says so and skips that pin;
do not regenerate `man/` from such a machine. Then:

```bash
R CMD build .
R CMD check --no-manual mesa_0.99.6.9000.tar.gz
Rscript -e 'BiocCheck::BiocCheck("mesa_0.99.6.9000.tar.gz")'
```

### The five things that actually bite on macOS

Each of these was hit and fixed; none is obvious from the error message.

| Symptom | Cause and fix |
|---|---|
| `R-4.6.0-arm64.pkg` 404s under `big-sur-arm64` | The 4.6 series moved to `sonoma-arm64`. rig handles this; a hand-built CRAN URL will not. |
| Installer fails: *"unexpected error while moving files to the final destination"* | Only `/Applications/R.app` failed — macOS blocks modifying an existing signed app bundle without App Management permission. **The framework installed fine and the ladder does not use R.app.** Check `R --version` before assuming the install failed. |
| `not available as a binary package` for `org.*.db`, `TxDb.*`, `BSgenome.*`, `MEDIPSData` | Bioconductor ships annotation and experiment data **source-only**. They are pure data, so `type = "source"` needs no compiler. A binary-only pass silently leaves them out. |
| Vignettes fail: `libXrender.1.dylib` not found under `/opt/X11` | P3M's `gdtools` binary links XQuartz's cairo. Rebuild from source against Homebrew's: needs `pkg-config`, `cairo`, and `LIBRARY_PATH` including `/opt/homebrew/opt/gettext/lib` (else `ld: library 'intl' not found`). No XQuartz required. |
| `symbol not found in flat namespace '___kmpc_barrier'` | An OpenMP package (`stringdist`, a BiocCheck dependency) needs the runtime. `brew install libomp` plus a `~/.R/Makevars` carrying **both** `-Xclang -fopenmp` (compile) and `-lomp` (link) — Apple clang adds neither on its own. |

Also: `R CMD build` needs `pandoc` (`brew install pandoc`). The Bioconductor docker images
bundle it, so CI never surfaces this.

### Reading a local result against CI

Local reproduces CI's findings exactly, plus **one extra NOTE**:

```
NOTE: Cannot determine whether maintainer is subscribed to the Bioc-Devel mailing list
```

That is an environment artifact — the check cannot reach the list locally. Do not treat it
as a regression, and do not "fix" it. Anything beyond it is a genuine divergence worth
investigating.

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
- **`Invalid package Version`** — BiocCheck rejects the four-part devel form
  `X.Y.Z.9000`. **Expected while a devel section is open**; it clears when the cycle is
  cut. Do not "fix" it by changing the version. See `bioc-release-cycle`.

Verify the tarball directly rather than assuming:

```bash
R CMD build . && tar tzf mesa_*.tar.gz | grep -vE '^mesa/(R|man|tests|data|vignettes|inst)/'
```

## Precedence

The `r-lib` plugin's `cran-extrachecks` skill covers CRAN's ad-hoc requirements. Useful
background, but **CRAN is not Bioconductor** — its advice does not substitute for
`BiocCheck`, and its release conventions do not apply here.
