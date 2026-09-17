# AGENTS.md

Canonical, tool-agnostic instructions for **any** AI agent working on **mesa**
(Claude Code, GitHub Copilot coding agent, Copilot Chat / agent mode, etc.).

> This file is the single source of truth for how agents work on this repo.
> `CLAUDE.md`, `.github/copilot-instructions.md`, `.github/instructions/*`, and the
> per-agent definitions in `.github/agents/*` and `.claude/agents/*` all defer to it.
> If they ever disagree, **AGENTS.md wins** — fix the other file.

---

## Project overview

**mesa** (Methylation Enrichment Sequencing Analysis) is an R/Bioconductor package for
analysing methylation enrichment sequencing data (MBD-seq, MeDIP-seq). It targets
Bioconductor submission and must conform to Bioconductor policies, coding standards, and
review guidelines at all times. Current version: pre-submission `0.99.x` series.

---

## Workflow contract (read this first)

Agents operate on a **dedicated branch as a sandbox**. The human owns when work becomes a
pull request and when anything merges.

- **Check the state of the repo first.** Run `git status` and `git log --oneline -5`
  before making changes, and say so if the working tree is dirty or HEAD is not where you
  expected — do not fold someone else's uncommitted work into your commit.
- **Branch off `main`.** Never commit directly to `main`. (`dev` is no longer the working
  branch — everything is cut from and PRs back to `main`.)
- Make the **smallest** set of changes needed. Do not refactor unrelated code.
- Commit in **atomic** steps, each a single logical change, using
  [Conventional Commits](#conventional-commit-messages).
- You **may push the feature branch** (with attribution — see below).
- You **may open a pull request, but only as a draft.** Never mark it ready for review,
  never request review, never merge.
- You may **never** push or force-push to `main`, delete branches or tags, or rewrite
  published history.
- When in doubt about an irreversible or outward-facing action, ask first.

These rules are also enforced mechanically by `.claude/hooks/guard-remote.py`. If the hook
and this file ever disagree, that is a bug — fix both.

---

## Attribution (required)

This repo is public and used for Bioconductor review, so **authorship must never be
misrepresented**. Every change an AI makes must be attributed:

- Add a co-author trailer to every commit the agent makes:
  - Claude Code: `Co-Authored-By: Claude <noreply@anthropic.com>`
  - Copilot coding agent: keep its default `Co-authored-by: Copilot` trailer.
- Clearly note AI authorship in any PR description, issue, or comment an agent writes.
- Never present AI-authored content as if written solely by the human.

---

## Branch strategy

```
main        ← primary branch. Branch all work off main; never commit directly.
feat/*      ← new features.            fix/*       ← bug fixes.
refactor/*  ← internal refactors.      docs/*      ← documentation only.
chore/*     ← maintenance (CI, deps).  test/*      ← test additions/fixes.
release/*   ← release stabilisation. Cut from main, merge to main + tag.
```

All branches are cut from `main` and PR back to `main`.

---

## Conventional Commit messages

```
<type>(<scope>): <short imperative description>

[optional body — wrap at 72 chars]

[optional footer: Fixes #<issue-number>]
Co-Authored-By: Claude <noreply@anthropic.com>
```

| Type | When to use |
|------|-------------|
| `feat` | New user-visible function or argument |
| `fix` | Bug fix |
| `docs` | Documentation only (roxygen, vignettes, README, NEWS) |
| `refactor` | Internal code change, no behaviour change |
| `test` | Adding or fixing tests |
| `chore` | CI, `data-raw` scripts, dependency bumps, tooling |
| `perf` | Performance improvement |

**Scope** (optional but recommended): the R file name without extension, e.g.
`makeDMRs`, `plotting`, `buildQset`, `DESCRIPTION`.

Reference issues in the footer with `Fixes #<n>` (closes on merge) or `Refs #<n>`.

---

## Release cycle

Every change set runs through three phases. Work PRs sit **between** phases 1 and 3 —
**never bundle a version bump into a work PR.**

Bioconductor version convention: `0.99.z` pre-submission, then `1.0.0` on acceptance.
Devel versions use an odd minor (e.g. `1.1.z`); release versions use an even minor
(e.g. `1.0.z`).

**Phase 1 — open the development section.** Branch `chore/bump-version-X.Y.Z.9000` off
`main`. Two commits, in order (reference: PR #96 = `59f870b`):

| commit message | change |
|---|---|
| `chore(DESCRIPTION): bump version to X.Y.Z.9000` | `DESCRIPTION` `Version:` only |
| `docs(NEWS): open X.Y.Z.9000 development section` | insert `# mesa X.Y.Z.9000` plus a blank line at the very top of `NEWS.md` |

**Phase 2 — the work.** Each feature/fix PR targets `main` and records its user-visible
changes under the `# mesa X.Y.Z.9000` heading, in the same PR that makes the change.

**Phase 3 — cut the release.** After the work PRs merge, branch
`chore/bump-version-X.Y.(Z+1)` off `main`. One commit
`chore(version): bump to X.Y.(Z+1)` (reference: PR #99 = `b5e80a0`): bump `DESCRIPTION`
and rename the `NEWS.md` heading.

**Stop after each phase.** Tell the human to push and open the PR. Tag format `vX.Y.Z`,
applied by the human after merge to `main`.

Full detail lives in the `bioc-release-cycle` skill.

---

## Issue workflow

Before creating an issue, **check whether a matching one already exists**:

```bash
gh issue list --repo cruk-mi/mesa --state open
```

If one does, report its number and title rather than opening a duplicate. If none does,
draft the title and body and show them to the human — say plainly that an agent wrote the
text, per the attribution rule above.

Reference issues from commits with `Fixes #<n>` (closes on merge) or `Refs #<n>`
(cross-reference only).

---

## Toolchain versions — single source of truth

**`DESCRIPTION` is the single source of truth; everything else derives from it.** If you
pin or bump a tool version, change it in the source-of-truth file and nowhere else. If a
value genuinely cannot be derived, it must appear in the inventory table in
[`.devcontainer/UPGRADE_GUIDE.md`](.devcontainer/UPGRADE_GUIDE.md) listing **every** place
it is duplicated.

- R version comes from `DESCRIPTION`'s `R (>= X.Y.Z)` line; the Bioconductor version and
  Docker base tag are derived from it by `.devcontainer/resolve_versions.sh`. CI and the
  image build both call that resolver.
- The dependency list comes from `DESCRIPTION`'s `Imports` / `Depends` / `Suggests`.
- **Never hardcode an R or Bioconductor version** in a workflow, Dockerfile or install
  script.

Read `UPGRADE_GUIDE.md` before changing any version. Detail lives in the `mesa-ci` skill.

---

## R / Bioconductor coding standards

Follow the [Bioconductor coding guidelines](https://contributions.bioconductor.org/r-code.html).

- 4-space indentation. No tabs.
- Line length ≤ 80 chars for R code; ≤ 100 chars for roxygen comments.
- Use `<-` for assignment, never `=` at top level.
- Never use `T` / `F` for `TRUE` / `FALSE`.
- Avoid `1:n`; use `seq_len(n)` or `seq_along(x)`.
- Avoid `:::` to reach unexported functions of other packages.
- No `library()` / `require()` inside package code — declare in `DESCRIPTION` `Imports`
  and call `pkg::fn()`.
- **Never** call `install.packages()` or `BiocManager::install()` from package source. A
  package must not install anything when it is loaded.
- Prefer vectorised ops and Bioconductor structures (`GRanges`, `DataFrame`, …) over base
  loops where natural.
- Keep functions short and single-purpose; split rather than nest deeply.
- No `<<-` (global assignment) without a documented reason.
- **Do not add new package dependencies without asking** — Bioconductor reviewers flag
  unnecessary imports.

### roxygen2 / documentation

- Every exported function needs `@title`, `@description`, `@param`, `@return`, and at least
  one runnable `@examples` block. Use `exampleMouse` / `exampleTumourNormal` data objects;
  avoid `\dontrun{}` unless an example genuinely needs network/files.
- Use `@seealso` to cross-reference related functions.
- **`man/*.Rd` and `NAMESPACE` are roxygen-generated — never hand-edit them.** Both carry
  roxygen2's `do not edit by hand` header. Regenerate with `roxygen2::roxygenise()`.
- **Before regenerating, check that your installed `roxygen2` matches `DESCRIPTION`'s
  `RoxygenNote`.** If it does not, **stop and tell the human** — regenerating with a
  different roxygen rewrites all 100+ man pages and bumps `RoxygenNote`, which must never
  ride along in a work PR. A roxygen upgrade is its own dedicated PR.

### Testing

- Tests live in `tests/testthat/`; use `testthat` (≥ 3rd edition).
- Every new function gets at least one test; every bug fix gets a regression test.
- Reuse the shared fixtures in `tests/testthat/helper-fixtures.R` rather than rebuilding
  example qsets. Guard heavy genome/annotation data with `skip_if_not_installed()`.
- BiocCheck wants ≥ 80 % coverage (`covr::package_coverage()`).
- Run `devtools::test()` before committing.

Detail lives in the `mesa-tests` skill.

---

## File-change behaviour

- Do not rename or move source files unless necessary.
- Do not rewrite working code for style preference. Keep diffs small.
- Preserve comments carrying biological/analytical context.
- `data-raw/` scripts generate `.rda` files in `data/` — never hand-edit generated `.rda`;
  re-run the script.
- Record every user-visible change in `NEWS.md` under the current devel heading.

---

## Local verification (run before declaring done)

```r
devtools::test()                 # unit tests
devtools::check()                # R CMD check
BiocCheck::BiocCheck(".")        # Bioconductor checks
covr::package_coverage()         # coverage (target ≥ 80%)
```

**Not every machine can run all of these.** `DESCRIPTION` requires a recent R, and
`BiocCheck` / `covr` are not part of a default install. Check what is actually available
before promising a result, and use the devcontainer or CI when it is not. The
`bioc-check-ladder` skill covers where each check can really run.

---

## Skills and plugins

Procedural detail lives in on-demand skills under `.claude/skills/`, so this file stays
short. Generic R practice comes from the installed
[`posit-dev/skills`](https://github.com/posit-dev/skills) plugins (`r-lib`, `github`,
`open-source`, `posit-dev`).

**Precedence: a mesa skill beats a plugin skill wherever they disagree.** The plugin
skills encode CRAN and tidyverse conventions, which are not Bioconductor's — notably
around release cadence, versioning and `BiocCheck`.

| Task | Skill |
|---|---|
| Version bumps, opening/closing a devel cycle | `bioc-release-cycle` |
| Running tests, `R CMD check`, `BiocCheck`, coverage | `bioc-check-ladder` |
| Writing or fixing tests | `mesa-tests` |
| Roxygen docs, `NEWS.md` entries | `mesa-docs-news` |
| CI workflows, devcontainer, toolchain versions | `mesa-ci` |
| Capturing a new procedure as a skill | `capture-skill` |

---

## Preferred decision rule

When multiple valid options exist, choose the one that:
1. keeps the human in control of PRs/merges,
2. keeps the diff smallest,
3. passes `R CMD check` and `BiocCheck` cleanly,
4. is easiest to maintain and review later.
