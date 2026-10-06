# AGENTS.md

Instructions for any AI agent working on **mesa**. `CLAUDE.md` imports this file, and
[`.claude/README.md`](.claude/README.md) explains every AI file in the repo.

## This package

mesa (Methylation Enrichment Sequencing Analysis) is an R/Bioconductor package for MBD-seq
and MeDIP-seq data. It is in Bioconductor pre-submission (`0.99.z`), so it must pass
`R CMD check` and `BiocCheck` and follow the
[Bioconductor guidelines](https://contributions.bioconductor.org/r-code.html).

## Key commands

```r
devtools::load_all()                 # run code against the working tree
devtools::test()                     # all tests
devtools::test(filter = "^{name}")   # tests/testthat/test-{name}.R
devtools::document()                 # regenerate man/ and NAMESPACE
devtools::check()                    # R CMD check
BiocCheck::BiocCheck(".")            # Bioconductor checks
covr::package_coverage()             # coverage, target >= 80 %
```

Run R with `Rscript -e "..."`. Not every machine can run the whole list, because
`DESCRIPTION` needs a recent R. Use the devcontainer or CI when a check can't run here, and
say which checks did not run. The `bioc-check-ladder` skill says where each one works.

## Code style

- 4-space indent, no tabs. At most 80 characters per line, roxygen included.
- `<-` for assignment. `TRUE`/`FALSE`, never `T`/`F`. `seq_len()`/`seq_along()`, never `1:n`.
  `vapply()`, never `sapply()`. camelCase names.
- [`.lintr`](.lintr) is the source of truth for these rules. `BiocCheck` in CI is the final check.
- No `library()`/`require()` in `R/`: declare the package in `DESCRIPTION` and call `pkg::fn()`.
  No `:::` into other packages. Never install packages from package code.
- Prefer Bioconductor classes (`GRanges`, `DataFrame`) and vectorised code. Keep functions
  short.
- **Ask before adding a dependency.** Bioconductor reviewers flag unneeded imports.
- Keep diffs small. Don't restyle working code, and keep comments that explain the biology.

## Tests

- `testthat` 3rd edition, in `tests/testthat/`. Every new function gets a test, and every
  bug fix gets a regression test.
- Reuse the fixtures in `tests/testthat/helper-fixtures.R`. Guard heavy annotation packages
  with `skip_if_not_installed()`. The `mesa-tests` skill has the details.

## Documentation

- Every exported function needs a title, a description, `@param` for each argument,
  `@return` and a runnable `@examples` that uses `exampleMouse` or `exampleTumourNormal`.
  The title and description are implicit (the first line and the paragraph after it), so
  don't add `@title`/`@description`. Add `@seealso` for related functions, and use
  `\dontrun{}` only when an example needs the network or files that aren't shipped.
- `man/*.Rd` and `NAMESPACE` are generated: never edit them by hand. Regenerate them with
  the roxygen2 version pinned in `DESCRIPTION` (`Config/roxygen2/version`). Another version
  rewrites every page, and CI fails on the drift.
- `data/*.rda` comes from `data-raw/`: re-run the script, never edit the file.

## NEWS.md

Every user-visible change gets a bullet under the top `# mesa X.Y.Z.9000` heading, in the
same PR as the change.

- File each bullet under a category heading (`## Bug fixes`, `## New features`,
  `## Documentation`, `## Testing`, `## Infrastructure`, `## Style`).
- Describe what changed for someone calling the function, not how it was done. Write
  function names as `` `makeQset()` `` and end the bullet with
  `([#NN](https://github.com/cruk-mi/mesa/issues/NN))`.
- Never bump `Version:` in a work PR. Version bumps have their own PRs, described in the
  `bioc-release-cycle` skill.

## Toolchain versions

`DESCRIPTION` is the single source of truth. The R version comes from its `R (>= …)` line,
and `.devcontainer/resolve_versions.sh` derives the Bioconductor version from it. Never
hardcode an R or Bioconductor version in a workflow, Dockerfile or script. Read
[`.devcontainer/UPGRADE_GUIDE.md`](.devcontainer/UPGRADE_GUIDE.md) before changing a
version.

## Git and pull requests

- Branch off `main` (`feat/`, `fix/`, `docs/`, `test/`, `refactor/`, `chore/`). Never commit
  or push to `main`, delete tags or rewrite published history. Never push to Bioconductor
  (`BiocStaging/mesa`): the maintainer does that.
- One issue per branch and PR. Put anything else you notice in the PR's **Parked** section
  instead of fixing it.
- [Conventional Commits](https://www.conventionalcommits.org):
  `type(scope): imperative summary`, with type `feat`, `fix`, `docs`, `test`, `refactor`,
  `perf` or `chore` and the R file name as the scope. Put `Fixes #n`/`Refs #n` in the
  footer. Keep each commit to one logical change.
- Open PRs as **drafts**. Only a human marks them ready. Merge (squash) only when a human
  explicitly asks for that merge; an approving review is not enough. Fill in
  [`.github/pull_request_template.md`](.github/pull_request_template.md), because
  `gh pr create --body-file` skips it. When the scope changes, rewrite the body to match
  the diff instead of appending an "Update:" note.
- Before opening an issue, search for one that already exists
  (`gh issue list --repo cruk-mi/mesa`).

## Attribution

This repo is public and under Bioconductor review, so authorship is never misrepresented.
Every AI-made commit carries a co-author trailer (Claude Code:
`Co-Authored-By: Claude <noreply@anthropic.com>`). Every PR, issue or comment an agent
writes says so on its first line.

## Skills

Claude Code loads these from `.claude/skills/` when the task matches. Other agents can read
them as plain Markdown.

| Task | Skill |
|---|---|
| Tests, `R CMD check`, `BiocCheck`, coverage, local R setup | `bioc-check-ladder` |
| Version bumps, devel cycles, releases, the Bioconductor push | `bioc-release-cycle` |
| Writing or fixing tests | `mesa-tests` |
