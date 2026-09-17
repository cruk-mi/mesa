---
applyTo: "R/**/*.R"
---

# R source code rules

Follow [`AGENTS.md`](../../AGENTS.md). For files under `R/`:

- 4-space indentation, no tabs. Line length ≤ 80 chars.
- `<-` for assignment, never `=` at top level. Never `T`/`F` for `TRUE`/`FALSE`.
- Use `seq_len()` / `seq_along()`, not `1:n`.
- No `library()` / `require()` / `install.packages()` / `BiocManager::install()` in package
  code. Declare deps in `DESCRIPTION` `Imports` and call `pkg::fn()`.
- Avoid `:::` and `<<-`. Prefer vectorised ops and Bioconductor structures (`GRanges`,
  `DataFrame`) over base loops.
- Keep functions short and single-purpose.
- **Do not add a new dependency without asking.**

## roxygen2 (every exported function)
- Needs `@title`, `@description`, `@param`, `@return`, and ≥ 1 runnable `@examples`.
- Examples use `exampleMouse` / `exampleTumourNormal`; avoid `\dontrun{}` unless network/files
  are genuinely required.
- Add `@seealso` cross-references.
- **`man/` and `NAMESPACE` are roxygen-generated — never hand-edit them.** Regenerate with
  `roxygen2::roxygenise()`, using the version pinned in `DESCRIPTION`'s
  `Config/roxygen2/version`. CI fails the build if they drift.

Any user-visible change must also be recorded in `NEWS.md` under the current devel heading
`# mesa X.Y.Z.9000`, and needs a test (see `tests` instructions).
