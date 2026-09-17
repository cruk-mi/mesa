---
applyTo: "man/**,vignettes/**,README.md,NEWS.md"
---

# Documentation rules

Follow [`AGENTS.md`](../../AGENTS.md). For docs:

- **`man/*.Rd` are generated** from roxygen2 comments in `R/`. Do not hand-edit `.Rd` files;
  change the roxygen above the function and rebuild man pages **manually**
  with `roxygen2::roxygenise()` (**never** `devtools::document()`, which rewrites more
  than expected). Check `packageVersion("roxygen2")` matches `DESCRIPTION`'s `RoxygenNote`
  first; if it does not, stop and say so.
- **`NAMESPACE` is hand-edited** — keep it in sync manually when exports change.
- Vignettes (`vignettes/*.Rmd`) must knit cleanly with `knitr`; keep runnable chunks light
  enough for the slim environment, or guard heavy data with `eval=FALSE` / availability checks.
- `NEWS.md`: record user-visible changes under the current devel heading
  `# mesa X.Y.Z.9000`. On release the `release-manager` agent renames it to `# mesa X.Y.Z`.
- `README.md`: keep examples consistent with the current API.
