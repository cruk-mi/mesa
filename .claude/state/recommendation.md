**Close #85 in favour of #102.** Both rewrite `calculateEnrichment` to drop MEDIPS; #102 is
the rebased successor and carries #85's fixes plus the four behaviour corrections. Keeping
both open invites merging the older one.

Two standing signals worth a look, neither urgent:

- **Coverage uploads have stopped.** Codecov last heard from us on 2026-06-17 (95 days), yet
  `main` has taken pushes since. `covr::codecov()` runs only when the push touches `R/`,
  `tests/`, `vignettes/`, `inst/`, `DESCRIPTION` or `NAMESPACE`, so some of that gap is
  expected — but 48.45% against the 80% target in `AGENTS.md` is the real gap.
- **The single BiocCheck ERROR is benign.** "Invalid package Version" is BiocCheck objecting
  to the four-part devel version `0.99.6.9000`; it clears when the cycle is cut to `0.99.7`.
  Do not "fix" it.
