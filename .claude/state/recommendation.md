**Review #102 next.** It is the only open PR that changes package behaviour (drops MEDIPS
from `calculateEnrichment`, fixes #81), and its CI is green. The tooling PRs can follow.

Two standing signals, neither urgent:

- **Coverage is well below target.** 48.45% against the 80% in `AGENTS.md`, and Codecov last
  heard from us on 2026-06-17. Uploads run only when a push touches `R/`, `tests/`,
  `vignettes/`, `inst/`, `DESCRIPTION` or `NAMESPACE`, so part of the gap is expected.
- **The single BiocCheck ERROR is benign.** "Invalid package Version" is BiocCheck objecting
  to the four-part devel version `0.99.6.9000`; it clears when the cycle is cut to `0.99.7`.
  Do not "fix" it.
