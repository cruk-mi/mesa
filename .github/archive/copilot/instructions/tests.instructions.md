---
applyTo: "tests/**"
---

# Test rules

Follow [`AGENTS.md`](../../AGENTS.md). For files under `tests/`:

- Use `testthat` (≥ 3rd edition); tests live in `tests/testthat/`.
- Every new function gets at least one test; **every bug fix gets a regression test** that
  fails before the fix and passes after.
- Use explicit expectations: `expect_equal()`, `expect_error()`, `expect_warning()`, etc.,
  with informative messages.
- Guard heavy genome/annotation data with `skip_if_not_installed()` (and the package's
  long-check skip helpers) so the slim CI/agent environment still passes.
- Use the bundled `exampleMouse` / `exampleTumourNormal` objects for fixtures.
- Target ≥ 80 % coverage (BiocCheck requirement); check with `covr::package_coverage()`.
- Run `devtools::test()` before committing.
