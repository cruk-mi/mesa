---
name: mesa-tests
description: mesa's testthat conventions — the cached qseaSet fixtures, the skip_long_checks option gate, and the skip guards for heavy genome/annotation packages. Use when writing, fixing or debugging tests in tests/testthat/.
---

# Testing mesa

`testthat` 3rd edition (`Config/testthat/edition: 3`). Tests live in `tests/testthat/`.

Every new function gets at least one test. Every bug fix gets a **regression test** that
fails before the fix and passes after — write it in that order and confirm it actually
fails first, or you have not tested the fix.

## Use the shared fixtures — do not rebuild example qsets

`tests/testthat/helper-fixtures.R` caches `qsea::getExampleQseaSet()` results across test
files within one R session, keyed by `(repl, expSamplingDepth)`:

```r
qs <- cachedExampleQset(repl = 8, expSamplingDepth = 1e5)
```

testthat sources `helper-*.R` once per run, before any `test_that()` block. Calling
`qsea::getExampleQseaSet()` directly re-runs the simulation and is the main cause of slow
test files. Reuse an existing `(repl, depth)` combination where you can — a new combination
costs a fresh simulation.

The package also ships `exampleMouse` and `exampleTumourNormal` data objects, which are what
roxygen `@examples` blocks should use.

## The long-check gate

`mesa:::skip_long_checks()` (defined in `R/utils.R`) skips when
`options(skip_long_checks = TRUE)` is set. `tests/testthat.R` sets that option, so **long
checks skip by default under `R CMD check`**, and around 15 tests currently depend on it.

```r
test_that("something slow", {
    skip_long_checks()
    # ...
})
```

CI runs the skipped long checks in a dedicated job so they are still verified — do not
assume a `skip` means untested.

## Skip guards for heavy data

Annotation and genome packages are `Suggests`, absent on the `slim` devcontainer and on
some CI legs. Guard them:

```r
skip_if_not_installed("org.Hs.eg.db")
skip_if_not_installed("BSgenome.Hsapiens.NCBI.GRCh38")
skip_if_not_installed("TxDb.Mmusculus.UCSC.mm10.knownGene")
```

Network-dependent tests (biomaRt) must also degrade gracefully — the existing pattern wraps
the call and skips on `try-error` rather than failing the suite when the service is down.

## Parallelism

`tests/testthat.R` registers `BiocParallel::MulticoreParam(workers = 1)`. Tests must not
assume more workers, and must not re-register a different backend without restoring it.

## Rules

- **Never change `R/` logic to make a test pass.** If a test and the code disagree, work out
  which is wrong and say so.
- Use explicit expectations — `expect_error(..., "message")` beats a bare `expect_error()`.
- BiocCheck wants ≥ 80 % coverage; check with `covr::package_coverage()` where available
  (see `bioc-check-ladder` — it is not installed everywhere).

## Precedence

The `r-lib` plugin's `testing-r-packages` skill is good general testthat 3e guidance and
applies here. This skill adds mesa's fixtures and gates on top; where they conflict, this
one wins.
