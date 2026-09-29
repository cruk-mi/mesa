---
name: Feature request
about: A new function, argument, or behaviour for mesa
title: "feat: <short description>"
labels: ["enhancement", "agent-ready"]
---

## Motivation
<!-- The analysis need this serves. -->

## Proposed behaviour
<!-- The function / argument and what it should do. Sketch the signature if you can. -->

```r
# e.g. makeDMRs(qset, minCpGs = 3)
```

## Scope / non-goals
<!-- Keep it small enough for one PR. Note anything explicitly out of scope. -->

## Acceptance criteria
- [ ] Implemented in `R/` following the Bioconductor coding standards
- [ ] roxygen complete (`@param`, `@return`, runnable `@examples`, `@seealso`); `man/` and
      `NAMESPACE` regenerated with `roxygen2::roxygenise()`, never hand-edited
- [ ] Test(s) added; `devtools::check()` + `BiocCheck::BiocCheck(".")` clean
- [ ] `NEWS.md` updated under the current `# mesa X.Y.Z.9000` heading
- [ ] No new dependency added without prior agreement
