---
name: Bug report
about: A reproducible bug in mesa
title: "fix: <short description>"
labels: ["bug", "agent-ready"]
---

## What happens
<!-- The incorrect behaviour. -->

## Expected behaviour
<!-- What should happen instead. -->

## Minimal reproducible example
```r
library(mesa)
# Use exampleMouse / exampleTumourNormal where possible so an agent can reproduce
# without external data.

```

## Affected function(s) / file(s)
<!-- e.g. R/makeDMRs.R — helps scope the fix. -->

## Acceptance criteria
- [ ] Root cause fixed with the smallest diff
- [ ] Regression test added (fails before, passes after)
- [ ] `devtools::check()` and `BiocCheck::BiocCheck(".")` clean
- [ ] `NEWS.md` updated under the current `# mesa X.Y.Z.9000` heading

## Session info
<details><summary><code>sessionInfo()</code></summary>

```

```
</details>
