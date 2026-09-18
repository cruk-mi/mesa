---
name: BiocCheck / R CMD check cleanup
about: A specific check NOTE/WARNING to clear before Bioconductor submission
title: "chore: clear <NOTE/WARNING> in <file>"
labels: ["bioccheck", "agent-ready"]
---

## Finding
<!-- Paste the exact NOTE/WARNING from BiocCheck or R CMD check. -->

```

```

## Where
<!-- File(s) / function(s) involved, e.g. R/plotting.R lines 40–55. -->

## Suggested fix
<!-- Optional: the deterministic change (line wrap, seq_len, T/F, etc.). -->

## Acceptance criteria
- [ ] Finding cleared with the smallest possible diff (no logic refactor)
- [ ] `BiocCheck::BiocCheck(".")` (and `devtools::check()` if relevant) shows it gone
- [ ] No new findings or test regressions introduced
