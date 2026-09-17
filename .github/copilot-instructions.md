# GitHub Copilot instructions for mesa

**Follow [`AGENTS.md`](../AGENTS.md) — it is the canonical source of truth** for coding
standards, branch strategy, commit conventions, the versioning cadence, and the workflow
contract. This file only adds Copilot-specific notes.

## Copilot-specific notes

- **mesa** is an R/Bioconductor package (MBD-seq / MeDIP-seq methylation analysis) targeting
  Bioconductor submission. Conform to Bioconductor policies at all times.
- **Sandbox workflow.** Work on a branch cut from `main`. Make **atomic** Conventional-Commit
  commits. **Keep your pull request in draft** — do not mark it ready for review, request
  review, or merge until a human explicitly asks. Never push or merge to `main`.
- **Attribution.** Keep the default `Co-authored-by: Copilot` trailer on commits and note
  AI authorship on PRs/issues. Authorship must never be misrepresented (public Bioc repo).
- **Smallest diff that works.** Do not refactor unrelated code. Do not add dependencies
  without asking. `man/` and `NAMESPACE` are roxygen-generated — never hand-edit them.
- **Path-specific rules** live in `.github/instructions/*.instructions.md` and apply
  automatically by file glob (R code, tests, docs, CI/devops).
- **Specialized agents** for recurring tasks live in `.github/agents/*` — use the one that
  matches the job: `bioc-reviewer` (read-only review), `test-author` (tests),
  `release-manager` (version bumps). Ordinary feature work is done directly, not delegated.
- **Validate** with `devtools::test()`, `devtools::check()`, and `BiocCheck::BiocCheck(".")`
  before declaring work done. The cloud environment is prepared by
  `.github/workflows/copilot-setup-steps.yml`.
