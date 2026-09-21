---
name: mesa-status
description: How mesa's project state is derived — STATUS.md, the status generator, the BiocCheck history and the dashboard. Use when reading or refreshing project status, when a figure in STATUS.md looks wrong, or when adding a signal to the status view.
---

# Project state

`STATUS.md` answers "what landed, what is in flight, what is next". It is **generated**, so
it cannot drift: every figure is read back from the places humans and agents both actually
write to.

```
git + gh + NEWS.md + gh-pages
            |
  .claude/scripts/mesa-status.py        read-only, deterministic, no model involved
            |
  .claude/state/status.json             gitignored intermediate
            |
    +-------+--------+
    |                |
STATUS.md      dashboard artifact
committed      private claude.ai page
```

## The rule

**Never hand-edit `STATUS.md`.** Never "correct" a number in it. A wrong figure means a
wrong probe — fix the script and re-run. Editing the output just means the next refresh
silently discards your correction.

The single hand-written input is `.claude/state/recommendation.md`: the "what to do next"
judgement a script cannot derive. The script renders it into `STATUS.md`. It is
**overwritten**, never appended to, or it decays into a log of stale advice. It is committed,
so regenerating on any machine reproduces `STATUS.md` byte for byte.

## Refreshing

```bash
python3 .claude/scripts/mesa-status.py              # facts only
python3 .claude/scripts/mesa-status.py --no-log     # skip the slow CI-log probe
python3 .claude/scripts/mesa-status.py --max-age N  # no-op if STATUS.md is under N seconds old
```

`/mesa-status` does the whole job: regenerate, rewrite the recommendation, republish the
dashboard. A `SessionStart` hook runs `--max-age 14400 --no-log --quiet`, so state is current
at the start of any session without paying a `gh` round-trip every time.

## Degradation is deliberate

There is no `set -e` equivalent here: every probe fails on its own and the run still
succeeds. A missing `gh`, an expired token or a dead network gives `unknown` fields and an
**Incomplete data** section listing why. Both outputs are written atomically, so an
interrupted run never leaves a half-written `STATUS.md`.

**When reading a degraded `STATUS.md`, say which fields are unknown before drawing any
conclusion from them.**

## The probes, and what is fragile about them

| Signal | Source | Fragility |
|---|---|---|
| Version, cycle phase | `DESCRIPTION` | none |
| Devel changes | first `# mesa X.Y.Z` block of `NEWS.md` | none |
| In flight, landed, next up | `gh pr list`, `gh issue list` | none |
| CI on main | `gh run list --workflow=R-CMD-check-bioc` | see below |
| Coverage | Codecov public API (no token) | endpoint shape |
| BiocCheck counts | parsed from the CI log | **log format** |
| Site freshness | `origin/gh-pages` commit subject | pkgdown's message format |
| Branch hygiene | branch names vs `gh pr list --state all` | none |

**BiocCheck counts** are the fragile one. They cannot be recomputed locally (see
`bioc-check-ladder`), so they are parsed out of the `Run BiocCheck` step of the latest
successful `main` run, matching the summary line `✖ 1 ERRORS | ⚠ 0 WARNINGS | ℹ 5 NOTES`.
Findings that wrap onto a second log line are joined back together. If that format changes,
the run degrades to "not parsed" rather than inventing numbers. Results are appended to
`.claude/state/bioccheck-history.jsonl` (committed, append-only) and keyed by commit SHA, so
the expensive log download happens once per checked commit and the trend survives re-clones.

**A BiocCheck `Invalid package Version` ERROR during devel is expected**, not a defect:
BiocCheck rejects the four-part `X.Y.Z.9000` form. It clears when the cycle is cut. Do not
try to fix it.

**CI on main is often behind main on purpose.** `check-bioc.yml` only triggers on pushes
touching `R/`, `tests/`, `vignettes/`, `inst/`, `DESCRIPTION` or `NAMESPACE`, so a docs- or
tooling-only commit leaves `main` legitimately un-run. The status view says how many commits
past the last run `main` is, so this reads as expected rather than as a failure.

## Branch hygiene never deletes

mesa squash-merges, so a landed branch's commits never appear on `main` and
`git branch --merged` cannot see that it landed. Classification therefore comes from PR
state, not from git. The output is a ready-to-paste `git branch -D` / `git push --delete`
block for the human to run: per `AGENTS.md`, agents do not delete branches, and
`guard-remote.py` blocks it mechanically.

## Adding a signal

Add a `collect_*()` that returns a plain dict and degrades to `None`/`unknown` on failure,
call it in `main()`, and render it in `render()`. Two constraints:

1. **Idempotent.** Two runs with no repo change must produce a byte-identical `STATUS.md`
   apart from the trailing timestamp, or every refresh churns a committed file. Provenance
   that varies between runs (cache vs live) belongs in `status.json`, not `STATUS.md`.
2. **Degrades quietly.** Call `degrade("why this is unknown")` and carry on. Never raise.

## Dashboard

`/mesa-status` publishes a private claude.ai artifact — the same data, laid out to scan. Its
URL lives in `.claude/state/artifact-url.txt` (gitignored, per-person) so refreshes update
the same page instead of spawning new ones.

It is a **snapshot with the data inlined**, not a live view: a published artifact can only
reach claude.ai connectors, and there is no GitHub connector, so it cannot query the repo
itself. It therefore shows its snapshot time and warns when over 24h old.
