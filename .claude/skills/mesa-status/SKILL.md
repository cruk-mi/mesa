---
name: mesa-status
description: How mesa's project state and roadmap are derived — STATUS.md, the status generator, the roadmap data in #124, the BiocCheck history, the Next steps page and the hooks that keep it current. Use when reading or refreshing project status, when a figure or a roadmap status looks wrong, when changing the plan, or when adding a signal to the status view.
---

# Project state

`STATUS.md` answers "what landed, what is in flight, what is next". It is **generated**, so
it cannot drift: every figure is read back from the places humans and agents both actually
write to.

```
git + gh + NEWS.md + gh-pages + #124 roadmap TOML
            |
  .claude/scripts/mesa-status.py        deterministic, no model involved; read-only except
                                        --sync-issue, a human-approved edit of #124
            |
  .claude/state/status.json             intermediate
            |
    +-------+--------+
    |                |
STATUS.md      dashboard.html
what Claude    what you look at, from
reads          disk or published
```

**All three outputs are gitignored.** They are derived, so committing one would put a
snapshot in the repo to go stale — the exact drift this design exists to avoid. Only the two
genuine inputs are committed: `.claude/state/recommendation.md` and
`.claude/state/bioccheck-history.jsonl`. A fresh clone has no `STATUS.md` until the script
runs, which the `SessionStart` hook does anyway.

## The rule

**Never hand-edit `STATUS.md`.** Never "correct" a number in it. A wrong figure means a
wrong probe — fix the script and re-run. Editing the output just means the next refresh
silently discards your correction.

The single hand-written input is `.claude/state/recommendation.md`: the "what to do next"
judgement a script cannot derive. The script renders it into `STATUS.md`. It is
**overwritten**, never appended to, or it decays into a log of stale advice. It is committed
— being an input, not an output — so regenerating on any machine reproduces the same
`STATUS.md`.

## Refreshing

```bash
python3 .claude/scripts/mesa-status.py              # facts only
python3 .claude/scripts/mesa-status.py --no-log     # skip the slow CI-log probe
python3 .claude/scripts/mesa-status.py --max-age N  # no-op if STATUS.md is under N seconds old
```

`/mesa-status` does the whole job: regenerate, show the diff #124 would get (and sync it
only when the human says yes), rewrite the recommendation and republish the page. A real
`--sync-issue` is kept out of the permission allowlist, so it always prompts. A `SessionStart` hook runs `--max-age 14400 --no-log --quiet`, so state
is current at the start of any session without paying a `gh` round-trip every time.

### Keeping the published page current

The page can only be republished by the model (the Artifact tool is not a script's), so
`.claude/hooks/refresh-flag.py` reminds the human at the right moment. The hooks never
have the model run `/mesa-status` itself, because it can edit the public #124:

| Hook | What it does |
|---|---|
| `PostToolUse` (Bash) | A command that changed GitHub (`git push`, `git tag`, `gh pr/issue create·edit·merge·close…`, `gh label`, `gh release`, a `gh api` write) raises `refresh-needed`. |
| `Stop` | If the flag is up, suggests `/mesa-status` to the human **once** (a `systemMessage`, never a block), lowering the flag to `refresh-pending` so it does not repeat every turn. |
| `SessionStart` | Tells the human about a refresh left pending by an earlier session. |

The flags and the page URL live in `<git common dir>/mesa-status/`, shared by every worktree
of the clone and never committed. `--mark-published` clears the flags. Changes made on
github.com raise no flag: the page's ">24h old" banner is the safety net for those.
`bash .claude/hooks/test-refresh-flag.sh` covers the hook and the roadmap and #124 sync logic
in `mesa-status.py`. `PYTHON=/usr/bin/python3 bash .claude/hooks/test-refresh-flag.sh` runs it
under macOS's Python 3.9, where the roadmap degrades (it needs 3.11's `tomllib`) and
everything else still works.

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

mesa squash-merges by default, so a landed branch's commits usually never appear on
`main` and `git branch --merged` cannot see that it landed. Classification therefore comes from PR
state, not from git — and a finished PR only vouches for the commit it ended on, so a
branch whose tip has moved past its PR's head is listed as unknown, not prunable. The output is a ready-to-paste `git branch -D` / `git push --delete`
block for the human to run. Per `AGENTS.md`, an agent deletes a branch only with the
human's explicit approval and only when it is safe. `guard-remote.py` enforces this: it
asks the human to confirm, and it refuses a branch that an open PR uses as its head or
base.

## The roadmap

The plan lives in the pinned roadmap issue **#124**, as a TOML block between
`<!-- roadmap:data:start -->` and `<!-- roadmap:data:end -->`. It holds only what GitHub
cannot know: the items, their order, wave, release, owner, dependencies and notes. The
comment at the top of the block documents every field. Maintainers edit it on github.com;
no PR is needed.

**No status is ever written by hand.** `collect_roadmap()` derives each item's status, and
`--sync-issue` rewrites the generated half of #124 (between the `roadmap:generated`
markers) only when a status changed:

| Status | Means | Derived when (first match wins) |
|---|---|---|
| **Done** | Merged, released, closed or decided | every PR ref merged, every issue ref closed and every `check` passes, or `done = true` (a PR closed without merging does not count) |
| **Later** | Deliberately after the Bioc release | the item's wave has `later = true` |
| **Waiting** | Blocked; the page says on what | `waiting_on` is set, or (below) an `after` item is not Done |
| **In review** | A PR is open | a ref is an open PR, or an open issue has an open PR that says `Fixes #N` |
| **Waiting** | … on another item | an `after` item is not Done |
| **To do** | Nobody has started it | none of the above |

`owner = "maintainer"` marks items a human has to act on (decisions, releases, the Bioc
push). It is a separate tag, not a status, so "delete branches" reads as **To do ·
Maintainer** until `git ls-remote` shows the branches gone, then **Done**.

Checks: `branch-gone:NAME` (the branch no longer exists on origin) and `tag:vX.Y.Z` (the
tag exists on origin). A release's row is "released" once its `vX.Y.Z` tag exists.

A TOML typo degrades only the roadmap: the page says so and the rest still renders. To try
a change before editing #124, save the block to a file and run
`mesa-status.py --roadmap-file FILE --sync-issue --dry-run`, which prints the diff #124
would get.

`--sync-issue` refuses, dry run included, when anything is incomplete:
- the roadmap or #124's current body could not be read;
- any probe degraded;
- any ref state or check is unknown.

#124 is public, and a single failed lookup would otherwise turn Done items back into To
do there. Fix the cause and run it again.

## Adding a signal

Add a `collect_*()` that returns a plain dict and degrades to `None`/`unknown` on failure,
call it in `main()`, and render it in `render()`. Two constraints:

1. **Idempotent.** Two runs with no repo change must produce a byte-identical `STATUS.md`
   apart from the trailing timestamp, so a refresh is a no-op when nothing moved. Provenance
   that varies between runs (cache vs live) belongs in `status.json`, not `STATUS.md`.
2. **Degrades quietly.** Call `degrade("why this is unknown")` and carry on. Never raise.

## The Next steps page

The page has two tabs: **Next steps** (the roadmap, "Do next", the recommendation, a
timeline, a status filter and the release calendar) and **Health** (CI, coverage, BiocCheck,
pkgdown, PRs in flight, NEWS, branch hygiene). `#health` in the link opens the second tab.

The page is **generated by the same script**, not written by a model each run:

```
.claude/scripts/dashboard-template.html     committed — the page: CSS, layout, rendering JS
            |                               carries one /*__STATUS_JSON__*/ marker
            |  mesa-status.py --html  (render_html: one string substitution)
            v
.claude/state/dashboard.html                gitignored — template + inlined status.json
            |
            v
      published to claude.ai
```

Same JSON in, byte-identical HTML out. The only thing that varies between real runs is the
`generated` timestamp the probes stamp into the state.

**To change how the page looks, edit the template.** Editing `.claude/state/dashboard.html`
is pointless — the next run overwrites it. The template renders from the `status.json` shape,
so a new probe needs a matching branch in the template's script to appear on the page.

You can open `.claude/state/dashboard.html` straight from disk in a browser; publishing is
only what makes it reachable from elsewhere.

`/mesa-status` publishes it to the URL that `mesa-status.py --artifact-url` prints (kept in
the git common dir, so each person's clone has its own page and all its worktrees share
it) rather than spawning a new artifact each time. The first run migrates a URL from the
old `.claude/state/artifact-url.txt`. Publishing to an artifact a conversation has not read
requires an `action: "read"` first.

It is a **snapshot with the data inlined**, not a live view: a published artifact can only
reach claude.ai connectors, and there is no GitHub connector, so it cannot query the repo
itself. It therefore shows its snapshot time and warns when over 24h old.

## The Parked section

`collect_parked()` reads the parking lot that `/park` writes (see "Session scope and the
parking lot" in `AGENTS.md`): one `- ` line per out-of-scope finding. The file lives in the
**main checkout's** `.claude/state/parking-lot.md`, resolved through
`git rev-parse --git-common-dir`, so a refresh run from any worktree sees the same items.
An empty or missing file renders as "Nothing parked".

Items are written by `mesa-status.py --park <issue>`, which reads the note from stdin. It
never goes through the shell, and newlines are collapsed so each finding stays one line. The
flag appends and exits; it does not refresh.

The lot is gitignored and **per-machine**. Sessions in the devcontainer, a Codespace or on
the web lose it with their checkout, so the PR description's **Parked** list is the durable
copy. It is not private: the items go into `status.json` and are inlined into
`dashboard.html`, which is published as a claude.ai artifact.
