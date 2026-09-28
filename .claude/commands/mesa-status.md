---
description: Refresh mesa's project state and roadmap, and republish the Next steps page
allowed-tools: Bash(python3 .claude/scripts/mesa-status.py:*), Bash(cat .claude/state/status.json), Read, Write, Edit, Artifact
---

Refresh mesa's project state and roadmap, then say what to do next.

## 1. Regenerate the facts

```bash
python3 .claude/scripts/mesa-status.py --html --sync-issue
```

This rewrites `.claude/state/status.json` and `STATUS.md` from `git`, `gh`, `NEWS.md` and the
`gh-pages` build commit. It also derives every roadmap item's status from the TOML block in
#124, and rewrites #124's generated checklist when a status changed. Everything in them is
derived: never hand-edit `STATUS.md` or the generated half of #124, and never "correct" a
status. If one looks wrong, either the probe is wrong (fix the script) or the plan is (edit
the **Roadmap data** block in #124).

Read `STATUS.md`. If it has an **Incomplete data** section, say which fields are unknown
before drawing any conclusion from them.

## 2. Write the judgement

The script derives facts; the recommendation is yours. Overwrite
`.claude/state/recommendation.md` with what should happen next — the thing a maintainer
returning after two weeks would otherwise miss. Then re-run the script, again with `--html`,
so the judgement renders into `STATUS.md` and into the page step 3 publishes:

```bash
python3 .claude/scripts/mesa-status.py --html --no-log
```

Keep it to one decision plus at most two standing observations. Ground each in a specific
number or issue from the status. Do not repeat the page's "Do next" list: say what it
cannot, such as which of those items to start first and why, or a risk to the deadline.
Overwrite it — never append, or it becomes a changelog of stale advice. Good: *"Close #85
in favour of #102 — both rewrite `calculateEnrichment`."* Useless: *"Continue working on
open issues."*

## 3. Republish the page

```bash
python3 .claude/scripts/mesa-status.py --artifact-url
```

prints the page's URL. It is kept in the git common dir, so every worktree of this clone
sees the same one.

- **A URL is printed:** `Artifact` with `action: "read"` and that `url` first (required before
  updating an artifact this conversation has not published), then publish
  `.claude/state/dashboard.html` to the same `url` so the link stays stable.
- **Nothing is printed:** publish `.claude/state/dashboard.html` as a new artifact, then
  record it with `python3 .claude/scripts/mesa-status.py --set-artifact-url <url>`. Only
  reach for the `artifact-design` skill if the template itself needs work.

The page inlines the JSON. It declares **no** `capabilities`: a snapshot needs none, and
there is no GitHub connector for it to read live anyway. It shows the snapshot time and warns
visibly when it is over 24h old, so it can never quietly mislead.

Then clear the flags that asked for this refresh:

```bash
python3 .claude/scripts/mesa-status.py --mark-published
```

## 4. Report

Lead with the roadmap items whose status changed since the last refresh (e.g. "#114 In
review → Done"), then the recommendation, in the labelled style `AGENTS.md` sets out. Give
the page link once. Do not re-print `STATUS.md`.
