---
description: Refresh mesa's project state (STATUS.md) and republish the dashboard
allowed-tools: Bash(python3 .claude/scripts/mesa-status.py:*), Bash(cat .claude/state/status.json), Bash(cat .claude/state/artifact-url.txt), Read, Write, Edit, Artifact
---

Refresh mesa's project state, then say what to do next.

## 1. Regenerate the facts

```bash
python3 .claude/scripts/mesa-status.py --html
```

This rewrites `.claude/state/status.json` and `STATUS.md` from `git`, `gh`, `NEWS.md` and the
`gh-pages` build commit. Everything in them is derived — never hand-edit `STATUS.md`, and
never "correct" a number in it. If a figure looks wrong, the probe is wrong: fix the script.

Read `STATUS.md`. If it has an **Incomplete data** section, say which fields are unknown
before drawing any conclusion from them.

## 2. Write the judgement

The script derives facts; the recommendation is yours. Overwrite
`.claude/state/recommendation.md` with what should happen next — the thing a maintainer
returning after two weeks would otherwise miss. Then re-run the script, again with `--html`,
so the judgement renders into `STATUS.md` and into the page step 3 publishes:

```bash
python3 .claude/scripts/mesa-status.py --html
```

Keep it to one decision plus at most two standing observations. Ground each in a specific
number or issue from the status. Overwrite it — never append, or it becomes a changelog of
stale advice. Good: *"Close #85 in favour of #102 — both rewrite `calculateEnrichment`."*
Useless: *"Continue working on open issues."*

## 3. Republish the dashboard

`.claude/state/artifact-url.txt` holds the dashboard's URL.

- **File exists:** `Artifact` with `action: "read"` and that `url` first (required before
  updating an artifact this conversation has not published), then publish to the same `url`
  so the link stays stable.
- **File missing:** publish `.claude/state/dashboard.html` — step 1 has just rebuilt it from
  the template — and write the returned URL to that file. It is gitignored: the URL is
  per-person. Only reach for the `artifact-design` skill if the template itself needs work.

The page inlines the JSON. It declares **no** `capabilities`: a snapshot needs none, and
there is no GitHub connector for it to read live anyway. Show the snapshot time and warn
visibly when it is over 24h old, so it can never quietly mislead.

## 4. Report

Lead with what changed since the last refresh, then the recommendation, in the labelled
style `AGENTS.md` sets out. Give the dashboard link once. Do not re-print `STATUS.md`.
