# CLAUDE.md

@AGENTS.md

**[`AGENTS.md`](AGENTS.md) is canonical** — coding standards, the branch/PR contract, the
release cycle, the toolchain version doctrine and the attribution rule all live there, and
it is imported above. This file adds only Claude-Code-specific notes. Where the two
disagree, `AGENTS.md` wins.

## How to answer

Lead with the result — no preamble, no restating the request, no closing recap. Use the
labelled parts in [`AGENTS.md`](AGENTS.md#output-format): one or two lines each, dropped
entirely when empty. A one-line fix gets a sentence and the commit subject, nothing more.
Do not narrate what you are about to do, and do not re-print tool output already on screen.
Say plainly what was verified and what could not run here.

## The contract, in one line

Branch off `main`, commit atomically, push the feature branch, open a **draft** PR if asked
— never mark ready, never merge, never touch `main`. `.claude/hooks/guard-remote.py`
enforces this; a blocked command means the contract said no, so do not work around it.

## Attribution

Every commit Claude makes carries `Co-Authored-By: Claude <noreply@anthropic.com>`, and any
PR, issue or comment Claude writes says plainly that Claude wrote it. This repo is public
and under Bioconductor review — authorship is never misrepresented.

## Where the detail lives

Load the skill rather than guessing; each one loads only when the task matches.

| Task | Skill |
|---|---|
| Version bumps, opening/closing a devel cycle | `bioc-release-cycle` |
| Tests, `R CMD check`, `BiocCheck`, coverage | `bioc-check-ladder` |
| Writing or fixing tests | `mesa-tests` |
| Roxygen docs, `NEWS.md` entries | `mesa-docs-news` |
| CI, devcontainer, toolchain versions | `mesa-ci` |
| Project state, `STATUS.md`, the status generator | `mesa-status` |
| Capturing a new procedure as a skill | `capture-skill` |

Generic R practice comes from the `posit-dev/skills` plugins (`r-lib`, `github`,
`open-source`, `posit-dev`). **A mesa skill beats a plugin skill** wherever they disagree —
the plugins know CRAN and tidyverse conventions, not Bioconductor's.

## Subagents

`bioc-reviewer` (read-only Bioconductor review of a diff), `test-author`, `release-manager`.
Use them for isolation on large or read-only jobs; do ordinary edits on the main thread.
