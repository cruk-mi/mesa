---
name: pr-review
description: Copilot-style review of a mesa pull request or branch — findings ranked 🔴 High / 🟠 Medium / 🟡 Low with file:line, a failure scenario and a fix. Replaces GitHub Copilot's automatic review. Use when asked to review a PR, a branch or a diff, and when the review-after-pr hook asks for one after a large `gh pr create`.
---

# PR review

A local, read-only review that stands in for Copilot's. It runs on the Claude plan, not on
Copilot quota. Large PRs get one automatically: after `gh pr create`,
`.claude/hooks/review-after-pr.py` sizes the PR with one `gh pr view` call and asks for a
review if it's large. Small PRs (version bumps, one-line docs) get one only when a human
asks: `/pr-review <PR number>`.

**Read-only.** Do not edit, commit or push while reviewing. Fixes are a separate step,
after the human has read the review.

## 1. Scope the real diff

```bash
gh pr view <N> --json number,title,body,baseRefName,headRefName,headRefOid,mergeable,commits
gh pr diff <N>
```

For a branch without a PR, use `git diff origin/main...HEAD`.

Check the base **before** reading the diff. Copilot misses this, and it matters most:

- **Stale stacked base.** If the PR says it is stacked on another PR, or its commits
  include another PR's, check whether that PR merged (`gh pr view <base-PR> --json
  state,mergedAt`). A squash-merged base leaves its pre-squash commits on this branch:
  the diff balloons, the PR shows `CONFLICTING`, and a merge replays the old history. Find
  the last commit that belongs to the base PR and compare
  `git diff --stat <that-commit> <head>` (the real change) with the PR's file list.
  That is a 🔴 High finding. The fix is `git rebase --onto origin/main <that-commit>`,
  and a human runs the force-push, because guard-remote.py blocks it for agents.
- Otherwise review only lines the PR adds or changes. Read the surrounding code for
  context, but do not flag code the PR leaves untouched.

## 2. What to look for

- **Correctness.** Wrong logic, off-by-one, unhandled `NULL`/`NA`/empty input, a
  `tryCatch` that swallows the error it should raise, a claim in the PR body that the
  code doesn't do.
- **PR description** (see *Pull request descriptions* in `AGENTS.md`). Compare the body
  with the diff and with `.github/pull_request_template.md`:
  - 🟠 A **NEWS** block that differs from the `NEWS.md` lines in the diff, a **What
    changes** bullet the diff doesn't make, a user-visible change it doesn't list, no
    `Fixes #N`/`Refs #N`, or a `Version:` checkbox the diff contradicts.
  - 🟡 No NEWS section, a missing Why/What changes/Verified, over about 40 lines, pasted
    review history or logs, or an agent-written body with no AI note on its first line.
    An extra section (e.g. Timings) is fine.
  - A release-cycle PR (`chore/bump-version-*`) bumps `Version:` and usually has no
    issue: skip the `Fixes`/`Version:` checks. Its NEWS block is the heading it adds or
    renames.
- **Bioconductor and mesa R rules.** Apply the checklist in
  [`.claude/agents/bioc-reviewer.md`](../../agents/bioc-reviewer.md): coding standards,
  no `library()` in `R/`, roxygen completeness, hand-edited `man/`/`NAMESPACE`, missing
  tests or regression tests, `NEWS.md`, and no `Version:` bump inside a work PR.
- **Hooks and scripts** (`.claude/hooks`, `.claude/scripts`, `.devcontainer`, workflows):
  - Quoting, and user text passed through the shell inside double quotes (`$`, backticks).
  - Shell state assumed to survive between separate agent Bash calls. It doesn't.
  - Commands a slash command runs but its `allowed-tools` doesn't cover.
  - `set -u` hazards, bash 3.2 vs 5 differences, and BSD vs GNU `sed`/`awk`.
  - Anything that should exit 0 in CI but doesn't.
- **Agent contract.** `AGENTS.md`, `CLAUDE.md`, `guard-remote.py` and `settings.json`
  must agree. If a rule is changed in one, check the others. A rule only Claude Code can
  follow (a slash command, a hook) needs a fallback for the other agents that read
  `AGENTS.md`.
- **Security.** Secrets in code, `.Renviron` contents, and new dependencies nobody agreed
  on.

Skip pure style nits and praise. When a PR is clean, say so: a review that invents
problems to look thorough is worse than no review.

## 3. Severity

| Tag | Meaning |
|---|---|
| 🔴 **High** | Bug, crash, data loss, security issue, or a PR that can't merge as it stands. Fix before merge. |
| 🟠 **Medium** | Likely defect, or a violation of a Bioconductor, mesa or contract rule. |
| 🟡 **Low** | Minor correctness or maintainability issue. |

Each finding says three things: what is wrong, a concrete failure scenario, and the fix.
Cite `file:line`.

## 4. Report

**By default, report locally only.** Give the findings in chat, most severe first,
grouped by severity. Say what was checked and found clean in one line.

**Posting to GitHub is outward-facing: only do it when the human asks** ("post it", or
`--post` in the arguments). Then:

1. Post each line-anchored finding as an inline comment. Its body starts with the
   severity tag:
   ```bash
   gh api repos/cruk-mi/mesa/pulls/<N>/comments --input - <<'EOF'
   {"body": "🟠 **Medium**: …", "commit_id": "<headRefOid>", "path": "<file>", "line": <n>, "side": "RIGHT"}
   EOF
   ```
   Build the JSON with a script (`json.dumps`) rather than by hand when a body contains
   quotes or code. Use `/pulls/<N>/comments`, never `/pulls/<N>/reviews`: guard-remote.py
   blocks the reviews endpoint because that endpoint can approve.
2. Post findings not tied to a line (a stale base, an outdated PR body) in one summary
   comment through `gh pr comment <N> --body-file <file>`. Head it `## Claude review` and
   give the counts per severity and a one-line verdict.
3. Every comment says Claude wrote it (see Attribution in `AGENTS.md`).

## Lighter alternative

The built-in `/code-review <N>` (add `--comment` to post) runs a generic correctness
review without mesa's rubric. It's enough for a quick second opinion.

## Turning Copilot off

No ruleset on cruk-mi/mesa requests Copilot reviews, so automatic Copilot reviews come
from each person's own setting: **github.com/settings/copilot → Automatic Copilot code
review**. Turn it off there.
