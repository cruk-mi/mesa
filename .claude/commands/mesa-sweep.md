---
description: Weekly repo sweep - roadmap, release readiness, parking lot, fixed issues, stale branches
allowed-tools: Bash(git status:*), Bash(git fetch:*), Bash(git worktree list), Bash(git cherry:*), Bash(git log:*), Bash(gh issue list:*), Bash(gh issue view:*), Bash(gh pr list:*), Bash(gh pr view:*), Bash(python3 .claude/scripts/mesa-status.py --html --sync-issue --dry-run), Bash(python3 .claude/scripts/mesa-status.py --html --no-log), Bash(python3 .claude/scripts/mesa-status.py --artifact-url), Bash(python3 .claude/scripts/mesa-status.py --mark-published), Read, Write, Edit, Artifact
---

Run mesa's weekly sweep, once a week and before every release cut. The steps are listed
in order.

**Every GitHub write needs the human's yes first.** That covers syncing or editing #124,
creating, labelling or closing an issue, and deleting a branch. Show exactly what will
change, then wait. Only the read-only commands are pre-approved, so each write also
prompts.

## 0. Start clean

Run `git status` and `git fetch --prune origin`. Stop and say so if the tree is dirty.

## 1. Refresh the facts

```bash
python3 .claude/scripts/mesa-status.py --html --sync-issue --dry-run
```

If the output has an **Incomplete data** section, run it once more, because `gh` sometimes
fails briefly. If the section is still there, name the unknown fields and do not sync.

Otherwise, list the roadmap items whose status would change. After a yes, run
`python3 .claude/scripts/mesa-status.py --html --sync-issue`.

## 2. What blocks the next release

Answer from `STATUS.md`, in at most five lines:

- the next release's roadmap items that are not **Done**;
- the parking lot count, which must be 0 before the cut (see AGENTS.md, Triage);
- BiocCheck, where only `Invalid package Version` is expected during devel;
- whether the `NEWS.md` devel section has entries.

## 3. Empty the parking lot

Follow **Triage** in `AGENTS.md`:

1. Look for a duplicate of each item.
2. Propose one plan per item: file it, fold it into an existing issue, or drop it.
3. Wait for a yes.
4. File each approved item with `gh issue create --label parked`.
5. Delete the item's line only once its issue exists.

Park anything new you find during the sweep in the same way.

## 4. Close issues that are already fixed

List the open issues named by a merged PR's `Fixes`, `Closes` or `Resolves`. A PR that
merged into a branch other than `main` closes nothing on its own.

```bash
gh pr list --state merged --limit 50 --json number,body \
  -q '.[] | "\(.number) \(.body // "")"' \
  | grep -owiE '(close[sd]?|fix(e[sd])?|resolve[sd]?) #[0-9]+' | sort -u
```

These are GitHub's own closing keywords. A looser pattern also matches phrases like "a
session fixing #115".

Check each one with `gh issue view N`. Propose a close only when the PR fixes the whole
issue. A partial fix or a `Refs` stays open. After a yes, run
`gh issue close N --comment "Fixed by #PR. Closed by Claude during /mesa-sweep, with the maintainer's approval."`.

## 5. Update the roadmap data

Do this when the plan in #124 has changed: new issues from step 3, or stale `waiting_on`
and notes.

1. Save the issue body with `gh issue view 124 --json body -q .body > FILE`.
2. Edit only the TOML between the `roadmap:data` markers. Write the TOML to its own file
   too, as the dry run in step 3 reads it.
3. Dry-run it with
   `python3 .claude/scripts/mesa-status.py --roadmap-file TOML --sync-issue --dry-run`.
4. Show the diff and wait for a yes.
5. Re-read #124. If it changed since step 1, start again.
6. Apply it with `gh issue edit 124 --body-file FILE`, then run
   `python3 .claude/scripts/mesa-status.py --html --sync-issue`.

## 6. Prune stale branches

`STATUS.md` → **Branch hygiene** lists the safe branches and the ones with no PR. For each
safe branch:

- `gh pr list --state open --base <b>` must be empty. Deleting a base branch closes the
  PRs stacked on it.
- If `git worktree list` shows it checked out, the human removes that worktree first with
  `git worktree remove <path>`. This refuses a dirty worktree, and that refusal is correct.

For each branch with **no PR**, run `git cherry origin/main origin/<b>` and report how many
of its commits are not in `main`. Delete only the branches the human names.

Deleting a branch goes through `guard-remote.py`, which asks for confirmation. In `auto`,
`bypassPermissions` or `dontAsk` mode, which show no prompt, it blocks the deletion
instead. Then print the `git push origin --delete -- …` and `git branch -D -- …` lines for
the human to run.

## 7. Recommend and republish

Do steps 2 and 3 of `/mesa-status`: overwrite `.claude/state/recommendation.md`,
regenerate with `--html --no-log`, publish the page and run `--mark-published`.
`recommendation.md` is committed, so the change goes on a branch, never on `main`.

## 8. Report

Report in the `AGENTS.md` output format:

- the roadmap statuses that moved;
- what blocks the next release;
- issues filed and closed;
- branches deleted, and branches left for the human with the commands to run;
- the page link.
