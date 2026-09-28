#!/usr/bin/env python3
"""PostToolUse hook: ask for a pr-review after `gh pr create` opens a large PR.

Small PRs (a version bump, a one-line docs fix) are not worth a review's time,
so the hook stays silent for them and a human runs `/pr-review <N>` if they
want one. Large PRs get the review automatically: the hook hands the session
an instruction to run the `pr-review` skill before it stops.

Size comes from GitHub, not from the session's checkout: the PR number and
repository are read from the URL `gh pr create` prints, then one
`gh pr view <N> --json baseRefName,files` call (15 s timeout) lists the files
and their added and deleted lines. The session's cwd, branch and the command's
`--base`/`--head` are never consulted, so `cd <worktree> && gh pr create` and a
stacked PR are measured as GitHub sees them. Files that are generated or
bookkeeping (man/*.Rd, NAMESPACE, DESCRIPTION, NEWS.md, .claude/state/) are left
out, so a version bump measures zero.

Never blocks and never fails the tool call: anything unexpected (no PR URL in
the output, `gh` missing, failing, slow or printing something unparsable) means
silence. Skipped in CI.

Exit code is always 0. Output, when there is any, is PostToolUse JSON with
`additionalContext`.
"""

import json
import os
import re
import shlex
import subprocess
import sys

# A PR is "large" when any of these is exceeded. Tune here.
MAX_LINES = 150       # counted added + deleted lines
MAX_FILES = 5         # counted files
HOT_PATHS = ("R/", ".claude/hooks/", ".github/workflows/")
MAX_HOT_LINES = 30    # counted lines under HOT_PATHS

IGNORED = re.compile(r"^(man/.*\.Rd|NAMESPACE|DESCRIPTION|NEWS\.md|\.claude/state/.*)$")
PR_URL = re.compile(r"https://github\.com/([^/\s]+/[^/\s]+)/pull/(\d+)")


def split_segments(command):
    """Split a shell command into separately-executed segments."""
    return [s for s in re.split(r"\|\||&&|;|\n|\||&", command) if s.strip()]


def tokenise(segment):
    try:
        return shlex.split(segment)
    except ValueError:
        return segment.split()


def pr_create_args(command):
    """The arguments after `gh pr create`, or None if no segment runs it."""
    for segment in split_segments(command):
        tokens = tokenise(segment)
        for i in range(len(tokens) - 2):
            if (os.path.basename(tokens[i]) == "gh"
                    and tokens[i + 1] == "pr" and tokens[i + 2] == "create"):
                return tokens[i + 3:]
    return None


def pr_size(repo, number):
    """(base, lines, files, hot_lines) for the PR as GitHub reports it, else None."""
    try:
        out = subprocess.run(
            ["gh", "pr", "view", number, "--repo", repo, "--json", "baseRefName,files"],
            capture_output=True, text=True, timeout=15)
        if out.returncode != 0:
            return None
        info = json.loads(out.stdout)
        entries = info["files"]
    except (OSError, subprocess.TimeoutExpired, ValueError, KeyError, TypeError):
        return None
    lines = files = hot = 0
    for entry in entries:
        path = entry.get("path") or ""
        if IGNORED.match(path):
            continue
        changed = int(entry.get("additions") or 0) + int(entry.get("deletions") or 0)
        files += 1
        lines += changed
        if path.startswith(HOT_PATHS):
            hot += changed
    return info.get("baseRefName") or "?", lines, files, hot


def main():
    if os.environ.get("CI") == "true" or os.environ.get("GITHUB_ACTIONS") == "true":
        return
    try:
        data = json.load(sys.stdin)
    except ValueError:
        return
    if data.get("tool_name") != "Bash":
        return
    if pr_create_args((data.get("tool_input") or {}).get("command") or "") is None:
        return
    # No PR URL in the output means the create failed: nothing to review.
    response = data.get("tool_response") or ""
    stdout = (response.get("stdout") or "") if isinstance(response, dict) else str(response)
    match = PR_URL.search(stdout)
    if not match:
        return
    repo, number = match.groups()

    size = pr_size(repo, number)
    if size is None:
        return
    base, lines, files, hot = size
    if lines <= MAX_LINES and files <= MAX_FILES and hot <= MAX_HOT_LINES:
        return

    context = (
        f"PR #{number} (into {base}) is large ({lines} changed lines across {files} files, "
        f"{hot} of them under {', '.join(HOT_PATHS)}; generated files not counted). "
        f"Before stopping, run the `pr-review` skill on PR #{number} and report the "
        "findings here. Do not post them to GitHub unless the human asks."
    )
    print(json.dumps({"hookSpecificOutput": {
        "hookEventName": "PostToolUse",
        "additionalContext": context,
    }}))


if __name__ == "__main__":
    try:
        main()
    except Exception:  # a review nudge must never break the session
        pass
    sys.exit(0)
