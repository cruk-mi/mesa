#!/usr/bin/env python3
"""PostToolUse hook: ask for a pr-review after `gh pr create` opens a large PR.

Small PRs (a version bump, a one-line docs fix) are not worth a review's time,
so the hook stays silent for them and a human runs `/pr-review <N>` if they
want one. Large PRs get the review automatically: the hook hands the session
an instruction to run the `pr-review` skill before it stops.

Size is measured locally against the merge base, with no network call. Files
that are generated or bookkeeping (man/*.Rd, NAMESPACE, DESCRIPTION, NEWS.md,
.claude/state/) are left out, so a version bump measures zero.

Never blocks and never fails the tool call: anything unexpected (no PR URL in
the output, an unresolvable base, git missing) means silence. Skipped in CI.

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
PR_URL = re.compile(r"https://github\.com/[^/\s]+/[^/\s]+/pull/(\d+)")


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


def option(args, names):
    """Value of the first of `names` (`--base x`, `-B x`, `--base=x`)."""
    for i, tok in enumerate(args):
        if tok in names and i + 1 < len(args):
            return args[i + 1]
        for name in names:
            if name.startswith("--") and tok.startswith(name + "="):
                return tok.split("=", 1)[1]
    return None


def git(cwd, *args):
    try:
        out = subprocess.run(["git", *args], cwd=cwd, capture_output=True,
                             text=True, timeout=15)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return out.stdout.strip() if out.returncode == 0 else None


def resolve(cwd, ref):
    """Prefer the remote-tracking ref, so a stale local branch is not used."""
    for candidate in (f"origin/{ref}", ref):
        if git(cwd, "rev-parse", "--verify", "--quiet", candidate + "^{commit}"):
            return candidate
    return None


def measure(cwd, base, head):
    """(lines, files, hot_lines) counted between merge-base(base, head) and head."""
    merge_base = git(cwd, "merge-base", base, head)
    if not merge_base:
        return None
    numstat = git(cwd, "diff", "--numstat", merge_base, head)
    if numstat is None:
        return None
    lines = files = hot = 0
    for row in numstat.splitlines():
        parts = row.split("\t")
        if len(parts) != 3:
            continue
        added, deleted, path = parts
        # A rename shows as "old => new"; judge the new path.
        path = re.sub(r"\{[^}]* => ([^}]*)\}", r"\1", path).split(" => ")[-1]
        if IGNORED.match(path):
            continue
        changed = (int(added) if added.isdigit() else 0) + (int(deleted) if deleted.isdigit() else 0)
        files += 1
        lines += changed
        if path.startswith(HOT_PATHS):
            hot += changed
    return lines, files, hot


def main():
    if os.environ.get("CI") == "true" or os.environ.get("GITHUB_ACTIONS") == "true":
        return
    try:
        data = json.load(sys.stdin)
    except ValueError:
        return
    if data.get("tool_name") != "Bash":
        return
    args = pr_create_args((data.get("tool_input") or {}).get("command") or "")
    if args is None:
        return
    # No PR URL in the output means the create failed: nothing to review.
    match = PR_URL.search(json.dumps(data.get("tool_response") or ""))
    if not match:
        return
    number = match.group(1)

    cwd = data.get("cwd") or os.getcwd()
    base = resolve(cwd, option(args, ("--base", "-B")) or "main")
    head_name = option(args, ("--head", "-H"))
    head = resolve(cwd, head_name) if head_name else "HEAD"
    if not base or not head:
        return
    size = measure(cwd, base, head)
    if size is None:
        return
    lines, files, hot = size
    if lines <= MAX_LINES and files <= MAX_FILES and hot <= MAX_HOT_LINES:
        return

    context = (
        f"PR #{number} is large ({lines} changed lines across {files} files, "
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
