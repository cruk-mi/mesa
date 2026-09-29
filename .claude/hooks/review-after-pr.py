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
# A bare PR URL on a line of its own, as `gh pr create` prints it. Not an
# `#issuecomment-…` / `#discussion_r…` link, and not a push's `/pull/new/<branch>`.
PR_URL = re.compile(r"^https://github\.com/([^/\s]+/[^/\s]+)/pull/(\d+)/?\s*$", re.M)
HEREDOC = re.compile(r"<<-?\s*(['\"]?)([A-Za-z_][A-Za-z0-9_]*)\1")
ASSIGNMENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*=")
SEPARATORS = set(";&|()\n")
GH_OPTS_WITH_VALUE = ("-R", "--repo")


def strip_heredocs(command):
    """Drop heredoc bodies: text fed to a command's stdin is not a command."""
    lines = command.split("\n")
    kept, i = [], 0
    while i < len(lines):
        line = lines[i]
        kept.append(line)
        i += 1
        for match in HEREDOC.finditer(line):
            while i < len(lines) and lines[i].strip() != match.group(2):
                i += 1
            i += 1  # the terminator line
    return "\n".join(kept)


def simple_commands(command):
    """The word list of each simple command. Quoted text stays one word, even
    across lines, so a body that mentions `gh pr create` is not a command."""
    lex = shlex.shlex(strip_heredocs(command), posix=True, punctuation_chars=";&|()\n")
    lex.whitespace = " \t\r"
    lex.whitespace_split = True
    lex.commenters = ""
    words = []
    for token in lex:
        if token and set(token) <= SEPARATORS:
            if words:
                yield words
            words = []
        else:
            words.append(token)
    if words:
        yield words


def runs_pr_create(words):
    """True if this simple command is `[env] [VAR=x …] gh [-R repo] pr create …`."""
    i = 0
    while i < len(words) and (ASSIGNMENT.match(words[i]) or os.path.basename(words[i]) == "env"):
        i += 1
    if i >= len(words) or os.path.basename(words[i]) != "gh":
        return False
    i += 1
    while i < len(words) and words[i].startswith("-"):
        i += 2 if words[i] in GH_OPTS_WITH_VALUE else 1
    return words[i:i + 2] == ["pr", "create"]


def is_pr_create(command):
    try:
        return any(runs_pr_create(words) for words in simple_commands(command))
    except ValueError:  # unbalanced quotes: not a command we can read, stay silent
        return False


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
    if not is_pr_create((data.get("tool_input") or {}).get("command") or ""):
        return
    # No PR URL in the output means the create failed: nothing to review.
    response = data.get("tool_response") or ""
    stdout = (response.get("stdout") or "") if isinstance(response, dict) else str(response)
    urls = PR_URL.findall(stdout)
    if not urls:
        return
    repo, number = urls[-1]

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
