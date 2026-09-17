#!/usr/bin/env python3
"""PreToolUse guard for mesa's agent workflow contract.

AGENTS.md lets an agent push a feature branch and open a DRAFT pull request.
Everything past that point belongs to a human. This hook enforces that
mechanically so the contract does not depend on the model remembering it.

Blocked:
  * any push targeting main or dev
  * force pushes (--force, --force-with-lease, +refs)
  * gh pr merge / ready / review --approve, and git merge into a protected branch
  * deleting branches or tags, locally or on the remote
  * history rewrites pushed to a remote

Exit codes: 0 = allow, 2 = block (stderr is shown to the agent).

Read-only inspection is never blocked. When this hook and AGENTS.md disagree,
that is a bug: fix both.
"""

import json
import re
import sys

PROTECTED = r"(?:main|dev|master)"

# (compiled pattern, message) — the message tells the agent what to do instead.
RULES = [
    (
        re.compile(
            r"\bgit\s+push\b(?=[^\n;&|]*?(?:"
            r"(?::(?:refs/heads/)?" + PROTECTED + r"\b)"      # HEAD:main
            r"|(?:\s" + PROTECTED + r"(?:\s|$))"              # push origin main
            r"))",
            re.IGNORECASE,
        ),
        "Pushing to a protected branch (main/dev) is not allowed. "
        "Push your feature branch instead and let a human merge.",
    ),
    (
        re.compile(
            r"\bgit\s+push\b[^\n;&|]*?(?:--force\b|--force-with-lease|\s-f\b|\s\+)",
            re.IGNORECASE,
        ),
        "Force pushing is not allowed — it can destroy published history. "
        "If a branch genuinely needs rewriting, ask the human to do it.",
    ),
    (
        re.compile(r"\bgh\s+pr\s+merge\b", re.IGNORECASE),
        "Merging pull requests is the human's decision, not the agent's.",
    ),
    (
        re.compile(r"\bgh\s+pr\s+(?:ready|edit\b[^\n;&|]*--ready)\b", re.IGNORECASE),
        "PRs opened by an agent stay in DRAFT. Only a human marks one ready for review.",
    ),
    (
        re.compile(r"\bgh\s+pr\s+review\b[^\n;&|]*--approve", re.IGNORECASE),
        "An agent does not approve pull requests on this repo.",
    ),
    (
        re.compile(r"\bgit\s+push\b[^\n;&|]*(?:--delete\b|--mirror\b|\s:\w)", re.IGNORECASE),
        "Deleting remote refs is not allowed.",
    ),
    (
        re.compile(r"\bgit\s+branch\b[^\n;&|]*\s-(?:D|-delete)\b", re.IGNORECASE),
        "Deleting branches is not allowed — they are the human's audit trail.",
    ),
    (
        re.compile(r"\bgit\s+tag\b[^\n;&|]*\s-(?:d|-delete)\b", re.IGNORECASE),
        "Deleting tags is not allowed. Tags are applied by the human after merge.",
    ),
    (
        re.compile(r"\bgit\s+(?:checkout|switch)\b[^\n;&|]*\s" + PROTECTED
                   + r"\b[^\n;&|]*&&[^\n;&|]*\bgit\s+(?:commit|merge)\b", re.IGNORECASE),
        "Committing or merging directly on a protected branch is not allowed. "
        "Branch off main and open a pull request.",
    ),
]


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        # A guard that crashes must not silently stop blocking.
        print("guard-remote: could not parse hook input; blocking to be safe.",
              file=sys.stderr)
        return 2

    if payload.get("tool_name") != "Bash":
        return 0

    command = payload.get("tool_input", {}).get("command", "")
    if not command:
        return 0

    # --dry-run inspects without mutating anything.
    if re.search(r"--dry-run\b", command):
        return 0

    for pattern, message in RULES:
        if pattern.search(command):
            print(
                f"Blocked by mesa's workflow contract (AGENTS.md):\n  {message}\n\n"
                f"Command: {command.strip()[:300]}\n\n"
                "Do not work around this — report it to the human instead.",
                file=sys.stderr,
            )
            return 2

    return 0


if __name__ == "__main__":
    sys.exit(main())
