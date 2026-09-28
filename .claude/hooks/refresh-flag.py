#!/usr/bin/env python3
"""Keep the published Next steps page from going stale after a session.

The page is a snapshot that only `/mesa-status` can republish (the Artifact
tool is the model's, not a script's), so the hooks cannot publish it
themselves. Instead they remind the human to run it, at the moment it
matters. They never have the model run it: /mesa-status can edit the public
#124, and that edit happens only when a human asks for it.

  post     PostToolUse(Bash). When a command changed GitHub state - a push,
           a tag, a PR or issue created, edited, merged or closed, a label or
           release changed, a `gh api` write - raise the `refresh-needed` flag.
  stop     Stop. If the flag is up, suggest /mesa-status to the human ONCE
           (a `systemMessage`, never a block, which would make the model act),
           and lower it to `refresh-pending` so the user is not told again on
           every later turn.
  session  SessionStart. Tell the human if a refresh was left pending.

The flags sit in the git common dir (see `shared_dir()` in mesa-status.py),
so every worktree of the clone shares them. `mesa-status.py --mark-published`
clears both once the page is republished.

Every mode exits 0 on anything unexpected: a broken reminder must never
block a tool call or trap the session.
"""

import importlib.util
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# The command parsing is guard-remote.py's, so the two hooks can never
# disagree about what a command does.
_spec = importlib.util.spec_from_file_location("guard_remote", os.path.join(HERE, "guard-remote.py"))
guard = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(guard)

GH_WRITES = {
    "pr": {"create", "merge", "ready", "close", "reopen", "edit"},
    "issue": {"create", "close", "reopen", "edit", "delete", "transfer", "pin", "unpin"},
    "label": {"create", "edit", "delete", "clone"},
    "release": {"create", "edit", "delete", "upload"},
}
GRAPHQL_MUTATION = re.compile(r"^\s*mutation\b|\bmutation\s*[({]", re.M)

STOP_MESSAGE = ("mesa-status: this session changed GitHub, so the Next steps page and the "
                "#124 checklist may be stale. Run /mesa-status when you want them refreshed.")
SESSION_MESSAGE = ("mesa-status: GitHub changed in an earlier session and the Next steps page was "
                   "not republished afterwards. Run /mesa-status before relying on it.")


def flag_dir():
    # Tests point this at a scratch directory, never at a real clone's flags.
    override = os.environ.get("MESA_STATUS_FLAG_DIR")
    if override:
        os.makedirs(override, exist_ok=True)
        return override
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--path-format=absolute", "--git-common-dir"],
            capture_output=True, text=True, timeout=10,
            cwd=os.environ.get("CLAUDE_PROJECT_DIR") or None,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if out.returncode != 0 or not out.stdout.strip():
        return None
    path = os.path.join(out.stdout.strip(), "mesa-status")
    os.makedirs(path, exist_ok=True)
    return path


def changes_github(segment):
    """True if this one shell segment writes to GitHub or the remote."""
    tokens = guard.tokenise(segment)
    sub, args = guard.git_subcommand(tokens)
    if sub == "push":
        return "--dry-run" not in args and "-n" not in args
    if sub == "tag":
        # Listing tags (`git tag`, `git tag -l ...`) is a read.
        return bool(args) and not any(a in ("-l", "--list", "-n", "--contains") for a in args)
    # Strip `env X=Y` and a path to gh, as git_subcommand does for git.
    i = 0
    while i < len(tokens) and "=" in tokens[i] and not tokens[i].startswith("-"):
        i += 1
    if i >= len(tokens) or os.path.basename(tokens[i]) != "gh":
        return False
    rest = tokens[i + 1:]
    if len(rest) >= 2 and rest[1] in GH_WRITES.get(rest[0], ()):
        return True
    if rest and rest[0] == "api":
        if guard.gh_api_endpoint(tokens) == "graphql":
            text = guard.payload_text(tokens)
            return text is None or bool(GRAPHQL_MUTATION.search(text))
        return guard.gh_api_method(tokens) != "GET"
    return False


def mode_post(payload, flags):
    if payload.get("tool_name") != "Bash":
        return 0
    command = (payload.get("tool_input") or {}).get("command") or ""
    if any(changes_github(seg) for seg in guard.split_segments(command)):
        with open(os.path.join(flags, "refresh-needed"), "w", encoding="utf-8") as handle:
            handle.write(command.strip()[:300] + "\n")
    return 0


def mode_stop(payload, flags):
    needed = os.path.join(flags, "refresh-needed")
    # stop_hook_active: another hook has already held this stop; stay quiet.
    if payload.get("stop_hook_active") or not os.path.exists(needed):
        return 0
    os.replace(needed, os.path.join(flags, "refresh-pending"))
    print(json.dumps({"systemMessage": STOP_MESSAGE}))
    return 0


def mode_session(payload, flags):
    if any(os.path.exists(os.path.join(flags, n)) for n in ("refresh-needed", "refresh-pending")):
        print(json.dumps({"systemMessage": SESSION_MESSAGE}))
    return 0


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        payload = {}
    handler = {"post": mode_post, "stop": mode_stop, "session": mode_session}.get(mode)
    if handler is None:
        return 0
    try:
        flags = flag_dir()
        if flags is None:
            return 0
        return handler(payload, flags)
    except Exception:  # a reminder that crashes must never block the session
        return 0


if __name__ == "__main__":
    sys.exit(main())
