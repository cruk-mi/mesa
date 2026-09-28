#!/usr/bin/env bash
# Regression tests for refresh-flag.py.
#
# A reminder that silently stops firing leaves the Next steps page stale with
# nobody noticing, so run this after editing the hook:
#   bash .claude/hooks/test-refresh-flag.sh
#
# Exit 0 = all cases behave as expected. PYTHON=/usr/bin/python3 runs them all
# under another interpreter (macOS ships 3.9).

set -uo pipefail
hook="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/refresh-flag.py"
MESA_STATUS_FLAG_DIR="$(mktemp -d)"
export MESA_STATUS_FLAG_DIR
flags="$MESA_STATUS_FLAG_DIR"
PY="${PYTHON:-python3}"
fails=0

post() { # post <command>: run the PostToolUse mode on one Bash command
    printf '%s' "$1" \
        | "$PY" -c 'import json,sys; print(json.dumps({"tool_name":"Bash","tool_input":{"command":sys.stdin.read()}}))' \
        | "$PY" "$hook" post >/dev/null 2>&1
}
fail() { printf 'FAIL: %s\n' "$1"; fails=$((fails + 1)); }

# --- must raise the flag ---------------------------------------------
for cmd in \
    'git push -u origin chore/x' \
    'git push' \
    'cd /repo && git push origin feat/y' \
    'git -C /repo push origin feat/y' \
    'git tag -a v0.99.8 -m "mesa 0.99.8"' \
    'gh pr create --draft --title "x" --body "y"' \
    'gh pr merge 135 --squash' \
    'gh pr edit 114 --add-label tests' \
    'gh pr close 110' \
    'gh issue create --title x --body y' \
    'gh issue close 34 --comment "decided"' \
    'gh issue edit 124 --body-file /tmp/b.md' \
    'gh label create "type: chore"' \
    'gh release create v0.99.8' \
    'gh api -X POST repos/cruk-mi/mesa/issues/34/comments -f body=x' \
    'gh api repos/cruk-mi/mesa/issues/34/labels -f "labels[]=parked"' \
    "gh api graphql -f query='mutation { closeIssue(input:{issueId:\"X\"}) { clientMutationId } }'"
do
    rm -f "$flags"/refresh-*
    post "$cmd"
    [ -f "$flags/refresh-needed" ] || fail "no flag for: $cmd"
done

# --- must NOT raise the flag -----------------------------------------
for cmd in \
    'git status' \
    'git push --dry-run origin feat/x' \
    'git tag' \
    'git tag -l "v0.99*"' \
    'git log --oneline -5' \
    'gh pr view 114' \
    'gh pr list --state open' \
    'gh pr checks 114' \
    'gh issue view 124 --json body' \
    'gh api repos/cruk-mi/mesa/pulls/106/comments' \
    "gh api graphql -f query='query { repository(owner:\"cruk-mi\", name:\"mesa\") { issue(number:124) { body } } }'" \
    'python3 .claude/scripts/mesa-status.py --html --sync-issue' \
    'echo "gh pr merge 1"'
do
    rm -f "$flags"/refresh-*
    post "$cmd"
    [ -f "$flags/refresh-needed" ] && fail "flag raised for read-only: $cmd"
done

# --- Stop blocks once, then lets go ------------------------------------
rm -f "$flags"/refresh-*
out="$(echo '{}' | "$PY" "$hook" stop)"
[ -z "$out" ] || fail "stop blocked with no flag raised"

post 'gh pr merge 1'
out="$(echo '{"stop_hook_active": false}' | "$PY" "$hook" stop)"
echo "$out" | grep -q '"decision": "block"' || fail "stop did not block with the flag raised"
[ -f "$flags/refresh-pending" ] || fail "stop did not lower the flag to pending"
out="$(echo '{"stop_hook_active": false}' | "$PY" "$hook" stop)"
[ -z "$out" ] || fail "stop blocked a second time for the same change"

post 'gh pr merge 2'
out="$(echo '{"stop_hook_active": true}' | "$PY" "$hook" stop)"
[ -z "$out" ] || fail "stop blocked while stop_hook_active (would loop)"

# --- SessionStart reports a pending refresh ----------------------------
out="$(echo '{}' | "$PY" "$hook" session)"
echo "$out" | grep -q 'not republished' || fail "session start did not report the pending refresh"
rm -f "$flags"/refresh-*
out="$(echo '{}' | "$PY" "$hook" session)"
[ -z "$out" ] || fail "session start reported a refresh with no flag"

# --- failure modes -------------------------------------------------------
[ "$(printf 'not json' | "$PY" "$hook" post >/dev/null 2>&1; echo $?)" = "0" ] \
    || fail "malformed input did not exit 0"
[ "$(echo '{}' | "$PY" "$hook" nonsense >/dev/null 2>&1; echo $?)" = "0" ] \
    || fail "unknown mode did not exit 0"

# --- mesa-status.py: the roadmap and the #124 sync -----------------------
# In-process, with run() and every GitHub fetch stubbed, so no case can reach
# the real gh or edit the real #124. A traceback counts as a failure.
status_script="$(dirname "$hook")/../scripts/mesa-status.py"
"$PY" "$status_script" --help >/dev/null 2>&1 \
    || fail "mesa-status.py --help fails under $("$PY" --version 2>&1)"
roadmap_out="$("$PY" - "$status_script" 2>&1 <<'PYEOF'
import importlib.util
import os
import sys
import tempfile

spec = importlib.util.spec_from_file_location("mesa_status", sys.argv[1])
ms = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ms)

CALLS = []


def fake_run(argv, timeout=30):
    CALLS.append(list(argv))
    return ""


ms.run = fake_run
ms.ROADMAP_MD = os.path.join(tempfile.mkdtemp(), "roadmap-issue.md")
try:
    import tomllib  # noqa: F401
    HAVE_TOML = True
except ImportError:
    HAVE_TOML = False

TOML = """title = "Roadmap"
[[wave]]
id = "w1"
name = "Wave 1"
[[item]]
wave = "w1"
title = "Fix A"
refs = [10]
[[item]]
wave = "w1"
title = "Fix B"
refs = [11]
after = ["10"]
check = ["branch-gone:old"]
"""
STATES = {
    10: {"__typename": "PullRequest", "state": "MERGED", "mergedAt": "2026-09-01T00:00:00Z",
         "url": "u10", "title": "A"},
    11: {"__typename": "Issue", "state": "OPEN", "url": "u11", "title": "B",
         "closedByPullRequestsReferences": {"nodes": []}},
}


def fail(message):
    print("FAIL: " + message)


def needs_toml(name):
    if not HAVE_TOML:
        print(f"SKIP: {name} (needs Python >= 3.11)")
    return HAVE_TOML


def data_block(toml):
    return f"{ms.DATA_START}\n<details><summary>data</summary>\n\n```toml\n{toml}```\n\n</details>\n{ms.DATA_END}"


def live_body(toml=TOML, gen=True, data=True):
    """#124 as the first sync leaves it: generated half, then the data block."""
    parts = [f"{ms.GEN_START}\n# Roadmap\nold checklist\n{ms.GEN_END}" if gen else "old checklist",
             "Hand-written note.",
             data_block(toml) if data else f"```toml\n{toml}```"]
    return "\n\n".join(parts) + "\n"


def roadmap_for(body, states=STATES, refs=({"main"}, set()), override=None):
    """collect_roadmap() against a stubbed GitHub. Clears DEGRADED and CALLS."""
    ms.DEGRADED.clear()
    CALLS.clear()
    ms.fetch_issue_body = lambda: body
    ms.fetch_ref_states = lambda numbers: None if states is None else dict(states)
    ms.remote_refs = lambda: refs
    return ms.collect_roadmap(override)


# tomllib is 3.11+: without it only the roadmap degrades.
saved = sys.modules.get("tomllib")
sys.modules["tomllib"] = None
ms.DEGRADED.clear()
try:
    if ms.parse_roadmap(TOML) is not None:
        fail("parse_roadmap without tomllib returned data")
    elif not any("3.11" in d for d in ms.DEGRADED):
        fail(f"no 'needs Python >= 3.11' reason without tomllib: {ms.DEGRADED}")
finally:
    if saved is None:
        sys.modules.pop("tomllib", None)
    else:
        sys.modules["tomllib"] = saved

# CRLF: github.com edits come back with \r\n line endings.
if needs_toml("CRLF body"):
    crlf = live_body().replace("\n", "\r\n")
    roadmap, body = roadmap_for(crlf)
    if roadmap is None:
        fail(f"CRLF body: roadmap unknown: {ms.DEGRADED}")
    else:
        current = ms.compose_issue_body(live_body(), ms.render_roadmap_issue(roadmap), None)
        verdict = ms.sync_issue(roadmap, current.replace("\n", "\r\n"), None, True, None)
        if "already current" not in verdict:
            fail(f"CRLF body that is current: {verdict!r}")
PYEOF
)"
while IFS= read -r line; do
    case "$line" in
        "") ;;
        SKIP:*) echo "$line" ;;
        FAIL:*) fail "${line#FAIL: }" ;;
        *) fail "mesa-status.py: $line" ;;
    esac
done <<<"$roadmap_out"

rm -rf "$flags"
if [ "$fails" -eq 0 ]; then
    echo "refresh-flag: all cases pass"
else
    echo "refresh-flag: $fails failing case(s)"
fi
exit $((fails > 0))
