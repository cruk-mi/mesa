#!/usr/bin/env bash
# Regression tests for refresh-flag.py.
#
# A reminder that silently stops firing leaves the Next steps page stale with
# nobody noticing, so run this after editing the hook:
#   bash .claude/hooks/test-refresh-flag.sh
#
# Exit 0 = all cases behave as expected.

set -uo pipefail
hook="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/refresh-flag.py"
MESA_STATUS_FLAG_DIR="$(mktemp -d)"
export MESA_STATUS_FLAG_DIR
flags="$MESA_STATUS_FLAG_DIR"
fails=0

post() { # post <command>: run the PostToolUse mode on one Bash command
    printf '%s' "$1" \
        | python3 -c 'import json,sys; print(json.dumps({"tool_name":"Bash","tool_input":{"command":sys.stdin.read()}}))' \
        | python3 "$hook" post >/dev/null 2>&1
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
out="$(echo '{}' | python3 "$hook" stop)"
[ -z "$out" ] || fail "stop blocked with no flag raised"

post 'gh pr merge 1'
out="$(echo '{"stop_hook_active": false}' | python3 "$hook" stop)"
echo "$out" | grep -q '"decision": "block"' || fail "stop did not block with the flag raised"
[ -f "$flags/refresh-pending" ] || fail "stop did not lower the flag to pending"
out="$(echo '{"stop_hook_active": false}' | python3 "$hook" stop)"
[ -z "$out" ] || fail "stop blocked a second time for the same change"

post 'gh pr merge 2'
out="$(echo '{"stop_hook_active": true}' | python3 "$hook" stop)"
[ -z "$out" ] || fail "stop blocked while stop_hook_active (would loop)"

# --- SessionStart reports a pending refresh ----------------------------
out="$(echo '{}' | python3 "$hook" session)"
echo "$out" | grep -q 'not republished' || fail "session start did not report the pending refresh"
rm -f "$flags"/refresh-*
out="$(echo '{}' | python3 "$hook" session)"
[ -z "$out" ] || fail "session start reported a refresh with no flag"

# --- failure modes -------------------------------------------------------
[ "$(printf 'not json' | python3 "$hook" post >/dev/null 2>&1; echo $?)" = "0" ] \
    || fail "malformed input did not exit 0"
[ "$(echo '{}' | python3 "$hook" nonsense >/dev/null 2>&1; echo $?)" = "0" ] \
    || fail "unknown mode did not exit 0"

rm -rf "$flags"
if [ "$fails" -eq 0 ]; then
    echo "refresh-flag: all cases pass"
else
    echo "refresh-flag: $fails failing case(s)"
fi
exit $((fails > 0))
