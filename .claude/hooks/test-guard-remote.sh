#!/usr/bin/env bash
# Regression tests for guard-remote.py.
#
# A guard that silently stops firing is worse than no guard, so run this after
# editing the hook:  bash .claude/hooks/test-guard-remote.sh
#
# Exit 0 = all cases behave as expected.

set -uo pipefail
hook="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/guard-remote.py"
fails=0

check() { # check <expected-exit> <command>
    local expected="$1" cmd="$2" actual
    actual="$(printf '%s' "$cmd" \
        | python3 -c 'import json,sys; print(json.dumps({"tool_name":"Bash","tool_input":{"command":sys.stdin.read()}}))' \
        | python3 "$hook" >/dev/null 2>&1; echo $?)"
    if [ "$actual" != "$expected" ]; then
        printf 'FAIL (expected %s, got %s): %s\n' "$expected" "$actual" "$cmd"
        fails=$((fails + 1))
    fi
}

# --- must be blocked -------------------------------------------------
for cmd in \
    'git push origin HEAD:main' \
    'git push origin main' \
    'git push -u origin main' \
    'git push origin dev' \
    'git push origin HEAD:refs/heads/main' \
    'git push --force origin HEAD' \
    'git push --force-with-lease origin my-branch' \
    'cd /somewhere && git push origin main' \
    'gh pr merge 102' \
    'gh pr merge 102 --squash' \
    'gh pr ready 102' \
    'gh pr review 102 --approve' \
    'git push origin --delete old-branch' \
    'git branch -D old-branch' \
    'git tag -d v0.99.6'
do check 2 "$cmd"; done

# --- must be allowed -------------------------------------------------
for cmd in \
    'git push origin chore/agent-setup' \
    'git push -u origin chore/agent-setup' \
    'git push origin chore/main-cleanup' \
    'git push --dry-run origin main' \
    'git status' \
    'git log --oneline -5' \
    'git diff main...HEAD' \
    'git fetch origin' \
    'git checkout main' \
    'gh pr create --draft --title "x" --body "y"' \
    'gh pr view 102' \
    'Rscript -e "devtools::test()"'
do check 0 "$cmd"; done

# --- failure modes ---------------------------------------------------
# Malformed input must fail closed (block), not fall open.
if [ "$(printf 'not json' | python3 "$hook" >/dev/null 2>&1; echo $?)" != "2" ]; then
    echo "FAIL: malformed input did not fail closed"
    fails=$((fails + 1))
fi
# Non-Bash tools are none of this hook's business.
if [ "$(echo '{"tool_name":"Read","tool_input":{"file_path":"/x"}}' | python3 "$hook" >/dev/null 2>&1; echo $?)" != "0" ]; then
    echo "FAIL: non-Bash tool was not passed through"
    fails=$((fails + 1))
fi

if [ "$fails" -eq 0 ]; then
    echo "guard-remote: all cases pass"
else
    echo "guard-remote: $fails failing case(s)"
fi
exit $((fails > 0))
