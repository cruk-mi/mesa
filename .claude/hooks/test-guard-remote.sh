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
    'gh pr edit 104 --ready' \
    'gh pr review 102 --approve' \
    'git push origin --delete old-branch' \
    'git branch -D old-branch' \
    'git tag -d v0.99.6' \
    'git -C /repo push origin main' \
    'git -c user.name=x push origin main' \
    'git push origin refs/heads/main' \
    'git push --dry-run origin main && git push origin main' \
    'git push origin main; echo done' \
    'git push --mirror origin' \
    'git push --all origin' \
    'git push origin :main' \
    'git push origin +main' \
    'gh api -X PUT repos/cruk-mi/mesa/pulls/106/merge' \
    'gh api repos/cruk-mi/mesa/pulls/106/merge -X PUT' \
    'gh api --method PATCH repos/cruk-mi/mesa/pulls/106 -f draft=false' \
    'gh api -X POST repos/cruk-mi/mesa/pulls/106/reviews -f event=APPROVE' \
    'gh api -X DELETE repos/cruk-mi/mesa/git/refs/heads/chore/status-tracking' \
    "gh api -X PUT 'repos/cruk-mi/mesa/pulls/106/merge?'" \
    "gh api -X PUT 'repos/cruk-mi/mesa/pulls/106/merge?merge_method=squash'" \
    "gh api -X PUT 'repos/cruk-mi/mesa/pulls/106/merge#x'" \
    "gh api -X PATCH 'repos/cruk-mi/mesa/pulls/106?' -f draft=false" \
    "gh api graphql -f query='mutation { mergePullRequest(input:{pullRequestId:\"X\"}) { clientMutationId } }'" \
    "gh api graphql -f query='mutation { enablePullRequestAutoMerge(input:{pullRequestId:\"X\"}) { clientMutationId } }'" \
    "gh api graphql -f query='mutation { markPullRequestReadyForReview(input:{pullRequestId:\"X\"}) { clientMutationId } }'" \
    "gh api graphql -f query='mutation { addPullRequestReview(input:{pullRequestId:\"X\", event: APPROVE}) { clientMutationId } }'" \
    "gh api graphql -f query='mutation { deleteRef(input:{refId:\"X\"}) { clientMutationId } }'" \
    "gh api graphql --raw-field query='mutation { mergePullRequest(input:{pullRequestId:\"X\"}) { clientMutationId } }'" \
    'gh api graphql --input -' \
    'gh api -X POST repos/cruk-mi/mesa/pulls/106/reviews/1/events -f event=APPROVE' \
    'gh api -X POST repos/cruk-mi/mesa/pulls/106/reviews --input -' \
    'gh api -X POST repos/cruk-mi/mesa/pulls/106/reviews/1/events --input -'
do check 2 "$cmd"; done

# A mutation hidden in a query file must be read and caught, not waved through.
qfile="$(mktemp)"
printf 'mutation { mergePullRequest(input:{pullRequestId:"X"}) { clientMutationId } }' >"$qfile"
check 2 "gh api graphql -F query=@$qfile"
check 2 "gh api graphql --input $qfile"
rm -f "$qfile"

# Same for a REST review payload: APPROVE in a file is caught, COMMENT is not.
rfile="$(mktemp)"
printf '{"event":"APPROVE"}' >"$rfile"
check 2 "gh api -X POST repos/cruk-mi/mesa/pulls/106/reviews --input $rfile"
check 2 "gh api -X POST repos/cruk-mi/mesa/pulls/106/reviews/1/events --input $rfile"
printf '{"event":"COMMENT","body":"x"}' >"$rfile"
check 0 "gh api -X POST repos/cruk-mi/mesa/pulls/106/reviews --input $rfile"
check 0 "gh api -X POST repos/cruk-mi/mesa/pulls/106/reviews/1/events -f event=COMMENT"
rm -f "$rfile"

# --- must be allowed -------------------------------------------------
for cmd in \
    'git push origin chore/agent-setup' \
    'git push -u origin chore/agent-setup' \
    'git push origin chore/main-cleanup' \
    'git push origin HEAD:refs/heads/fix/81' \
    'git push --dry-run origin main' \
    'git status' \
    'git log --oneline -5' \
    'git diff main...HEAD' \
    'git fetch origin' \
    'git checkout main' \
    'git branch --show-current' \
    'gh pr create --draft --title "x" --body "y"' \
    'gh pr view 102' \
    'gh pr list --state open' \
    'gh api repos/cruk-mi/mesa/pulls/106/comments' \
    'gh api repos/cruk-mi/mesa/pulls/106 --jq .state' \
    'gh api -X POST repos/cruk-mi/mesa/pulls/106/comments/4063116787/replies -f body=fixed' \
    'gh api -X POST repos/cruk-mi/mesa/pulls/106/comments/1/replies -f body="see /pulls/106/merge"' \
    "gh api graphql -f query='mutation { resolveReviewThread(input:{threadId:\"X\"}) { thread { isResolved } } }'" \
    "gh api graphql -f query='query { repository(owner:\"cruk-mi\", name:\"mesa\") { pullRequest(number:106) { state } } }'" \
    "gh api graphql -f query='mutation { addPullRequestReview(input:{pullRequestId:\"X\", event: COMMENT, body:\"x\"}) { clientMutationId } }'" \
    "gh api 'repos/cruk-mi/mesa/pulls/106/comments?per_page=100'" \
    "gh api repos/cruk-mi/mesa/pulls/106/reviews --jq '.[] | select(.state == \"APPROVED\")'" \
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

# --- stateful rule ---------------------------------------------------
# `git commit` / `git merge` are blocked by inspecting the CURRENT branch, not
# the command string, so the rule survives being split across separate tool
# calls (checkout in one, commit in the next). Exercised in a throwaway repo
# because it depends on real repo state.
tmp="$(mktemp -d)"
(
    cd "$tmp" || exit 1
    git init -q -b main . && git commit -q --allow-empty -m init
) >/dev/null 2>&1
state_check() { # state_check <expected> <cmd>
    local expected="$1" cmd="$2" actual
    actual="$(cd "$tmp" && printf '%s' "$cmd" \
        | python3 -c 'import json,sys; print(json.dumps({"tool_name":"Bash","tool_input":{"command":sys.stdin.read()}}))' \
        | python3 "$hook" >/dev/null 2>&1; echo $?)"
    if [ "$actual" != "$expected" ]; then
        printf 'FAIL (stateful, expected %s, got %s): %s\n' "$expected" "$actual" "$cmd"
        fails=$((fails + 1))
    fi
}
state_check 2 'git commit -m "x"'          # on main -> blocked
state_check 2 'git merge some-branch'      # on main -> blocked
state_check 2 'git push'                   # bare push on main -> blocked
(cd "$tmp" && git checkout -q -b feat/x) >/dev/null 2>&1
state_check 0 'git commit -m "x"'          # on a feature branch -> allowed
state_check 0 'git push'                   # bare push on a feature branch -> allowed
rm -rf "$tmp"

if [ "$fails" -eq 0 ]; then
    echo "guard-remote: all cases pass"
else
    echo "guard-remote: $fails failing case(s)"
fi
exit $((fails > 0))
