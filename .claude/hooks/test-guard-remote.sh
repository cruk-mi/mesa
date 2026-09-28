#!/usr/bin/env bash
# Regression tests for guard-remote.py.
#
# A guard that silently stops firing is worse than no guard, so run this after
# editing the hook:  bash .claude/hooks/test-guard-remote.sh
#
# Exit 0 = all cases behave as expected.

set -uo pipefail
hook="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/guard-remote.py"
# The interpreter under test: PY=/usr/bin/python3 checks the macOS system 3.9.
PY="${PY:-python3}"
fails=0

# The hook asks `gh` about open PRs and tags. A fake `gh` first on PATH answers
# from a fixed set, so the refusal paths run offline and never depend on which
# PRs happen to be open. #110 is stacked on #108's branch `stacked-base`.
stub_dir="$(mktemp -d)"
cat >"$stub_dir/gh" <<'STUB'
#!/usr/bin/env bash
[ "${GH_STUB_FAIL:-}" = 1 ] && exit 1
case "$1 $2" in
"pr list")
    echo '[{"number":110,"baseRefName":"stacked-base","headRefName":"feat/on-top"},
           {"number":108,"baseRefName":"main","headRefName":"stacked-base"},
           {"number":114,"baseRefName":"main","headRefName":"chore/open-head"},
           {"number":102,"baseRefName":"main","headRefName":"feat/clean"}]' ;;
"pr view")
    case "$3" in
    108) echo '{"number":108,"headRefName":"stacked-base"}' ;;
    110) echo '{"number":110,"headRefName":"feat/on-top"}' ;;
    102|feat/clean|-*) echo '{"number":102,"headRefName":"feat/clean"}' ;;
    *) exit 1 ;;
    esac ;;
"api repos/{owner}/{repo}/git/matching-refs/tags/"*)
    case "${2##*/tags/}" in
    v*) echo '[{"ref":"refs/tags/v0.99.6"}]' ;;
    *) echo '[]' ;;
    esac ;;
*) exit 1 ;;
esac
STUB
chmod +x "$stub_dir/gh"
export PATH="$stub_dir:$PATH"
trap 'rm -rf "$stub_dir"' EXIT

# Hook input as Claude Code sends it; MODE sets permission_mode.
payload() {
    MODE="${MODE:-default}" "$PY" -c 'import json,os,sys; print(json.dumps({"tool_name":"Bash","permission_mode":os.environ["MODE"],"tool_input":{"command":sys.stdin.read()}}))'
}

check() { # check <expected-exit> <command>
    local expected="$1" cmd="$2" actual
    actual="$(printf '%s' "$cmd" | payload | "$PY" "$hook" >/dev/null 2>&1; echo $?)"
    if [ "$actual" != "$expected" ]; then
        printf 'FAIL (expected %s, got %s): %s\n' "$expected" "$actual" "$cmd"
        fails=$((fails + 1))
    fi
}

ask_check() { # ask_check <command> [text]: allowed only via a permission prompt
    local cmd="$1" text="${2:-}" out code
    out="$(printf '%s' "$cmd" | payload | "$PY" "$hook" 2>/dev/null)"
    code=$?
    if [ "$code" != 0 ] || ! printf '%s' "$out" | grep -q '"permissionDecision": "ask"'; then
        printf 'FAIL (expected ask, got exit %s): %s\n' "$code" "$cmd"
        fails=$((fails + 1))
    elif [ -n "$text" ] && ! printf '%s' "$out" | grep -qF -- "$text"; then
        printf 'FAIL (ask does not say "%s"): %s\n  %s\n' "$text" "$cmd" "$out"
        fails=$((fails + 1))
    fi
}

# --- must ask the human (merging needs explicit approval) ------------
for cmd in \
    'gh pr merge 102' \
    'gh pr merge 102 --squash' \
    'gh api -X PUT repos/cruk-mi/mesa/pulls/106/merge' \
    'gh api repos/cruk-mi/mesa/pulls/106/merge -X PUT' \
    "gh api -X PUT 'repos/cruk-mi/mesa/pulls/106/merge?'" \
    "gh api -X PUT 'repos/cruk-mi/mesa/pulls/106/merge?merge_method=squash'" \
    "gh api -X PUT 'repos/cruk-mi/mesa/pulls/106/merge#x'" \
    "gh api graphql -f query='mutation { mergePullRequest(input:{pullRequestId:\"X\"}) { clientMutationId } }'" \
    "gh api graphql -f query='mutation { enablePullRequestAutoMerge(input:{pullRequestId:\"X\"}) { clientMutationId } }'" \
    "gh api graphql --raw-field query='mutation { mergePullRequest(input:{pullRequestId:\"X\"}) { clientMutationId } }'" \
    'gh pr view 102 && gh pr merge 102 --squash' \
    'gh pr merge 102 --squash --delete-branch' \
    'gh pr merge 102 -d' \
    'gh pr merge -R cruk-mi/mesa 102 --squash --delete-branch' \
    'git push origin --delete old-branch' \
    'git push -d origin old-branch other-branch' \
    'git push origin :old-branch' \
    'git branch -D old-branch' \
    'git branch -d old-branch' \
    'gh api -X DELETE repos/cruk-mi/mesa/git/refs/heads/chore/status-tracking' \
    "gh api graphql -f query='mutation { deleteRef(input:{refId:\"X\"}) { clientMutationId } }'"
do ask_check "$cmd"; done

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
    'gh pr ready 102' \
    'git push origin --delete main' \
    'git push -d origin dev' \
    'git branch -D main' \
    'gh api -X DELETE repos/cruk-mi/mesa/git/refs/heads/main' \
    'gh api -X DELETE repos/cruk-mi/mesa/git/refs/tags/v0.99.6' \
    'git push origin --delete refs/tags/v0.99.6' \
    'gh pr merge 102 -d && git push origin main' \
    'gh pr merge 102 --squash && git push origin main' \
    'gh pr edit 104 --ready' \
    'gh pr review 102 --approve' \
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
    'gh api --method PATCH repos/cruk-mi/mesa/pulls/106 -f draft=false' \
    'gh api -X POST repos/cruk-mi/mesa/pulls/106/reviews -f event=APPROVE' \
    "gh api -X PATCH 'repos/cruk-mi/mesa/pulls/106?' -f draft=false" \
    "gh api graphql -f query='mutation { markPullRequestReadyForReview(input:{pullRequestId:\"X\"}) { clientMutationId } }'" \
    "gh api graphql -f query='mutation { addPullRequestReview(input:{pullRequestId:\"X\", event: APPROVE}) { clientMutationId } }'" \
    'gh api graphql --input -' \
    'gh api -X POST repos/cruk-mi/mesa/pulls/106/reviews/1/events -f event=APPROVE' \
    'gh api -X POST repos/cruk-mi/mesa/pulls/106/reviews --input -' \
    'gh api -X POST repos/cruk-mi/mesa/pulls/106/reviews/1/events --input -'
do check 2 "$cmd"; done

# A mutation hidden in a query file must be read and caught, not waved through.
qfile="$(mktemp)"
printf 'mutation { mergePullRequest(input:{pullRequestId:"X"}) { clientMutationId } }' >"$qfile"
ask_check "gh api graphql -F query=@$qfile"
ask_check "gh api graphql --input $qfile"
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

# --- stacked PRs (answered by the fake gh) --------------------------
# Deleting a branch an open PR is based on closes that PR: refused.
check 2 'git push origin --delete stacked-base'
check 2 'git push origin :stacked-base'
check 2 'gh api -X DELETE repos/cruk-mi/mesa/git/refs/heads/stacked-base'
check 2 'gh pr merge 108 --delete-branch'
check 2 'gh pr merge 108 -d'
# Nothing is built on feat/clean: asked, and the prompt says it was checked.
ask_check 'git push origin --delete feat/clean' 'no open PR'
ask_check 'gh pr merge 102 --squash --delete-branch' 'no open PR'

# gh (pflag) also spells --delete-branch as --delete-branch=true, and bundles
# -d with other shorthand. Each must reach the same check, and the prompt must
# name the PR and the branch it deletes.
check 2 'gh pr merge 108 --delete-branch=true'
check 2 'gh pr merge 108 -sd'
check 2 'gh pr merge 108 -ds'
check 2 'gh pr close 108 -d'
check 2 'gh pr close 108 --delete-branch'
check 2 'gh pr merge 999 -d'                # PR not found: cannot tell what goes
ask_check 'gh pr merge 102 --delete-branch=true' "PR #102 and delete its branch 'feat/clean'"
ask_check 'gh pr merge 102 -sd' "merge PR #102 and delete its branch 'feat/clean'"
ask_check 'gh pr merge 102 -ds' "delete its branch 'feat/clean'"
ask_check 'gh pr merge 102 -d' "delete its branch 'feat/clean'"
ask_check 'gh pr close 102 -d' "close PR #102 and delete its branch 'feat/clean'"
check 0 'gh pr close 102'

# --- failure modes ---------------------------------------------------
# Malformed input must fail closed (block), not fall open.
if [ "$(printf 'not json' | "$PY" "$hook" >/dev/null 2>&1; echo $?)" != "2" ]; then
    echo "FAIL: malformed input did not fail closed"
    fails=$((fails + 1))
fi
# Non-Bash tools are none of this hook's business.
if [ "$(echo '{"tool_name":"Read","tool_input":{"file_path":"/x"}}' | "$PY" "$hook" >/dev/null 2>&1; echo $?)" != "0" ]; then
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
    actual="$(cd "$tmp" && printf '%s' "$cmd" | payload | "$PY" "$hook" >/dev/null 2>&1; echo $?)"
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
