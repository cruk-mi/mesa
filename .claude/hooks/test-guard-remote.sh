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
# GH_STUB_FAIL=1 fails every call; =list fails only the open-PR listing.
case "${GH_STUB_FAIL:-}" in
1) exit 1 ;;
list) [ "$1 $2" = "pr list" ] && exit 1 ;;
esac
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

# Hook input as Claude Code sends it; MODE sets permission_mode (- leaves it out).
payload() {
    MODE="${MODE:-default}" "$PY" -c '
import json, os, sys
p = {"tool_name": "Bash", "tool_input": {"command": sys.stdin.read()}}
if os.environ["MODE"] != "-":
    p["permission_mode"] = os.environ["MODE"]
print(json.dumps(p))'
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
    'gh api -X DELETE repos/cruk-mi/mesa/git/refs/heads/chore/status-tracking'
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
    'gh api -X POST repos/cruk-mi/mesa/pulls/106/reviews/1/events --input -' \
    "gh api graphql -f query='mutation { deleteRef(input:{refId:\"X\"}) { clientMutationId } }'" \
    "gh api graphql -f query='mutation { updateRefs(input:{repositoryId:\"X\", refUpdates:[]}) { clientMutationId } }'" \
    "gh api graphql -f query='mutation { mergeBranch(input:{repositoryId:\"X\", base:\"main\", head:\"x\"}) { clientMutationId } }'" \
    'gh api -X POST repos/cruk-mi/mesa/merges -f base=main -f head=feat/x' \
    'gh api repos/cruk-mi/mesa/merges -f base=feat/y -f head=feat/x' \
    'gh api -X PATCH repos/cruk-mi/mesa/git/refs/heads/main -f sha=abc' \
    'gh api -X PATCH repos/cruk-mi/mesa/git/refs/heads/feat/x -f sha=abc -F force=true' \
    'gh api -X PATCH repos/cruk-mi/mesa/git/refs/tags/v0.99.6 -f sha=abc' \
    'gh api -X POST repos/cruk-mi/mesa/git/refs -f ref=refs/tags/v0.99.9 -f sha=abc' \
    'gh api -X POST repos/cruk-mi/mesa/git/refs -f ref=refs/heads/main -f sha=abc' \
    'gh api -XDELETE repos/cruk-mi/mesa/git/refs/tags/v0.99.6' \
    'gh api -X=DELETE repos/cruk-mi/mesa/git/refs/heads/main' \
    'gh api -XPATCH repos/cruk-mi/mesa/git/refs/heads/feat/x -Fforce=true -fsha=abc' \
    'gh api repos/cruk-mi/mesa/merges -fbase=main -fhead=feat/x' \
    "gh api graphql -fquery='mutation { deleteRef(input:{refId:\"X\"}) { clientMutationId } }'" \
    'gh api repos/cruk-mi/mesa/pulls/106/reviews -fevent=APPROVE' \
    "gh api https://api.github.com/graphql -f query='mutation { deleteRef(input:{refId:\"X\"}) { clientMutationId } }'"
do check 2 "$cmd"; done
# The same merge spelled with attached values or a full URL is still asked.
ask_check "gh api graphql -fquery='mutation { mergePullRequest(input:{pullRequestId:\"X\"}) { clientMutationId } }'"
ask_check "gh api https://api.github.com/graphql -f query='mutation { mergePullRequest(input:{pullRequestId:\"X\"}) { clientMutationId } }'"
ask_check 'gh api -XPUT repos/cruk-mi/mesa/pulls/106/merge'
check 0 'gh api -XPOST repos/cruk-mi/mesa/pulls/106/comments/1/replies -fbody=fixed'
# -i (--include) may lead a bundle with the value flag after it.
ask_check 'gh api -iX PUT repos/cruk-mi/mesa/pulls/1/merge'
ask_check "gh api graphql -ifquery='mutation { mergePullRequest(input:{pullRequestId:\"X\"}) { clientMutationId } }'"
for cmd in \
    'gh api -iXDELETE repos/cruk-mi/mesa/git/refs/tags/v0.99.6' \
    'gh api -iX DELETE repos/cruk-mi/mesa/git/refs/heads/stacked-base' \
    'gh api -iXDELETE repos/cruk-mi/mesa/git/refs/heads/main' \
    "gh api graphql -if query='mutation { deleteRef(input:{refId:\"X\"}) { clientMutationId } }'" \
    'gh api repos/cruk-mi/mesa/pulls/1/reviews -iFevent=APPROVE'
do check 2 "$cmd"; done
check 0 'gh api -i repos/cruk-mi/mesa/pulls/106'

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
    'gh api repos/cruk-mi/mesa/git/refs/heads/main' \
    'gh api -X POST repos/cruk-mi/mesa/git/refs -f ref=refs/heads/feat/x -f sha=abc' \
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
# Nothing uses feat/merged (its PR is closed): asked, and the prompt says so.
ask_check 'git push origin --delete feat/merged' 'no open PR'
ask_check 'gh pr merge 102 --squash --delete-branch' 'no open PR'

# gh (pflag) also spells --delete-branch as --delete-branch=true, and bundles
# -d with other shorthand. Each must reach the same check, and the prompt must
# name the PR and the branch it deletes.
check 2 'gh pr merge 108 --delete-branch=true'
check 2 'gh pr merge 108 -sd'
check 2 'gh pr merge 108 -ds'
check 2 'gh pr close 108 -d'
check 2 'gh pr close 108 --delete-branch'
check 2 'gh pr merge 108 -d=true'
check 2 'gh pr merge 108 -sd=true'
check 2 'gh pr close 108 -d=1'
check 2 'gh pr merge 999 -d'              # PR not found: cannot tell what goes
ask_check 'gh pr merge 102 --delete-branch=true' "PR #102 and delete its branch 'feat/clean'"
ask_check 'gh pr merge 102 -sd' "merge PR #102 and delete its branch 'feat/clean'"
ask_check 'gh pr merge 102 -ds' "delete its branch 'feat/clean'"
ask_check 'gh pr merge 102 -d' "delete its branch 'feat/clean'"
ask_check 'gh pr close 102 -d' "close PR #102 and delete its branch 'feat/clean'"
ask_check 'gh pr merge 102 -d=true' "merge PR #102 and delete its branch 'feat/clean'"
check 0 'gh pr close 102'

# Deleting the head branch of an open PR closes that PR: refused, except for
# the PR that `gh pr merge|close --delete-branch` itself acts on.
check 2 'git push origin --delete chore/open-head'
check 2 'gh api -X DELETE repos/cruk-mi/mesa/git/refs/heads/chore/open-head'
ask_check 'gh pr merge 110 --delete-branch' "PR #110 and delete its branch 'feat/on-top'"
ask_check 'git branch -D chore/open-head'   # local only: no PR is touched
# If gh cannot answer, nothing is verified: refused, not asked.
export GH_STUB_FAIL=1
check 2 'git push origin --delete feat/merged'
check 2 'gh api -X DELETE repos/cruk-mi/mesa/git/refs/heads/feat/merged'
check 2 'gh pr merge 102 --delete-branch'
ask_check 'gh pr merge 102'                 # no deletion, nothing to look up
export GH_STUB_FAIL=list
check 2 'git push origin --delete feat/merged'
unset GH_STUB_FAIL

# --- other spellings of a protected branch or a tag --------------------
# The remote resolves heads/main to refs/heads/main, and a short name to a tag
# when one exists (the fake gh knows v0.99.6), so each must be refused.
for cmd in \
    'git push origin :heads/main' \
    'git push origin :refs/heads/main' \
    'git push origin --delete heads/dev' \
    'git push origin HEAD:heads/main' \
    'git branch -D heads/main' \
    'git push origin --delete v0.99.6' \
    'git push origin :v0.99.6' \
    'git push origin --delete tags/v0.99.6' \
    'git push origin :tags/v0.99.6' \
    'git push origin :refs/tags/v0.99.6' \
    'git push origin --delete tag v0.99.6' \
    'git push origin --delete refs/remotes/origin/feat/merged' \
    'git push -dq origin v0.99.6' \
    'git push -uf origin my-branch' \
    'git push -o ci.skip -f origin my-branch' \
    'git push --prune origin refs/heads/*:refs/heads/*' \
    "git push origin 'refs/heads/*'" \
    'git tag --delete v0.99.6' \
    'git update-ref -d refs/tags/v0.99.6'
do check 2 "$cmd"; done
# git takes any unambiguous prefix of a long option: `--del` is `--delete`.
for cmd in \
    'git push origin --del stacked-base' \
    'git push origin --del v0.99.6' \
    'git push --force-w origin feat/x' \
    'git push --force-w=feat/x origin feat/x' \
    'git push origin --mirr' \
    'git push origin --pru' \
    'git push origin --al' \
    'git tag --del v0.99.6' \
    'git branch --del main'
do check 2 "$cmd"; done
ask_check 'git push origin --del feat/merged' "'feat/merged'"
ask_check 'git branch --del old-branch'
check 0 'git push --dry origin main'
check 0 'git push --set-up origin chore/agent-setup'
ask_check 'git push -dq origin feat/merged' "'feat/merged'"
ask_check 'git push -o ci.skip origin --delete feat/merged' "'feat/merged'"
check 0 'git push -o ci.skip origin chore/agent-setup'

# --- wrappers, prefixes and nested shells -------------------------------
# The rules follow the command wherever it actually runs: behind VAR=x, env,
# sudo, timeout or xargs; with gh's -R before or after `pr`; inside bash -c,
# sh -c or eval; and in $( ), backticks or a subshell.
for cmd in \
    'FOO=1 gh pr merge 102' \
    'env gh pr merge 102' \
    'env -i FOO=1 gh pr merge 102' \
    'command gh pr merge 102' \
    'timeout 30 gh pr merge 102' \
    'echo 102 | xargs gh pr merge' \
    '/opt/homebrew/bin/gh pr merge 102' \
    'gh -R cruk-mi/mesa pr merge 102' \
    'gh pr -R cruk-mi/mesa merge 102' \
    'gh --repo=cruk-mi/mesa pr merge 102' \
    'bash -c "gh pr merge 102"' \
    "sh -c 'gh pr merge 102'" \
    'bash -lc "gh pr merge 102"' \
    'eval gh pr merge 102' \
    "eval 'gh pr merge 102'" \
    'echo $(gh pr merge 102)' \
    'echo `gh pr merge 102`' \
    '(gh pr merge 102)'
do ask_check "$cmd"; done
for cmd in \
    'env git push origin main' \
    'sudo git push origin main' \
    'nohup git push --force origin x' \
    'bash -c "git push origin main"' \
    'sh -c "gh pr ready 102"' \
    'zsh -c "gh pr merge 108 -d"' \
    'eval "git push --force origin x"' \
    'echo $(git push origin main)' \
    'x=`git push origin main`' \
    'cat <(git push origin :main)' \
    'bash -c "bash -c \"git push origin main\""' \
    'bash -c "gh pr merge 102" && git push origin main' \
    'gh -R cruk-mi/mesa pr ready 102' \
    'gh -R cruk-mi/mesa pr merge 108 --delete-branch' \
    'gh pr -R cruk-mi/mesa review 102 --approve' \
    'FOO=1 gh pr edit 104 --ready'
do check 2 "$cmd"; done
for cmd in \
    'bash -c "git status"' \
    'echo $(git rev-parse HEAD)' \
    'env FOO=1 git push origin chore/agent-setup' \
    'gh -R cruk-mi/mesa pr view 102'
do check 0 "$cmd"; done
# --- permission modes -------------------------------------------------
# An "ask" is only a safeguard if a human sees it. bypassPermissions approves
# it automatically, auto may, and dontAsk denies it silently; a payload with
# no mode is not trusted either. There a merge or deletion is blocked outright.
mode_check() { # mode_check <mode> <expected-exit> <command>
    local actual err
    err="$(printf '%s' "$3" | MODE="$1" payload | "$PY" "$hook" 2>&1 >/dev/null)"
    actual=$?
    if [ "$actual" != "$2" ]; then
        printf 'FAIL (mode %s, expected %s, got %s): %s\n' "$1" "$2" "$actual" "$3"
        fails=$((fails + 1))
    elif [ "$2" = 2 ] && ! printf '%s' "$err" | grep -q 'run it themselves'; then
        printf 'FAIL (mode %s, block does not tell the human to run it): %s\n' "$1" "$3"
        fails=$((fails + 1))
    fi
}
for mode in bypassPermissions auto dontAsk -; do
    mode_check "$mode" 2 'gh pr merge 102'
    mode_check "$mode" 2 'gh pr merge 102 --squash --delete-branch'
    mode_check "$mode" 2 'git push origin --delete feat/merged'
    mode_check "$mode" 2 'git branch -D old-branch'
    mode_check "$mode" 0 'git status'
    mode_check "$mode" 0 'git push origin chore/agent-setup'
done
mode_check bypassPermissions 2 'bash -c "gh pr merge 102"'
for mode in default acceptEdits plan; do
    ask=$(printf '%s' 'gh pr merge 102' | MODE="$mode" payload | "$PY" "$hook" 2>/dev/null)
    printf '%s' "$ask" | grep -q '"ask"' || { echo "FAIL (mode $mode did not ask)"; fails=$((fails + 1)); }
done

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
