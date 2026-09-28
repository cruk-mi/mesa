#!/usr/bin/env bash
# Regression tests for review-after-pr.py.
#
# Run after editing the hook:  bash .claude/hooks/test-review-after-pr.sh
# Under another interpreter:   PYTHON=/usr/bin/python3 bash .claude/hooks/test-review-after-pr.sh
# `gh` is a fake executable first on PATH, so nothing reaches GitHub. Throwaway
# git repos serve as a decoy cwd: the hook must size the PR from `gh`, not them.
#
# Exit 0 = all cases behave as expected.

set -uo pipefail
hook="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/review-after-pr.py"
PYTHON="${PYTHON:-python3}"
fails=0
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
url='https://github.com/cruk-mi/mesa/pull/999'

# Fake gh: logs its arguments, prints $FAKE_GH_JSON, or fails when FAKE_GH_FAIL=1.
mkdir -p "$tmp/bin" "$tmp/nogh"
cat > "$tmp/bin/gh" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$*" >> "$GH_LOG"
[ "${FAKE_GH_FAIL-}" = 1 ] && exit 1
printf '%s\n' "$FAKE_GH_JSON"
EOF
chmod +x "$tmp/bin/gh"
ln -s "$(command -v "$PYTHON")" "$tmp/nogh/python3"
export GH_LOG="$tmp/gh.log"
HOOK_PATH="$tmp/bin:$PATH"

# fixture <path:added:deleted>...  -> what the fake `gh pr view --json baseRefName,files` prints
fixture() {
    FAKE_GH_JSON="$("$PYTHON" -c 'import json,sys
files = [dict(zip(("path","additions","deletions"), a.rsplit(":", 2))) for a in sys.argv[1:]]
for f in files: f["additions"], f["deletions"] = int(f["additions"]), int(f["deletions"])
print(json.dumps({"baseRefName": "main", "files": files}))' "$@")"
    export FAKE_GH_JSON
}

# run_hook <cwd> <command> [tool stdout]  -> prints the hook's stdout, and its exit code last
run_hook() {
    local dir="$1" cmd="$2" out="${3-$url}"
    "$PYTHON" -c 'import json,sys; print(json.dumps({"tool_name":"Bash","cwd":sys.argv[1],
        "tool_input":{"command":sys.argv[2]},
        "tool_response":{"stdout":sys.argv[3],"stderr":"","interrupted":False}}))' \
        "$dir" "$cmd" "$out" | PATH="$HOOK_PATH" "$PYTHON" "$hook"
    echo "exit=$?"
}

# check <reviews|silent> <label> <cwd> <command> [tool stdout]
check() {
    local expected="$1" label="$2" result actual
    shift 2
    result="$(run_hook "$@")"
    case "$result" in
        *additionalContext*) actual=reviews ;;
        *) actual=silent ;;
    esac
    if [ "$actual" != "$expected" ] || [ "${result##*exit=}" != "0" ]; then
        printf 'FAIL (expected %s, got %s): %s\n' "$expected" "$actual" "$label"
        fails=$((fails + 1))
    fi
}

# new_repo <name>: a repo whose main has one commit, checked out on branch feat
new_repo() {
    local dir="$tmp/$1"
    mkdir -p "$dir"
    (
        cd "$dir" &&
        git init -q -b main &&
        git config user.email t@t && git config user.name t &&
        mkdir -p vignettes && printf '# mesa\n' > NEWS.md &&
        git add -A && git commit -qm init &&
        git checkout -qb feat
    ) >/dev/null 2>&1
    printf '%s' "$dir"
}

create='gh pr create --draft --title t --body b'
none="$tmp"   # not a git repo: the hook must not need one

# --- large: must ask for a review ------------------------------------
fixture vignettes/a.Rmd:200:0
check reviews '200 lines of vignette'            "$none" "$create"
check reviews 'compound push && create'          "$none" "git push -u origin feat && $create --base main"
check reviews 'absolute path to gh'              "$none" "/opt/homebrew/bin/gh pr create --draft"
fixture R/new.R:40:0
check reviews '40 lines under R/'                "$none" "$create"
fixture vignettes/f1.md:1:0 vignettes/f2.md:1:0 vignettes/f3.md:1:0 \
        vignettes/f4.md:1:0 vignettes/f5.md:1:0 vignettes/f6.md:1:0
check reviews '6 small files'                    "$none" "$create"

# --- small or not applicable: must stay silent -----------------------
fixture DESCRIPTION:1:1 NEWS.md:60:0
check silent 'version bump (DESCRIPTION + NEWS)' "$none" "$create"
fixture vignettes/a.Rmd:20:0
check silent '20-line docs change'               "$none" "$create"
fixture man/f1.Rd:100:0 man/f2.Rd:100:0 man/f3.Rd:100:0 man/f4.Rd:100:0 \
        man/f5.Rd:100:0 man/f6.Rd:100:0 NAMESPACE:1:0
check silent 'regenerated man/ and NAMESPACE'    "$none" "$create"

fixture vignettes/a.Rmd:200:0
check silent 'not gh pr create'                  "$none" 'git status'
check silent 'gh pr view, not create'            "$none" 'gh pr view 999'
check silent 'create in a quoted body only'      "$none" "git commit -m 'run gh pr create later'"
check silent 'create failed (no PR URL)'         "$none" "$create" 'pull request create failed'
FAKE_GH_FAIL=1 check silent 'gh pr view fails'   "$none" "$create"
FAKE_GH_JSON='not json' check silent 'gh prints bad JSON' "$none" "$create"
HOOK_PATH="$tmp/nogh" check silent 'gh not on PATH' "$none" "$create"

# --- the PR, not the session's checkout, decides (review #138) --------
ahead="$(new_repo ahead)"   # feat is 200 lines ahead of main
seq 1 200 > "$ahead/vignettes/a.Rmd"
(cd "$ahead" && git add -A && git commit -qm change) >/dev/null 2>&1
level="$(new_repo level)"   # feat is level with main

fixture R/new.R:200:0
check reviews 'cd <worktree> && gh pr create'    "$level" "cd '$ahead' && $create"
: > "$GH_LOG"
check reviews 'sized by one gh pr view call'     "$level" "$create"
if [ "$(cat "$GH_LOG")" != "pr view 999 --repo cruk-mi/mesa --json baseRefName,files" ]; then
    printf 'FAIL (gh called as: %s): sized by one gh pr view call\n' "$(cat "$GH_LOG")"
    fails=$((fails + 1))
fi

fixture R/stacked.R:1:0
heredoc="gh pr create --title t --body \"\$(cat <<'EOF'
Stacked on feat.
EOF
)\" --base feat"
check silent '--base after a heredoc body'       "$ahead" "$heredoc"

# --- only a real `gh … pr create` command, only a bare PR URL (review #138)
fixture R/new.R:200:0
check reviews 'gh -R owner/repo pr create'       "$none" "gh -R cruk-mi/mesa pr create --draft"
check reviews 'gh --repo=owner/repo pr create'   "$none" "gh --repo=cruk-mi/mesa pr create --draft"
check reviews 'env VAR= prefix'                  "$none" "env GH_PROMPT_DISABLED=1 gh pr create --draft"
check reviews 'VAR= prefix'                      "$none" "GH_REPO=cruk-mi/mesa gh pr create --draft"
check reviews 'progress line before the URL'     "$none" "$create" "$(printf 'Creating draft pull request for feat into main\n\n%s\n' "$url")"
check reviews 'create inside ( … )'              "$none" "(cd /x && gh pr create --draft)"
check silent 'echo gh pr create'                 "$none" "echo gh pr create"
comment="gh pr comment 125 --body \"\$(cat <<'EOF'
Next step: gh pr create --draft
EOF
)\""
check silent 'gh pr comment body mentions create' "$none" "$comment" 'https://github.com/cruk-mi/mesa/pull/125#issuecomment-123'
unquoted="cat > notes.md <<EOF
gh pr create --draft
EOF"
check silent 'create inside an unquoted heredoc' "$none" "$unquoted"
check silent 'issuecomment link, not a PR URL'   "$none" "$create" 'https://github.com/cruk-mi/mesa/pull/999#issuecomment-1'
check silent 'discussion link, not a PR URL'     "$none" "$create" 'https://github.com/cruk-mi/mesa/pull/999#discussion_r1'
check silent 'push hint /pull/new/<branch>'      "$none" "git push -u origin feat" 'remote:   https://github.com/cruk-mi/mesa/pull/new/feat'

result="$(CI=true run_hook "$none" "$create")"
case "$result" in
    *additionalContext*) printf 'FAIL (expected silent in CI)\n'; fails=$((fails + 1)) ;;
esac

if [ "$fails" -eq 0 ]; then
    echo "review-after-pr: all cases pass"
else
    echo "review-after-pr: $fails failing case(s)"
fi
exit $((fails > 0))
