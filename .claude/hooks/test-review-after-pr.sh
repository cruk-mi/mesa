#!/usr/bin/env bash
# Regression tests for review-after-pr.py.
#
# Run after editing the hook:  bash .claude/hooks/test-review-after-pr.sh
# Builds throwaway repos in a temp dir; touches nothing in this checkout.
#
# Exit 0 = all cases behave as expected.

set -uo pipefail
hook="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/review-after-pr.py"
fails=0
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
url='https://github.com/cruk-mi/mesa/pull/999'

# run_hook <dir> <command> [tool stdout]  -> prints the hook's stdout, and its exit code last
run_hook() {
    local dir="$1" cmd="$2" out="${3-$url}"
    python3 -c 'import json,sys; print(json.dumps({"tool_name":"Bash","cwd":sys.argv[1],
        "tool_input":{"command":sys.argv[2]},
        "tool_response":{"stdout":sys.argv[3],"stderr":"","interrupted":False}}))' \
        "$dir" "$cmd" "$out" | python3 "$hook"
    echo "exit=$?"
}

# check <reviews|silent> <label> <dir> <command> [tool stdout]
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
        mkdir -p R man vignettes &&
        printf 'Version: 0.99.7\n' > DESCRIPTION && printf '# mesa\n' > NEWS.md &&
        git add -A && git commit -qm init &&
        git checkout -qb feat
    ) >/dev/null 2>&1
    printf '%s' "$dir"
}

lines() { seq 1 "$1" | sed 's/^/x <- /'; }
commit_all() { (cd "$1" && git add -A && git commit -qm change) >/dev/null 2>&1; }

create='gh pr create --draft --title t --body b'

# --- large: must ask for a review ------------------------------------
big="$(new_repo big)"
lines 200 > "$big/vignettes/a.Rmd"; commit_all "$big"
check reviews '200 lines of vignette'            "$big" "$create"
check reviews 'compound push && create'          "$big" "git push -u origin feat && $create --base main"
check reviews '--base=main form'                 "$big" "$create --base=main"
check reviews 'absolute path to gh'              "$big" "/opt/homebrew/bin/gh pr create --draft"

hot="$(new_repo hot)"
lines 40 > "$hot/R/new.R"; commit_all "$hot"
check reviews '40 lines under R/'                "$hot" "$create"

many="$(new_repo many)"
for i in 1 2 3 4 5 6; do echo "# $i" > "$many/vignettes/f$i.md"; done; commit_all "$many"
check reviews '6 small files'                    "$many" "$create"

# --- small or not applicable: must stay silent -----------------------
bump="$(new_repo bump)"
printf 'Version: 0.99.7.9000\n' > "$bump/DESCRIPTION"; lines 60 >> "$bump/NEWS.md"
commit_all "$bump"
check silent 'version bump (DESCRIPTION + NEWS)' "$bump" "$create"

docs="$(new_repo docs)"
lines 20 > "$docs/vignettes/a.Rmd"; commit_all "$docs"
check silent '20-line docs change'               "$docs" "$create"

gen="$(new_repo gen)"
for i in 1 2 3 4 5 6 7 8; do lines 100 > "$gen/man/f$i.Rd"; done
printf 'export(x)\n' > "$gen/NAMESPACE"; commit_all "$gen"
check silent 'regenerated man/ and NAMESPACE'    "$gen" "$create"

check silent 'not gh pr create'                  "$big" 'git status'
check silent 'gh pr view, not create'            "$big" 'gh pr view 999'
check silent 'create in a quoted body only'      "$big" "git commit -m 'run gh pr create later'"
check silent 'create failed (no PR URL)'         "$big" "$create" 'pull request create failed'
check silent 'unknown base branch'               "$big" "$create --base no-such-branch"
check silent 'not a git repo'                    "$tmp" "$create"

result="$(CI=true run_hook "$big" "$create")"
case "$result" in
    *additionalContext*) printf 'FAIL (expected silent in CI)\n'; fails=$((fails + 1)) ;;
esac

if [ "$fails" -eq 0 ]; then
    echo "review-after-pr: all cases pass"
else
    echo "review-after-pr: $fails failing case(s)"
fi
exit $((fails > 0))
