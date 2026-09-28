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

# --- Stop suggests /mesa-status once, and never drives it -------------
# Hooks must never lead to an unprompted edit of the public #124: Stop only
# tells the human, and it never blocks (a block makes the model act).
rm -f "$flags"/refresh-*
out="$(echo '{}' | "$PY" "$hook" stop)"
[ -z "$out" ] || fail "stop spoke with no flag raised"

post 'gh pr merge 1'
out="$(echo '{"stop_hook_active": false}' | "$PY" "$hook" stop)"
echo "$out" | grep -q '"systemMessage"' || fail "stop did not suggest /mesa-status with the flag raised"
echo "$out" | grep -q '/mesa-status' || fail "stop's suggestion does not name /mesa-status"
echo "$out" | grep -q '"decision"' && fail "stop blocked, which has the model act instead of the human"
echo "$out" | grep -q 'sync-issue' && fail "stop's message asks for --sync-issue"
[ -f "$flags/refresh-pending" ] || fail "stop did not lower the flag to pending"
out="$(echo '{"stop_hook_active": false}' | "$PY" "$hook" stop)"
[ -z "$out" ] || fail "stop spoke a second time for the same change"

post 'gh pr merge 2'
out="$(echo '{"stop_hook_active": true}' | "$PY" "$hook" stop)"
[ -z "$out" ] || fail "stop spoke while stop_hook_active"

# --- SessionStart tells the human about a pending refresh ---------------
out="$(echo '{}' | "$PY" "$hook" session)"
echo "$out" | grep -q '"systemMessage"' || fail "session start did not address the human"
echo "$out" | grep -q 'not republished' || fail "session start did not report the pending refresh"
echo "$out" | grep -q 'sync-issue' && fail "session start's message asks for --sync-issue"
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
# Output goes to a file: bash 3.2 misparses a heredoc inside $(...).
"$PY" - "$status_script" >"$flags/roadmap-cases.out" 2>&1 <<'PYEOF'
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
        verdict = ms.sync_issue(roadmap, current.replace("\n", "\r\n"), None, True,
                                roadmap.pop("toml"), False)
        if "already current" not in verdict:
            fail(f"CRLF body that is current: {verdict!r}")

# TOML that parses but has the wrong shape degrades; it never raises.
BAD_SHAPES = {
    "refs = 5": TOML.replace("refs = [10]", "refs = 5"),
    'wave = "w1" at the top': 'wave = "w1"\n[[item]]\nwave = "w1"\ntitle = "x"\n',
    "a single [wave] table": '[wave]\nid = "w1"\n[[item]]\nwave = "w1"\ntitle = "x"\n',
    "wave without id, item without wave": '[[wave]]\nname = "W"\n[[item]]\ntitle = "x"\n',
    'after = "x"': TOML.replace('after = ["10"]', 'after = "10"'),
    'check = "tag:v1"': TOML.replace('check = ["branch-gone:old"]', 'check = "tag:v1"'),
    'release = "0.99.8"': 'release = "0.99.8"\n' + TOML,
    "title = 5": TOML.replace('title = "Fix A"', "title = 5"),
}
for name, text in BAD_SHAPES.items():
    if not needs_toml(f"bad shape: {name}"):
        break
    ms.DEGRADED.clear()
    try:
        if ms.parse_roadmap(text) is not None:
            fail(f"bad shape passed validation: {name}")
        elif not ms.DEGRADED:
            fail(f"bad shape gave no reason: {name}")
    except Exception as exc:  # the bug: a traceback aborts the whole refresh
        fail(f"bad shape raised {type(exc).__name__}: {name}")

# `after` must point backwards: forward, self and cyclic waits are reported.
TWO = '[[wave]]\nid = "w1"\n[[item]]\nkey = "a"\nwave = "w1"\ntitle = "A"\n{a}\n' \
      '[[item]]\nkey = "b"\nwave = "w1"\ntitle = "B"\ndone = true\n{b}\n'
BAD_AFTER = {
    "forward": TWO.format(a='after = ["b"]', b=""),
    "self": TWO.format(a='after = ["a"]', b=""),
    "cycle": TWO.format(a='after = ["b"]', b='after = ["a"]'),
}
for name, text in BAD_AFTER.items():
    if not needs_toml(f"after: {name}"):
        break
    ms.DEGRADED.clear()
    if ms.parse_roadmap(text) is not None:
        fail(f"{name} `after` passed validation")
    elif not any("not listed before it" in d for d in ms.DEGRADED):
        fail(f"{name} `after` gave the wrong reason: {ms.DEGRADED}")

# A PR closed without merging landed nothing: its item is not Done.
closed = {"__typename": "PullRequest", "state": "CLOSED", "closedAt": "2026-09-02T00:00:00Z",
          "url": "u10", "title": "A"}
if ms.describe_ref(10, closed)["resolved"]:
    fail("a PR closed without merging counts as resolved")
if needs_toml("closed PR status"):
    roadmap, _ = roadmap_for(live_body(), states={**STATES, 10: closed})
    if roadmap is None or roadmap["items"][0]["status"] == "done":
        fail("an item whose PR was closed without merging shows Done")


def synced(roadmap, body, toml_text, from_file=False):
    """Run a real (stubbed) sync as main() does; (verdict, body written or None)."""
    if os.path.exists(ms.ROADMAP_MD):
        os.unlink(ms.ROADMAP_MD)
    verdict = ms.sync_issue(roadmap, body, None, False, toml_text, from_file)
    if not os.path.exists(ms.ROADMAP_MD):
        return verdict, None
    with open(ms.ROADMAP_MD, encoding="utf-8") as handle:
        return verdict, handle.read()


# A bare ```toml fence (data markers deleted on github.com): no KeyError, and
# the new marked block replaces the fence instead of sitting next to it.
if needs_toml("bare fence, generated markers present"):
    roadmap, body = roadmap_for(live_body(data=False))
    toml_text = roadmap.pop("toml")  # as main() does, before the sync
    try:
        verdict, new = synced(roadmap, body, toml_text)
    except Exception as exc:
        fail(f"bare fence: sync raised {type(exc).__name__}: {exc}")
    else:
        if new is None:
            fail(f"bare fence: nothing written ({verdict!r})")
        elif new.count(ms.DATA_START) != 1 or new.count("```toml") != 1:
            fail("bare fence: #124 would carry the TOML twice")
        elif "Hand-written note." not in new:
            fail("bare fence: hand-written text dropped")

# Data markers kept but generated markers tidied away: no TypeError, one data
# block, and the hand-written text survives.
if needs_toml("data markers without generated markers"):
    roadmap, body = roadmap_for(live_body(gen=False))
    toml_text = roadmap.pop("toml")
    try:
        verdict, new = synced(roadmap, body, toml_text)
    except Exception as exc:
        fail(f"no generated markers: sync raised {type(exc).__name__}: {exc}")
    else:
        if new is None:
            fail(f"no generated markers: nothing written ({verdict!r})")
        elif new.count(ms.DATA_START) != 1 or new.count("```toml") != 1:
            fail("no generated markers: the data block is duplicated or lost")
        elif new.count(ms.GEN_START) != 1:
            fail("no generated markers: the generated half was not put back")
        elif "Hand-written note." not in new:
            fail("no generated markers: hand-written text dropped")

# --sync-issue never writes #124 from incomplete data, dry run or not.
def edits():
    return [c for c in CALLS if c[:3] == ["gh", "issue", "edit"]]


def expect_refusal(name, roadmap, body, toml_text, from_file=False):
    for dry_run in (True, False):
        CALLS.clear()
        verdict = ms.sync_issue(roadmap, body, None, dry_run, toml_text, from_file)
        if "left alone" not in verdict or edits():
            fail(f"{name}: synced #124 anyway (dry_run={dry_run}, {verdict!r})")


expect_refusal("roadmap unknown", None, live_body(), None)
if needs_toml("refusals on degraded data"):
    roadmap, body = roadmap_for(live_body(), states=None)
    expect_refusal("issue/PR states not fetched", roadmap, body, roadmap.pop("toml"))
    roadmap, body = roadmap_for(live_body(), refs=(None, None))
    expect_refusal("ls-remote failed", roadmap, body, roadmap.pop("toml"))
    roadmap, body = roadmap_for(live_body(), states={10: STATES[10]})
    if ms.DEGRADED:
        fail(f"a ref missing from the answer degraded by itself: {ms.DEGRADED}")
    expect_refusal("one ref state unknown", roadmap, body, roadmap.pop("toml"))
    seed = os.path.join(tempfile.mkdtemp(), "seed.toml")
    with open(seed, "w", encoding="utf-8") as handle:
        handle.write(TOML)
    roadmap, body = roadmap_for(None, override=seed)
    ms.DEGRADED.clear()  # isolate the missing body from the fetch failure
    expect_refusal("#124 body unknown with --roadmap-file", roadmap, body,
                   roadmap.pop("toml"), from_file=True)
    # Positive control: complete data does reach `gh issue edit`, exactly once.
    roadmap, body = roadmap_for(live_body())
    verdict, new = synced(roadmap, body, roadmap.pop("toml"))
    if len(edits()) != 1 or new is None:
        fail(f"complete data did not sync: {verdict!r}, {len(edits())} edit(s)")
PYEOF
while IFS= read -r line; do
    case "$line" in
        "") ;;
        SKIP:*) echo "$line" ;;
        FAIL:*) fail "${line#FAIL: }" ;;
        *) fail "mesa-status.py: $line" ;;
    esac
done <"$flags/roadmap-cases.out"

rm -rf "$flags"
if [ "$fails" -eq 0 ]; then
    echo "refresh-flag: all cases pass"
else
    echo "refresh-flag: $fails failing case(s)"
fi
exit $((fails > 0))
