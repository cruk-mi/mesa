#!/usr/bin/env python3
"""PreToolUse guard for mesa's agent workflow contract.

AGENTS.md lets an agent push a feature branch and open a DRAFT pull request.
Everything past that point belongs to a human. This hook enforces that
mechanically so the contract does not depend on the model remembering it.

Blocked:
  * any push that targets a protected branch (main / dev / master)
  * force pushes, --mirror and --all
  * committing or merging while HEAD is on a protected branch
  * gh pr ready / review --approve
  * the same actions reached through `gh api`
  * deleting a protected branch, or a branch an open PR uses as its base
    (deleting it would close that PR instead of retargeting it)
  * deleting tags

Asked, not blocked:
  * merging a pull request (`gh pr merge`, or the same through `gh api`).
  * deleting any other branch (`git branch -d/-D`, `git push --delete`,
    `gh pr merge --delete-branch`, or the same through `gh api`).
  AGENTS.md allows these only when the human explicitly approves or asks, so
  the hook returns a permission "ask": Claude Code shows the command and the
  human confirms each one.

The command is split on shell separators and tokenised, so each segment is
judged on its own. That matters: `git push --dry-run origin main && git push
origin main` must not be waved through because the first half is harmless, and
`git -C /repo push origin main` must not slip past because `git` and `push` are
not adjacent.

Exit codes: 0 = allow, 2 = block (stderr is shown to the agent). A merge or a
safe branch deletion exits 0 with a PreToolUse "ask" decision on stdout.
Read-only inspection is never blocked. When this hook and AGENTS.md disagree,
that is a bug: fix both.
"""

import functools
import json
import os
import re
import shlex
import subprocess
import sys

PROTECTED = {"main", "dev", "master"}

# Git's own options that swallow the following token, so the subcommand can be
# located without mistaking an option's argument for it.
GIT_OPTS_WITH_VALUE = {"-C", "-c", "--git-dir", "--work-tree", "--namespace", "--exec-path"}

# `gh api` reaches every endpoint `gh pr` does, so the same four rules have to
# hold there or the contract is one flag away from being bypassed. Matched on
# the endpoint rather than the method: the path is what names the action.
GH_API_OPTS_WITH_VALUE = {
    "-X", "--method", "-f", "--field", "-F", "--raw-field", "-H", "--header",
    "-q", "--jq", "-t", "--template", "--input", "--cache", "-p", "--preview",
    "--hostname",
}
GH_API_MERGE = re.compile(r"/pulls/\d+/merge/?$")
GH_API_PULL = re.compile(r"/pulls/\d+/?$")
# Also the nested submit of a pending review: .../reviews/{id}/events.
GH_API_REVIEWS = re.compile(r"/pulls/\d+/reviews(/\d+/events)?/?$")
GH_API_REF = re.compile(r"/git/refs?(/|$)")

# `gh api graphql` carries the action in the query body, not the endpoint, so
# the path patterns above never see it. These are the mutations that do what
# the four rules forbid; everything else (resolveReviewThread, comments, ...)
# stays allowed.
GH_GRAPHQL_MERGE = re.compile(r"\b(mergePullRequest|enablePullRequestAutoMerge)\b")
GH_GRAPHQL_READY = re.compile(r"\bmarkPullRequestReadyForReview\b")
GH_GRAPHQL_REVIEW = re.compile(r"\b(addPullRequestReview|submitPullRequestReview)\b")
GH_GRAPHQL_DELETE_REF = re.compile(r"\bdeleteRef\b")



class Ask:
    """Not refused outright: needs the human's explicit approval, so main()
    turns it into a permission prompt instead of a block."""

    def __init__(self, what):
        self.what = what


MERGE = Ask("merge a pull request")

def gh_json(args):
    """Parsed JSON from a read-only `gh` call; None if it fails in any way.

    The tests put a fake `gh` first on PATH, so the lookups run offline.
    """
    try:
        out = subprocess.run(["gh", *args], capture_output=True, text=True, timeout=15)
        return json.loads(out.stdout) if out.returncode == 0 else None
    except (OSError, subprocess.SubprocessError, ValueError):
        return None


@functools.lru_cache(maxsize=None)
def open_prs():
    """Every open PR as {number, baseRefName, headRefName}; None if unknown."""
    prs = gh_json(["pr", "list", "--state", "open", "--limit", "500",
                   "--json", "number,baseRefName,headRefName"])
    if not isinstance(prs, list) or not all(isinstance(p, dict) for p in prs):
        return None
    return prs


def open_prs_based_on(branch):
    """Numbers of open PRs whose base is `branch`; None if it cannot be checked."""
    prs = open_prs()
    if prs is None:
        return None
    return [str(p.get("number")) for p in prs if p.get("baseRefName") == branch]


def pr_number_and_head(pr_args):
    """(number, head branch) of the PR a `gh pr` command targets; None if unknown."""
    pr = gh_json(["pr", "view", *pr_args, "--json", "number,headRefName"])
    if not isinstance(pr, dict):
        return None
    number, head = pr.get("number"), pr.get("headRefName")
    if not isinstance(number, int) or not isinstance(head, str) or not head:
        return None
    return number, head


def check_branch_deletion(branches, remote, what="delete branch"):
    """Block unsafe deletions; ask for the rest.

    A protected branch is never deleted. A remote branch that an open PR uses
    as its base is refused too: GitHub closes such a PR rather than
    retargeting it (#114 and #125 were closed that way), so it has to be
    retargeted first.
    """
    names = [re.sub(r"^refs/heads/", "", b) for b in branches if b]
    if not names:
        return None
    for name in names:
        if name in PROTECTED:
            return f"Deleting the protected branch '{name}' is not allowed."
    unchecked = []
    if remote:
        for name in names:
            stacked = open_prs_based_on(name)
            if stacked is None:
                unchecked.append(name)
            elif stacked:
                prs = ", ".join("#" + n for n in stacked)
                return (f"Open PR(s) {prs} use '{name}' as their base; deleting it "
                        f"would close them. Retarget first: gh pr edit <N> --base main.")
    what += " " + ", ".join(f"'{n}'" for n in names)
    if remote and not unchecked:
        what += " (checked: no open PR is based on it)"
    if unchecked:
        what += (" (could not check whether an open PR is based on "
                 + ", ".join(unchecked) + ")")
    return Ask(what)


def split_segments(command):
    """Split a shell command into separately-executed segments."""
    return [s for s in re.split(r"\|\||&&|;|\n|\||&", command) if s.strip()]


def tokenise(segment):
    try:
        return shlex.split(segment)
    except ValueError:
        # Unbalanced quotes — fall back to whitespace so we still inspect it.
        return segment.split()


def git_subcommand(tokens):
    """Return (subcommand, remaining_args) for a git invocation, else (None, [])."""
    if not tokens:
        return None, []
    # Strip a leading `env FOO=bar` or absolute path to git.
    i = 0
    while i < len(tokens) and ("=" in tokens[i] and not tokens[i].startswith("-")):
        i += 1
    if i >= len(tokens) or os.path.basename(tokens[i]) != "git":
        return None, []
    i += 1
    while i < len(tokens):
        tok = tokens[i]
        if tok in GIT_OPTS_WITH_VALUE:
            i += 2
            continue
        if tok.startswith("-"):
            i += 1
            continue
        return tok, tokens[i + 1:]
    return None, []


def targets_protected_ref(args):
    """True if any positional refspec of a push resolves to a protected branch."""
    for arg in args:
        if arg.startswith("-"):
            continue
        ref = arg.lstrip("+")
        # A refspec's destination is what matters: src:dst -> dst.
        if ":" in ref:
            ref = ref.split(":", 1)[1]
        ref = re.sub(r"^refs/heads/", "", ref)
        if ref in PROTECTED:
            return True
    return False


def gh_api_endpoint(tokens):
    """The endpoint argument of `gh api`, skipping flags and their values.

    Parsed rather than grepped so a reply body that merely mentions a path
    (`-f body='... /pulls/1/merge ...'`) is not mistaken for one.
    """
    try:
        i = next(n for n, t in enumerate(tokens) if not t.startswith("-")
                 and os.path.basename(t) != "gh") + 1
    except StopIteration:
        return None
    while i < len(tokens):
        tok = tokens[i]
        if tok in GH_API_OPTS_WITH_VALUE:
            i += 2
            continue
        if tok.startswith("-"):
            i += 1
            continue
        # Drop any query string or fragment: the patterns below are anchored
        # at the end of the path, so `.../merge?` would otherwise slip past.
        return re.split(r"[?#]", tok, maxsplit=1)[0]
    return None


def gh_api_method(tokens):
    """The HTTP method a `gh api` call will use."""
    for i, tok in enumerate(tokens):
        if tok in ("-X", "--method") and i + 1 < len(tokens):
            return tokens[i + 1].upper()
        if tok.startswith("--method="):
            return tok.split("=", 1)[1].upper()
    # gh switches to POST as soon as a field is supplied.
    if any(t in ("-f", "--field", "-F", "--raw-field") or
           t.startswith(("--field=", "--raw-field=")) for t in tokens):
        return "POST"
    return "GET"


def payload_text(tokens):
    """Everything a `gh api` call could send as its request body.

    A field can arrive inline (`-f query=...`), from a file (`-F query=@f`,
    `--input f`) or on stdin (`--input -`). Files are read so a mutation or an
    APPROVE event cannot hide in one; stdin cannot be inspected, so it is
    reported as None.
    """
    parts = []
    for i, tok in enumerate(tokens):
        value = None
        if tok in ("-f", "--field", "-F", "--raw-field", "--input") and i + 1 < len(tokens):
            value = tokens[i + 1]
        elif tok.startswith(("--field=", "--raw-field=", "--input=")):
            value = tok.split("=", 1)[1]
        if value is None:
            continue
        path = None
        if tok.startswith("--input"):
            path = value
        elif "=@" in value:
            path = value.split("=@", 1)[1]
        if path == "-":
            return None
        if path:
            try:
                with open(os.path.expanduser(path), encoding="utf-8") as fh:
                    parts.append(fh.read())
            except OSError:
                return None
        else:
            parts.append(value)
    return "\n".join(parts)


def check_gh_graphql(tokens):
    text = payload_text(tokens)
    if text is None:
        return ("A `gh api graphql` query read from stdin or an unreadable file "
                "cannot be checked, so it is refused. Pass it with -f query=...")
    if GH_GRAPHQL_MERGE.search(text):
        return MERGE
    if GH_GRAPHQL_READY.search(text):
        return ("PRs opened by an agent stay in DRAFT. Only a human marks one "
                "ready for review.")
    if GH_GRAPHQL_REVIEW.search(text) and "APPROVE" in text.upper():
        return "An agent does not approve pull requests on this repo."
    if GH_GRAPHQL_DELETE_REF.search(text):
        return Ask("delete a ref through the deleteRef mutation (the branch and "
                   "any PR stacked on it cannot be checked; prefer git push --delete)")
    return None


def check_gh_api(tokens):
    """`gh api` is not a read-only escape hatch from the `gh pr` rules.

    Reads stay allowed, and so do the writes an agent is meant to make —
    posting a review reply is `POST .../pulls/N/comments/ID/replies`, which
    none of these patterns match.
    """
    endpoint = gh_api_endpoint(tokens)
    if not endpoint:
        return None
    if endpoint.strip("/") == "graphql":
        return check_gh_graphql(tokens)
    method = gh_api_method(tokens)
    if GH_API_MERGE.search(endpoint):
        return MERGE
    # Only a write can submit a review; a GET that filters on "APPROVED" is a read.
    if GH_API_REVIEWS.search(endpoint) and method != "GET":
        body = payload_text(tokens)
        if body is None:
            return ("A review payload read from stdin or an unreadable file cannot "
                    "be checked, so it is refused. Pass it with -f event=...")
        if "APPROVE" in (" ".join(tokens) + "\n" + body).upper():
            return "An agent does not approve pull requests on this repo."
    if GH_API_PULL.search(endpoint) and method != "GET":
        return ("PRs opened by an agent stay in DRAFT. Only a human marks one "
                "ready for review.")
    if GH_API_REF.search(endpoint) and method == "DELETE":
        ref = re.split(r"/git/refs?/", endpoint, maxsplit=1)[-1].strip("/")
        if not ref.startswith("heads/"):
            return "Deleting tags or other refs is not allowed."
        return check_branch_deletion([ref[len("heads/"):]], remote=True)
    return None


GH_PR_OPTS_WITH_VALUE = {"-R", "--repo", "-t", "--subject", "-b", "--body",
                         "-F", "--body-file", "--match-head-commit", "-A",
                         "--author-email"}


def gh_pr_target(tokens):
    """Arguments that make `gh pr view` look at the PR `gh pr merge` targets."""
    args, i, seen = [], 0, []
    while i < len(tokens):
        tok = tokens[i]
        if tok in GH_PR_OPTS_WITH_VALUE and i + 1 < len(tokens):
            if tok in ("-R", "--repo"):
                args += [tok, tokens[i + 1]]
            i += 2
            continue
        if not tok.startswith("-"):
            seen.append(tok)
        i += 1
    # seen: gh, pr, merge|close, [number | url | branch]
    return (seen[3:4] if len(seen) > 3 else []) + args


def deletes_branch(tokens):
    """True if `gh pr merge` / `gh pr close` also deletes the head branch.

    gh (pflag) accepts `--delete-branch=true` and bundled shorthand such as
    `-sd`, so matching the exact tokens is not enough. A `d` bundled with any
    other letter counts; the cost of a false match is only an extra check.
    """
    return any(t == "--delete-branch" or t.startswith("--delete-branch=")
               or (re.fullmatch(r"-[A-Za-z]+", t) and "d" in t) for t in tokens)


def check_pr_branch_deletion(verb, tokens):
    """`gh pr merge|close --delete-branch`: name the PR and the branch it deletes."""
    pr = pr_number_and_head(gh_pr_target(tokens))
    if pr is None:
        return (f"Cannot tell which branch `gh pr {verb} --delete-branch` would delete, "
                "so it is refused. Ask the human to run it.")
    number, head = pr
    return check_branch_deletion([head], remote=True,
                                 what=f"{verb} PR #{number} and delete its branch")


def current_branch():
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            capture_output=True, text=True, timeout=5,
        )
        return out.stdout.strip() if out.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        return None


def check_segment(segment):
    """Return a refusal message for this segment, or None to allow it."""
    tokens = tokenise(segment)
    if not tokens:
        return None

    # --- gh ---------------------------------------------------------------
    if os.path.basename(tokens[0]) == "gh":
        rest = [t for t in tokens[1:] if not t.startswith("-")]
        if rest[:1] == ["api"]:
            return check_gh_api(tokens)
        if rest[:2] in (["pr", "merge"], ["pr", "close"]):
            if deletes_branch(tokens):
                return check_pr_branch_deletion(rest[1], tokens)
            return MERGE if rest[1] == "merge" else None
        if rest[:2] == ["pr", "ready"]:
            return ("PRs opened by an agent stay in DRAFT. Only a human marks one "
                    "ready for review.")
        if rest[:2] == ["pr", "edit"] and "--ready" in tokens:
            return "Only a human marks a pull request ready for review."
        if rest[:2] == ["pr", "review"] and "--approve" in tokens:
            return "An agent does not approve pull requests on this repo."
        return None

    sub, args = git_subcommand(tokens)
    if sub is None:
        return None

    # --- git push ---------------------------------------------------------
    if sub == "push":
        # A dry run mutates nothing, but only this segment is exempt.
        if "--dry-run" in args or "-n" in args:
            return None
        if any(a in ("--force", "-f", "--force-with-lease") or
               a.startswith("--force-with-lease=") for a in args):
            return ("Force pushing is not allowed — it can destroy published history. "
                    "If a branch genuinely needs rewriting, ask the human to do it.")
        if any(a.lstrip("+").startswith("+") or a.startswith("+") for a in args
               if not a.startswith("-")):
            return "Force pushing via a '+refspec' is not allowed."
        if "--mirror" in args:
            return "`git push --mirror` can delete remote refs and is not allowed."
        if "--all" in args:
            return "`git push --all` would push protected branches too."
        positional = [a for a in args if not a.startswith("-")]
        deleting = []
        if "--delete" in args or "-d" in args:
            deleting = positional[1:]
        deleting += [a[1:] for a in positional[1:] if a.startswith(":")]
        if deleting:
            if any(d.startswith("refs/tags/") for d in deleting):
                return "Deleting tags is not allowed. Tags are applied by the human after merge."
            return check_branch_deletion(deleting, remote=True)
        if targets_protected_ref(args):
            return ("Pushing to a protected branch (main/dev) is not allowed. "
                    "Push your feature branch instead and open a pull request.")
        # No refspec: git pushes the current branch to its upstream.
        if not [a for a in args if not a.startswith("-")][1:]:
            branch = current_branch()
            if branch in PROTECTED:
                return (f"HEAD is on '{branch}', so a bare `git push` would publish a "
                        "protected branch. Move your work to a feature branch.")
        return None

    # --- committing or merging on a protected branch ----------------------
    # Checked against live repo state, so it holds across separate tool calls,
    # not just within one `checkout && commit` string.
    if sub in ("commit", "merge", "rebase", "cherry-pick", "revert", "am"):
        branch = current_branch()
        if branch in PROTECTED:
            return (f"HEAD is on protected branch '{branch}'. Branch off main and open a "
                    "pull request instead of committing here.")
        return None

    # --- deletions --------------------------------------------------------
    if sub == "branch" and any(a in ("-D", "-d", "--delete") for a in args):
        return check_branch_deletion([a for a in args if not a.startswith("-")],
                                     remote=False)
    if sub == "tag" and any(a in ("-d", "--delete") for a in args):
        return "Deleting tags is not allowed. Tags are applied by the human after merge."

    return None


def main():
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        # A guard that crashes must not silently stop blocking.
        print("guard-remote: could not parse hook input; blocking to be safe.",
              file=sys.stderr)
        return 2

    if payload.get("tool_name") != "Bash":
        return 0

    command = payload.get("tool_input", {}).get("command", "")
    if not command:
        return 0

    asks = []
    for segment in split_segments(command):
        message = check_segment(segment)
        if isinstance(message, Ask):
            asks.append(message.what)
        elif message:
            print(
                f"Blocked by mesa's workflow contract (AGENTS.md):\n  {message}\n\n"
                f"Command: {command.strip()[:300]}\n\n"
                "Do not work around this — report it to the human instead.",
                file=sys.stderr,
            )
            return 2
    if asks:
        print(json.dumps({"hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "ask",
            "permissionDecisionReason": (
                "mesa (AGENTS.md) needs your explicit approval to "
                + "; ".join(dict.fromkeys(asks))
                + ". Allow only if you asked for it."),
        }}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
