#!/usr/bin/env python3
"""PreToolUse guard for mesa's agent workflow contract.

AGENTS.md lets an agent push a feature branch and open a DRAFT pull request.
Everything past that point belongs to a human. This hook enforces that
mechanically so the contract does not depend on the model remembering it.

Blocked:
  * any push that targets a protected branch (main / dev / master)
  * force pushes, --mirror, --all, and remote ref deletion
  * committing or merging while HEAD is on a protected branch
  * gh pr merge / ready / review --approve
  * the same four actions reached through `gh api`
  * deleting branches or tags

The command is split on shell separators and tokenised, so each segment is
judged on its own. That matters: `git push --dry-run origin main && git push
origin main` must not be waved through because the first half is harmless, and
`git -C /repo push origin main` must not slip past because `git` and `push` are
not adjacent.

Exit codes: 0 = allow, 2 = block (stderr is shown to the agent).
Read-only inspection is never blocked. When this hook and AGENTS.md disagree,
that is a bug: fix both.
"""

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
        return "Merging pull requests is the human's decision, not the agent's."
    if GH_GRAPHQL_READY.search(text):
        return ("PRs opened by an agent stay in DRAFT. Only a human marks one "
                "ready for review.")
    if GH_GRAPHQL_REVIEW.search(text) and "APPROVE" in text.upper():
        return "An agent does not approve pull requests on this repo."
    if GH_GRAPHQL_DELETE_REF.search(text):
        return "Deleting branches is not allowed — they are the human's audit trail."
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
        return "Merging pull requests is the human's decision, not the agent's."
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
        return "Deleting branches is not allowed — they are the human's audit trail."
    return None


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
        if rest[:2] == ["pr", "merge"]:
            return "Merging pull requests is the human's decision, not the agent's."
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
        if "--delete" in args or "-d" in args:
            return "Deleting remote refs is not allowed."
        if any(a.startswith(":") for a in args if not a.startswith("-")):
            return "Deleting a remote ref via ':ref' is not allowed."
        if targets_protected_ref(args):
            return ("Pushing to a protected branch (main/dev) is not allowed. "
                    "Push your feature branch instead and let a human merge.")
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
        return "Deleting branches is not allowed — they are the human's audit trail."
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

    for segment in split_segments(command):
        message = check_segment(segment)
        if message:
            print(
                f"Blocked by mesa's workflow contract (AGENTS.md):\n  {message}\n\n"
                f"Command: {command.strip()[:300]}\n\n"
                "Do not work around this — report it to the human instead.",
                file=sys.stderr,
            )
            return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
