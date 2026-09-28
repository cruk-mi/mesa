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
import urllib.parse

PROTECTED = {"main", "dev", "master"}

# Permission modes in which Claude Code shows a hook's "ask" to the human.
# bypassPermissions approves it automatically, auto may (a hook ask forces a
# prompt there only from v2.1.211), and dontAsk denies it without a word. In
# those, or when the mode is missing, an ask would not reach a human, so a
# merge or deletion is blocked outright instead.
PROMPTING_MODES = {"default", "acceptEdits", "plan"}

# Git's own options that swallow the following token, so the subcommand can be
# located without mistaking an option's argument for it.
GIT_OPTS_WITH_VALUE = {"-C", "-c", "--git-dir", "--work-tree", "--namespace", "--exec-path"}
# Same for `git push`, so an option's value is not read as the remote or a ref.
PUSH_OPTS_WITH_VALUE = {"-o", "--push-option", "--repo", "--receive-pack", "--exec"}

# Commands that run the command after them. The rules follow gh or git behind
# them (`env FOO=1 gh pr merge`, `sudo git push`, `xargs gh pr merge`), and
# into the string a shell -c or eval runs.
WRAPPERS = {"env", "command", "builtin", "exec", "nohup", "time", "nice", "timeout",
            "sudo", "xargs", "stdbuf", "caffeinate"}
SHELLS = {"bash", "sh", "zsh", "dash", "ksh"}
# Beyond this many levels of bash -c / eval / $( ) the command is refused.
MAX_DEPTH = 5

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
# POST .../merges merges one branch into another, main included, with no PR.
GH_API_MERGES = re.compile(r"/merges/?$")

# `gh api graphql` carries the action in the query body, not the endpoint, so
# the path patterns above never see it. These are the mutations that do what
# the four rules forbid; everything else (resolveReviewThread, comments, ...)
# stays allowed.
GH_GRAPHQL_MERGE = re.compile(r"\b(mergePullRequest|enablePullRequestAutoMerge)\b")
GH_GRAPHQL_READY = re.compile(r"\bmarkPullRequestReadyForReview\b")
GH_GRAPHQL_REVIEW = re.compile(r"\b(addPullRequestReview|submitPullRequestReview)\b")
# These take an opaque ref ID or a batch, so which ref they delete, move or
# merge into cannot be checked: main or a tag looks the same as a feature
# branch. git push covers every legitimate use.
GH_GRAPHQL_REF_WRITE = re.compile(r"\b(deleteRef|updateRefs?|createRef|mergeBranch)\b")



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


def pr_number_and_head(pr_args):
    """(number, head branch) of the PR a `gh pr` command targets; None if unknown."""
    pr = gh_json(["pr", "view", *pr_args, "--json", "number,headRefName"])
    if not isinstance(pr, dict):
        return None
    number, head = pr.get("number"), pr.get("headRefName")
    if not isinstance(number, int) or not isinstance(head, str) or not head:
        return None
    return number, head


def branch_name(ref):
    """The branch `ref` names, or None if it names a tag or any other ref.

    The remote resolves `heads/main` and `refs/heads/main` alike to main, so
    both prefixes are dropped; anything else under refs/ or tags/ is not a
    branch.
    """
    ref = re.sub(r"^(refs/)?heads/", "", ref)
    return None if ref.startswith(("refs/", "tags/")) else ref


def is_remote_tag(name):
    """True if `name` is a tag on GitHub; None if that cannot be checked.

    A short name in `git push --delete` resolves to a tag when one exists.
    """
    refs = gh_json(["api", "repos/{owner}/{repo}/git/matching-refs/tags/"
                    + urllib.parse.quote(name, safe="/")])
    if not isinstance(refs, list):
        return None
    return any(isinstance(r, dict) and r.get("ref") == "refs/tags/" + name for r in refs)


def check_branch_deletion(branches, remote, what="delete branch", exempt=None):
    """Block unsafe deletions; ask for the rest.

    A protected branch or a tag is never deleted, however it is spelled. A
    remote branch is refused while an open PR uses it: as its head, because
    deleting it closes that PR (AGENTS.md: "its PR is merged or closed"), or
    as its base, because GitHub then closes the stacked PR rather than
    retargeting it (#114 and #125 were closed that way). `exempt` is the PR
    that `gh pr merge|close --delete-branch` acts on. If any of this cannot
    be checked, the deletion is refused rather than asked.
    """
    names = []
    for ref in branches:
        if not ref:
            continue
        name = branch_name(ref)
        if name is None:
            return f"Deleting '{ref}' is not allowed: only branches may be deleted, never tags."
        names.append(name)
    if not names:
        return None
    for name in names:
        if name in PROTECTED:
            return f"Deleting the protected branch '{name}' is not allowed."
    if remote:
        prs = open_prs()
        if prs is None:
            return ("Could not list the open PRs to check that deleting "
                    + ", ".join(f"'{n}'" for n in names)
                    + " closes none of them, so it is refused. Ask the human to run it.")
        for name in names:
            tag = is_remote_tag(name)
            if tag is None:
                return (f"Could not check whether '{name}' is a tag on GitHub, so deleting "
                        "it is refused. Ask the human to run it.")
            if tag:
                return f"'{name}' is a tag on GitHub. Deleting tags is not allowed."
        for name in names:
            heads = [p.get("number") for p in prs
                     if p.get("headRefName") == name and p.get("number") != exempt]
            if heads:
                return (f"'{name}' is the head branch of open PR(s) "
                        + ", ".join(f"#{n}" for n in heads)
                        + "; deleting it would close them. Merge or close them first.")
            stacked = [p.get("number") for p in prs if p.get("baseRefName") == name]
            if stacked:
                return (f"Open PR(s) " + ", ".join(f"#{n}" for n in stacked)
                        + f" use '{name}' as their base; deleting it would close them. "
                        "Retarget first: gh pr edit <N> --base main.")
    what += " " + ", ".join(f"'{n}'" for n in names)
    if remote:
        what += " (checked: no open PR uses it as its head or base)"
    return Ask(what)


def substitutions(command):
    """Commands run inside $( ), <( ), >( ) or backticks, at every depth.

    Read on the raw text, even inside quotes: a false match costs only a
    check of harmless text, a missed one lets a merge through unseen.
    """
    found = re.findall(r"`([^`]*)`", command)
    for start in [m.end() for m in re.finditer(r"[$<>]\(", command)]:
        depth = 1
        for end in range(start, len(command)):
            depth += {"(": 1, ")": -1}.get(command[end], 0)
            if depth == 0:
                found.append(command[start:end])
                break
        else:
            found.append(command[start:])
    return found


def split_segments(command):
    """Split a shell command into separately-executed segments."""
    return [s for s in re.split(r"\|\||&&|;|\n|\||&", command) if s.strip()]


def tokenise(segment):
    try:
        return shlex.split(segment)
    except ValueError:
        # Unbalanced quotes — fall back to whitespace so we still inspect it.
        return segment.split()


def command_tokens(tokens):
    """The tokens from the command that actually runs.

    Drops `VAR=value` prefixes, a subshell's brackets, and wrappers such as
    env, sudo, timeout or xargs (with their options), so `FOO=1 gh pr merge`
    is judged as `gh pr merge`.
    """
    tokens = list(tokens)
    if tokens and tokens[0][:1] in ("(", "{", "!"):
        tokens[0] = tokens[0].lstrip("({! ")
        if tokens[-1].endswith((")", "}")):
            tokens[-1] = tokens[-1].rstrip(")} ;")
        tokens = [t for t in tokens if t]
    known = WRAPPERS | SHELLS | {"gh", "git", "eval"}
    while tokens:
        name = os.path.basename(tokens[0])
        if re.match(r"[A-Za-z_][A-Za-z0-9_]*=", tokens[0]):
            tokens = tokens[1:]
        elif name in WRAPPERS:
            # A wrapper's own options and values come first; the command is
            # the first known program after it.
            nxt = next((n for n, t in enumerate(tokens[1:], 1)
                        if os.path.basename(t) in known), None)
            if nxt is None:
                return tokens
            tokens = tokens[nxt:]
        else:
            return tokens
    return tokens


def gh_positionals(tokens):
    """gh's subcommand words, skipping `-R/--repo <repo>` wherever it sits."""
    out, i = [], 1
    while i < len(tokens):
        if tokens[i] in ("-R", "--repo"):
            i += 2
            continue
        if not tokens[i].startswith("-"):
            out.append(tokens[i])
        i += 1
    return out


def git_subcommand(tokens):
    """Return (subcommand, remaining_args) for a git invocation, else (None, [])."""
    if not tokens or os.path.basename(tokens[0]) != "git":
        return None, []
    i = 1
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


def short_flags(args):
    """Letters of every bundled short option (`-uf` -> {"u", "f"})."""
    return {c for a in args if re.fullmatch(r"-[A-Za-z0-9]+", a) for c in a[1:]}


def push_positionals(args):
    """The remote and refspecs of a `git push`, skipping options and their values."""
    out, i = [], 0
    while i < len(args):
        if args[i] in PUSH_OPTS_WITH_VALUE:
            i += 2
            continue
        if not args[i].startswith("-"):
            out.append(args[i])
        i += 1
    return out


def targets_protected_ref(positional):
    """True if any refspec of a push resolves to a protected branch."""
    for arg in positional[1:]:
        ref = arg.lstrip("+")
        # A refspec's destination is what matters: src:dst -> dst.
        if ":" in ref:
            ref = ref.split(":", 1)[1]
        if branch_name(ref) in PROTECTED:
            return True
    return False


def gh_api_endpoint(tokens):
    """The endpoint argument of `gh api`, skipping flags and their values.

    Parsed rather than grepped so a reply body that merely mentions a path
    (`-f body='... /pulls/1/merge ...'`) is not mistaken for one.
    """
    try:
        i = tokens.index("api") + 1
    except ValueError:
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
    ref_write = GH_GRAPHQL_REF_WRITE.search(text)
    if ref_write:
        return (f"The {ref_write.group(1)} mutation is not allowed: which ref it changes "
                "cannot be checked. Use git push (or git push origin --delete <branch>).")
    return None


def check_ref_write(endpoint, method, tokens):
    """PATCH or POST on .../git/refs: moving or creating main or a tag, or
    forcing any ref, is the push the rules forbid, made through the API."""
    body = payload_text(tokens)
    if body is None:
        return ("A ref update read from stdin or an unreadable file cannot be checked, "
                "so it is refused.")
    if re.search(r'force"?\s*[=:]\s*"?true', body, re.IGNORECASE):
        return ("Force-updating a ref is not allowed — it can destroy published history. "
                "If a branch genuinely needs rewriting, ask the human to do it.")
    if method == "PATCH":
        ref = re.split(r"/git/refs?/", endpoint, maxsplit=1)[-1].strip("/")
    else:
        found = re.search(r'(?:^|[\s"{,])ref"?\s*[=:]\s*"?([^"\s,}]+)', body)
        ref = found.group(1) if found else ""
    name = branch_name(ref) if ref else None
    if name is None:
        return f"Writing the ref '{ref or '?'}' through the API is not allowed: only branches."
    if name in PROTECTED:
        return ("Pushing to a protected branch (main/dev) is not allowed. "
                "Push your feature branch instead and open a pull request.")
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
    if GH_API_MERGES.search(endpoint) and method != "GET":
        return ("Merging branches through the merges endpoint is not allowed; "
                "merge a pull request instead.")
    if GH_API_REF.search(endpoint) and method in ("PATCH", "POST", "PUT"):
        return check_ref_write(endpoint, method, tokens)
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
        if tok.startswith("--repo="):
            args.append(tok)
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
    return check_branch_deletion([head], remote=True, exempt=number,
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


def check_nested(command, depth):
    """Judge a command that bash -c or eval runs, as one segment's verdict."""
    block, asks = evaluate(command, depth + 1)
    if block:
        return block
    return Ask("; ".join(dict.fromkeys(asks))) if asks else None


def check_segment(segment, depth=0):
    """Return a refusal message for this segment, an Ask, or None to allow it."""
    tokens = command_tokens(tokenise(segment))
    if not tokens:
        return None
    program = os.path.basename(tokens[0])

    # --- a shell running a string: judge the string ------------------------
    if program == "eval":
        return check_nested(" ".join(tokens[1:]), depth)
    if program in SHELLS:
        for n, tok in enumerate(tokens[1:-1], 1):
            if re.fullmatch(r"-[A-Za-z]*c[A-Za-z]*", tok):
                return check_nested(tokens[n + 1], depth)
        return None

    # --- gh ---------------------------------------------------------------
    if program == "gh":
        rest = gh_positionals(tokens)
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
        # Short options bundle (`-uf`, `-dq`), so flags are read letter by letter.
        flags = short_flags(args)
        positional = push_positionals(args)
        # A dry run mutates nothing, but only this segment is exempt.
        if "--dry-run" in args or "n" in flags:
            return None
        if "f" in flags or any(a in ("--force", "--force-with-lease") or
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
        if "--prune" in args:
            return "`git push --prune` deletes remote branches and is not allowed."
        if any("*" in a for a in positional[1:]):
            return "A wildcard refspec can reach protected branches and is not allowed."
        deleting = []
        if "--delete" in args or "d" in flags:
            if "tag" in positional[1:]:
                return "Deleting tags is not allowed. Tags are applied by the human after merge."
            deleting = positional[1:]
        deleting += [a[1:] for a in positional[1:] if a.startswith(":")]
        if deleting:
            return check_branch_deletion(deleting, remote=True)
        if targets_protected_ref(positional):
            return ("Pushing to a protected branch (main/dev) is not allowed. "
                    "Push your feature branch instead and open a pull request.")
        # No refspec: git pushes the current branch to its upstream.
        if not positional[1:]:
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
    if sub == "branch" and ("--delete" in args or {"d", "D"} & short_flags(args)):
        return check_branch_deletion([a for a in args if not a.startswith("-")],
                                     remote=False)
    if sub == "tag" and ("--delete" in args or "d" in short_flags(args)):
        return "Deleting tags is not allowed. Tags are applied by the human after merge."
    if sub == "update-ref" and "-d" in args and any("tags/" in a for a in args):
        return "Deleting tags is not allowed. Tags are applied by the human after merge."

    return None


def evaluate(command, depth=0):
    """(refusal, asks) for a whole command line, nested commands included."""
    if depth > MAX_DEPTH:
        return "Commands nested this deeply cannot be checked, so this is refused.", []
    asks = []
    for segment in split_segments(command):
        message = check_segment(segment, depth)
        if isinstance(message, Ask):
            asks.append(message.what)
        elif message:
            return message, []
    for inner in substitutions(command):
        message, more = evaluate(inner, depth + 1)
        if message:
            return message, []
        asks += more
    return None, asks


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

    message, asks = evaluate(command)
    if message:
        print(
            f"Blocked by mesa's workflow contract (AGENTS.md):\n  {message}\n\n"
            f"Command: {command.strip()[:300]}\n\n"
            "Do not work around this — report it to the human instead.",
            file=sys.stderr,
        )
        return 2
    if asks and payload.get("permission_mode") not in PROMPTING_MODES:
        mode = payload.get("permission_mode") or "unknown"
        print(
            "Blocked by mesa's workflow contract (AGENTS.md):\n  This needs the human's "
            "explicit approval to " + "; ".join(dict.fromkeys(asks)) + f", but the "
            f"'{mode}' permission mode would not show them a confirmation.\n\n"
            f"Command: {command.strip()[:300]}\n\n"
            "Do not work around this. Tell the human to run it themselves, or to "
            "switch to a permission mode that prompts and ask again.",
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
