#!/usr/bin/env python3
"""Derive mesa's project state from git, gh, NEWS.md and the gh-pages site.

Writes two files:

  .claude/state/status.json   machine-readable state (gitignored)
  STATUS.md                   the human/agent summary (gitignored: it is derived)

Nothing here is hand-maintained, so nothing can go stale: every fact is read
back from the place humans and agents both actually write to. A change made on
github.com by a co-maintainer shows up on the next run with no further work.

Every probe degrades on its own. A missing `gh`, a revoked token or a dead
network yields "unknown" fields and a `degraded` list, never a traceback and
never a half-written file (both outputs are written atomically at the end).

The roadmap is the one input that is judgement rather than fact: which items,
in which order, for which release. It lives as a TOML block in the pinned
roadmap issue (#124), editable on github.com. Every item's *status* is still
derived - from the state of the issues and PRs it names, remote branches and
release tags - and the issue's readable checklist is regenerated from it.

Usage:
  mesa-status.py                  refresh unconditionally
  mesa-status.py --max-age 14400  refresh only if STATUS.md is older than 4h
  mesa-status.py --no-log         skip the BiocCheck CI-log fetch (the slow probe)
  mesa-status.py --json-only      write status.json but not STATUS.md
  mesa-status.py --html           also build the page /mesa-status publishes
  mesa-status.py --sync-issue     rewrite #124's generated checklist if it changed
  mesa-status.py --artifact-url   print the published page's URL, if known
  mesa-status.py --mark-published clear the "page needs a refresh" flags
  mesa-status.py --park 115 <<'EOF'
  <note>
  EOF                             append a finding to the parking lot, then stop
"""

import argparse
import difflib
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from datetime import date, datetime, timedelta, timezone

REPO = "cruk-mi/mesa"
WORKFLOW = "R-CMD-check-bioc"
PROTECTED = {"main", "dev", "master", "gh-pages", "HEAD"}
COVERAGE_TARGET = 80  # AGENTS.md "Local verification"
# Only labels that exist on GitHub: a name here that nobody can apply is noise.
FIX_LABELS = ("ERROR_fix", "NOTE_fix", "Bioc_feedback_fix")
ROADMAP_ISSUE = 124

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
STATE_DIR = os.path.join(ROOT, ".claude", "state")
JSON_PATH = os.path.join(STATE_DIR, "status.json")
MD_PATH = os.path.join(ROOT, "STATUS.md")
BIOC_HISTORY = os.path.join(STATE_DIR, "bioccheck-history.jsonl")
RECOMMENDATION = os.path.join(STATE_DIR, "recommendation.md")
TEMPLATE = os.path.join(ROOT, ".claude", "scripts", "dashboard-template.html")
HTML_OUT = os.path.join(STATE_DIR, "dashboard.html")
ROADMAP_MD = os.path.join(STATE_DIR, "roadmap-issue.md")

DEGRADED = []


def shared_dir():
    """Per-clone state every worktree must see: the page URL, the refresh flags.

    `.claude/state` is per worktree, so a URL saved there by one checkout is
    invisible to the next, and /mesa-status run from a fresh worktree would
    publish a second page instead of updating the first. The git common dir
    is the one place all worktrees of a clone share, and git never commits it.
    """
    common = git("rev-parse", "--path-format=absolute", "--git-common-dir")
    path = os.path.join(common, "mesa-status") if common else STATE_DIR
    os.makedirs(path, exist_ok=True)
    return path


def degrade(message):
    """Record why a field is unknown, so STATUS.md can say so out loud."""
    if message not in DEGRADED:
        DEGRADED.append(message)


def fetch_refs():
    """Update the remote-tracking refs every probe below reads.

    Without this, `origin/main` and `origin/gh-pages` are only as fresh as the
    last time a human fetched, so a refresh could report a stale SHA and an
    "N commits behind" count that contradicts GitHub - the exact drift this
    script exists to prevent. Degrades like every other probe: on failure the
    comparisons still run, against whatever is already local.
    """
    # The configured refspec (every branch) plus --prune, not just main and
    # gh-pages: the branch report prints paste-ready `git push origin --delete`
    # lines, and those must not be built from remote-tracking refs for
    # branches that no longer exist.
    if git("fetch", "--quiet", "--prune", "origin", timeout=60) is None:
        degrade("`git fetch origin` failed, so main/gh-pages comparisons use "
                "possibly stale local refs")


def run(argv, timeout=30):
    """Run a command in the repo. Returns stdout, or None on any failure."""
    try:
        proc = subprocess.run(
            argv, cwd=ROOT, capture_output=True, text=True, timeout=timeout
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if proc.returncode != 0:
        return None
    return proc.stdout


def git(*args, timeout=30):
    out = run(["git", *args], timeout=timeout)
    return out.strip() if out is not None else None


def gh_json(args, timeout=60):
    """Run a `gh ... --json` command and parse it. None on any failure."""
    out = run(["gh", *args], timeout=timeout)
    if out is None:
        degrade("gh " + " ".join(args[:3]) + " failed (not installed, not authenticated, or offline)")
        return None
    try:
        return json.loads(out)
    except (json.JSONDecodeError, ValueError):
        degrade("gh " + " ".join(args[:3]) + " returned output that is not JSON")
        return None


def now():
    return datetime.now(timezone.utc)


def iso(dt):
    return dt.strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_ts(text):
    """Parse a GitHub ISO-8601 timestamp. None if unparseable."""
    if not text:
        return None
    try:
        return datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return None


def day(text):
    ts = parse_ts(text)
    return ts.strftime("%Y-%m-%d") if ts else "?"


def age_days(text):
    ts = parse_ts(text)
    return None if ts is None else (now() - ts).days


# --------------------------------------------------------------------------
# probes
# --------------------------------------------------------------------------

def collect_release():
    """Version from DESCRIPTION, and which phase of the cycle that means.

    Per the bioc-release-cycle skill: a four-part x.y.z.9000 version means a
    devel section is open and work PRs land against it; a three-part version
    means the cycle is closed and the next step is opening the next devel.
    """
    version, phase = None, "unknown"
    try:
        with open(os.path.join(ROOT, "DESCRIPTION"), encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("Version:"):
                    version = line.split(":", 1)[1].strip()
                    break
    except OSError:
        degrade("DESCRIPTION could not be read, so the version is unknown")
    if version:
        parts = version.split(".")
        if len(parts) == 4 and parts[3] == "9000":
            phase = "devel open - work PRs land here"
        elif len(parts) == 3:
            phase = "released - next step is opening the next devel section"
        else:
            phase = "unrecognised version shape"
    return {"version": version, "phase": phase}


def collect_news():
    """The first NEWS.md section, verbatim. That is 'what landed recently'."""
    try:
        with open(os.path.join(ROOT, "NEWS.md"), encoding="utf-8") as handle:
            lines = handle.read().splitlines()
    except OSError:
        degrade("NEWS.md could not be read")
        return {"heading": None, "body": []}
    start = next((i for i, line in enumerate(lines) if line.startswith("# mesa ")), None)
    if start is None:
        degrade("NEWS.md has no '# mesa X.Y.Z' heading")
        return {"heading": None, "body": []}
    end = next(
        (i for i in range(start + 1, len(lines)) if lines[i].startswith("# mesa ")),
        len(lines),
    )
    body = lines[start + 1:end]
    while body and not body[0].strip():
        body.pop(0)
    while body and not body[-1].strip():
        body.pop()
    return {"heading": lines[start].lstrip("# ").strip(), "body": body}


def collect_worktree():
    """Where this checkout actually is, relative to origin/main."""
    branch = git("rev-parse", "--abbrev-ref", "HEAD")
    porcelain = git("status", "--porcelain")
    dirty = [line for line in (porcelain or "").splitlines() if line.strip()]
    ahead = behind = None
    counts = git("rev-list", "--left-right", "--count", "origin/main...HEAD")
    if counts:
        try:
            behind, ahead = (int(part) for part in counts.split())
        except ValueError:
            pass
    recent = (git("log", "--oneline", "-5") or "").splitlines()
    return {
        "branch": branch,
        "dirty": dirty,
        "ahead_of_main": ahead,
        "behind_main": behind,
        "recent_commits": recent,
    }


def _rollup(rollup):
    """Reduce a PR's statusCheckRollup to one word."""
    if not rollup:
        return "none"
    states = []
    for check in rollup:
        state = check.get("conclusion") or check.get("state") or check.get("status")
        if state:
            states.append(str(state).lower())
    if not states:
        return "none"
    if any(s in ("failure", "timed_out", "cancelled", "action_required", "error") for s in states):
        return "failing"
    if any(s in ("pending", "in_progress", "queued", "waiting", "expected") for s in states):
        return "running"
    if all(s in ("success", "neutral", "skipped") for s in states):
        return "green"
    return "mixed"


def collect_in_flight():
    """Open PRs, with their CI verdict pulled from the same call."""
    data = gh_json([
        "pr", "list", "--repo", REPO, "--state", "open", "--limit", "50",
        "--json", "number,title,headRefName,isDraft,updatedAt,statusCheckRollup",
    ])
    if data is None:
        return None
    return [
        {
            "number": pr.get("number"),
            "title": pr.get("title"),
            "branch": pr.get("headRefName"),
            "draft": bool(pr.get("isDraft")),
            "updated": day(pr.get("updatedAt")),
            "ci": _rollup(pr.get("statusCheckRollup")),
        }
        for pr in data
    ]


def collect_landed():
    data = gh_json([
        "pr", "list", "--repo", REPO, "--state", "merged", "--limit", "10",
        "--json", "number,title,mergedAt",
    ])
    if data is None:
        return None
    return [
        {"number": pr.get("number"), "title": pr.get("title"), "merged": day(pr.get("mergedAt"))}
        for pr in data
    ]


def collect_next_up():
    """Open issues, grouped on the repo's existing labels."""
    data = gh_json([
        "issue", "list", "--repo", REPO, "--state", "open", "--limit", "100",
        "--json", "number,title,labels,updatedAt",
    ])
    if data is None:
        return None
    groups = {}
    for issue in data:
        # A `status:` label says how far along an issue is, not what it is
        # about, so it must never become the group an issue is filed under.
        names = [label.get("name") for label in issue.get("labels") or []
                 if not str(label.get("name", "")).startswith("status:")]
        key = next((name for name in FIX_LABELS if name in names), None)
        if key is None:
            key = names[0] if names else "unlabelled"
        groups.setdefault(key, []).append({
            "number": issue.get("number"),
            "title": issue.get("title"),
            "updated": day(issue.get("updatedAt")),
        })
    for issues in groups.values():
        issues.sort(key=lambda i: i["number"] or 0)
    return groups


def latest_main_run(successful=False):
    """The newest main run, or the newest that actually finished green.

    CI wants the newest run whatever its verdict - that is the current state.
    BiocCheck wants the newest *successful* one: a failed or in-progress run
    may never have reached the BiocCheck step, and parsing its log would report
    "not parsed" while pointing at a run that was never checked.
    """
    argv = [
        "run", "list", "--repo", REPO, "--workflow", WORKFLOW, "--branch", "main",
        "--limit", "1", "--json", "databaseId,conclusion,headSha,createdAt,url",
    ]
    if successful:
        argv += ["--status", "success"]
    data = gh_json(argv)
    if not data:
        return None
    return data[0]


def collect_ci(run_info):
    """CI on main, plus whether main has moved past the last run.

    check-bioc.yml only triggers on pushes touching R/, tests/, vignettes/,
    inst/, DESCRIPTION or NAMESPACE, so a tooling-only commit legitimately
    leaves main un-run. That is worth saying rather than reading as a failure.
    """
    if run_info is None:
        return {"conclusion": "unknown", "sha": None, "when": None, "url": None,
                "main_sha": git("rev-parse", "--short", "origin/main"), "behind": None}
    main_sha = git("rev-parse", "--short", "origin/main")
    run_sha = (run_info.get("headSha") or "")[:7]
    behind = None
    if main_sha and run_sha:
        count = git("rev-list", "--count", f"{run_sha}..origin/main")
        try:
            behind = int(count) if count is not None else None
        except ValueError:
            behind = None
    return {
        "conclusion": run_info.get("conclusion") or "unknown",
        "sha": run_sha,
        "when": day(run_info.get("createdAt")),
        "url": run_info.get("url"),
        "main_sha": main_sha,
        "behind": behind,
    }


def collect_coverage():
    """Coverage from Codecov's public API (no token needed for a public repo)."""
    url = "https://api.codecov.io/api/v2/github/cruk-mi/repos/mesa/"
    try:
        with urllib.request.urlopen(url, timeout=20) as response:
            payload = json.load(response)
    except (urllib.error.URLError, OSError, ValueError, json.JSONDecodeError) as exc:
        degrade(f"Codecov API unreachable ({type(exc).__name__}), so coverage is unknown")
        return {"pct": None, "updated": None, "stale_days": None, "target": COVERAGE_TARGET}
    totals = payload.get("totals") or {}
    updated = payload.get("updatestamp")
    return {
        "pct": totals.get("coverage"),
        "updated": day(updated),
        "stale_days": age_days(updated),
        "target": COVERAGE_TARGET,
    }


BIOC_SUMMARY = re.compile(
    r"(\d+)\s+ERRORS?\s*\|\s*\D*(\d+)\s+WARNINGS?\s*\|\s*\D*(\d+)\s+NOTES?", re.I
)
# `gh run view --log` prefixes every line with "<job>\t<step>\t<timestamp> ".
# The step column is not usable here - gh reports "UNKNOWN STEP" for every line
# of this workflow's logs - so the BiocCheck section is delimited by BiocCheck's
# own banners instead. Without that bound, an "ERROR: " from dependency install
# or R CMD check would be collected and rendered as a BiocCheck finding.
LOG_PREFIX = re.compile(r"^.*?\d{4}-\d{2}-\d{2}T[\d:.]+Z\s?")
BIOC_START = re.compile(r"Running BiocCheck on ")
BIOC_END = re.compile(r"BiocCheck v[\d.]+ results")
# BiocCheck prefixes its own findings with a glyph; a bare "ERROR: " anywhere
# else in the log is some other tool's.
BIOC_ERROR = re.compile(r"^[*\u2716]\s*ERROR:\s*(.+)$")
# A BiocCheck finding starts with one of these; anything else continues the
# previous one, because long messages wrap onto their own lines.
FINDING_START = re.compile(r"^[*\u2716\u2139\u26a0\u2714]")


def log_message(line):
    """Strip gh's job/step/timestamp prefix, leaving the tool's own output."""
    return LOG_PREFIX.sub("", line).strip()


def read_bioc_history():
    entries = []
    try:
        with open(BIOC_HISTORY, encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                try:
                    entries.append(json.loads(line))
                except (json.JSONDecodeError, ValueError):
                    continue
    except OSError:
        pass
    return entries


def collect_bioccheck(run_info, allow_log):
    """BiocCheck ERROR/WARNING/NOTE counts for the latest checked main commit.

    These cannot be recomputed locally (see the bioc-check-ladder skill), so
    they are parsed out of CI. Downloading a run log is the slowest probe here,
    so it only happens when the run has not already been recorded: main pushes
    are rare, making this a cache hit almost every time.
    """
    history = read_bioc_history()
    latest = history[-1] if history else None
    if run_info is None:
        return {"latest": latest, "history": history, "source": "cache"}
    sha = (run_info.get("headSha") or "")[:7]
    # Keyed on the sha anywhere in history, not on the last row: a re-run of an
    # older commit would otherwise append it after newer rows and render a
    # regression that never happened, durably, since the file is append-only.
    if any(entry.get("sha") == sha for entry in history):
        return {"latest": latest, "history": history, "source": "cache"}
    if not allow_log:
        return {"latest": latest, "history": history, "source": "cache (log fetch skipped)"}

    log = run(["gh", "run", "view", str(run_info.get("databaseId")),
               "--repo", REPO, "--log"], timeout=180)
    if log is None:
        degrade("CI log for the latest main run could not be fetched, so BiocCheck counts are the cached ones")
        return {"latest": latest, "history": history, "source": "cache"}

    counts, errors, pending, inside = None, [], None, False
    for line in log.splitlines():
        message = log_message(line)
        if not inside:
            if BIOC_START.search(message):
                inside = True
            continue
        match = BIOC_SUMMARY.search(message)
        if match:
            counts = {
                "error": int(match.group(1)),
                "warning": int(match.group(2)),
                "note": int(match.group(3)),
            }
            break
        if BIOC_END.search(message):
            # The banner precedes the summary; keep reading for one more line.
            pending = None
            continue
        found = BIOC_ERROR.match(message)
        if found:
            pending = found.group(1).strip()
            errors.append(pending)
        elif pending is not None:
            if message and not FINDING_START.match(message):
                errors[-1] = (errors[-1].rstrip(", ") + " " + message).strip()
            else:
                pending = None
    errors = list(dict.fromkeys(error for error in errors if error))
    if not inside:
        degrade("BiocCheck output not found in the CI log for this run")
        return {"latest": latest, "history": history, "source": "not parsed"}
    if counts is None:
        degrade("BiocCheck summary line not found in the CI log (output format may have changed)")
        return {"latest": latest, "history": history, "source": "not parsed"}

    entry = {
        "sha": sha,
        "date": day(run_info.get("createdAt")),
        **counts,
        "errors": errors[:5],
    }
    history.append(entry)
    try:
        os.makedirs(STATE_DIR, exist_ok=True)
        with open(BIOC_HISTORY, "a", encoding="utf-8") as handle:
            handle.write(json.dumps(entry, sort_keys=True) + "\n")
    except OSError:
        degrade("BiocCheck history file could not be appended to")
    return {"latest": entry, "history": history, "source": "CI log"}


def collect_recommendation():
    """The one thing a script cannot derive: what to actually do next.

    /mesa-status writes this file, then re-runs the generator, so the judgement
    is rendered into STATUS.md rather than pasted on top of it and lost at the
    next refresh. It is committed, so regenerating anywhere reproduces the same
    file byte for byte.
    """
    try:
        with open(RECOMMENDATION, encoding="utf-8") as handle:
            text = handle.read().strip()
    except OSError:
        return None
    return text or None


def parking_lot_path():
    """The parking lot lives in the main checkout, never in a worktree.

    `/park` is used from per-issue worktrees, and each worktree has its own
    `.claude/state/`. Resolving through the common git dir sends every session
    to the same gitignored file, so nothing is split across worktrees.
    """
    common = git("rev-parse", "--path-format=absolute", "--git-common-dir")
    if not common:
        return os.path.join(STATE_DIR, "parking-lot.md")
    return os.path.join(os.path.dirname(common.rstrip("/")), ".claude", "state", "parking-lot.md")


def collect_parked():
    """Out-of-scope findings waiting for triage (AGENTS.md "Session scope").

    One `- ` line per item, written by /park. The file is gitignored, so this is
    per-machine: it shows what this checkout has parked, not the whole team's.
    """
    try:
        with open(parking_lot_path(), encoding="utf-8") as handle:
            lines = handle.read().splitlines()
    except OSError:
        return []
    return [line[2:].strip() for line in lines if line.startswith("- ")]


def park(issue, note):
    """Append one finding to the parking lot. Returns the line written, or None.

    The note arrives on stdin, not as an argument, so it never passes through
    the shell: `$`, backticks and quotes in R code survive as written. Newlines
    are collapsed so each finding stays one `- ` line, which is all
    collect_parked() reads.
    """
    note = " ".join(note.split())
    if not note:
        return None
    issue = issue.strip()
    if issue not in ("", "-") and not issue.startswith("#"):
        issue = "#" + issue
    branch = git("branch", "--show-current") or "-"
    today = datetime.now(timezone.utc).astimezone().strftime("%Y-%m-%d")
    line = f"- {today} ({branch}, {issue or '-'}) {note}"
    path = parking_lot_path()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a", encoding="utf-8") as handle:
        handle.write(line + "\n")
    return line


def collect_pkgdown():
    """Is the published site built from current main?

    pkgdown::deploy_to_branch() writes 'Built site for mesa@<version>: <sha>'.
    """
    subject = git("log", "-1", "--format=%s", "origin/gh-pages")
    main_sha = git("rev-parse", "--short", "origin/main")
    site_sha = None
    if subject:
        match = re.search(r":\s*([0-9a-f]{7,40})\s*$", subject)
        if match:
            site_sha = match.group(1)[:7]
    else:
        degrade("origin/gh-pages not found locally; run `git fetch origin` for site freshness")
    drift = None
    if site_sha and main_sha:
        count = git("rev-list", "--count", f"{site_sha}..origin/main")
        try:
            drift = int(count) if count is not None else None
        except ValueError:
            drift = None
    return {
        "site_sha": site_sha,
        "main_sha": main_sha,
        "commits_behind": drift,
        "has_config": os.path.exists(os.path.join(ROOT, "_pkgdown.yml")),
    }


def collect_branches():
    """Classify every branch against PR history.

    mesa squash-merges, so `git branch --merged` cannot see that a branch
    landed - its commits never appear on main. Only PR state knows.
    """
    local = [b for b in (git("for-each-ref", "--format=%(refname:short)", "refs/heads") or "").splitlines() if b]
    remote = [
        b.split("/", 1)[1]
        for b in (git("for-each-ref", "--format=%(refname:short)", "refs/remotes/origin") or "").splitlines()
        if b and "/" in b
    ]
    prs = gh_json([
        "pr", "list", "--repo", REPO, "--state", "all", "--limit", "300",
        "--json", "number,headRefName,headRefOid,state,isCrossRepository",
    ])
    if prs is None:
        return None
    # An OPEN PR outranks everything: a branch backing live work must never be
    # classified as landed, because that is what fills the prune block below.
    RANK = {"OPEN": 3, "MERGED": 2, "CLOSED": 1}
    state_of = {}
    for pr in prs:
        ref = pr.get("headRefName")
        if not ref:
            continue
        # headRefName is only the short name, so a fork's branch can collide
        # with an unrelated ref of ours. The branches listed above are all this
        # repository's, so a fork PR can say nothing about them.
        if pr.get("isCrossRepository"):
            continue
        rank = RANK.get(pr.get("state"), 0)
        held = state_of.get(ref)
        held_rank = RANK.get(held["state"], 0) if held else 0
        if rank > held_rank:
            state_of[ref] = {"state": pr.get("state"), "number": pr.get("number"),
                             "oids": {pr.get("headRefOid")}}
        elif rank == held_rank:
            held["oids"].add(pr.get("headRefOid"))

    current = git("rev-parse", "--abbrev-ref", "HEAD")
    buckets = {"landed": [], "abandoned": [], "open": [], "unknown": []}
    for name in sorted(set(local) | set(remote)):
        # Protected branches, and the branch in hand, are never prune candidates.
        if name in PROTECTED or name == current:
            continue
        info = state_of.get(name)
        where = "local+remote" if name in local and name in remote else ("local" if name in local else "remote")
        row = {"branch": name, "where": where, "pr": info.get("number") if info else None}
        # A finished PR only speaks for the commit it ended on. A reused name
        # that has moved on since is live work, so it is not a prune candidate.
        tips = [git("rev-parse", "--verify", "-q", ref) for ref in
                ([f"refs/heads/{name}"] if name in local else []) +
                ([f"refs/remotes/origin/{name}"] if name in remote else [])]
        if info is None or (info["state"] != "OPEN" and
                            not all(t and t in info["oids"] for t in tips)):
            buckets["unknown"].append(row)
        elif info["state"] == "MERGED":
            buckets["landed"].append(row)
        elif info["state"] == "OPEN":
            buckets["open"].append(row)
        else:
            buckets["abandoned"].append(row)
    return buckets


# --------------------------------------------------------------------------
# roadmap
# --------------------------------------------------------------------------

GEN_START, GEN_END = "<!-- roadmap:generated:start -->", "<!-- roadmap:generated:end -->"
DATA_START, DATA_END = "<!-- roadmap:data:start -->", "<!-- roadmap:data:end -->"
TOML_FENCE = re.compile(r"```toml[ \t]*\n(.*?)\n```", re.S)

# The five statuses, in the order the page lists them. Each names where an
# item *is*, never an action: "delete" could mean done or to-do, "In review"
# cannot. Who has to act is the separate `owner` field.
STATUSES = ("todo", "waiting", "review", "done", "later")
STATUS_LABEL = {"todo": "To do", "waiting": "Waiting", "review": "In review",
                "done": "Done", "later": "Later"}
STATUS_MARK = {"todo": "⬜", "waiting": "⏳", "review": "\U0001f50d",
               "done": "✅", "later": "\U0001f4a4"}
ROLLUP_WORD = {"SUCCESS": "green", "FAILURE": "failing", "ERROR": "failing",
               "PENDING": "running", "EXPECTED": "running"}


def _iso_date(value):
    return value.isoformat() if isinstance(value, (date, datetime)) else value


def _fmt_day(text):
    """'2026-09-25T14:24:16Z' -> '25 Sep', for the one-line ref summaries."""
    ts = parse_ts(text)
    return f"{ts.day} {ts.strftime('%b')}" if ts else "?"


def fetch_issue_body():
    data = gh_json(["issue", "view", str(ROADMAP_ISSUE), "--repo", REPO, "--json", "body"])
    return None if data is None else (data.get("body") or "")


def read_roadmap_source(override=None):
    """(issue body, TOML text) - from a local file when testing, else #124."""
    body = fetch_issue_body()
    if override:
        try:
            with open(override, encoding="utf-8") as handle:
                return body, handle.read()
        except OSError:
            degrade(f"roadmap file {override} could not be read")
            return body, None
    if body is None:
        return None, None
    section = body
    if DATA_START in body and DATA_END in body:
        section = body.split(DATA_START, 1)[1].split(DATA_END, 1)[0]
    match = TOML_FENCE.search(section)
    if not match:
        degrade(f"#{ROADMAP_ISSUE} has no ```toml roadmap block, so the roadmap is unknown")
        return body, None
    return body, match.group(1)


def parse_roadmap(text):
    """Parse and check the roadmap TOML. None, with a reason, if it is unusable.

    Maintainers edit this on github.com, so a typo is expected, not
    exceptional: it must cost the roadmap section, never the whole refresh.
    """
    try:
        import tomllib  # 3.11+; macOS's /usr/bin/python3 is 3.9
    except ImportError:
        degrade("the roadmap needs Python >= 3.11 (tomllib); this is Python "
                + ".".join(map(str, sys.version_info[:2])))
        return None
    try:
        data = tomllib.loads(text)
    except tomllib.TOMLDecodeError as exc:
        degrade(f"roadmap TOML in #{ROADMAP_ISSUE} does not parse ({exc})")
        return None
    waves = data.get("wave") or []
    items = data.get("item") or []
    wave_ids = [w.get("id") for w in waves]
    problems = []
    keys = set()
    for index, item in enumerate(items):
        refs = item.get("refs") or []
        # A key only needs writing when two items would otherwise share one,
        # or when the item has no ref to borrow it from.
        item.setdefault("key", str(refs[0]) if refs else f"item{index + 1}")
        item["key"] = str(item["key"])
        if item["key"] in keys:
            problems.append(f"duplicate item key {item['key']!r}")
        keys.add(item["key"])
        if item.get("wave") not in wave_ids:
            problems.append(f"item {item['key']!r} names unknown wave {item.get('wave')!r}")
        if not item.get("title"):
            problems.append(f"item {item['key']!r} has no title")
        if not all(isinstance(ref, int) for ref in refs):
            problems.append(f"item {item['key']!r} has a non-numeric ref")
    for item in items:
        for dep in item.get("after") or []:
            if str(dep) not in keys:
                problems.append(f"item {item['key']!r} waits on unknown key {dep!r}")
    if problems:
        degrade("roadmap TOML is inconsistent: " + "; ".join(problems[:4]))
        return None
    return data


def fetch_ref_states(numbers):
    """State of every issue/PR the roadmap names, in one GraphQL round trip."""
    if not numbers:
        return {}
    owner, name = REPO.split("/")
    rollup = "commits(last: 1) { nodes { commit { statusCheckRollup { state } } } }"
    fields = []
    for number in sorted(numbers):
        fields.append(
            f"n{number}: issueOrPullRequest(number: {number}) {{ __typename "
            f"... on Issue {{ title state closedAt url "
            f"closedByPullRequestsReferences(first: 5, includeClosedPrs: false) "
            f"{{ nodes {{ number state isDraft url {rollup} }} }} }} "
            f"... on PullRequest {{ title state isDraft mergedAt closedAt url {rollup} }} }}"
        )
    query = f'query {{ repository(owner: "{owner}", name: "{name}") {{ {" ".join(fields)} }} }}'
    data = gh_json(["api", "graphql", "-f", f"query={query}"])
    if data is None:
        return None
    repo = (data.get("data") or {}).get("repository") or {}
    return {int(key[1:]): value for key, value in repo.items() if value}


def _ci_word(node):
    commits = ((node.get("commits") or {}).get("nodes")) or []
    state = ((commits[0].get("commit") or {}).get("statusCheckRollup") or {}).get("state") if commits else None
    return ROLLUP_WORD.get(state, "no CI")


def remote_refs():
    """(branch names, tag names) on origin right now - not the local cache."""
    out = git("ls-remote", "--heads", "--tags", "origin", timeout=60)
    if out is None:
        degrade("`git ls-remote origin` failed, so branch and tag checks are unknown")
        return None, None
    heads, tags = set(), set()
    for line in out.splitlines():
        ref = line.split("\t", 1)[-1]
        if ref.startswith("refs/heads/"):
            heads.add(ref[len("refs/heads/"):])
        elif ref.startswith("refs/tags/"):
            tags.add(ref[len("refs/tags/"):].removesuffix("^{}"))
    return heads, tags


def describe_ref(number, info):
    """One ref as the page shows it: resolved or not, plus a short summary."""
    url = f"https://github.com/{REPO}/issues/{number}"
    if info is None:
        return {"number": number, "kind": "?", "resolved": False, "open_prs": [],
                "text": "state unknown", "url": url}
    if info.get("__typename") == "PullRequest":
        state = info.get("state")
        if state == "MERGED":
            text = f"merged {_fmt_day(info.get('mergedAt'))}"
        elif state == "CLOSED":
            text = "closed without merging"
        else:
            text = ("draft" if info.get("isDraft") else "ready") + " · CI " + _ci_word(info)
        return {"number": number, "kind": "pr", "resolved": state != "OPEN",
                "open_prs": [number] if state == "OPEN" else [],
                "text": text, "url": info.get("url") or url, "title": info.get("title")}
    if info.get("state") == "CLOSED":
        return {"number": number, "kind": "issue", "resolved": True, "open_prs": [],
                "text": f"closed {_fmt_day(info.get('closedAt'))}",
                "url": info.get("url") or url, "title": info.get("title")}
    # An open issue with an open PR that says "Fixes #N" is already in review,
    # whether or not the roadmap lists the PR: GitHub knows the link.
    linked = [pr for pr in ((info.get("closedByPullRequestsReferences") or {}).get("nodes") or [])
              if pr.get("state") == "OPEN"]
    text = "open"
    if linked:
        text += " → " + ", ".join(
            f"PR #{pr['number']} ({'draft' if pr.get('isDraft') else 'ready'}, CI {_ci_word(pr)})"
            for pr in linked)
    return {"number": number, "kind": "issue", "resolved": False,
            "open_prs": [pr["number"] for pr in linked],
            "text": text, "url": info.get("url") or url, "title": info.get("title")}


def run_check(check, heads, tags):
    """A `kind:value` check the page can verify. (passed or None, label)."""
    kind, _, value = str(check).partition(":")
    if kind == "branch-gone":
        return (None if heads is None else value not in heads), f"branch {value} deleted"
    if kind == "tag":
        return (None if tags is None else value in tags), f"tag {value} pushed"
    degrade(f"roadmap check {check!r} is not a known kind (branch-gone:, tag:)")
    return None, str(check)


def derive_status(item, wave, refs, checks, done_keys, key_label):
    """(status, reason) for one item. The table in the mesa-status skill.

    Order matters: finished beats everything, "Later" beats open work, an
    explicit `waiting_on` beats an open PR (e.g. a PR waiting on re-review),
    and an open PR beats a pending dependency (someone is already on it).
    """
    signals = bool(refs or checks)
    finished = (all(r["resolved"] for r in refs)
                and all(passed is True for passed, _ in checks))
    if item.get("done") or (signals and finished):
        return "done", ""
    if wave.get("later"):
        return "later", ""
    if item.get("waiting_on"):
        return "waiting", "on " + item["waiting_on"]
    open_prs = [n for r in refs for n in r["open_prs"]]
    if open_prs:
        return "review", ""
    pending = [str(dep) for dep in item.get("after") or [] if str(dep) not in done_keys]
    if pending:
        return "waiting", "on " + ", ".join(key_label[dep] for dep in pending)
    return "todo", ""


def collect_roadmap(override=None):
    body, text = read_roadmap_source(override)
    if text is None:
        return None, body
    data = parse_roadmap(text)
    if data is None:
        return None, body
    items = data.get("item") or []
    numbers = {ref for item in items for ref in item.get("refs") or []}
    states = fetch_ref_states(numbers)
    if states is None:
        degrade("issue/PR states for the roadmap could not be fetched")
        states = {}
    needs_remote = any(item.get("check") for item in items)
    heads, tags = remote_refs() if needs_remote or data.get("release") else (set(), set())

    waves = {w["id"]: w for w in data.get("wave") or []}
    key_label = {item["key"]: item.get("short") or
                 (f"#{item['key']}" if item["key"].isdigit() else item["title"])
                 for item in items}
    # Dependencies point backwards through the file (a plan lists what comes
    # first, first), so one ordered pass sees each dependency settled before
    # the items waiting on it.
    done_keys, out_items = set(), []
    for position, item in enumerate(items, 1):
        wave = waves[item["wave"]]
        refs = [describe_ref(n, states.get(n)) for n in item.get("refs") or []]
        checks = [run_check(c, heads, tags) for c in item.get("check") or []]
        status, why = derive_status(item, wave, refs, checks, done_keys, key_label)
        if status == "done":
            done_keys.add(item["key"])
        out_items.append({
            "key": item["key"],
            "order": position,
            "wave": item["wave"],
            "title": item["title"],
            "note": item.get("note"),
            "release": item.get("release"),
            "owner": item.get("owner", "claude"),
            "status": status,
            "why": why,
            "refs": refs,
            "checks": [{"passed": p, "label": label} for p, label in checks],
        })

    releases = []
    for rel in data.get("release") or []:
        version = str(rel.get("version"))
        # The item that cuts the release (its tag check) is the release, not
        # something it ships.
        shipped = [i for i in out_items if i["release"] == version
                   and not any(c["label"] == f"tag v{version} pushed" for c in i["checks"])]
        releases.append({
            "version": version,
            "date": _iso_date(rel.get("date")),
            "approx": bool(rel.get("approx")),
            "note": rel.get("note"),
            "released": (None if tags is None else f"v{version}" in tags),
            "items": [{"key": i["key"], "title": i["title"], "status": i["status"]} for i in shipped],
        })

    counts = {s: sum(1 for i in out_items if i["status"] == s) for s in STATUSES}
    return {
        "issue": ROADMAP_ISSUE,
        "title": data.get("title") or "Roadmap",
        "target": data.get("target"),
        "deadline": _iso_date(data.get("deadline")),
        "deadline_note": data.get("deadline_note"),
        "waves": [{"id": w["id"], "name": w.get("name"), "start": _iso_date(w.get("start")),
                   "end": _iso_date(w.get("end")), "later": bool(w.get("later")),
                   "note": w.get("note")}
                  for w in data.get("wave") or []],
        "items": out_items,
        "releases": releases,
        "counts": counts,
        "next": [i["key"] for i in out_items if i["status"] == "todo"][:3],
        "toml": text,
    }, body


def render_roadmap_issue(roadmap, page_url=None):
    """The readable half of #124. Regenerated, so it carries no checkboxes:
    a box ticked by hand would be overwritten on the next refresh."""
    out = [GEN_START, f"# {roadmap['title']}", ""]
    lead = []
    if page_url:
        lead.append(f"**Readable view:** {page_url}")
    if roadmap.get("target"):
        lead.append(f"**Target:** {roadmap['target']}")
    if roadmap.get("deadline"):
        lead.append(f"**Acceptance deadline ≈ {roadmap['deadline']}**"
                    + (f" ({roadmap['deadline_note']})" if roadmap.get("deadline_note") else ""))
    out += [" · ".join(lead), ""] if lead else []
    out += [
        "> **Generated - do not edit this part.** Each status is derived from GitHub "
        "(issue/PR state, branches, tags) by `.claude/scripts/mesa-status.py`. To change "
        "the plan (add, drop, reorder, re-note an item), edit the **Roadmap data** block "
        "below, then run `/mesa-status`.",
        "",
        "Status: " + " · ".join(f"{STATUS_MARK[s]} {STATUS_LABEL[s]}" for s in STATUSES)
        + " · \U0001f464 = a maintainer has to act",
        "",
    ]
    for wave in roadmap["waves"]:
        span = f" ({wave['start']} → {wave['end']})" if wave.get("start") else ""
        out += [f"## {wave['name']}{span}", ""]
        for item in (i for i in roadmap["items"] if i["wave"] == wave["id"]):
            refs = ", ".join(f"#{r['number']}" for r in item["refs"])
            line = f"- {STATUS_MARK[item['status']]} **{item['title']}**"
            if refs:
                line += f" ({refs})"
            if item["owner"] == "maintainer":
                line += " \U0001f464"
            if item["why"]:
                line += f" - _waiting {item['why']}_"
            out.append(line)
        out.append("")
    if roadmap["releases"]:
        out += ["## Releases", "", "| Version | When | Ships |", "|---|---|---|"]
        for rel in roadmap["releases"]:
            when = ("released" if rel["released"] else
                    (("~" if rel["approx"] else "") + (rel["date"] or "as needed")))
            ships = ", ".join(f"{STATUS_MARK[i['status']]} {i['title']}" for i in rel["items"])
            out.append(f"| {rel['version']} | {when} | {ships or rel.get('note') or ''} |")
        out.append("")
    out += [f"_Generated by Claude (AI) via `/mesa-status`; statuses as of {iso(now())}._", GEN_END]
    return "\n".join(out)


def compose_issue_body(old_body, generated, toml_text):
    """Swap the generated half into #124 and leave everything else as it is.

    The first sync (no markers yet) replaces the hand-written checklist, and
    carries over only its <details> blocks - the review summary.
    """
    data = None if toml_text is None else (
        f"{DATA_START}\n<details><summary><b>Roadmap data</b> - edit this to change "
        f"the plan, then run <code>/mesa-status</code></summary>\n\n```toml\n"
        f"{toml_text.strip()}\n```\n\n</details>\n{DATA_END}")
    if old_body and GEN_START in old_body and GEN_END in old_body:
        head, rest = old_body.split(GEN_START, 1)
        tail = rest.split(GEN_END, 1)[1]
        body = head + generated + tail
        if data is None:
            return body
        if DATA_START not in body:
            body = body.rstrip() + "\n\n" + data + "\n"
        else:
            before, after = body.split(DATA_START, 1)
            body = before + data + after.split(DATA_END, 1)[1]
        return body
    kept = re.findall(r"<details>.*?</details>", old_body or "", re.S)
    return "\n\n".join([generated, data] + kept) + "\n"


def sync_issue(roadmap, old_body, page_url, dry_run, seed_toml):
    """Rewrite #124 when its generated half is out of date. Returns a verdict."""
    if roadmap is None:
        return "roadmap unknown - #124 left alone"
    toml_text = seed_toml if seed_toml is not None else (
        None if old_body and DATA_START in old_body else roadmap["toml"])
    generated = render_roadmap_issue(roadmap, page_url)
    new_body = compose_issue_body(old_body, generated, toml_text)

    def comparable(text):  # the timestamp alone is not a change
        return re.sub(r"statuses as of \S+", "", text or "").strip()

    if comparable(new_body) == comparable(old_body):
        return f"#{ROADMAP_ISSUE} already current"
    if dry_run:
        sys.stdout.writelines(difflib.unified_diff(
            (old_body or "").splitlines(True), new_body.splitlines(True),
            f"#{ROADMAP_ISSUE} (now)", f"#{ROADMAP_ISSUE} (after sync)"))
        return f"#{ROADMAP_ISSUE} would change (dry run)"
    write_atomic(ROADMAP_MD, new_body)
    if run(["gh", "issue", "edit", str(ROADMAP_ISSUE), "--repo", REPO,
            "--body-file", ROADMAP_MD], timeout=60) is None:
        degrade(f"#{ROADMAP_ISSUE} could not be updated (gh issue edit failed)")
        return f"#{ROADMAP_ISSUE} update FAILED"
    return f"#{ROADMAP_ISSUE} updated"


# --------------------------------------------------------------------------
# render
# --------------------------------------------------------------------------

BANNER = "<!-- GENERATED by .claude/scripts/mesa-status.py - do not edit by hand. Refresh with /mesa-status. -->"


def render(state):
    out = [BANNER, "", "# mesa project status", ""]

    release, ci, cov = state["release"], state["ci"], state["coverage"]
    bioc, pkg = state["bioccheck"], state["pkgdown"]

    cov_text = "unknown"
    if cov["pct"] is not None:
        cov_text = f"{cov['pct']}% (target {cov['target']}%)"
        if cov["stale_days"] and cov["stale_days"] > 30:
            cov_text += f", last upload {cov['updated']} - {cov['stale_days']}d ago"
    ci_text = "unknown"
    if ci["sha"]:
        ci_text = f"{ci['conclusion']} at `{ci['sha']}` ({ci['when']})"
        # None is "could not compare", not "zero behind": same three-way split
        # as the pkgdown tile, so an unknown distance is never read as current.
        if ci["behind"] is None:
            ci_text += " - how far main has moved past this could not be determined"
        elif ci["behind"]:
            # check-bioc.yml only fires on pushes touching R/, tests/, vignettes/,
            # inst/, DESCRIPTION or NAMESPACE, so a docs- or tooling-only commit
            # leaves main un-run. That is expected, not a failure.
            ci_text += f" - main is {ci['behind']} commit(s) past this"

    latest_bioc = bioc.get("latest")
    bioc_text = "unknown"
    if latest_bioc:
        bioc_text = (f"{latest_bioc.get('error')} ERROR / {latest_bioc.get('warning')} WARNING"
                     f" / {latest_bioc.get('note')} NOTE at {latest_bioc.get('sha')}")

    out += [
        "| | |",
        "|---|---|",
        f"| Version | `{release['version']}` |",
        f"| Phase | {release['phase']} |",
        f"| CI on main | {ci_text} |",
        f"| Coverage | {cov_text} |",
        f"| BiocCheck | {bioc_text} |",
        "",
    ]

    if state["recommendation"]:
        out += ["## Recommended next step", "", state["recommendation"], ""]

    # --- roadmap ----------------------------------------------------------
    out += [f"## Roadmap (#{ROADMAP_ISSUE})", ""]
    roadmap = state.get("roadmap")
    if roadmap is None:
        out += ["_unknown - see Incomplete data._", ""]
    else:
        out.append(" · ".join(f"{roadmap['counts'][s]} {STATUS_LABEL[s]}" for s in STATUSES))
        out.append("")
        # Everything still open, in plan order: the list to work down.
        for item in roadmap["items"]:
            if item["status"] in ("done", "later"):
                continue
            refs = " ".join(f"#{r['number']}" for r in item["refs"])
            owner = " (maintainer)" if item["owner"] == "maintainer" else ""
            why = f" - {item['why']}" if item["why"] else ""
            out.append(f"- **{STATUS_LABEL[item['status']]}**{owner} {item['title']}"
                       f"{' ' + refs if refs else ''}{why}")
        out.append("")

    # --- in flight --------------------------------------------------------
    out += ["## In flight", ""]
    in_flight = state["in_flight"]
    if in_flight is None:
        out += ["_unknown - could not reach GitHub._", ""]
    elif not in_flight:
        out += ["Nothing open.", ""]
    else:
        out += ["| PR | Title | Branch | CI | Updated |", "|---|---|---|---|---|"]
        for pr in in_flight:
            draft = " _(draft)_" if pr["draft"] else ""
            out.append(f"| #{pr['number']} | {pr['title']}{draft} | `{pr['branch']}` | {pr['ci']} | {pr['updated']} |")
        out.append("")

    # --- next up ----------------------------------------------------------
    out += ["## Next up", ""]
    groups = state["next_up"]
    if groups is None:
        out += ["_unknown - could not reach GitHub._", ""]
    elif not groups:
        out += ["No open issues.", ""]
    else:
        order = [label for label in FIX_LABELS if label in groups]
        order += sorted(label for label in groups if label not in FIX_LABELS)
        for label in order:
            out.append(f"**{label}**")
            out.append("")
            for issue in groups[label]:
                out.append(f"- #{issue['number']} {issue['title']}")
            out.append("")

    # --- parked -------------------------------------------------------------
    parked = state.get("parked") or []
    out += [f"## Parked ({len(parked)})", ""]
    if parked:
        out += ["Out-of-scope findings from `/park`, waiting for triage:", ""]
        out += [f"- {item}" for item in parked]
    else:
        out.append("Nothing parked.")
    out.append("")

    # --- bioc readiness ---------------------------------------------------
    out += ["## Bioconductor readiness", ""]
    if latest_bioc:
        # Provenance (CI log vs cache) stays in status.json: it changes between
        # otherwise identical runs and would churn this committed file.
        out.append(f"Latest checked commit `{latest_bioc.get('sha')}` ({latest_bioc.get('date')}).")
        out.append("")
        for detail in latest_bioc.get("errors") or []:
            out.append(f"- ERROR: {detail}")
        if latest_bioc.get("errors"):
            out.append("")
        history = bioc.get("history") or []
        if len(history) > 1:
            out += ["| Commit | Date | E | W | N |", "|---|---|---|---|---|"]
            for entry in history[-6:]:
                out.append(f"| `{entry.get('sha')}` | {entry.get('date')} | {entry.get('error')} "
                           f"| {entry.get('warning')} | {entry.get('note')} |")
            out.append("")
    else:
        out += ["_No BiocCheck run recorded yet._", ""]

    # --- published site ---------------------------------------------------
    out += ["## Published site", ""]
    if pkg["site_sha"]:
        # None means the comparison could not be made at all. Folding that into
        # "current with main" would be a freshness claim precisely when nothing
        # is known, so the three cases stay distinct.
        if pkg["commits_behind"] is None:
            out.append(f"pkgdown site built from `{pkg['site_sha']}`; how far main has moved "
                       "past it could not be determined.")
        elif pkg["commits_behind"]:
            out.append(f"pkgdown site built from `{pkg['site_sha']}`; main is `{pkg['main_sha']}` "
                       f"- **{pkg['commits_behind']} commit(s) behind**.")
        else:
            out.append(f"pkgdown site built from `{pkg['site_sha']}` - current with main.")
    else:
        out.append("_Site build commit could not be determined._")
    if not pkg["has_config"]:
        out.append("")
        out.append("No `_pkgdown.yml` - the site builds on pkgdown defaults.")
    out.append("")

    # --- devel news -------------------------------------------------------
    news = state["news"]
    out += [f"## {news['heading'] or 'Devel changes'} (from NEWS.md)", ""]
    out += news["body"] or [
        "_No entries yet - nothing user-visible has landed in this cycle._"
    ]
    out.append("")

    # --- landed -----------------------------------------------------------
    out += ["## Recently landed", ""]
    landed = state["landed"]
    if landed is None:
        out += ["_unknown - could not reach GitHub._", ""]
    else:
        for pr in landed:
            out.append(f"- #{pr['number']} {pr['title']} ({pr['merged']})")
        out.append("")

    # Which branch this checkout is on, and what is uncommitted in it, is
    # per-machine: rendering it into a committed file would mean STATUS.md
    # differed for every clone and changed every time the tree got dirty.
    # It stays in status.json, which the dashboard reads.

    # --- branches ---------------------------------------------------------
    out += ["## Branch hygiene", ""]
    buckets = state["branches"]
    if buckets is None:
        out += ["_unknown - could not reach GitHub._", ""]
    else:
        out.append(f"{len(buckets['landed'])} landed, {len(buckets['abandoned'])} abandoned, "
                   f"{len(buckets['open'])} open, {len(buckets['unknown'])} with no PR.")
        out.append("")
        prunable = buckets["landed"] + buckets["abandoned"]
        if prunable:
            out += ["Safe to prune (landed or closed without merging). **Run these yourself** - "
                    "agents do not delete branches.", "", "```bash"]
            # A git ref may contain ; $ ` & | ( ), all of which a shell acts on.
            # This block is meant to be pasted into one, so every name is quoted
            # and `--` ends the option list.
            local = [shlex.quote(row["branch"]) for row in prunable if "local" in row["where"]]
            remote = [shlex.quote(row["branch"]) for row in prunable if "remote" in row["where"]]
            if local:
                out.append("git branch -D -- " + " ".join(local))
            if remote:
                out.append("git push origin --delete -- " + " ".join(remote))
            out += ["```", ""]
        if buckets["unknown"]:
            out += ["No PR found for these - check before deleting:", ""]
            for row in buckets["unknown"]:
                out.append(f"- `{row['branch']}` ({row['where']})")
            out.append("")

    if DEGRADED:
        out += ["## Incomplete data", ""]
        for message in DEGRADED:
            out.append(f"- {message}")
        out.append("")

    out += ["---", "", f"_generated {state['generated']}_"]
    return "\n".join(out) + "\n"


def render_html(state):
    """Inline the state into the dashboard template, ready to publish.

    The page is a snapshot: a published artifact can only reach claude.ai
    connectors and there is no GitHub connector, so it cannot query the repo
    itself. Emitting it from the same dict as STATUS.md is what keeps the two
    from disagreeing.
    """
    try:
        with open(TEMPLATE, encoding="utf-8") as handle:
            template = handle.read()
    except OSError:
        degrade("dashboard template missing, so no HTML was written")
        return None
    marker = "/*__STATUS_JSON__*/"
    if marker not in template:
        degrade("dashboard template has no /*__STATUS_JSON__*/ marker")
        return None
    # The payload sits in a <script type="application/json">, so the only
    # sequence that could break out of it is a literal "</script>".
    payload = json.dumps(state, sort_keys=True).replace("</", "<\\/")
    return template.replace(marker, payload, 1)


def write_atomic(path, text):
    """Replace `path` in one step, without a temp name anyone else can share.

    The SessionStart refresh and a manual `/mesa-status` can overlap, and a
    fixed `<path>.tmp` would have them writing the same file. A unique name in
    the same directory keeps `os.replace` atomic, and the `finally` means an
    interrupted run leaves nothing behind - `STATUS.md.tmp` sits in the package
    root, where a stray file is not just untidy.
    """
    directory = os.path.dirname(path) or "."
    # Prefixed with the target's own name, not hidden with a leading dot, so
    # `^STATUS\.md` in .Rbuildignore covers a leftover as well as the real file.
    fd, tmp = tempfile.mkstemp(dir=directory, prefix=os.path.basename(path) + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(text)
        # mkstemp is 0600 by design; these are ordinary derived files, so give
        # them the mode a plain open() would have.
        umask = os.umask(0)
        os.umask(umask)
        os.chmod(tmp, 0o666 & ~umask)
        os.replace(tmp, path)
        tmp = None
    finally:
        if tmp is not None:
            try:
                os.unlink(tmp)
            except OSError:
                pass


def main():
    parser = argparse.ArgumentParser(description="Derive mesa's project state.")
    parser.add_argument("--max-age", type=int, default=None,
                        help="skip entirely if STATUS.md is newer than this many seconds")
    parser.add_argument("--no-log", action="store_true",
                        help="skip the BiocCheck CI-log fetch (the slow probe)")
    parser.add_argument("--json-only", action="store_true", help="do not write STATUS.md")
    parser.add_argument("--html", action="store_true",
                        help="also write .claude/state/dashboard.html, ready to publish")
    parser.add_argument("--quiet", action="store_true", help="print nothing on success")
    parser.add_argument("--sync-issue", action="store_true",
                        help=f"rewrite #{ROADMAP_ISSUE}'s generated checklist if it is out of date")
    parser.add_argument("--dry-run", action="store_true",
                        help="with --sync-issue: print the diff instead of editing the issue")
    parser.add_argument("--roadmap-file", metavar="PATH",
                        help=f"read the roadmap TOML from PATH instead of #{ROADMAP_ISSUE}; with "
                             "--sync-issue, also write it into the issue's data block")
    parser.add_argument("--artifact-url", action="store_true",
                        help="print the published page's URL (empty if none yet) and exit")
    parser.add_argument("--set-artifact-url", metavar="URL",
                        help="record the published page's URL and exit")
    parser.add_argument("--mark-published", action="store_true",
                        help="clear the refresh flags the hooks set, and exit")
    parser.add_argument("--park", metavar="ISSUE",
                        help="append the note on stdin to the parking lot for ISSUE "
                             "(#N, N or -), then exit without refreshing")
    args = parser.parse_args()

    if args.park is not None:
        line = park(args.park, sys.stdin.read())
        if line is None:
            print("Nothing parked: the note on stdin is empty.", file=sys.stderr)
            return 1
        print(f"Parked in {parking_lot_path()}: {line[2:]}")
        return 0

    shared = shared_dir()
    url_file = os.path.join(shared, "artifact-url.txt")
    if args.set_artifact_url:
        write_atomic(url_file, args.set_artifact_url.strip() + "\n")
        return 0
    if args.artifact_url:
        print(read_artifact_url(url_file) or "")
        return 0
    if args.mark_published:
        for name in ("refresh-needed", "refresh-pending"):
            try:
                os.unlink(os.path.join(shared, name))
            except FileNotFoundError:
                pass
        return 0

    if args.max_age is not None and os.path.exists(MD_PATH):
        if time.time() - os.path.getmtime(MD_PATH) < args.max_age:
            if not args.quiet:
                print(f"STATUS.md is under {args.max_age}s old; not refreshing.")
            return 0

    fetch_refs()
    run_info = latest_main_run()
    bioc_run = run_info
    if run_info is not None and run_info.get("conclusion") != "success":
        bioc_run = latest_main_run(successful=True)
    state = {
        "generated": iso(now()),
        "repo": REPO,
        "release": collect_release(),
        "news": collect_news(),
        "recommendation": collect_recommendation(),
        "worktree": collect_worktree(),
        "in_flight": collect_in_flight(),
        "landed": collect_landed(),
        "next_up": collect_next_up(),
        "parked": collect_parked(),
        "ci": collect_ci(run_info),
        "coverage": collect_coverage(),
        "bioccheck": collect_bioccheck(bioc_run, allow_log=not args.no_log),
        "pkgdown": collect_pkgdown(),
        "branches": collect_branches(),
    }
    roadmap, issue_body = collect_roadmap(args.roadmap_file)
    seed_toml = roadmap.pop("toml") if roadmap else None
    state["roadmap"] = roadmap
    sync_verdict = None
    if args.sync_issue:
        sync_verdict = sync_issue(roadmap, issue_body, read_artifact_url(url_file),
                                  args.dry_run, seed_toml if args.roadmap_file else None)
    state["degraded"] = DEGRADED

    # Render first: a template problem is recorded in DEGRADED, and it must
    # reach status.json and STATUS.md like every other incomplete field.
    html = render_html(state) if args.html else None

    os.makedirs(STATE_DIR, exist_ok=True)
    write_atomic(JSON_PATH, json.dumps(state, indent=2, sort_keys=True) + "\n")
    written = ["status.json"]
    if not args.json_only:
        write_atomic(MD_PATH, render(state))
        written.append("STATUS.md")
    if html is not None:
        write_atomic(HTML_OUT, html)
        written.append("dashboard.html")
    elif args.html and os.path.exists(HTML_OUT):
        # Never leave an old page behind for /mesa-status to republish as new.
        os.unlink(HTML_OUT)
    if not args.quiet:
        print("Wrote " + ", ".join(written) + "."
              + (f" {sync_verdict}." if sync_verdict else "")
              + (f" {len(DEGRADED)} field(s) incomplete." if DEGRADED else ""))
    return 0


def read_artifact_url(url_file):
    """The page's URL, migrating it from the old per-worktree location once."""
    try:
        with open(url_file, encoding="utf-8") as handle:
            return handle.read().strip() or None
    except OSError:
        pass
    # Before the shared dir, the URL sat in .claude/state of whichever checkout
    # first published - usually the main one, listed first by `worktree list`.
    main_tree = (git("worktree", "list", "--porcelain") or "").split("\n", 1)[0]
    candidates = [os.path.join(STATE_DIR, "artifact-url.txt")]
    if main_tree.startswith("worktree "):
        candidates.append(os.path.join(main_tree[len("worktree "):], ".claude", "state", "artifact-url.txt"))
    for legacy in candidates:
        try:
            with open(legacy, encoding="utf-8") as handle:
                url = handle.read().strip()
        except OSError:
            continue
        if url:
            write_atomic(url_file, url + "\n")
            return url
    return None


if __name__ == "__main__":
    sys.exit(main())
