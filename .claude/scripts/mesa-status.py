#!/usr/bin/env python3
"""Derive mesa's project state from git, gh, NEWS.md and the gh-pages site.

Writes two files:

  .claude/state/status.json   machine-readable state (gitignored)
  STATUS.md                   the committed human/agent summary

Nothing here is hand-maintained, so nothing can go stale: every fact is read
back from the place humans and agents both actually write to. A change made on
github.com by a co-maintainer shows up on the next run with no further work.

Every probe degrades on its own. A missing `gh`, a revoked token or a dead
network yields "unknown" fields and a `degraded` list, never a traceback and
never a half-written file (both outputs are written atomically at the end).

Usage:
  mesa-status.py                  refresh unconditionally
  mesa-status.py --max-age 14400  refresh only if STATUS.md is older than 4h
  mesa-status.py --no-log         skip the BiocCheck CI-log fetch (the slow probe)
  mesa-status.py --json-only      write status.json but not STATUS.md
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timedelta, timezone

REPO = "cruk-mi/mesa"
WORKFLOW = "R-CMD-check-bioc"
PROTECTED = {"main", "dev", "master", "gh-pages", "HEAD"}
COVERAGE_TARGET = 80  # AGENTS.md "Local verification"
FIX_LABELS = ("ERROR_fix", "WARNING_fix", "NOTE_fix", "Bioc_feedback_fix")

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
STATE_DIR = os.path.join(ROOT, ".claude", "state")
JSON_PATH = os.path.join(STATE_DIR, "status.json")
MD_PATH = os.path.join(ROOT, "STATUS.md")
BIOC_HISTORY = os.path.join(STATE_DIR, "bioccheck-history.jsonl")
RECOMMENDATION = os.path.join(STATE_DIR, "recommendation.md")
TEMPLATE = os.path.join(ROOT, ".claude", "scripts", "dashboard-template.html")
HTML_OUT = os.path.join(STATE_DIR, "dashboard.html")

DEGRADED = []


def degrade(message):
    """Record why a field is unknown, so STATUS.md can say so out loud."""
    if message not in DEGRADED:
        DEGRADED.append(message)


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
        names = [label.get("name") for label in issue.get("labels") or []]
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


def latest_main_run():
    data = gh_json([
        "run", "list", "--repo", REPO, "--workflow", WORKFLOW, "--branch", "main",
        "--limit", "1", "--json", "databaseId,conclusion,headSha,createdAt,url",
    ])
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
LOG_PREFIX = re.compile(r"^.*?\d{4}-\d{2}-\d{2}T[\d:.]+Z\s?")
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
    if latest and latest.get("sha") == sha:
        return {"latest": latest, "history": history, "source": "cache"}
    if not allow_log:
        return {"latest": latest, "history": history, "source": "cache (log fetch skipped)"}

    log = run(["gh", "run", "view", str(run_info.get("databaseId")),
               "--repo", REPO, "--log"], timeout=180)
    if log is None:
        degrade("CI log for the latest main run could not be fetched, so BiocCheck counts are the cached ones")
        return {"latest": latest, "history": history, "source": "cache"}

    counts, errors, pending = None, [], None
    for line in log.splitlines():
        message = log_message(line)
        match = BIOC_SUMMARY.search(message)
        if match:
            counts = {
                "error": int(match.group(1)),
                "warning": int(match.group(2)),
                "note": int(match.group(3)),
            }
            pending = None
            continue
        if "ERROR: " in message:
            pending = message.split("ERROR: ", 1)[1].strip()
            errors.append(pending)
        elif pending is not None:
            if message and not FINDING_START.match(message):
                errors[-1] = (errors[-1].rstrip(", ") + " " + message).strip()
            else:
                pending = None
    errors = list(dict.fromkeys(error for error in errors if error))
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
        "--json", "number,headRefName,state",
    ])
    if prs is None:
        return None
    state_of = {}
    for pr in prs:
        ref = pr.get("headRefName")
        if not ref:
            continue
        # A branch reused across PRs: the most advanced state wins.
        rank = {"MERGED": 3, "OPEN": 2, "CLOSED": 1}
        if rank.get(pr.get("state"), 0) >= rank.get(state_of.get(ref, {}).get("state"), 0):
            state_of[ref] = {"state": pr.get("state"), "number": pr.get("number")}

    current = git("rev-parse", "--abbrev-ref", "HEAD")
    buckets = {"landed": [], "abandoned": [], "open": [], "unknown": []}
    for name in sorted(set(local) | set(remote)):
        # Protected branches, and the branch in hand, are never prune candidates.
        if name in PROTECTED or name == current:
            continue
        info = state_of.get(name)
        where = "local+remote" if name in local and name in remote else ("local" if name in local else "remote")
        row = {"branch": name, "where": where, "pr": info.get("number") if info else None}
        if info is None:
            buckets["unknown"].append(row)
        elif info["state"] == "MERGED":
            buckets["landed"].append(row)
        elif info["state"] == "OPEN":
            buckets["open"].append(row)
        else:
            buckets["abandoned"].append(row)
    return buckets


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
        if ci["behind"]:
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
        if pkg["commits_behind"]:
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
            local = [row["branch"] for row in prunable if "local" in row["where"]]
            remote = [row["branch"] for row in prunable if "remote" in row["where"]]
            if local:
                out.append("git branch -D " + " ".join(local))
            if remote:
                out.append("git push origin --delete " + " ".join(remote))
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
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        handle.write(text)
    os.replace(tmp, path)


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
    args = parser.parse_args()

    if args.max_age is not None and os.path.exists(MD_PATH):
        if time.time() - os.path.getmtime(MD_PATH) < args.max_age:
            if not args.quiet:
                print(f"STATUS.md is under {args.max_age}s old; not refreshing.")
            return 0

    run_info = latest_main_run()
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
        "ci": collect_ci(run_info),
        "coverage": collect_coverage(),
        "bioccheck": collect_bioccheck(run_info, allow_log=not args.no_log),
        "pkgdown": collect_pkgdown(),
        "branches": collect_branches(),
    }
    state["degraded"] = DEGRADED

    os.makedirs(STATE_DIR, exist_ok=True)
    write_atomic(JSON_PATH, json.dumps(state, indent=2, sort_keys=True) + "\n")
    written = ["status.json"]
    if not args.json_only:
        write_atomic(MD_PATH, render(state))
        written.append("STATUS.md")
    if args.html:
        html = render_html(state)
        if html is not None:
            write_atomic(HTML_OUT, html)
            written.append("dashboard.html")
    if not args.quiet:
        print("Wrote " + ", ".join(written) + "."
              + (f" {len(DEGRADED)} field(s) incomplete." if DEGRADED else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
