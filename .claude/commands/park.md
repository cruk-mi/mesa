---
description: Park an out-of-scope finding for later triage instead of fixing it now
argument-hint: <what you found, and where>
allowed-tools: Bash(git rev-parse:*), Bash(git branch --show-current), Bash(mkdir -p:*), Bash(printf:*)
---

Park this finding: $ARGUMENTS

It is out of scope for the current session (see "Session scope and the parking lot" in
`AGENTS.md`). **Do not fix it, do not open an issue for it, and do not touch any file for it
other than the parking lot.**

1. Resolve the parking lot in the **main checkout**, so every worktree writes to the same
   file:

   ```bash
   LOT="$(dirname "$(git rev-parse --path-format=absolute --git-common-dir)")/.claude/state/parking-lot.md"
   mkdir -p "$(dirname "$LOT")"
   ```

2. Append exactly one line, where `<issue>` is the issue this session is working on
   (`#N`, or `-` if none) and `<note>` is the finding rewritten as one self-contained
   sentence that names the file (and line, if known):

   ```bash
   printf -- '- %s (%s, %s) %s\n' "$(date +%F)" "$(git branch --show-current)" "<issue>" "<note>" >> "$LOT"
   ```

3. Reply with one line, `Parked: <note>`, then carry on with the session's issue.

At the end of the session, list everything parked from it in the PR description under
**Parked**.
