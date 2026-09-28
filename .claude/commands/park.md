---
description: Park an out-of-scope finding for later triage instead of fixing it now
argument-hint: <what you found, and where>
allowed-tools: Bash(python3 .claude/scripts/mesa-status.py --park:*)
---

Park this finding: $ARGUMENTS

It is out of scope for the current session (see "Session scope and the parking lot" in
`AGENTS.md`). **Do not fix it, do not open an issue for it, and do not touch any file for it
other than the parking lot.**

1. Rewrite the finding as `<note>`: one self-contained sentence that names the file (and
   line, if known).

2. Append it with **one** command, where `<issue>` is the issue this session is working on
   (`115`, or `-` if none):

   ```bash
   python3 .claude/scripts/mesa-status.py --park <issue> <<'EOF'
   <note>
   EOF
   ```

   The note goes in on stdin through a quoted heredoc, so `$`, backticks and quotes in it
   reach the file as written. The script adds the date and branch, collapses the note to
   one line, and writes to the **main checkout's** `.claude/state/parking-lot.md`, so
   every worktree shares one lot.

3. Reply with one line, `Parked: <note>`, then carry on with the session's issue.

At the end of the session, list everything parked from it in the PR description under
**Parked**. That list is the durable copy: the lot itself is per-machine.
