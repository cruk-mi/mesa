---
name: capture-skill
description: Capture a newly worked-out mesa procedure as a project skill, so it is not re-derived next session. Use after solving something non-obvious and repeatable, or when the human says to remember or write down how something is done.
---

# Capturing a procedure as a skill

When you work out a mesa procedure that was not written down and will come up again, write
it down here rather than letting it evaporate.

## Is it worth a skill?

Write a skill when **all** of these hold:

- It will recur. A one-off investigation is not a skill.
- It is **specific to mesa or Bioconductor.** Generic R practice is already covered by the
  `posit-dev/skills` plugins — do not re-document what `r-package-development` or
  `testing-r-packages` already says.
- It was non-obvious. If the next agent would get it right from the code alone, skip it.
- It is procedural. A single invariant belongs in `AGENTS.md`; a multi-step procedure
  belongs in a skill.

If it is a short, always-relevant rule → add a line to `AGENTS.md` instead.
If it belongs to an existing skill → **extend that skill**, do not create a near-duplicate.

## How

Use the `skill-creator` skill for the authoring mechanics. mesa's conventions on top:

- Location: `.claude/skills/<kebab-name>/SKILL.md`.
- Frontmatter needs `name` and `description`. **The description is the part that is always
  in context** — write it so the next agent knows when to load the skill. Name the
  triggering situations ("Use when…"), not just the topic.
- Keep the body focused. Link to files in the repo rather than copying their contents,
  which goes stale.
- State precedence against any overlapping plugin skill.
- Record concrete anchors — file paths, function names, reference PR/commit SHAs. Those are
  what make a skill worth more than a general description.

## After writing it

1. Add it to the skill tables in **both** `AGENTS.md` and `CLAUDE.md`.
2. Confirm `.Rbuildignore` still excludes `^\.claude$` — nothing in this directory may ever
   reach the package tarball.
3. Commit as `chore(agents): add <name> skill` with the AI co-author trailer.
