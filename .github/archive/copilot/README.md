# Archived GitHub Copilot configuration

**Status: inactive.** Nothing in this directory is loaded by any tool. Copilot is not
currently used on mesa, but this is kept rather than deleted so it can be revived without
being rewritten from scratch.

`AGENTS.md` at the repo root is the canonical, tool-agnostic instruction file. Copilot Chat
and most current agents read `AGENTS.md` directly, so a lot of what is here is only needed
for the **Copilot coding agent** and for path-scoped rules.

## What is here

| Archived path | Live path if reactivated | What it does |
|---|---|---|
| `copilot-instructions.md` | `.github/copilot-instructions.md` | Repo-wide Copilot instructions; defers to `AGENTS.md`. |
| `instructions/*.instructions.md` | `.github/instructions/` | Path-scoped rules applied by glob (`R/**`, `tests/**`, docs, CI). |
| `agents/*.md` | `.github/agents/` | Copilot-side role definitions mirroring `.claude/agents/`. |

## To reactivate

1. `git mv` each directory back to the live path in the table above.
2. **Re-read them against `AGENTS.md` first.** They were accurate when archived and will
   have drifted since — `AGENTS.md` is canonical and wins on any disagreement.
3. The `.claude/agents/*` files are self-contained and no longer reference
   `.github/agents/*`. If you reactivate the Copilot roles, the two sets must be kept in
   sync by hand; that duplication is why this was archived.
4. If the Copilot **coding agent** (not just Chat) is being used, it also needs a
   `copilot-setup-steps.yml` workflow to preinstall the Bioconductor dependency stack in
   its cloud runner. That file was never merged; see the `chore/ai-agent-setup` branch
   history (commit `74fd4f1`) for a starting point.

Archived 2026-09-17, when the setup was consolidated around `AGENTS.md` and `.claude/`.
