---
applyTo: ".github/**,.devcontainer/**"
---

# CI / DevOps rules

Follow [`AGENTS.md`](../../AGENTS.md). For CI and container config:

- CI lives in `.github/workflows/` — `check-bioc.yml` (R CMD check + BiocCheck) and
  `build-image.yml` (devcontainer image). Keep workflow changes minimal and reviewable.
- The Copilot coding agent's environment is prepared by `copilot-setup-steps.yml`; its job
  **must** be named `copilot-setup-steps`. Reuse the prebuilt image
  `ghcr.io/cruk-mi/mesa/devcontainer:codespaces-slim` so the Bioc dependency stack is
  already present.
- The devcontainer (`.devcontainer/`) is the source of truth for the dev/agent environment:
  `Dockerfile`, `install.R` (installs deps from `DESCRIPTION`), `versions.env`,
  `resolve_versions.sh`. `DESCRIPTION` is the single source of truth for *which* packages —
  do not hardcode dependency lists elsewhere.
- Do not hardcode R or Bioconductor version numbers in scripts; read them from the env vars
  the resolver sets (`BIOC_VERSION`, etc.).
- Prefer Posit Package Manager (RSPM) binaries before source installs for speed.
- Pin GitHub-only deps to a specific SHA (see `install.R`).
