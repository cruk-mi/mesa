---
name: bioc-release-cycle
description: mesa's three-phase version cadence — open a X.Y.Z.9000 devel section, land work PRs, cut X.Y.(Z+1). Use when bumping DESCRIPTION Version, opening or closing a development cycle, renaming a NEWS.md devel heading, or preparing a release or tag.
---

# mesa release cycle

Every change set runs through three phases. Work PRs sit **between** phases 1 and 3.

> **Never bundle a version bump into a work PR.** A bump is always its own branch and its
> own PR. This is the single rule that matters most here.

## Version convention

Bioconductor pre-submission uses `0.99.z`, becoming `1.0.0` on acceptance. After that,
**odd minor = devel** (`1.1.z`), **even minor = release** (`1.0.z`). The `.9000` suffix
marks an open development section.

Current series: `0.99.x` (pre-submission).

## Phase 1 — open the development section

Branch `chore/bump-version-X.Y.Z.9000` off `main`. **Two commits, in this order.**
Reference: PR #96 = `59f870b` (commits `d4b0413` then `cad49bd`).

| commit message | change |
|---|---|
| `chore(DESCRIPTION): bump version to X.Y.Z.9000` | `DESCRIPTION` `Version:` only |
| `docs(NEWS): open X.Y.Z.9000 development section` | insert `# mesa X.Y.Z.9000` plus a blank line at the very top of `NEWS.md` |

The new `NEWS.md` section starts empty. Work PRs fill it.

## Phase 2 — the work

Each feature/fix PR targets `main` and records its user-visible changes under the
`# mesa X.Y.Z.9000` heading, **in the same PR that makes the change**. See `mesa-docs-news`
for the entry format. `DESCRIPTION` `Version:` is not touched during this phase.

## Phase 3 — cut the release

After the work PRs merge, branch `chore/bump-version-X.Y.(Z+1)` off `main`. **One commit**,
`chore(version): bump to X.Y.(Z+1)`. Reference: PR #99 = `b5e80a0` (also `d13bbeb`,
`61c890a`).

- `DESCRIPTION`: `X.Y.Z.9000` → `X.Y.(Z+1)`
- `NEWS.md`: rename the heading `# mesa X.Y.Z.9000` → `# mesa X.Y.(Z+1)`
- Commit body: *"Bump DESCRIPTION Version and rename the NEWS.md devel heading from
  X.Y.Z.9000 to the X.Y.(Z+1) release, mirroring PR #NN."*

Before committing phase 3, run the full check ladder (see `bioc-check-ladder`) and fix any
errors or warnings.

## After each phase

**Stop.** Tell the human the branch is ready. Tags (`vX.Y.Z`) are applied by the human
after merge to `main` — never by an agent.

## Precedence

The `open-source` plugin's `create-release-checklist` skill describes generic semver
release hygiene. **This skill wins.** mesa's cadence is Bioconductor's, not semver's.

## Related

Open issue [#97](https://github.com/cruk-mi/mesa/issues/97) tracks automating these two
bump PRs. If you are working on that issue, this file is the specification.
