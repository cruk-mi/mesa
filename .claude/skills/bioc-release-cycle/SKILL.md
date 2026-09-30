---
name: bioc-release-cycle
description: mesa's three-phase version cadence — open a X.Y.Z.9000 devel section, land work PRs, cut X.Y.(Z+1) — and the post-merge tag and push to Bioconductor (BiocStaging). Use when bumping DESCRIPTION Version, opening or closing a development cycle, renaming a NEWS.md devel heading, or preparing a release, tag or Bioconductor push.
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

**The parking lot must be empty first.** If `/mesa-status` lists anything under Parked,
promote it as "Triage" in `AGENTS.md` describes, and cut only when nothing is left.

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

## After phase 3 merges — tag and push to Bioconductor

A release is not out until Bioconductor has it. Once the phase-3 PR merges, hand the human
this checklist, and do not open phase 1 until they confirm step 2 is done:

1. **Tag on GitHub.** `origin/main` is now the phase-3 squash commit.

   ```bash
   git fetch origin
   git tag -a vX.Y.(Z+1) origin/main -m "mesa X.Y.(Z+1)"
   git push origin vX.Y.(Z+1)
   ```

   Pushing the tag starts the `release` workflow (`.github/workflows/release.yml`, about
   two hours). It checks the tag against `DESCRIPTION` and `NEWS.md`, runs the r-universe
   Bioconductor build, and publishes the GitHub Release. **Check the run succeeded before
   step 2**, because that build is the gate Bioconductor's own build would otherwise be:

   ```bash
   gh run list --repo cruk-mi/mesa --workflow release.yml --branch vX.Y.(Z+1) --limit 1   # note the run ID
   gh run watch <run-id> --repo cruk-mi/mesa --exit-status
   gh release view vX.Y.(Z+1) --repo cruk-mi/mesa                     # tarball attached
   ```

   A failed run has published nothing. If a job failed for a transient reason (a runner or
   network error), re-run it: `gh run rerun <run-id> --repo cruk-mi/mesa --failed`.
   Otherwise re-tag. If `verify` fails, the tag is on the wrong commit or the bump is
   incomplete. If `runiverse` or `build` fails, fix it on `main` first. Then move the tag
   to the right commit (`git tag -f -a vX.Y.(Z+1) <commit> -m …` and
   `git push -f origin vX.Y.(Z+1)`), which starts a fresh run. Never move a tag while its
   run is still going: that run refuses to publish, and only the new one releases.

2. **Push to Bioconductor.** During review, Bioconductor builds from the `devel` branch of
   [`BiocStaging/mesa`](https://github.com/BiocStaging/mesa), not from
   `git.bioconductor.org` (that copy lags behind and is not what the reviewers build). Its
   `devel` is an ancestor of our `main`, so this is a plain fast-forward. There is no local
   `devel` branch, so push `main`'s commit to it:

   ```bash
   # the remote persists in the clone; this adds it only if it is missing (fresh clone)
   git remote get-url biocstaging >/dev/null 2>&1 ||
     git remote add biocstaging https://github.com/BiocStaging/mesa.git
   git fetch biocstaging
   git merge-base --is-ancestor biocstaging/devel origin/main && echo fast-forward
   git push biocstaging origin/main:devel
   ```

   HTTPS works with an existing `gh` login that has push access. The version bump triggers
   a new build, and its report appears on the package's submission issue.

3. **Then open the next cycle:** phase 1, `X.Y.(Z+1).9000`.

4. **Refresh the roadmap:** `/mesa-status`. The release's `tag:vX.Y.(Z+1)` check turns its
   "Cut …" item Done and marks the release as released on the Next steps page and in #124.

**Pushing to Bioconductor is outward-facing and starts a public build. An agent never does
it without asking.** Ask the human for explicit permission every time, even when an earlier
push was approved, and show the exact command and the commits it will send
(`git log --oneline biocstaging/devel..origin/main`). By default the human runs it.

Before handing over the checklist, check that the previous release reached Bioconductor
(`git ls-remote https://github.com/BiocStaging/mesa refs/heads/devel` should be at the
last `vX.Y.Z` tag), that the tag exists (`git ls-remote --tags origin`), and that it has a
GitHub Release (`gh release view vX.Y.Z --repo cruk-mi/mesa`; only from v0.99.8 on, since
earlier tags predate the `release` workflow). If any of these is
missing, say so first. 0.99.6 went out without its GitHub tag, which was only added when
0.99.7 was cut.

## Precedence

The `open-source` plugin's `create-release-checklist` skill describes generic semver
release hygiene. **This skill wins.** mesa's cadence is Bioconductor's, not semver's.

## Related

Open issue [#97](https://github.com/cruk-mi/mesa/issues/97) tracks automating these two
bump PRs. If you are working on that issue, this file is the specification.
