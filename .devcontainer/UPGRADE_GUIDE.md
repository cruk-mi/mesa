# MESA R Version Upgrade Guide

> **Who this is for:** anyone upgrading MESA to a new R / Bioconductor release.
> **Time required:** ~2–3 hours (most of which is `devtools::check()` running unattended).
> **When to do it:** when a new Bioconductor release is announced. Check dates at
> https://bioconductor.org/about/release-announcements/

---

## TL;DR — what you actually edit

```
DESCRIPTION   ← bump the one line:  R (>= X.Y.0)
```

**That's it.** Everything else derives from that single line:

- `.devcontainer/resolve_versions.sh` parses `R (>= X.Y.0)` from `DESCRIPTION`,
  then looks up the latest *released* Bioconductor for that R version from
  Bioconductor's published `config.yaml`.
- **CI** (`check-bioc.yml`) and the **devcontainer image build**
  (`build-image.yml`) both call that resolver, so they pick the right R, the
  right Bioc, and the right `bioconductor/bioconductor_docker:RELEASE_X_Y` base
  image automatically.
- Both devcontainer variants (`slim`, `full`) install mesa's dependencies
  straight from `DESCRIPTION`, so the package list is never edited by hand.

### Optional: pinning a specific Bioconductor version

Each R version covers **two** Bioc releases (e.g. R 4.6 ↔ Bioc 3.23 and 3.24).
By default the resolver picks the **newest released** one. To stay on an older
Bioc while a newer one exists, uncomment a line in `.devcontainer/versions.env`:

```bash
BIOC_VERSION=3.23   # pin; otherwise auto-derived from DESCRIPTION
```

---

## Version map (reference only — you do not edit this)

| R version | Bioconductor    | Status              |
|-----------|-----------------|---------------------|
| 4.3       | 3.17 / 3.18     | superseded          |
| 4.4       | 3.19 / 3.20     | superseded          |
| 4.5       | 3.21 / 3.22     | superseded          |
| 4.6       | 3.23 / 3.24     | current             |

**Rules:**
- Each R version covers two Bioc releases. The resolver uses the latest released
  one unless you pin (see above).
- Bioconductor releases in April and October. R major releases happen once per year.
- Do not upgrade mid-project. Coordinate with the team before starting.

---

## Step 1 - Check compatibility (15 min)

Run this **before** touching any file. It queries the target Bioc repo and checks
every DESCRIPTION dependency is available.

```bash
Rscript .devcontainer/check_bioc_compat.R <new_bioc_version>
```

### Reading the output

| Result | Meaning | Action |
|--------|---------|--------|
| `✅ All packages available` | Safe to proceed | Continue to Step 2 |
| `✗ BLOCKING (Imports/Depends)` | Package missing from new Bioc | See troubleshooting below |
| `✗ non-blocking (Suggests)` | Optional package missing | Install manually later, not urgent |

### If a blocking package is missing

```bash
# 1. Check if it was renamed - scan the Bioc release notes
#    https://bioconductor.org/news/

# 2. If it's genuinely dropped, consider moving it from Imports → Suggests
#    in DESCRIPTION (only if it's not called in every workflow)
```

---

## Step 2 - Create a branch (2 min)

Always work on a dedicated branch. Never upgrade directly on `dev` or `main`.

```bash
git checkout dev
git pull origin dev
git checkout -b update_to_R_<version>
```

---

## Step 3 - Update DESCRIPTION (2 min)

One line change only:

```bash
# Replace the old R version constraint with the new one (example):
sed -i 's/R (>= 4.5.0)/R (>= 4.6.0)/' DESCRIPTION

# Verify
grep "R (>=" DESCRIPTION
```

Then confirm the whole toolchain resolves correctly from that one edit:

```bash
bash .devcontainer/resolve_versions.sh
# Expected (example):
#   R_VERSION=4.6
#   BIOC_VERSION=3.23
#   BIOC_RELEASE=RELEASE_3_23
```

> **CI and the image build need no edits.** `check-bioc.yml` and
> `build-image.yml` both call `resolve_versions.sh`, so they pick up the new
> versions automatically once DESCRIPTION is updated.
>
> **RSPM Ubuntu codename:** the only value still hardcoded is the `rspm` URL
> codename (`noble`) in `check-bioc.yml` and `install.R`. It tracks the Bioc
> Docker base image's Ubuntu version, not the R/Bioc version. Only change it if
> the base image changed Ubuntu (focal = 20.04, jammy = 22.04, noble = 24.04).

---

## Step 4 - Update your R library on HPC (20 min)

On the HPC with the new R loaded:

```r
# Set Bioc to the new version (use the version resolve_versions.sh reported)
BiocManager::install(version = "3.23", ask = FALSE)

# Update all installed packages to their new Bioc-compatible versions
BiocManager::install(ask = FALSE, update = TRUE)
```

If R itself is not yet available on the HPC, ask your sysadmin — do not install R
yourself on a shared system.

---

## Step 5 - Test MESA on HPC (30–60 min)

```r
devtools::load_all()
devtools::check()
```

Look at the **final summary only**:

```
0 errors ✔ | 0 warnings ✔ | N notes   ← good, proceed
X errors   | X warnings              ← must fix, see below
```

### Triage guide for common errors

**1. Package stopped re-exporting a function**
```
Error: 'mutate' is not an exported object from 'namespace:plyranges'
```
```bash
grep -rn "plyranges::mutate" R/
sed -i 's/plyranges::mutate(/dplyr::mutate(/g' R/*.R
```

**2. Function renamed or removed**
```
Error: could not find function "X"
```
```bash
# Check the package changelog at https://bioconductor.org/packages/<PKG>/
grep -rn "PKG::X" R/
```

**3. Broken @seealso cross-reference**
```
Warning: @seealso refers to unavailable topic GenomeInfoDb::genome
```
```bash
grep -rn "seealso.*GenomeInfoDb::genome" R/
```

**4. S4 class missing package anchor**
```
Note: Rd \link{} targets missing package anchors: GRanges-class
```
```bash
sed -i 's/\\linkS4class{GRanges}/\\link[GenomicRanges]{GRanges-class}/g' R/*.R
```

**5. Vignette fails in headless environment (no display)**
```
Error: dev.control() called without an open graphics device
```
```r
# Add to the first knitr::opts_chunk$set() in the failing vignette:
knitr::opts_chunk$set(dev = "png")
```

---

## Step 6 - Test the devcontainer in Codespaces (30 min)

This verifies the Docker image built correctly with the new R version.

### 6.1 - Trigger the image build

`build-image.yml` builds and pushes the `slim` and `full` images to ghcr.io. It
triggers automatically on push to `dev`/`main` when `.devcontainer/**` or
`DESCRIPTION` changes. For a feature branch, trigger it manually:

```
github.com/cruk-mi/mesa
→ Actions → Build Devcontainer Images
→ Run workflow → select your branch → Run workflow
```

Wait for both `build (slim)` and `build (full)` to complete (~20–40 min).

> **Note:** the "Run workflow" button only appears when `build-image.yml` exists
> on the default branch. If you don't see it, temporarily add your branch to the
> `branches:` list in `build-image.yml`, push, then revert after testing.

### 6.2 - Open a new Codespace on your branch

```
github.com/cruk-mi/mesa
→ Code (green button) → Codespaces tab
→ "..." → New with options
→ Branch: your branch
→ Machine type: 4-core (minimum for R)
→ Create codespace
```

> **Important:** use "New with options" to select your branch — a plain "New
> codespace" defaults to `main`.
>
> **Existing Codespace?** It uses a cached old image. Force a rebuild:
> `Ctrl+Shift+P` → "Codespaces: Rebuild Container".

### 6.3 - Verify R, Bioc and mesa

The container is the Bioconductor base image — R is on the default `PATH`
(`/usr/local/bin/R`); there is **no conda environment to activate**.

```bash
# Confirm the resolved versions
R --version | head -1            # Expected: R version 4.6.x (matches DESCRIPTION)
bash .devcontainer/resolve_versions.sh
```

```r
# Start R, then:
BiocManager::version()           # should match the resolved BIOC_VERSION
packageVersion("mesa")
library(mesa)

# Confirm slim vs full variant:
# FALSE on slim (heavy data package absent), TRUE on full
requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38", quietly = TRUE)

# Run the test suite (data-dependent tests skip on slim)
devtools::test()
```

#### Expected output for the slim variant

```
R version 4.6.x
BiocManager: 3.23
requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38"): FALSE  ← correct for slim
devtools::test(): [ PASS x | FAIL 0 | WARN 0 | SKIP x ]
```

---

## Step 7 - Commit and open PR (5 min)

```bash
git add DESCRIPTION \
        R/*.R \       # only if code fixes were needed in Step 5
        man/*.Rd      # only if docs were regenerated

git commit -m "fix: update compatibility for R X.Y / Bioconductor A.B

- Bump minimum R requirement to X.Y.0 in DESCRIPTION
- [any code fixes from Step 5]"

git push origin update_to_R_<version>

# Open PR against dev (not main)
```

---

## Checklist

Copy this into your PR description:

```
- [ ] check_bioc_compat.R shows all green for target Bioc version
- [ ] DESCRIPTION: R (>= X.Y.0) updated
- [ ] resolve_versions.sh reports the expected R / Bioc / RELEASE tag
- [ ] BiocManager updated to new Bioc version on HPC
- [ ] devtools::check() shows 0 errors, 0 warnings
- [ ] Codespaces: both images built; correct R version, mesa loads, tests pass
- [ ] PR opened against dev
```

---

## Quick reference - files changed per upgrade

| File | What to change | Command |
|------|---------------|---------|
| `DESCRIPTION` | `R (>= X.Y.0)` | `sed -i 's/R (>= 4.5.0)/R (>= 4.6.0)/' DESCRIPTION` |
| `R/*.R` | Only if breaking changes found in Step 5 | varies |
| `man/*.Rd` | Auto-regenerated by `devtools::document()` | do not edit manually |
| `.devcontainer/versions.env` | Optional — only to pin a non-latest Bioc | uncomment `BIOC_VERSION=` |
