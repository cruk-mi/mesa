# MESA R Version Upgrade Guide

> **Who this is for:** anyone upgrading MESA to a new R / Bioconductor release.  
> **Time required:** ~2–3 hours (most of which is `devtools::check()` running unattended).  
> **When to do it:** when a new Bioconductor release is announced. Check dates at
> https://bioconductor.org/about/release-announcements/

---

## Version map

| R version | Bioconductor    | Status              |
|-----------|-----------------|---------------------|
| 4.3       | 3.17 / 3.18     | superseded          |
| 4.4       | 3.19 / 3.20     | superseded          |
| 4.5       | 3.21 / **3.22** | **current**         |
| 4.6       | 3.23            | ~April 29 2026      |

**Rules:**
- Each R version covers two Bioc releases. Always use the latest Bioc for a given R.
- Bioconductor releases in April and October. R major releases happen once per year.
- Do not upgrade mid-project. Coordinate with the team before starting.

---

## Files you will touch

```
DESCRIPTION                              ← bump R (>= X.X.0)
.devcontainer/versions.env               ← bump R_VERSION and BIOC_VERSION
.github/workflows/check-bioc.yml         ← bump R version and Bioc Docker image tag
```

Everything else (Dockerfile, package lists, install.R, build-image.yml) reads from
`versions.env` automatically - you should not need to edit them.

---

## Step 1 - Check compatibility (15 min)

Run this **before** touching any file. It queries the target Bioc repo and checks
every DESCRIPTION dependency is available.

```bash
Rscript .devcontainer/check_bioc_compat.R <new_bioc_version>

# Example:
Rscript .devcontainer/check_bioc_compat.R 3.23
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
#    https://bioconductor.org/news/bioc_X.XX_release/

# 2. Check if it's available on bioconda under a different name
mamba search -c bioconda <package_name>

# 3. If it's genuinely dropped, consider moving it from Imports → Suggests
#    in DESCRIPTION (only if it's not called in every workflow)
```

---

## Step 2 - Verify R is available on your system (5 min)

### HPC

```bash
module avail R

# Load the new version
module load R/4.6.0       # adjust as needed
R --version               # confirm version
```

If the module is not available yet, ask your sysadmin. Do not proceed until
the new R is properly installed - do not try to install R yourself on a shared HPC.

### Codespaces / local

The container will be rebuilt automatically in Step 7 once you update `versions.env`.
You do not need to install R manually.

---

## Step 3 - Create a branch (2 min)

Always work on a dedicated branch. Never upgrade directly on `dev` or `main`.

```bash
git checkout dev
git pull origin dev
git checkout -b update_to_R_<version>

# Example:
git checkout -b update_to_R_4.6
```

---

## Step 4 - Update DESCRIPTION (2 min)

One line change only:

```bash
# Replace old R version constraint with new one
sed -i 's/R (>= 4.5.0)/R (>= 4.6.0)/' DESCRIPTION

# Verify
grep "R (>=" DESCRIPTION
```

Expected output:
```
    R (>= 4.6.0),
```

---

## Step 5 - Update versions.env (2 min)

This is the single source of truth. Changing this file triggers the Docker
image rebuild in CI.

```bash
sed -i 's/R_VERSION=4.5/R_VERSION=4.6/' .devcontainer/versions.env
sed -i 's/BIOC_VERSION=3.22/BIOC_VERSION=3.23/' .devcontainer/versions.env

# Verify
cat .devcontainer/versions.env
```

Expected output:
```
R_VERSION=4.6
BIOC_VERSION=3.23
IMAGE_TAG_SLIM=codespaces-slim
IMAGE_TAG_FULL=full
```

---

## Step 6 - Update CI workflow (2 min)

`check-bioc.yml` has the R version and Bioc Docker image tag hardcoded
in the matrix config. Update both:

```bash
# Update R version
sed -i "s/r: '4.5.1'/r: '4.6.0'/" .github/workflows/check-bioc.yml

# Update Bioc Docker image tag
sed -i 's/RELEASE_3_22/RELEASE_3_23/' .github/workflows/check-bioc.yml

# Update RSPM URL Ubuntu codename if needed (check Bioc Docker image release notes)
# focal = Ubuntu 20.04, jammy = Ubuntu 22.04, noble = Ubuntu 24.04
# Only change this if the Bioc Docker base image changed Ubuntu version
sed -i 's/jammy/jammy/' .github/workflows/check-bioc.yml   # no-op if unchanged

# Verify
grep -E "r:|bioc:|cont:|rspm:" .github/workflows/check-bioc.yml
```

---

## Step 7 - Update your R library on HPC (20 min)

On the HPC with the new R loaded:

```r
# Set Bioc to the new version
BiocManager::install(version = "3.23", ask = FALSE)

# Update all installed packages to their new Bioc-compatible versions
# Choose "a" (all) when prompted
BiocManager::install(ask = FALSE, update = TRUE)

# Install any Suggests packages that are not auto-installed
BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3", ask = FALSE)
```

---

## Step 8 - Test MESA on HPC (30–60 min)

```r
devtools::load_all()
devtools::check()
```

Most of the time is `devtools::check()` running examples and vignettes unattended.
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
# Find all calls using the old namespace
grep -rn "plyranges::mutate" R/

# Fix: replace with the true owner
sed -i 's/plyranges::mutate(/dplyr::mutate(/g' R/*.R
```

**2. Function renamed or removed**
```
Error: could not find function "X"
```
```bash
# Check the package changelog
browseURL("https://bioconductor.org/packages/<PKG>/news/")

# Find all usages in MESA
grep -rn "PKG::X" R/
```

**3. Broken @seealso cross-reference**
```
Warning: @seealso refers to unavailable topic GenomeInfoDb::genome
```
```bash
# Find and remove or fix the broken link
grep -rn "seealso.*GenomeInfoDb::genome" R/

# Fix: remove the broken entry, or update to the correct function name
```

**4. S4 class missing package anchor**
```
Note: Rd \link{} targets missing package anchors: GRanges-class
```
```bash
# Fix all instances across R/ files
sed -i 's/\\linkS4class{GRanges}/\\link[GenomicRanges]{GRanges-class}/g' R/*.R
```

**5. Vignette fails in headless environment (no display)**
```
Error: dev.control() called without an open graphics device
```
```r
# Add to the first knitr::opts_chunk$set() in the failing vignette:
knitr::opts_chunk$set(dev = "png")
# Remove any dev.args = list(type = "cairo") - requires a display
```

---

## Step 9 - Test the devcontainer in Codespaces (30 min)

This step verifies the Docker image built correctly with the new R version
before merging.

### 9.1 - Trigger the image build

The `build-image.yml` workflow builds and pushes the Docker image to ghcr.io.
It triggers automatically on push to `dev` and `main` when `.devcontainer/**`
changes. For a feature branch, trigger it manually:

```
github.com/cruk-mi/mesa
→ Actions → Build Devcontainer Images
→ Run workflow → select your branch → Run workflow
```

Wait for both `build (slim)` and `build (full)` jobs to complete (~20–40 min).

> **Note:** The "Run workflow" button only appears when `build-image.yml` exists
> on the default branch (`main`). If you don't see it, temporarily add your
> branch to the `branches:` list in `build-image.yml`, push, then revert after
> testing.

### 9.2 - Open a new Codespace on your branch

```
github.com/cruk-mi/mesa
→ Code (green button) → Codespaces tab
→ "..." → New with options
→ Branch: update_to_R_4.6
→ Machine type: 4-core (minimum for R)
→ Create codespace
```

> **Important:** Use "New with options" to select your branch.
> A plain "New codespace" defaults to `main`.

> **If you have an existing Codespace:** it will use a cached old image.
> Force a rebuild: `Ctrl+Shift+P` → "Codespaces: Rebuild Container".

### 9.3 - Activate the R environment

The `mesa_env` conda environment is not activated by default in new terminal
sessions. Run:

```bash
source /opt/conda/etc/profile.d/conda.sh
conda activate mesa_env
```

To make this permanent for all future terminal sessions in this Codespace:

```bash
echo "source /opt/conda/etc/profile.d/conda.sh && conda activate mesa_env" >> ~/.bashrc
source ~/.bashrc
```

### 9.4 - Verify R version and environment

```bash
# Confirm R version matches versions.env
R --version | head -1
# Expected: R version 4.6.x

# Confirm versions.env was used correctly
cat .devcontainer/versions.env

# Confirm which image variant is running (slim vs full)
cat .devcontainer/devcontainer.json | grep image
```

### 9.5 - Verify mesa is installed inside R

```r
# Start R
R

# Check R and Bioc versions
R.version$version.string
BiocManager::version()     # should match BIOC_VERSION in versions.env

# Check mesa is installed
packageVersion("mesa")
library(mesa)

# Confirm slim vs full variant
# Returns FALSE for slim (full-only package), TRUE for full
requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38", quietly = TRUE)

# Run the test suite
devtools::test()
```

### 9.6 - Expected output for slim variant

```
R version 4.6.x
BiocManager: 3.23
packageVersion("mesa"): 0.99.0 (or current version)
requireNamespace("BSgenome.Hsapiens.NCBI.GRCh38"): FALSE  ← correct for slim
devtools::test(): [ PASS x | FAIL 0 | WARN 0 | SKIP x ]
```

---

## Step 10 - Snapshot installed packages (5 min)

Save a record of what's installed for future comparison and for the PR description.

```r
pkgs <- as.data.frame(
  installed.packages(lib.loc = .libPaths()[1])[, c("Package", "Version")],
  stringsAsFactors = FALSE
)
write.csv(pkgs,
  sprintf("installed_packages_R%s_bioc%s.csv",
    paste0(R.version$major, ".", substr(R.version$minor, 1, 1)),
    BiocManager::version()),
  row.names = FALSE)
```

To compare against the previous snapshot:

```r
old <- read.csv("installed_packages_R4.5_bioc3.22.csv")
new <- read.csv("installed_packages_R4.6_bioc3.23.csv")

dplyr::full_join(old, new, by = "Package", suffix = c("_old", "_new")) |>
  dplyr::mutate(status = dplyr::case_when(
    is.na(Version_old) ~ "new",
    is.na(Version_new) ~ "removed",
    Version_old == Version_new ~ "unchanged",
    TRUE ~ "updated"
  )) |>
  dplyr::filter(status != "unchanged") |>
  dplyr::arrange(status, Package) |>
  print(n = Inf)
```

---

## Step 11 - Commit and open PR (5 min)

```bash
# Stage all changed files
git add DESCRIPTION \
        .devcontainer/versions.env \
        .github/workflows/check-bioc.yml \
        R/*.R \       # only if code fixes were needed in Step 8
        man/*.Rd      # only if docs were regenerated

# Commit - adjust the bullet list to match what actually changed
git commit -m "fix: update compatibility for R 4.6 / Bioconductor 3.23

- Bump minimum R requirement to 4.6.0 in DESCRIPTION
- Update versions.env to R_VERSION=4.6, BIOC_VERSION=3.23
- Update check-bioc.yml matrix to R 4.6.0 / RELEASE_3_23
- [any code fixes from Step 8]"

git push origin update_to_R_4.6

# Open PR against dev (not main)
```

---

## Checklist

Copy this into your PR description:

```
- [ ] check_bioc_compat.R shows all green for target Bioc version
- [ ] DESCRIPTION: R (>= X.X.0) updated
- [ ] versions.env: R_VERSION and BIOC_VERSION updated
- [ ] check-bioc.yml: R version and RELEASE tag updated
- [ ] BiocManager updated to new Bioc version on HPC
- [ ] All packages updated via BiocManager::install(update=TRUE)
- [ ] devtools::check() shows 0 errors, 0 warnings
- [ ] Codespaces: image built and verified (correct R version, mesa loads, tests pass)
- [ ] Package snapshot CSV saved (installed_packages_RX.X_biocX.XX.csv)
- [ ] PR opened against dev
```

---

## Quick reference - files changed per upgrade

| File | What to change | Command |
|------|---------------|---------|
| `DESCRIPTION` | `R (>= X.X.0)` | `sed -i 's/R (>= 4.5.0)/R (>= 4.6.0)/' DESCRIPTION` |
| `versions.env` | `R_VERSION`, `BIOC_VERSION` | `sed -i 's/R_VERSION=4.5/R_VERSION=4.6/'` |
| `check-bioc.yml` | `r:`, `cont: RELEASE_X_XX` | `sed -i "s/r: '4.5.1'/r: '4.6.0'/"` |
| `R/*.R` | Only if breaking changes found in Step 8 | varies |
| `man/*.Rd` | Auto-regenerated by `devtools::document()` | do not edit manually |