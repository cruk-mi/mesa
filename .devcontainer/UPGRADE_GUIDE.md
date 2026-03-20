# MESA R Version Upgrade Guide

## R / Bioconductor version map

| R version | Bioconductor | Expected release  |
|-----------|-------------|-------------------|
| 4.3       | 3.18        | released          |
| 4.4       | 3.19 / 3.20 | released          |
| 4.5       | 3.21 / 3.22 | released ← current|
| 4.6       | 3.23        | ~April 29 2026    |

💡 Bioconductor releases twice per year. Each R version covers two
   Bioc releases. Always use the latest Bioc for a given R version.

---

## Step 1 — Check compatibility BEFORE upgrading (15 min)

Run this the day a new Bioc version is announced:

    Rscript .devcontainer/check_bioc_compat.R <new_bioc_version>

    # Example for R 4.6 / Bioc 3.23:
    Rscript .devcontainer/check_bioc_compat.R 3.23

### What to look for:
- ✅ All packages found          → safe to proceed
- ✗  BLOCKING (Imports/Depends)  → must fix before upgrading
- ✗  non-blocking (Suggests)     → install manually, low priority

If any Imports/Depends packages are missing:
  - Check if the package was renamed (search Bioconductor release notes)
  - Check if it's still on bioconda: mamba search -c bioconda <package>
  - Consider moving it to Suggests if truly optional

---

## Step 2 — Upgrade your R environment

On the HPC, check if the new R module is available:

    module avail R

If not available yet, ask your sysadmin or wait — don't proceed
until the new R version is properly installed.

Load the new R and verify:

    module load R/4.6.0   # adjust version as needed
    R --version
    R -e "BiocManager::version()"   # should match target Bioc version

---

## Step 3 — Update DESCRIPTION (5 min)

One line change only:

    # In DESCRIPTION, update the R version constraint:
    R (>= 4.5.0)  →  R (>= 4.6.0)

    # Quick command:
    sed -i 's/R (>= 4.5.0)/R (>= 4.6.0)/' DESCRIPTION
    grep "R (>=" DESCRIPTION   # verify

---

## Step 4 — Update your R library (20 min)

In R with the new version loaded:

    # Set correct Bioc version
    BiocManager::install(version = "3.23", ask = FALSE)

    # Update all packages
    BiocManager::install(ask = FALSE, update = TRUE)
    # → choose "a" for all when prompted

    # Install any missing Suggests packages
    BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3", ask = FALSE)

---

## Step 5 — Test MESA (30-60 min)

    # Load and check
    R -e "devtools::load_all()"
    R -e "devtools::check()"

### Triage the output:

    0 errors, 0 warnings    → done, go to Step 6
    errors or warnings      → see common issues below

### Common issues and fixes:

1. Package stopped re-exporting a function
   Error: 'X' is not an exported object from 'namespace:Y'
   Fix:   Replace Y::X() with the true owner pkg::X()
          grep -rn "Y::X" R/

2. Function renamed or removed in updated package
   Error: could not find function "X"
   Fix:   Check package changelog, find replacement function
          browseURL("https://bioconductor.org/packages/PKG")  

3. Broken Rd cross-reference
   Warning: @seealso refers to unavailable topic X::Y
   Fix:   Remove link or update to correct function name
          grep -rn "seealso.*X::Y" R/

4. S4 class link missing package anchor
   Note: Rd \link{} targets missing package anchors: GRanges-class
   Fix:  Replace \linkS4class{X} with \link[PKG]{X-class}
         sed -i 's/\\linkS4class{GRanges}/\\link[GenomicRanges]{GRanges-class}/g' R/*.R

---

## Step 6 — Compare old vs new packages (5 min)

Useful to document what changed for the PR description:

    R -e "
    pkgs_new <- as.data.frame(
      installed.packages(lib.loc = .libPaths()[1])[, c('Package','Version')])
    write.csv(pkgs_new, 'installed_packages_R4.6_bioc3.23.csv', row.names=FALSE)
    "

    # Diff against previous snapshot
    # (we saved installed_packages_R4.5_bioc3.22.csv during this upgrade)

---

## Step 7 — Commit and PR

    git checkout -b update_to_R_4.6
    git add DESCRIPTION
    git commit -m "fix: update compatibility for R 4.6 / Bioconductor 3.23

    - Bump minimum R requirement to 4.6.0 in DESCRIPTION
    - [list any code fixes made in Step 5]"

    git push origin update_to_R_4.6
    # Open PR against dev

---

## Checklist

- [ ] check_bioc_compat.R shows all green for new Bioc version
- [ ] R version constraint updated in DESCRIPTION
- [ ] BiocManager updated to new Bioc version
- [ ] All packages updated via BiocManager
- [ ] devtools::check() shows 0 errors, 0 warnings
- [ ] package_comparison CSV saved for reference
- [ ] PR opened against dev
