# mesa

## Overview
mesa is an R package for Methylation Enrichment Sequencing Analysis, to investigate sequencing data that has been enriched for regions of CG methylation. 
This includes both MEDIP-Seq (Methylated DNA immunoprecipitation sequencing) or MBD-seq (methyl-CpG binding domain protein enriched sequencing).

This package builds off the [qsea package](https://github.com/MatthiasLienhard/qsea), using the `qseaSet` object type defined in that package as the object that holds the data. One of the most useful features of the package is in defining tidyverse verbs to act on a qseaSet, such as `filter`, `mutate` and `select`, but there also includes a range of functionality for other tasks, such as plotting and quality control.

This package has been used in two published papers, in [Nature Cancer](https://www.nature.com/articles/s43018-022-00415-9) and [Nature Communications](https://www.nature.com/articles/s41467-024-47195-7).

Development was undertaken internally between September 2022 and March 2024; these commits are now incorporated into this repository. This has involved an extensive rewriting of the package functions, with many function names changing, see [the NEWS.md file for all the details](NEWS.md).

The package is under active development with plans to submit to BioConductor, and a full vignette is under construction. If you are interested in using the package, we're happy to hear from you.


## Install the package

mesa requires R ≥ 4.5.0 and Bioconductor 3.22. The easiest way to install
mesa and all its dependencies in one step is with the provided installer script:
```r
source("https://raw.githubusercontent.com/cruk-mi/mesa/main/inst/scripts/install_mesa.R")
```

This script automatically sets the correct Bioconductor version for your R
installation, installs mesa with all dependencies, and validates that your
environment is consistent.

**Manual installation:**

If you prefer to install manually:
```r
# 1. Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 2. Set Bioconductor to the correct version for your R
BiocManager::install(ask = FALSE)

# 3. Install mesa
BiocManager::install("cruk-mi/mesa", dependencies = TRUE)

# 4. Verify everything is consistent
BiocManager::valid()
```

**Development version:**

To install the latest development version directly from GitHub:
```r
install.packages("devtools")
devtools::install_github("cruk-mi/mesa")
```

> ⚠️ The development version may be unstable. Use the installer script above
> for a validated installation.
