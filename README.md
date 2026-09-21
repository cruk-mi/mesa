# mesa

<!-- badges: start -->
[![R-CMD-check-bioc](https://github.com/cruk-mi/mesa/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/cruk-mi/mesa/actions/workflows/check-bioc.yml)
[![Codecov test coverage](https://codecov.io/gh/cruk-mi/mesa/branch/main/graph/badge.svg)](https://app.codecov.io/gh/cruk-mi/mesa)
[![License: GPL (>= 2)](https://img.shields.io/badge/license-GPL%20(%3E%3D%202)-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

<!--
Bioconductor shields. These 404 until mesa is accepted and appears in the
Bioconductor repository at 1.0, so they stay commented out until then.
Tracked in issue #109.

[![BioC release](https://bioconductor.org/shields/build/release/bioc/mesa.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/mesa)
[![BioC devel](https://bioconductor.org/shields/build/devel/bioc/mesa.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/mesa)
[![Platforms](https://bioconductor.org/shields/availability/devel/mesa.svg)](https://bioconductor.org/packages/devel/bioc/html/mesa.html)
[![Downloads](https://bioconductor.org/shields/downloads/release/mesa.svg)](https://bioconductor.org/packages/stats/bioc/mesa/)
[![Support posts](https://bioconductor.org/shields/posts/mesa.svg)](https://support.bioconductor.org/tag/mesa)
[![Years in BioC](https://bioconductor.org/shields/years-in-bioc/mesa.svg)](https://bioconductor.org/packages/release/bioc/html/mesa.html)
-->

## Overview
mesa is an R package for Methylation Enrichment Sequencing Analysis, to investigate sequencing data that has been enriched for regions of CG methylation. 
This includes both MEDIP-Seq (Methylated DNA immunoprecipitation sequencing) or MBD-seq (methyl-CpG binding domain protein enriched sequencing).

This package builds off the [qsea package](https://github.com/MatthiasLienhard/qsea), using the `qseaSet` object type defined in that package as the object that holds the data. One of the most useful features of the package is in defining tidyverse verbs to act on a qseaSet, such as `filter`, `mutate` and `select`, but there also includes a range of functionality for other tasks, such as plotting and quality control.

This package has been used in two published papers, in [Nature Cancer](https://www.nature.com/articles/s43018-022-00415-9) and [Nature Communications](https://www.nature.com/articles/s41467-024-47195-7).

Development was undertaken internally between September 2022 and March 2024; these commits are now incorporated into this repository. This has involved an extensive rewriting of the package functions, with many function names changing, see [the NEWS.md file for all the details](NEWS.md).

The package is under active development, with plans to submit to Bioconductor.

## Install the package

**Manual installation:**

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

## Quick start

mesa ships an example `qseaSet`, `exampleTumourNormal`: 10 samples (5 tumours and
5 matched adjacent normals) over a small part of chromosome 7. Everything below
runs without any data of your own.

```r
library(mesa)

# A qseaSet holds the sample metadata, counts, windows, CNV and library
# information for an experiment in a single object.
exampleTumourNormal

# The sample table is ordinary sample metadata, one row per sample.
exampleTumourNormal |>
    getSampleTable()
```

The distinguishing feature of mesa is that tidyverse verbs act on the `qseaSet`
itself, filtering or annotating the sample table while reducing every other part
of the object to match, behind the scenes:

```r
lungOnly <- exampleTumourNormal |>
    filter(tissue == "Lung") |>
    mutate(label = paste(patient, type, sep = "_"))
```

Methylation values come out as tables, either as normalised reads per million or
as beta values scaled between 0 (unmethylated) and 1 (methylated):

```r
exampleTumourNormal |>
    getNRPMTable()

exampleTumourNormal |>
    getBetaTable()
```

Dimensionality reduction is two steps, calculate then plot:

```r
exampleTumourNormal |>
    getPCA(nPC = 3) |>
    plotDimRed(colour = "group")
```

See the articles below for generating a `qseaSet` from BAM files, quality
control, and calling differentially methylated regions.

## Documentation

Full function reference and articles: <https://cruk-mi.github.io/mesa/>

| Article | Covers |
|---|---|
| [Introduction to Mesa](https://cruk-mi.github.io/mesa/articles/introduction.html) | The `qseaSet` object, its parts, and the tidyverse verbs that act on it |
| [Generating qseaSets](https://cruk-mi.github.io/mesa/articles/generation.html) | Building a `qseaSet` from aligned, de-duplicated BAM files |
| [Data Tables and Quality Control](https://cruk-mi.github.io/mesa/articles/data-and-qc.html) | Count, NRPM and beta tables, and assessing whether enrichment worked |
| [Dimensionality Reduction - PCAs and UMAPs](https://cruk-mi.github.io/mesa/articles/pca.html) | `getPCA()`, `getUMAP()` and plotting the results |
| [Differentially Methylated Regions and Annotation](https://cruk-mi.github.io/mesa/articles/differentially-methylated-regions.html) | Calling DMRs between groups and annotating them to genomic features |

The same articles are available from R once the package is installed, via
`vignette(package = "mesa")`.

## Citation

mesa underpins the analysis in two published papers. If you use it, please cite
whichever is applicable to your work:

> Chemi F, Pearce SP, Clipson A, *et al.* (2022). cfDNA methylome profiling for
> detection and subtyping of small cell lung cancers. *Nature Cancer* **3**(10),
> 1260-1270. <https://doi.org/10.1038/s43018-022-00415-9>

> Conway A-M, Pearce SP, Clipson A, *et al.* (2024). A cfDNA methylation-based
> tissue-of-origin classifier for cancers of unknown primary. *Nature
> Communications* **15**, 3292.
> <https://doi.org/10.1038/s41467-024-47195-7>

## Contributing and support

If you are interested in using the package, we're happy to hear from you.

- Questions, bug reports and feature requests: [open an issue](https://github.com/cruk-mi/mesa/issues)
- Contributors are expected to follow the [Code of Conduct](.github/CODE_OF_CONDUCT.md)

mesa is developed at the [Cancer Research UK National Biomarker Centre](https://cruknbc.org/).
