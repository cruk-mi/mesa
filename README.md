# mesa

## Overview
mesa is an R package for Methylation Enrichment Sequencing Analysis, to investigate sequencing data that has been enriched for regions of CG methylation. 
This includes both MEDIP-Seq (Methylated DNA immunoprecipitation sequencing) or MBD-seq (methyl-CpG binding domain protein enriched sequencing).

This package builds off the [qsea package](https://github.com/MatthiasLienhard/qsea), using the `qseaSet` object type defined in that package as the object that holds the data. One of the most useful features of the package is in defining tidyverse verbs to act on a qseaSet, such as `filter`, `mutate` and `select`, but there also includes a range of functionality for other tasks, such as plotting and quality control.

This package has been used in a Nature Cancer paper, [cfDNA methylome profiling for detection and subtyping of small cell lung cancers] [https://www.nature.com/articles/s43018-022-00415-9), as well as an in [CUPiD: A cfDNA methylation-based tissue-of-origin classifier for Cancers of Unknown Primary](https://assets.researchsquare.com/files/rs-3758456/v1/6f505510-a351-49f3-bba4-882e111be410.pdf).

Development was undertaken internally between September 2022 and March 2024; these commits are now incorporated into this repository. This has involved an extensive rewriting of the package functions, with many function names changing, see [the NEWS.md file for all the details](NEWS.md).

The package is under active development with plans to submit to BioConductor, and a full vignette is under construction. If you are interested in using the package, we're happy to hear from you.

## Install the package
Currently, the easiest way to install the package is with devtools:
```{r}
install.packages("devtools")
devtools::install_github("cruk-mi/mesa")
```
