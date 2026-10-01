#' mesa: Methylation Enrichment Sequencing Analysis
#'
#' mesa analyses DNA methylation measured by enrichment sequencing, such as
#' MeDIP-seq (methylated DNA immunoprecipitation) and MBD-seq (methyl-CpG
#' binding domain capture), including low-input samples such as cell-free DNA.
#' These assays pull down methylated fragments instead of converting
#' unmethylated cytosines, so the signal is a read count per genomic window
#' whose relation to methylation depends on the local CpG density. mesa builds
#' on the `qseaSet` class from \pkg{qsea}, which models that relation, and adds
#' functions to build, quality-control, explore and compare `qseaSet` objects
#' in a tidyverse style.
#'
#' @section Workflow:
#' A typical analysis follows these steps:
#'
#' 1. **Build** a `qseaSet` from BAM files with [makeQset()], which tiles the
#'    genome into windows, counts fragments and estimates copy number (see
#'    [addHMMcopyCNV()]). [addNormalisation()] fits the CpG-density-dependent
#'    enrichment used to estimate beta values.
#' 2. **Check quality** with [getSampleQCSummary()] and
#'    [plotCorrelationMatrix()], and subset samples or windows with
#'    [filter()][dplyr::filter], [filterWindows()] or [subsetQset()].
#' 3. **Extract** counts, normalised reads per million (nrpm) or beta values
#'    per window with [getDataTable()].
#' 4. **Explore** sample structure with principal components ([getPCA()]) or
#'    UMAP ([getUMAP()]).
#' 5. **Compare** groups of samples with [calculateDMRs()], which fits a
#'    negative binomial model per window and returns differentially methylated
#'    regions (DMRs), and summarise them with [summariseDMRsByContrast()].
#' 6. **Annotate** windows or DMRs with nearby genes and genomic features with
#'    [annotateWindows()] and [plotGenomicFeatureDistribution()].
#'
#' Two small `qseaSet` objects are shipped for examples:
#' [exampleTumourNormal] (human, paired tumour/normal samples) and
#' [exampleMouse] (mouse cell line).
#'
#' @section Vignettes:
#' * `vignette("introduction", package = "mesa")`: overview and quick start.
#' * `vignette("generation", package = "mesa")`: building `qseaSet` objects.
#' * `vignette("data-and-qc", package = "mesa")`: data tables and quality
#'   control.
#' * `vignette("pca", package = "mesa")`: PCA and UMAP.
#' * `vignette("differentially-methylated-regions", package = "mesa")`: DMRs
#'   and their annotation.
#'
#' @examples
#' # A qseaSet of 10 samples (5 tumour/normal pairs) over windows on chr7
#' exampleTumourNormal
#' getSampleTable(exampleTumourNormal)
#'
#' # Beta values (methylation estimates) per window and sample
#' exampleTumourNormal %>%
#'     getDataTable(normMethod = "beta") %>%
#'     head()
#'
#' # Differentially methylated regions between tumour and normal samples
#' exampleTumourNormal %>%
#'     calculateDMRs(variable = "tumour", contrasts = "first") %>%
#'     head()
#' @import qsea
"_PACKAGE"
