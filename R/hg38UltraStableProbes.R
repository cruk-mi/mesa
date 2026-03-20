#' Ultra-stable methylated regions (GRCh38)
#'
#' Genomic loci that are consistently methylated (“ultra-stable”) across many human tissues
#' and diseases, mapped to GRCh38. Regions correspond to Illumina array probe loci
#' (cg identifiers) derived as in Edgar et al. (2014), and are intended for QC and summary
#' metrics (e.g., [addHyperStableFraction()]).
#'
#' @format
#' A `GRanges` object with 974 ranges and 1 metadata column:
#' \describe{
#'   \item{Probe_ID}{Character cg identifier for the probe/locus (e.g., `"cg11733071"`).}
#' }
#' 
#' @details
#' - **Seqnames:** numeric chromosomes `"1"`–`"22"` (no `"chr"` prefix; UCSC style
#'   can be set if needed).  
#' - **Strand:** `"*"` for all ranges.  
#' - **Widths:** 3 bp windows (probe CpG ±1 bp).  
#' - **Genome metadata:** currently unset; set to `"hg38"` if you need an explicit
#'   tag.  
#'
#' When combining with other `GRanges`, harmonise styles/metadata as needed, e.g.
#' use `GenomeInfoDb::seqlevelsStyle()` to convert to `"UCSC"` (`"chr1"`, …) and
#' `GenomeInfoDb::genome()` to set `"hg38"`.
#'
#' @source
#' Processed from public array resources following Edgar et al. (2014).
#'
#' @references
#' Edgar R, Tan PPC, Portales-Casamar E, Pavlidis P (2014).
#' *Meta-analysis of human methylomes reveals stably methylated sequences surrounding CpG islands associated with high gene expression*.
#' <https://pubmed.ncbi.nlm.nih.gov/25493099/>
#'
#' @seealso
#' [addHyperStableFraction()]
#'
#' @examples
#' data(hg38UltraStableProbes)
#'
#' # Basic overview
#' hg38UltraStableProbes
#' length(hg38UltraStableProbes)             # 974
#' all(GenomicRanges::width(hg38UltraStableProbes) == 3)  # TRUE (3-bp windows)
#'
#' # The object uses numeric chromosomes and an unset genome tag
#' GenomeInfoDb::seqlevels(hg38UltraStableProbes)[1:5]
#' GenomeInfoDb::genome(hg38UltraStableProbes)  # empty
#'
#' # Harmonise style/metadata if needed
#' gr <- hg38UltraStableProbes
#' GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"  # or "NCBI" / "Ensembl" for 1,2... instead of chr1, chr2...
#' GenomeInfoDb::genome(gr) <- "hg38"  # or "GRCh38"
#'
#' # Count overlaps with a toy region
#' toy <- GenomicRanges::GRanges("1", IRanges::IRanges(951160, 951170))
#' GenomicRanges::countOverlaps(toy, hg38UltraStableProbes)
#'
#' # With UCSC-style example
#' toy_ucsc <- GenomicRanges::GRanges("chr1", IRanges::IRanges(951160, 951170))
#' GenomicRanges::countOverlaps(toy_ucsc, gr)   # after style harmonisation above
"hg38UltraStableProbes"

