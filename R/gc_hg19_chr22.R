#' GC content of chromosome 22 (hg19), in 100 kb bins
#'
#' Fixed-size bins with the GC fraction of each bin, for the GC correction in
#' [addHMMcopyCNV()]. It matches the chr22 reads shipped with the MEDIPSData
#' package, so that the function has a small example that runs on hg19.
#'
#' @docType data
#' @name gc_hg19_chr22_100kb
#' @usage data(gc_hg19_chr22_100kb)
#' @format A `data.frame` with 513 rows and the columns:
#'   \describe{
#'     \item{chr}{Chromosome, `"chr22"`.}
#'     \item{start}{1-based start of the bin.}
#'     \item{end}{1-based end of the bin.}
#'     \item{gc}{GC fraction of the bin (0–1), or `-1` for a bin that is all
#'       `N`.}
#'   }
#' @source Computed from `BSgenome.Hsapiens.UCSC.hg19` by
#'   `data-raw/gc_hg19_chr22_100kb.R`.
#' @seealso [map_hg19_chr22_100kb], [addHMMcopyCNV()]
#' @keywords datasets genomics
"gc_hg19_chr22_100kb"
