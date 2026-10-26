#' Mappability of chromosome 22 (hg19), in 100 kb bins
#'
#' Fixed-size bins with the mean mappability of each bin, for the mappability
#' correction in [addHMMcopyCNV()]. It matches the chr22 reads shipped with the
#' MEDIPSData package, so that the function has a small example that runs on
#' hg19.
#'
#' @docType data
#' @name map_hg19_chr22_100kb
#' @usage data(map_hg19_chr22_100kb)
#' @format A `data.frame` with 513 rows and the columns:
#'   \describe{
#'     \item{chr}{Chromosome, `"chr22"`.}
#'     \item{start}{1-based start of the bin.}
#'     \item{end}{1-based end of the bin.}
#'     \item{map}{Mean mappability of the bin (0–1).}
#'   }
#' @source The ENCODE/CRG 100mer alignability track
#'   `wgEncodeCrgMapabilityAlign100mer.bigWig` from the UCSC hg19 downloads,
#'   averaged by `data-raw/map_hg19_chr22_100kb.R`.
#' @seealso [gc_hg19_chr22_100kb], [addHMMcopyCNV()]
#' @keywords datasets genomics
"map_hg19_chr22_100kb"
