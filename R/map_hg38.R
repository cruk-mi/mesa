#' Mappability across the human genome (hg38/GRCh38)
#'
#' Fixed-size bins with average mappability per bin. Bins are ordered by
#' chromosome size (not strictly numeric).
#'
#' @docType data
#' @name map_hg38_1000kb
#' @aliases map_hg38_50kb map_hg38_500kb
#' @usage data(map_hg38_1000kb); data(map_hg38_500kb); data(map_hg38_50kb)
#' @format One of:
#' \enumerate{
#' \item A \code{data.table}/\code{data.frame} with columns:
#'   \describe{
#'     \item{chrom}{Chromosome (e.g., "chr1", ..., "chrX", "chrY").}
#'     \item{start}{0-based start of the bin.}
#'     \item{end}{1-based end of the bin.}
#'     \item{map}{Average mappability in the bin (0–1).}
#'   }
#' \item A \code{GRanges} with genomic bins and \code{mcols(.)$map} storing
#'   the average mappability (0–1) per bin.
#' }
#' @return \code{map_hg38_1000kb}, \code{map_hg38_500kb}, and
#' \code{map_hg38_50kb} return objects as described in \code{@format} with
#' bin sizes of 1000 kb, 500 kb, and 50 kb, respectively.
#' @source <https://github.com/broadinstitute/ichorCNA>
#' @keywords datasets genomics
"map_hg38_1000kb"

#' @rdname map_hg38_1000kb
"map_hg38_50kb"

#' @rdname map_hg38_1000kb
"map_hg38_500kb"
