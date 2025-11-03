#' GC content across the human genome (hg38/GRCh38)
#'
#' Fixed-size bins with average GC content per bin. Bins are ordered by
#' chromosome size (not strictly numeric).
#'
#' @docType data
#' @name gc_hg38_1000kb
#' @aliases gc_hg38_50kb gc_hg38_500kb
#' @usage data(gc_hg38_1000kb); data(gc_hg38_500kb); data(gc_hg38_50kb)
#' @format One of:
#' \enumerate{
#' \item A \code{data.table}/\code{data.frame} with columns:
#'   \describe{
#'     \item{chrom}{Chromosome (e.g., "chr1", ..., "chrX", "chrY").}
#'     \item{start}{0-based start of the bin.}
#'     \item{end}{1-based end of the bin.}
#'     \item{gc}{Average GC fraction in the bin (0–1).}
#'   }
#' \item A \code{GRanges} with genomic bins and \code{mcols(.)$gc} storing
#'   the average GC fraction (0–1) per bin.
#' }
#' @return \code{gc_hg38_1000kb}, \code{gc_hg38_500kb}, and \code{gc_hg38_50kb}
#' return objects as described in \code{@format} with bin sizes of 1000 kb,
#' 500 kb, and 50 kb, respectively.
#' @source <https://github.com/broadinstitute/ichorCNA>
#' @keywords datasets genomics
"gc_hg38_1000kb"

#' @rdname gc_hg38_1000kb
"gc_hg38_50kb"

#' @rdname gc_hg38_1000kb
"gc_hg38_500kb"
