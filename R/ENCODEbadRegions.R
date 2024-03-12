#' 2019 ENCODE list of poorly mapped regions in hg38
#'
#' The 2019 ENCODE list of regions that are poorly mapped (2019 https://doi.org/10.1038/s41598-019-45839-z).
#' List downloaded from https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz and edited to remove "chr".
#'
#' @format A Granges object with 636 ranges:
#' \describe{
#'   \item{seqnames}{chromosome}
#'   \item{ranges}{position on the chromosome}
#'   ...
#' }
#' @source \url{https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz}
"ENCODEbadRegions"
