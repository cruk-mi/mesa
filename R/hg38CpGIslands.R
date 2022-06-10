#' CpG islands for hg38.
#'
#' A GRanges object containing CpG islands in hg38 coordinates. Generated from http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz
#'
#' @format A GRanges object with 31144 ranges and 7 metadata columns:
#' \describe{
#'   \item{length}{Number of bases in each island}
#'   \item{cpgNum}{Number of CGs in each island}
#'   \item{gcNum}{Number of Cs or Gs in the island}
#'   \item{perCpg}{Percent CGs in the island}
#'   \item{perGc}{Percent of Cs or Gs in the island}
#'   \item{obsExp}{Unclear what this column means}
#'   ...
#' }
"hg38CpGIslands"
