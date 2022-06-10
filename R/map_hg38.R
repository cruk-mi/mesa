#' GC content across the human genome build hg38/GRCh38
#'
#' A set of data.table/GRanges object containing mapabillity values for fixed size windows across the genome
#' Downloaded from ichorCNA,  https://github.com/broadinstitute/ichorCNA, formatted for hmmcopy.
#' Note chromosomes are in size order, not numerical order.
#' Each file contains the chromosome, start and end of the bin and an average GC content metric for that bin.
#'
#' @source <https://github.com/broadinstitute/ichorCNA>
#' @format NULL
"map_hg38_1000kb"

#' @rdname map_hg38_1000kb
#' @format NULL
"map_hg38_10kb"

#' @rdname map_hg38_1000kb
#' @format NULL
"map_hg38_50kb"

#' @rdname map_hg38_1000kb
#' @format NULL
"map_hg38_500kb"
