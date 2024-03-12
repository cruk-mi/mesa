#' GC content across the human genome build hg38/GRCh38
#'
#' A set of data.table/GRanges object containing GC content for fixed size windows across the genome
#' Downloaded from ichorCNA,  https://github.com/broadinstitute/ichorCNA, formatted for hmmcopy.
#' Note chromosomes are in size order, not numerical order.
#' Each file contains the chromosome, start and end of the bin and an average GC content metric for that bin.
#'
#' @source <https://github.com/broadinstitute/ichorCNA>
#' @format NULL
"gc_hg38_1000kb"

#' @rdname gc_hg38_1000kb
#' @format NULL
"gc_hg38_50kb"

#' @rdname gc_hg38_1000kb
#' @format NULL
"gc_hg38_500kb"
