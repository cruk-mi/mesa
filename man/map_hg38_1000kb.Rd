% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_hg38.R
\docType{data}
\name{map_hg38_1000kb}
\alias{map_hg38_1000kb}
\alias{map_hg38_50kb}
\alias{map_hg38_500kb}
\title{GC content across the human genome build hg38/GRCh38}
\source{
\url{https://github.com/broadinstitute/ichorCNA}
}
\usage{
map_hg38_1000kb

map_hg38_50kb

map_hg38_500kb
}
\description{
A set of data.table/GRanges object containing mapabillity values for fixed size windows across the genome
Downloaded from ichorCNA,  https://github.com/broadinstitute/ichorCNA, formatted for hmmcopy.
Note chromosomes are in size order, not numerical order.
Each file contains the chromosome, start and end of the bin and an average GC content metric for that bin.
}
\keyword{datasets}
