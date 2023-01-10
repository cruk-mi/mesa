## code to prepare `ENCODEbadRegions` dataset goes here

ENCODEbadRegions <- regioneR::toGRanges("/data/cep/Methylation/refData/hg38-blacklist.v2.noChr.noDesc.bed")

usethis::use_data(ENCODEbadRegions, overwrite = TRUE, compress = "xz")
