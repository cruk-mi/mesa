## code to prepare `encodeBlacklist` dataset goes here

encodeBlacklist <- regioneR::toGRanges("/data/cep/Methylation/refData/hg38-blacklist.v2.noChr.noDesc.bed")


usethis::use_data(encodeBlacklist, overwrite = TRUE, internal = TRUE)
