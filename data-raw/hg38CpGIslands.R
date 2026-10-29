## code to prepare `hg38CpGIslands` dataset goes here

## UCSC "CpG Islands" track (cpgIslandExt) for hg38, downloaded on 2026-10-29.
## The table has a bin column first and writes the name as "CpG: n". Drop the
## bin, close up the name, remove "chr" and sort by chromosome then start (as
## bedtools sort did when this was first built from the same table).

url <- paste0("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/",
    "cpgIslandExt.txt.gz")
path <- tempfile(fileext = ".txt.gz")
download.file(url, path, mode = "wb")

cpgIslandExt <- read.table(path, sep = "", quote = "")
cpgIslandExt <- data.frame(
    chrom = sub("chr", "", cpgIslandExt$V2),
    start = cpgIslandExt$V3,
    end = cpgIslandExt$V4,
    name = paste0(cpgIslandExt$V5, cpgIslandExt$V6),
    cpgIslandExt[, 7:12]
)
colnames(cpgIslandExt)[5:10] <- c("length", "cpgNum", "gcNum", "perCpg",
                                  "perGc", "obsExp")
cpgIslandExt <- cpgIslandExt[order(cpgIslandExt$chrom, cpgIslandExt$start,
                                   method = "radix"), ]
rownames(cpgIslandExt) <- NULL

hg38CpGIslands <- regioneR::toGRanges(cpgIslandExt) %>%
  dplyr::select(-name)

usethis::use_data(hg38CpGIslands, overwrite = TRUE, internal = TRUE, compress = "xz")
