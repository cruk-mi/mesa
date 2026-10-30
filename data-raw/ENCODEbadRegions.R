## code to prepare `ENCODEbadRegions` dataset goes here

## ENCODE Blacklist v2 for hg38 (Amemiya, Kundaje and Boyle 2019,
## doi:10.1038/s41598-019-45839-z), from the Boyle-Lab/Blacklist repository.
## Remove "chr" and the description column.

url <- paste0("https://github.com/Boyle-Lab/Blacklist/raw/",
    "d0527cdea356bf57d30535b48783804ca8d21031/lists/",
    "hg38-blacklist.v2.bed.gz")
path <- tempfile(fileext = ".bed.gz")
download.file(url, path, mode = "wb")

bedFile <- tempfile(fileext = ".bed")
blacklist <- read.table(path, sep = "\t", quote = "", comment.char = "")
blacklist <- data.frame(sub("chr", "", blacklist$V1), blacklist$V2,
                        blacklist$V3)
write.table(blacklist, bedFile, sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = FALSE)

ENCODEbadRegions <- regioneR::toGRanges(bedFile)

usethis::use_data(ENCODEbadRegions, overwrite = TRUE, compress = "xz")
