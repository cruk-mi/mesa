## code to prepare `map_hg19_chr22_100kb` dataset goes here
##
## Mean mappability of 100 kb bins on chr22 of hg19, from the public
## ENCODE/CRG 100mer alignability track on the UCSC download server.

binSize <- 100000L
trackUrl <- paste0(
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/",
    "wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig"
)
chrLength <- GenomeInfoDb::seqlengths(
    BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
)[["chr22"]]

# the same tiling as addHMMcopyCNV(), which drops the last partial bin
starts <- seq(1L, chrLength - binSize + 1L, by = binSize)
track <- rtracklayer::import(
    rtracklayer::BigWigFile(trackUrl),
    which = GenomicRanges::GRanges("chr22", IRanges::IRanges(1, chrLength)),
    as = "RleList"
)[["chr22"]]

map_hg19_chr22_100kb <- data.frame(
    chr = "chr22",
    start = starts,
    end = starts + binSize - 1L,
    map = as.numeric(IRanges::viewMeans(
        IRanges::Views(track, starts, starts + binSize - 1L)
    ))
)

usethis::use_data(map_hg19_chr22_100kb, overwrite = TRUE, compress = "xz")
