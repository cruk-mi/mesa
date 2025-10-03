#!/usr/bin/env Rscript
if (!requireNamespace("desc", quietly = TRUE)) {
  install.packages("desc", repos = "https://cloud.r-project.org")
}
library(desc)

d <- desc::desc(file = "DESCRIPTION")
deps <- d$get_deps()
deps <- deps[deps$package != "R", ]

# Base R stdlib that must NOT go into conda list
base_r <- c("base","compiler","datasets","graphics","grDevices",
            "grid","methods","parallel","splines","stats",
            "stats4","tcltk","tools","translations","utils")

deps <- deps[!(deps$package %in% base_r), ]

# Bioconductor packages to map to bioconductor-*
bioc_pkgs <- c(
  "GenomicRanges","IRanges","BSgenome","ChIPseeker","qsea",
  "BiocGenerics","Biostrings","ComplexHeatmap","GenomeInfoDb",
  "GenomicAlignments","Rsamtools","BiocParallel","biomaRt","HMMcopy",
  "MEDIPS","rtracklayer","workflows","MEDIPSData",
  "org.Hs.eg.db","org.Mm.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene","TxDb.Mmusculus.UCSC.mm10.knownGene",
  "BSgenome.Hsapiens.NCBI.GRCh38","BSgenome.Hsapiens.UCSC.hg19"
)

# Map to conda names
deps$conda_name <- ifelse(deps$package %in% bioc_pkgs,
                          paste0("bioconductor-", tolower(deps$package)),
                          paste0("r-", tolower(deps$package)))

# Finalize: pin R, sort, unique
deps_sorted <- sort(unique(deps$conda_name))
out <- c("r-base=4.3", deps_sorted)
writeLines(out, "package_list.txt")

cat("✅ package_list.txt written with", length(out), "lines.\n")
