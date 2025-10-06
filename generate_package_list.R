#!/usr/bin/env Rscript
if (!requireNamespace("desc", quietly = TRUE)) {
  install.packages("desc", repos = "https://cloud.r-project.org")
}
library(desc)

d <- desc::desc(file = "DESCRIPTION")
deps <- d$get_deps()
deps <- subset(deps, package != "R")

# Base R libs to drop
base_r <- c("base","compiler","datasets","graphics","grDevices","grid","methods",
            "parallel","splines","stats","stats4","tcltk","tools","translations","utils")
deps <- subset(deps, !(package %in% base_r))

# Bioconductor packages (add here as needed)
bioc_pkgs <- c(
  "limma","plyranges","GenomicRanges","IRanges","BSgenome","ChIPseeker","qsea",
  "BiocGenerics","Biostrings","ComplexHeatmap","GenomeInfoDb","GenomicAlignments",
  "Rsamtools","BiocParallel","biomaRt","HMMcopy","MEDIPS","rtracklayer","MEDIPSData",
  "org.Hs.eg.db","org.Mm.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene","TxDb.Mmusculus.UCSC.mm10.knownGene",
  "BSgenome.Hsapiens.NCBI.GRCh38","BSgenome.Hsapiens.UCSC.hg19"
)

# Map (note: workflows is CRAN!)
deps$conda_name <- ifelse(deps$package %in% bioc_pkgs,
                          paste0("bioconductor-", tolower(deps$package)),
                          paste0("r-", tolower(deps$package)))

# Optional: pick your R version. If you need ggplot2 >= 4, prefer R 4.4
r_line <- "r-base=4."  # use "r-base=4.3" and add "r-ggplot2<4" if you must stay on 4.3

# Sort + unique
pkg_lines <- sort(unique(deps$conda_name))
out <- c(r_line, pkg_lines)

# If you must stay on R 4.3, uncomment:
out <- c("r-base=4.3", pkg_lines, "r-ggplot2<4")

writeLines(out, "package_list.txt")
cat("✅ package_list.txt written with", length(out), "lines.\n")
