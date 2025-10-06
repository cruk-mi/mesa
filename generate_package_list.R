#!/usr/bin/env Rscript
# generate_package_list.R
# - Keep R 4.3
# - Force ggplot2 < 4 (compatible with R 4.3)
# - Remove any preexisting ggplot2 entries

if (!requireNamespace("desc", quietly = TRUE)) {
  install.packages("desc", repos = "https://cloud.r-project.org")
}
library(desc)

# Read DESCRIPTION deps
d <- desc::desc(file = "DESCRIPTION")
deps <- d$get_deps()
deps <- subset(deps, package != "R")

# Drop base R stdlib
base_r <- c("base","compiler","datasets","graphics","grDevices","grid","methods",
            "parallel","splines","stats","stats4","tcltk","tools","translations","utils")
deps <- subset(deps, !(package %in% base_r))

# Bioconductor packages to map to bioconductor-*
bioc_pkgs <- c(
  "limma","plyranges","GenomicRanges","IRanges","BSgenome","ChIPseeker","qsea",
  "BiocGenerics","Biostrings","ComplexHeatmap","GenomeInfoDb","GenomicAlignments",
  "Rsamtools","BiocParallel","biomaRt","HMMcopy","MEDIPS","rtracklayer","MEDIPSData",
  "org.Hs.eg.db","org.Mm.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene","TxDb.Mmusculus.UCSC.mm10.knownGene",
  "BSgenome.Hsapiens.NCBI.GRCh38","BSgenome.Hsapiens.UCSC.hg19"
)

# Map to conda names
deps$conda_name <- ifelse(deps$package %in% bioc_pkgs,
                          paste0("bioconductor-", tolower(deps$package)),
                          paste0("r-", tolower(deps$package)))

# Dedup/sort
pkg_lines <- sort(unique(deps$conda_name))

# ---- Enforce ggplot2 < 4 for R 4.3 ----
# Remove any ggplot2 entries (plain or pinned)
pkg_lines <- pkg_lines[!grepl("^r-ggplot2(\\b|[<>=]).*$", pkg_lines)]
# Add the required pin
pkg_lines <- sort(unique(c(pkg_lines, "r-ggplot2<4")))

# Pin R to 4.3
out <- c("r-base=4.3", pkg_lines)

writeLines(out, "package_list.txt")
cat("✅ package_list.txt written with", length(out), "lines (R 4.3, ggplot2<4).\n")
