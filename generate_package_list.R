#!/usr/bin/env Rscript

# Script: generate_package_list.R
# Purpose: Parse DESCRIPTION and generate a conda-compatible package_list.txt
# Usage:   Rscript generate_package_list.R

if (!requireNamespace("desc", quietly = TRUE)) {
  install.packages("desc", repos = "https://cloud.r-project.org")
}
library(desc)

# Load DESCRIPTION
d <- desc::desc(file = "DESCRIPTION")

# Extract all dependencies
deps <- d$get_deps()

# Drop R itself from the list (will add separately)
deps <- deps[deps$package != "R", ]

# Known Bioconductor packages (expand if needed)
bioc_pkgs <- c(
  "GenomicRanges", "IRanges", "BSgenome", "ChIPseeker", "qsea",
  "BiocGenerics", "Biostrings", "ComplexHeatmap", "GenomeInfoDb",
  "GenomicAlignments", "Rsamtools", "BiocParallel", "biomaRt", "HMMcopy",
  "MEDIPS", "rtracklayer", "workflows", "MEDIPSData",
  "org.Hs.eg.db", "org.Mm.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "BSgenome.Hsapiens.NCBI.GRCh38", "BSgenome.Hsapiens.UCSC.hg19"
)

# Map to conda package names
deps$conda_name <- ifelse(deps$package %in% bioc_pkgs,
                          paste0("bioconductor-", tolower(deps$package)),
                          paste0("r-", tolower(deps$package)))

# Deduplicate & sort alphabetically (keep r-base pinned first)
deps_sorted <- sort(unique(deps$conda_name))
deps_final <- c("r-base=4.3", deps_sorted)

# Write to package_list.txt
writeLines(deps_final, "package_list.txt")

cat("✅ package_list.txt generated with", length(deps_final), "packages.\n")
