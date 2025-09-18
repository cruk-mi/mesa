# Ensure user library path exists
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(user_lib, .libPaths()))

options(repos = c(CRAN = "https://cloud.r-project.org"))

# --- Bootstrap remotes ---
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", lib = user_lib, repos = "https://cloud.r-project.org")
}

# --- Pin exact versions to avoid fgsea/ggtree issues ---
remotes::install_version(
  "BH",
  version = "1.81.0-1",
  lib = user_lib,
  repos = "https://cloud.r-project.org",
  type = "source"
)

remotes::install_version(
  "ggplot2",
  version = "3.4.4",
  lib = user_lib,
  repos = "https://cloud.r-project.org",
  type = "source"
)

# --- Install BiocManager if missing ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = user_lib, repos = "https://cloud.r-project.org")
}

# --- Install ggtree from Bioconductor 3.18, don't upgrade pinned packages ---
BiocManager::install(
  "ggtree",
  version = "3.18",
  ask = FALSE,
  update = FALSE,   # 👈 prevents BH/ggplot2 from being overwritten
  lib = user_lib
)

# --- Commented for now (extend later when stable) ---
# install.packages("devtools", lib = user_lib, repos = "https://cloud.r-project.org")
# BiocManager::install(c("qsea", "MEDIPS", "clusterProfiler"), lib = user_lib, ask = FALSE)

message("✅ Test environment ready with pinned BH 1.81.0-1 + ggplot2 3.4.4 + ggtree 3.10.1")

# # Ensure user library path exists
# user_lib <- Sys.getenv("R_LIBS_USER")
# if (!dir.exists(user_lib)) {
#   dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
# }
# .libPaths(c(user_lib, .libPaths()))
# 
# options(repos = c(CRAN = "https://cloud.r-project.org"))
# 
# # Bootstrap remotes
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes", lib = user_lib, repos = "https://cloud.r-project.org")
# }
# 
# # --- Pin critical versions to avoid build failures ---
# remotes::install_version("BH", version = "1.81.0-1", lib = user_lib, repos = "https://cloud.r-project.org")
# remotes::install_version("ggplot2", version = "3.4.4", lib = user_lib, repos = "https://cloud.r-project.org")
# 
# 
# # Install devtools (needed for devtools::test)
# install.packages("devtools", lib = user_lib, repos = "https://cloud.r-project.org")
# 
# # Core CRAN dependencies (avoid full tidyverse to skip Google API deps)
# install.packages(c(
#   "httr", "png", "RCurl", "igraph", "ggraph",
#   "rlang", "dplyr", "tibble", "tidyr", "readr",
#   "purrr", "stringr", "forcats",
#   "Rcpp", "RcppAnnoy", "RcppProgress",
#   "uwot"
# ), lib = user_lib)
# 
# # Bioconductor version 3.18
# BiocManager::install(version = "3.18", ask = FALSE)
# 
# # Core Bioconductor packages
# BiocManager::install(c(
#   "BiocGenerics", "GenomeInfoDb", "IRanges", "S4Vectors",
#   "Biostrings", "GenomicRanges", "Rhtslib", "Rsamtools",
#   "SummarizedExperiment"
# ), lib = user_lib, ask = FALSE, update = TRUE)
# 
# # Annotation + genomics
# BiocManager::install(c(
#   "AnnotationDbi", "GenomicFeatures", "GenomicAlignments", "BSgenome",
#   "TxDb.Hsapiens.UCSC.hg38.knownGene",
#   "TxDb.Mmusculus.UCSC.mm10.knownGene",
#   "org.Mm.eg.db"
# ), lib = user_lib, ask = FALSE, update = TRUE)
# 
# # High-level analysis packages
# BiocManager::install(c(
#   "qsea", "MEDIPS", "ChIPseeker", "clusterProfiler",
#   "DOSE", "GOSemSim", "enrichplot", "ReactomePA",
#   "ggtree", "treeio", "tidytree"
# ), lib = user_lib, ask = FALSE, update = TRUE)
# 
# message("✅ R dev environment ready with user library at: ", user_lib)
