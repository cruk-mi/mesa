# ============================================================
# install.R
# Purpose: Install R packages not managed by conda inside the
#          Codespaces container. Called once by postCreateCommand.
#
# DO NOT hardcode R or Bioconductor version numbers here.
# Both are read from environment variables set by the conda
# activation script (written by the Dockerfile from versions.env).
# ============================================================

# --- User library setup -------------------------------------
user_lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(user_lib)) {
  # Fallback: construct path dynamically if env var not set
  user_lib <- file.path(
    Sys.getenv("HOME"), "R",
    paste0(version[["arch"]], "-conda-linux-gnu-library"),
    paste0(version[["major"]], ".", substr(version[["minor"]], 1, 1))
  )
}
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(user_lib, .libPaths()))

# --- Read versions from environment -------------------------
# Set by Dockerfile via conda activation script
r_ver   <- Sys.getenv("R_VERSION",   unset = paste0(version$major, ".", substr(version$minor, 1, 1)))
bioc_ver <- Sys.getenv("BIOC_VERSION", unset = "3.22")

message(sprintf("── Environment: R %s  |  Bioc %s ──────────────────", r_ver, bioc_ver))
message(sprintf("── Installing to: %s", user_lib))

# --- Package repositories -----------------------------------
# Posit Package Manager provides pre-built Linux binaries —
# critical for fast builds (avoids compiling from source)
options(
  repos = c(
    RSPM = "https://packagemanager.posit.co/cran/__linux__/jammy/latest",
    CRAN = "https://cloud.r-project.org"
  ),
  install.packages.check.source = "no"
)
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

# --- Bootstrap remotes + BiocManager ------------------------
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", lib = user_lib)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = user_lib)
}

# --- Set Bioconductor version --------------------------------
# Version read from env var — no hardcoding needed
BiocManager::install(version = bioc_ver, ask = FALSE)
message(sprintf("✅ Bioconductor version set to: %s", BiocManager::version()))

# --- Install devtools ----------------------------------------
install.packages("devtools", lib = user_lib,
                 dependencies = TRUE,
                 INSTALL_opts = "--no-multiarch")

# --- Core CRAN packages -------------------------------------
install.packages(c(
  "httr", "png", "RCurl", "igraph", "ggraph",
  "rlang", "dplyr", "tibble", "tidyr", "readr",
  "purrr", "stringr", "forcats",
  "Rcpp", "RcppAnnoy", "RcppProgress", "uwot"
), lib = user_lib, update = FALSE)

# --- Critical Bioc packages ---------------------------------
BiocManager::install("fgsea", lib = user_lib, ask = FALSE, update = FALSE)
BiocManager::install("qsea",  lib = user_lib, ask = FALSE, update = FALSE)

# Core infrastructure
BiocManager::install(c(
  "BiocGenerics", "GenomeInfoDb", "IRanges", "S4Vectors",
  "Biostrings", "GenomicRanges", "Rhtslib", "Rsamtools",
  "SummarizedExperiment"
), lib = user_lib, ask = FALSE, update = FALSE)

# Annotation + genomics
BiocManager::install(c(
  "AnnotationDbi", "GenomicFeatures", "GenomicAlignments", "BSgenome",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "org.Mm.eg.db"
), lib = user_lib, ask = FALSE, update = FALSE)

# High-level analysis
BiocManager::install(c(
  "qsea", "MEDIPS", "ChIPseeker", "clusterProfiler",
  "DOSE", "GOSemSim", "enrichplot", "ReactomePA",
  "treeio", "tidytree"
), lib = user_lib, ask = FALSE, update = FALSE)

# --- Additional packages ------------------------------------
install.packages(c(
  "pheatmap", "janitor", "hues"
), lib = user_lib)

BiocManager::install(c(
  "plyranges", "ComplexHeatmap", "biomaRt", "circlize"
), lib = user_lib, ask = FALSE, update = FALSE)

# --- ggtree (dev version) -----------------------------------
# Requires ggplot2 >= 4.0.0 — confirmed available in Bioc 3.22
remotes::install_github("YuLab-SMU/ggtree", lib = user_lib, upgrade = "never")

message(sprintf("✅ Full R dev environment ready in: %s", user_lib))