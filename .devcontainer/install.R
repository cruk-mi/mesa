# Ensure user library path exists
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(user_lib, .libPaths()))

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  install.packages.check.source = "no" # don’t “update” to source if binary differs
)
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

# --- Bootstrap remotes + BiocManager ---
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", lib = user_lib, repos = "https://cloud.r-project.org")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = user_lib, repos = "https://cloud.r-project.org")
}

# --- Pin CRAN critical versions ---
remotes::install_version("BH", version = "1.81.0-1", lib = user_lib,
                         repos = "https://cloud.r-project.org", type = "source", upgrade = "never")
remotes::install_version("ggplot2", version = "4.0.0", lib = user_lib,
                         repos = "https://cloud.r-project.org", type = "source", upgrade = "never")
remotes::install_version("ggiraph", version = "0.9.1", lib = user_lib,
                         repos = "https://cloud.r-project.org", upgrade = "never")

# --- Install ggtree (dev version, requires pinned deps) ---
remotes::install_github("YuLab-SMU/ggtree", lib = user_lib, upgrade = "never")

# --- Install devtools (needed for devtools::test) ---
install.packages("devtools", lib = user_lib, repos = "https://cloud.r-project.org",
                 dependencies = TRUE, INSTALL_opts = "--no-multiarch")

message("✅ Core graphics stack ready: BH 1.81.0-1, ggplot2 4.0.0, ggiraph 0.9.1, ggtree dev")

# --- CRAN core deps (without full tidyverse) ---
install.packages(c(
  "httr", "png", "RCurl", "igraph", "ggraph",
  "rlang", "dplyr", "tibble", "tidyr", "readr",
  "purrr", "stringr", "forcats",
  "Rcpp", "RcppAnnoy", "RcppProgress", "uwot"
), lib = user_lib, update = FALSE)

# --- Sanity check before Bioconductor installs ---
ip <- installed.packages(lib.loc = user_lib)

required <- list(
  BH       = "1.81.0-1", # pinned for fgsea
  ggplot2  = "4.0.0",    # required for ggtree >= 3.99.0
  ggiraph  = "0.9.1",    # required for ggtree dev
  ggtree   = NA,         # dev version, skip version check
  devtools = NA          # latest ok
)

for (pkg in names(required)) {
  if (!pkg %in% rownames(ip)) {
    stop("❌ Package ", pkg, " is not installed in ", user_lib)
  }
  target_ver <- required[[pkg]]
  if (!is.na(target_ver)) {
    inst_ver <- ip[pkg, "Version"]
    if (inst_ver != target_ver) {
      stop("❌ Package ", pkg, " version mismatch: installed ", inst_ver,
           " but expected ", target_ver)
    }
  }
}


# --- Bioconductor 3.18 (set version only) ---
BiocManager::install(version = "3.18", ask = FALSE)

# --- Enforce BH before fgsea ---
ip <- installed.packages(lib.loc = user_lib)
if (!"BH" %in% rownames(ip) || ip["BH", "Version"] != "1.81.0-1") {
  message("⚠️ Reinstalling BH 1.81.0-1 to satisfy fgsea build...")
  remove.packages("BH", lib = user_lib)  # clean out any wrong version
  remotes::install_version("BH",
                           version = "1.81.0-1",
                           lib = user_lib,
                           repos = "https://cloud.r-project.org",
                           type = "source",
                           upgrade = "never"
  )
}

# --- Critical Bioc packages ---
BiocManager::install("fgsea", lib = user_lib, ask = FALSE, update = FALSE)
BiocManager::install("qsea",  lib = user_lib, ask = FALSE, update = FALSE)

# Core Bioconductor packages
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

# High-level analysis packages
BiocManager::install(c(
  "qsea", "MEDIPS", "ChIPseeker", "clusterProfiler",
  "DOSE", "GOSemSim", "enrichplot", "ReactomePA",
  "treeio", "tidytree"
), lib = user_lib, ask = FALSE, update = FALSE)

# --- Additional test/runtime dependencies ---
install.packages(c(
  "pheatmap", "janitor", "hues"
), lib = user_lib, repos = "https://cloud.r-project.org")

BiocManager::install(c(
  "plyranges", "ComplexHeatmap", "biomaRt", "circlize"
), lib = user_lib, ask = FALSE, update = FALSE)


message("✅ Full R dev environment ready in: ", user_lib)
