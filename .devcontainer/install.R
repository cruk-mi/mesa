# ============================================================
# install.R
# Purpose: Install mesa's dependency stack into the devcontainer
#          image, on top of the Bioconductor base image.
#
# The base image already provides R, the Bioc core stack, pandoc and
# all system libraries, so this script only installs:
#   1. mesa's declared dependencies, read from DESCRIPTION
#      (DESCRIPTION is the single source of truth for *which*
#      packages; already-present base packages are no-ops).
#   2. A small set of genuine extras the base image / DESCRIPTION
#      do not cover (dev versions + IDE tooling).
#
# The slim/full split is one explicit list: `full_only` below.
#
# DO NOT hardcode R or Bioconductor version numbers here — the Bioc
# version is read from the BIOC_VERSION env var (set by the Dockerfile
# from the resolver). mesa itself is NOT installed here; it is
# installed from the live workspace by devcontainer.json.
# ============================================================

# --- Read context from environment --------------------------
bioc_ver <- Sys.getenv("BIOC_VERSION")
variant  <- Sys.getenv("MESA_VARIANT", unset = "slim")
pkg_dir  <- Sys.getenv("MESA_PKG_DIR", unset = ".")

message(sprintf("── Environment: Bioc %s  |  Variant %s ──", bioc_ver, variant))

# --- Bootstrap remotes + BiocManager (usually present in base) ----
options(repos = c(CRAN = "https://cloud.r-project.org"),
        install.packages.check.source = "no")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Re-assert the Bioc version (the base image already pins it; this is a
# safe no-op when they match, and honours an override if one was set).
if (nzchar(bioc_ver)) {
  BiocManager::install(version = bioc_ver, ask = FALSE, update = FALSE)
}
message(sprintf("✅ Bioconductor version: %s", BiocManager::version()))

# --- Package repositories -----------------------------------
# Resolve both CRAN and Bioconductor packages. Posit Package Manager
# (RSPM) provides pre-built Linux binaries first — critical for fast
# builds (avoids compiling from source) — then the Bioc + CRAN repos
# for the pinned Bioc version so every declared dependency is found.
options(repos = c(
  RSPM = "https://packagemanager.posit.co/cran/__linux__/jammy/latest",
  BiocManager::repositories()
))

# --- Dependencies, derived from DESCRIPTION -----------------
deps <- remotes::dev_package_deps(pkg_dir, dependencies = TRUE)$package

# The ONLY difference between slim and full: these heavy genome /
# annotation data packages (declared in DESCRIPTION Suggests). slim
# omits them; the test suite skips the cases that need them via
# skip_if_not_installed() / skip_long_checks().
full_only <- c(
  "BSgenome.Hsapiens.NCBI.GRCh38",
  "BSgenome.Hsapiens.UCSC.hg19",
  "BSgenome.Scerevisiae.UCSC.sacCer3",
  "MEDIPSData",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "TxDb.Mmusculus.UCSC.mm10.knownGene"
)
if (!identical(variant, "full")) {
  deps <- setdiff(deps, full_only)
}

message(sprintf("── Installing %d declared dependencies ──", length(deps)))
BiocManager::install(deps, ask = FALSE, update = FALSE)

# --- Genuine extras (not in DESCRIPTION) --------------------
# IDE tooling and dev-only / GitHub-only packages.
for (pkg in c("languageserver", "imsig")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
for (pkg in c("devtools", "roxygen2", "rcmdcheck")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

# ggtree dev version (needs ggplot2 >= 4.0.0) and immunedeconv from GitHub.
remotes::install_github("YuLab-SMU/ggtree", upgrade = "never")
remotes::install_github("omnideconv/immunedeconv", upgrade = "never")

message("✅ mesa dependency stack ready")
