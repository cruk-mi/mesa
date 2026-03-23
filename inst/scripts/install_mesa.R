# inst/scripts/install_mesa.R
# ============================================================
# One-shot installer for MESA and all its dependencies.
#
# This script is version-agnostic — it automatically selects
# the correct Bioconductor version for your R installation.
# You do not need to update this script when R or Bioc releases.
#
# Usage (run once before first use):
#
#   From GitHub (recommended):
#   source("https://raw.githubusercontent.com/cruk-mi/mesa/main/inst/scripts/install_mesa.R")
#
#   Locally (if mesa is already installed):
#   source(system.file("scripts/install_mesa.R", package = "mesa"))
#
# Requirements:
#   R >= 4.3.0  (see https://cran.r-project.org for upgrade instructions)
# ============================================================

# --- Check minimum R version --------------------------------
# This is the oldest R version mesa has ever supported.
# It is intentionally conservative — users on newer R are fine.
minimum_r <- "4.3.0"

if (getRversion() < minimum_r) {
  stop(sprintf(
    "mesa requires R >= %s. You are running %s.\n%s",
    minimum_r,
    R.version$version.string,
    "Please upgrade R: https://cran.r-project.org"
  ), call. = FALSE)
}

message("── mesa installer ───────────────────────────────────────────")
message(sprintf("   R version : %s", R.version$version.string))

# --- Bootstrap BiocManager ----------------------------------
# BiocManager is the only package we need to install manually.
# Everything else flows from it.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("   Installing BiocManager...")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# --- Set correct Bioc version for this R --------------------
# BiocManager knows which Bioc release matches the user's R version.
# No version number needed here — it always picks the right one.
message(sprintf(
  "   Bioconductor: setting to latest release for R %s.%s...",
  R.version$major,
  substr(R.version$minor, 1, 1)
))
BiocManager::install(ask = FALSE)
message(sprintf("   Bioconductor: %s ✓", BiocManager::version()))

# --- Install remotes ----------------------------------------
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

# --- Install MESA -------------------------------------------
message("── Installing mesa from GitHub ─────────────────────────────")
remotes::install_github(
  "cruk-mi/mesa",
  dependencies = TRUE,  # installs Imports + Suggests automatically
  upgrade      = "never"
)

# --- Validate environment consistency -----------------------
# BiocManager::valid() checks that all installed packages are
# the correct versions for the active Bioc release.
# Out-of-date packages can cause subtle bugs — always fix them.
message("── Validating environment ──────────────────────────────────")
check <- BiocManager::valid(checkBuilt = TRUE)

if (isTRUE(check)) {
  message(sprintf(
    "✅ Environment consistent with Bioconductor %s",
    BiocManager::version()
  ))
} else {
  out_of_date <- names(check$out_of_date)
  too_new     <- names(check$too_new)

  if (length(out_of_date) > 0) {
    message(sprintf("   Updating %d out-of-date package(s)...", length(out_of_date)))
    BiocManager::install(out_of_date, update = TRUE, ask = FALSE)
  }

  if (length(too_new) > 0) {
    # "too new" means the user has a newer version than Bioc expects.
    # This usually happens when CRAN releases a patch between Bioc cycles.
    # It's rarely a real problem — just report it.
    message(sprintf(
      "   ℹ️  %d package(s) newer than Bioc %s expects: %s",
      length(too_new),
      BiocManager::version(),
      paste(too_new, collapse = ", ")
    ))
    message("   This is usually harmless. If mesa fails, try:")
    message(sprintf(
      "   BiocManager::install(c(%s), force = TRUE)",
      paste0('"', too_new, '"', collapse = ", ")
    ))
  }

  message("✅ Done. Run BiocManager::valid() to confirm.")
}

message(sprintf(
  "\n✅ mesa is ready. Load it with: library(mesa)\n"
))
