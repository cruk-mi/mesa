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
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("   Installing BiocManager...")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# --- Set correct Bioc version for this R --------------------
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
  dependencies = TRUE,
  upgrade      = "never"
)

# --- Validate and fix environment ---------------------------
message("── Validating environment ──────────────────────────────────")

fix_environment <- function(max_attempts = 3) {
  
  for (attempt in seq_len(max_attempts)) {
    
    check <- BiocManager::valid(checkBuilt = TRUE)
    
    # All good
    if (isTRUE(check)) {
      message(sprintf(
        "✅ Environment consistent with Bioconductor %s",
        BiocManager::version()
      ))
      return(invisible(TRUE))
    }
    
    out_of_date <- names(check$out_of_date)
    too_new     <- names(check$too_new)
    
    # ── Handle "too new" ──────────────────────────────────
    # Happens when CRAN releases a patch between Bioc cycles.
    # BiocManager itself is a common false positive here —
    # CRAN version is often ahead of what Bioc expects.
    # This is cosmetic and never affects functionality.
    real_too_new <- setdiff(too_new, "BiocManager")
    if (length(real_too_new) > 0) {
      message(sprintf(
        "   ℹ️  %d package(s) newer than Bioc %s expects (usually harmless): %s",
        length(real_too_new),
        BiocManager::version(),
        paste(real_too_new, collapse = ", ")
      ))
    }
    
    # ── Handle out-of-date ────────────────────────────────
    if (length(out_of_date) == 0) {
      # Only BiocManager too_new was the issue — treat as clean
      message(sprintf(
        "✅ Environment consistent with Bioconductor %s",
        BiocManager::version()
      ))
      return(invisible(TRUE))
    }
    
    # ── Special case: BiocManager updating itself ─────────
    # BiocManager can install a newer version to disk but the
    # running session keeps the old version in memory.
    # We install it via install.packages() and tell the user
    # to restart — we cannot do more within the same session.
    biocmanager_only <- identical(out_of_date, "BiocManager")
    if (biocmanager_only) {
      message("   Updating BiocManager...")
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
      message(paste(
        "   ℹ️  BiocManager was updated on disk but the current session",
        "still holds the old version in memory.\n",
        "  Please restart R and run BiocManager::valid() to confirm.",
        "  Everything else is correctly installed."
      ))
      return(invisible(TRUE))
    }
    
    # ── General out-of-date packages ──────────────────────
    message(sprintf(
      "   Attempt %d/%d: updating %d out-of-date package(s): %s",
      attempt, max_attempts,
      length(out_of_date),
      paste(out_of_date, collapse = ", ")
    ))
    BiocManager::install(out_of_date, update = TRUE, ask = FALSE)
  }
  
  # If we get here, max attempts reached with packages still out of date
  remaining_check <- BiocManager::valid(checkBuilt = TRUE)
  if (!isTRUE(remaining_check)) {
    remaining <- setdiff(names(remaining_check$out_of_date), "BiocManager")
    if (length(remaining) > 0) {
      message(sprintf(
        "⚠️  %d package(s) still out of date after %d attempts: %s\n%s",
        length(remaining),
        max_attempts,
        paste(remaining, collapse = ", "),
        "   Run BiocManager::install(ask=FALSE, update=TRUE) to fix manually."
      ))
    }
  }
}

fix_environment()

message(sprintf("\n✅ mesa is ready. Load it with: library(mesa)\n"))