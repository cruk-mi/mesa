# inst/scripts/install_mesa.R
# ============================================================
# One-shot installer for MESA and all its dependencies.
#
# This script is version-agnostic — it automatically selects
# the correct Bioconductor version for your R installation.
# You do not need to update this script when R or Bioc releases.
#
# Works on: personal laptops, HPC shared systems, Docker containers.
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
suppressWarnings(BiocManager::install(ask = FALSE))
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

# Helper: identify which library paths are writable by the current user.
# Works universally — on a laptop all paths are typically writable;
# on HPC the system library is read-only and correctly excluded.
get_writable_libs <- function() {
  Filter(function(p) file.access(p, mode = 2) == 0, .libPaths())
}

# Helper: given a named list from BiocManager::valid()$out_of_date or
# $too_new, split packages into those the user can update (writable lib)
# and those they cannot (read-only system lib).
split_by_writability <- function(pkg_list) {
  writable <- get_writable_libs()
  can_update <- vapply(names(pkg_list), function(pkg) {
    pkg_lib <- pkg_list[[pkg]]["LibPath"]
    any(vapply(writable, function(lib) {
      grepl(lib, pkg_lib, fixed = TRUE)
    }, logical(1)))
  }, logical(1))
  list(
    user   = names(pkg_list)[can_update],
    system = names(pkg_list)[!can_update]
  )
}

fix_environment <- function(max_attempts = 3) {
  
  for (attempt in seq_len(max_attempts)) {
    
    # Suppress BiocManager's own warning — we report results ourselves
    check <- suppressWarnings(BiocManager::valid(checkBuilt = TRUE))
    
    # All good
    if (isTRUE(check)) {
      message(sprintf(
        "✅ Environment consistent with Bioconductor %s",
        BiocManager::version()
      ))
      return(invisible(TRUE))
    }
    
    # Split out-of-date packages by whether the user can update them
    ood_split <- split_by_writability(check$out_of_date)
    user_ood   <- ood_split$user
    system_ood <- ood_split$system
    
    # Report read-only packages clearly — not a user problem
    # This handles HPC system libraries gracefully without
    # confusing laptop users (who will never hit this case)
    if (length(system_ood) > 0) {
      message(sprintf(
        "   ℹ️  %d package(s) out of date in a read-only library (cannot update): %s\n      This is expected on shared systems (HPC, managed desktops).",
        length(system_ood),
        paste(system_ood, collapse = ", ")
      ))
    }
    
    # Report "too new" packages — usually cosmetic
    # BiocManager itself is a common false positive: CRAN releases
    # patch versions faster than Bioc cycles, harmless to ignore
    too_new_split <- split_by_writability(check$too_new)
    real_too_new  <- setdiff(too_new_split$user, "BiocManager")
    if (length(real_too_new) > 0) {
      message(sprintf(
        "   ℹ️  %d package(s) newer than Bioc %s expects (usually harmless): %s",
        length(real_too_new),
        BiocManager::version(),
        paste(real_too_new, collapse = ", ")
      ))
    }
    
    # Nothing the user can act on — environment is effectively clean
    if (length(user_ood) == 0) {
      message(sprintf(
        "✅ Environment consistent with Bioconductor %s",
        BiocManager::version()
      ))
      return(invisible(TRUE))
    }
    
    # ── Special case: BiocManager updating itself ─────────────────
    # BiocManager installs a newer version to disk but the running
    # session keeps the old version in memory until R is restarted.
    if (identical(user_ood, "BiocManager")) {
      message("   Updating BiocManager...")
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
      message(
        "   ℹ️  BiocManager updated on disk. Please restart R and run\n",
        "      BiocManager::valid() to confirm. Everything else is installed."
      )
      return(invisible(TRUE))
    }
    
    # ── General out-of-date packages in user library ───────────────
    message(sprintf(
      "   Attempt %d/%d: updating %d out-of-date package(s): %s",
      attempt, max_attempts,
      length(user_ood),
      paste(user_ood, collapse = ", ")
    ))
    suppressWarnings(
      BiocManager::install(user_ood, update = TRUE, ask = FALSE)
    )
  }
  
  # Max attempts reached — report anything still outstanding
  remaining_check <- suppressWarnings(BiocManager::valid(checkBuilt = TRUE))
  if (!isTRUE(remaining_check)) {
    remaining_ood <- split_by_writability(remaining_check$out_of_date)$user
    remaining     <- setdiff(remaining_ood, "BiocManager")
    if (length(remaining) > 0) {
      message(sprintf(
        "⚠️  %d package(s) still out of date after %d attempts: %s\n%s",
        length(remaining),
        max_attempts,
        paste(remaining, collapse = ", "),
        "   Run BiocManager::install(ask=FALSE, update=TRUE) to fix manually."
      ))
    } else {
      # Only system packages or BiocManager left — effectively clean
      message(sprintf(
        "✅ Environment consistent with Bioconductor %s",
        BiocManager::version()
      ))
    }
  }
}

fix_environment()

message(sprintf("\n✅ mesa is ready. Load it with: library(mesa)\n"))