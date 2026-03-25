#!/usr/bin/env Rscript
# ============================================================
# check_bioc_compat.R
# Purpose: Before upgrading R, check all DESCRIPTION packages
#          are available in the target Bioconductor version.
#
# Usage:
#   Rscript .devcontainer/check_bioc_compat.R 3.22
#   Rscript .devcontainer/check_bioc_compat.R 3.23
# ============================================================

# --- Setup --------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: check_bioc_compat.R <bioc_version>  e.g. 3.21")
target_bioc <- args[1]

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("desc", quietly = TRUE)) {
  install.packages("desc", repos = "https://cloud.r-project.org")
}

# --- Read all packages from DESCRIPTION ---------------------
# desc package parses DESCRIPTION cleanly — no regex needed
d <- desc::desc("DESCRIPTION")

# Get packages from each field, strip version constraints like (>= 1.0)
get_pkgs <- function(field) {
  pkgs <- d$get_deps()
  pkgs <- pkgs[pkgs$type == field, "package"]
  pkgs <- pkgs[!pkgs %in% c("R")]   # exclude the R version constraint itself
  sort(unique(pkgs))
}

# Base and recommended packages ship with R — never in repos
# Get them dynamically rather than hardcoding a list
base_pkgs <- rownames(installed.packages(priority = c("base", "recommended")))

imports  <- setdiff(get_pkgs("Imports"),  base_pkgs)
depends  <- setdiff(get_pkgs("Depends"),  base_pkgs)
suggests <- setdiff(get_pkgs("Suggests"), base_pkgs)

all_pkgs <- list(
  Depends  = depends,
  Imports  = imports,
  Suggests = suggests
)

cat(sprintf("\n══ Checking DESCRIPTION packages for Bioc %s ══\n\n", target_bioc))
cat(sprintf("  Depends  : %d packages\n", length(depends)))
cat(sprintf("  Imports  : %d packages\n", length(imports)))
cat(sprintf("  Suggests : %d packages\n", length(suggests)))
cat(sprintf("  Total    : %d packages\n\n", length(unlist(all_pkgs))))

# --- Get available packages for target Bioc version ---------
# BiocManager knows which CRAN + Bioc repos correspond to each version
repos <- BiocManager::repositories(version = target_bioc)

# Download the master package list from each repo
# This is what install.packages() does internally — we're just reading it
available <- tryCatch({
  ap <- available.packages(repos = repos, fields = c("Package", "Version"))
  as.data.frame(ap[, c("Package", "Version")], stringsAsFactors = FALSE)
}, error = function(e) {
  stop("Could not fetch package list. Check your internet connection.\n", e$message)
})

cat(sprintf("  Fetched %d available packages from Bioc %s repos\n\n",
            nrow(available), target_bioc))

# --- Check each field ---------------------------------------
results <- list()

for (field in names(all_pkgs)) {
  pkgs <- all_pkgs[[field]]
  if (length(pkgs) == 0) next
  
  for (pkg in pkgs) {
    found   <- pkg %in% available$Package
    version <- if (found) available$Version[available$Package == pkg] else NA
    
    results[[length(results) + 1]] <- data.frame(
      field    = field,
      package  = pkg,
      found    = found,
      version  = ifelse(found, version, "NOT FOUND"),
      stringsAsFactors = FALSE
    )
  }
}

results_df <- do.call(rbind, results)

# --- Report -------------------------------------------------
missing <- results_df[!results_df$found, ]
found   <- results_df[results_df$found, ]

if (nrow(missing) == 0) {
  cat("✅ All packages available for Bioc", target_bioc, "\n\n")
} else {
  cat(sprintf("⚠️  %d package(s) NOT available for Bioc %s:\n\n", 
              nrow(missing), target_bioc))
  
  # Group by field so you know how critical each missing package is
  for (field in c("Depends", "Imports", "Suggests")) {
    m <- missing[missing$field == field, ]
    if (nrow(m) == 0) next
    
    # Severity label — Depends/Imports are blocking, Suggests are not
    severity <- if (field %in% c("Depends", "Imports")) "BLOCKING" else "non-blocking"
    cat(sprintf("  [%s — %s]\n", field, severity))
    for (pkg in m$package) cat(sprintf("    ✗ %s\n", pkg))
    cat("\n")
  }
}

# --- Summary table of what IS found (useful to sanity check versions) ---
cat(sprintf("── Found packages (%d) ─────────────────────────────────────\n",
            nrow(found)))
for (field in c("Depends", "Imports", "Suggests")) {
  f <- found[found$field == field, ]
  if (nrow(f) == 0) next
  cat(sprintf("\n  %s:\n", field))
  for (i in seq_len(nrow(f))) {
    cat(sprintf("    ✓ %-40s %s\n", f$package[i], f$version[i]))
  }
}

# --- Write CSV report for reference -------------------------
out_file <- sprintf("bioc_%s_compat_report.csv", gsub("\\.", "_", target_bioc))
write.csv(results_df, out_file, row.names = FALSE)
cat(sprintf("\n\nFull report written to: %s\n\n", out_file))

# --- Exit code — lets CI/CD catch failures ------------------
# Non-zero exit if any Depends/Imports packages are missing (blocking)
blocking_missing <- missing[missing$field %in% c("Depends", "Imports"), ]
if (nrow(blocking_missing) > 0) {
  quit(status = 1)  # signals failure to shell/CI
}
