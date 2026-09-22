#!/usr/bin/env Rscript
# One-shot local setup of mesa's check ladder: Bioconductor 3.23 on R 4.6.0.
#
# Assumes R 4.6.0 is already the default (`rig default 4.6.0`) — installing the
# framework itself needs sudo and is the one step this cannot do.
#
# Everything installs as binaries into the user library, so nothing here needs
# a compiler. Run:  Rscript .claude/scripts/setup-r-toolchain.R

BIOC <- "3.23"
R_WANT <- "4.6.0"

say <- function(...) cat("==>", ..., "\n")

# --- 0. right R? ---------------------------------------------------------
have <- paste(R.version$major, R.version$minor, sep = ".")
say("R", have, "at", R.home())
if (have != R_WANT) {
    stop("expected R ", R_WANT, " but this is ", have,
         ". Run `rig default ", R_WANT, "` first.", call. = FALSE)
}

# --- 1. a writable user library, chosen rather than prompted for ---------
lib <- Sys.getenv("R_LIBS_USER")
lib <- strsplit(lib, .Platform$path.sep)[[1]][1]
if (!nzchar(lib) || is.na(lib)) {
    lib <- file.path("~", "Library", "R", R.version$arch,
                     paste(R.version$major, substr(R.version$minor, 1, 1), sep = "."),
                     "library")
}
lib <- path.expand(lib)
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
say("user library:", lib)

# --- 2. binary repos: the same CRAN mirror CI uses -----------------------
options(
    repos = c(CRAN = "https://p3m.dev/cran/latest"),
    pkgType = "binary",
    install.packages.check.source = "no",
    Ncpus = max(1L, parallel::detectCores() - 1L)
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    say("installing BiocManager")
    install.packages("BiocManager", lib = lib)
}

if (!identical(as.character(BiocManager::version()), BIOC)) {
    say("switching Bioconductor to", BIOC)
    BiocManager::install(version = BIOC, ask = FALSE, update = FALSE, lib = lib)
}
stopifnot(identical(as.character(BiocManager::version()), BIOC))
say("Bioconductor", BIOC)
options(repos = BiocManager::repositories())

# --- 3. mesa's dependencies ---------------------------------------------
desc <- read.dcf("DESCRIPTION")
fields <- c("Depends", "Imports", "Suggests")
deps <- unlist(lapply(fields, function(f) {
    if (!f %in% colnames(desc)) return(character())
    parts <- strsplit(desc[1, f], ",")[[1]]
    trimws(sub("\\(.*", "", parts))
}))
deps <- setdiff(unique(deps[nzchar(deps)]), c("R", rownames(installed.packages())))
say(length(deps), "dependencies to install")

# One batch call: resolving the graph 40 times over is the slow way, and the
# annotation packages are large enough that download order matters more than
# per-package error handling. Anything that did not land is reported after.
if (length(deps)) {
    try(BiocManager::install(deps, ask = FALSE, update = FALSE, lib = lib))
}
failed <- setdiff(deps, rownames(installed.packages()))

# Bioconductor ships annotation and experiment-data packages source-only, so a
# binary-only pass leaves them behind with "not available as a binary package".
# They are pure data - a source install of these needs no compiler, unlike a
# source install of R code - so retry exactly those as source.
if (length(failed)) {
    say("retrying", length(failed), "as source (annotation data ships source-only):",
        paste(failed, collapse = ", "))
    try(BiocManager::install(failed, type = "source", ask = FALSE, update = FALSE, lib = lib))
    failed <- setdiff(deps, rownames(installed.packages()))
}

# --- 4. the check tooling CI uses ---------------------------------------
for (pkg in c("rcmdcheck", "BiocCheck", "remotes")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg, ask = FALSE, update = FALSE, lib = lib)
    }
}

# --- 5. report -----------------------------------------------------------
say("BiocManager::valid():")
print(BiocManager::valid())
if (length(failed)) {
    say("COULD NOT INSTALL:", paste(failed, collapse = ", "))
    say("These are Suggests-or-Imports gaps; note them, do not paper over them.")
} else {
    say("every dependency installed")
}
say("next: R CMD build . && R CMD check --no-manual mesa_*.tar.gz")
