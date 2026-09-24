#!/usr/bin/env Rscript
# One-shot local setup of mesa's check ladder.
#
# No version is chosen here. `.devcontainer/resolve_versions.sh` is the one
# place that maps DESCRIPTION's `R (>= X.Y.Z)` to a Bioconductor release, and
# DESCRIPTION also records the roxygen2 that generated the committed man/
# pages. This script asks the resolver and installs what it says — hardcoding
# a version here would be exactly the duplicate that the toolchain doctrine in
# AGENTS.md exists to prevent.
#
# Assumes the resolved R is already the default (`rig default <version>`) —
# installing the framework itself needs sudo and is the one step this cannot do.
#
# Everything installs as binaries into the user library, so nothing here needs
# a compiler. Run from the repo root:
#   Rscript .claude/scripts/setup-r-toolchain.R

say <- function(...) cat("==>", ..., "\n")

# --- 0. versions, from the single source of truth ------------------------
if (!file.exists("DESCRIPTION")) {
    stop("run this from the repo root — no DESCRIPTION here.", call. = FALSE)
}
resolver <- file.path(".devcontainer", "resolve_versions.sh")
if (!file.exists(resolver)) stop("missing ", resolver, call. = FALSE)

resolved <- system2("bash", shQuote(resolver), stdout = TRUE)
status <- attr(resolved, "status")
if (!is.null(status) && status != 0L) {
    stop("resolve_versions.sh failed — see its message above.", call. = FALSE)
}
pairs <- strsplit(resolved, "=", fixed = TRUE)
vers <- setNames(
    vapply(pairs, function(p) paste(p[-1], collapse = "="), character(1)),
    vapply(pairs, `[`, character(1), 1L)
)
need <- function(key) {
    if (!key %in% names(vers) || !nzchar(vers[[key]])) {
        stop("resolver did not report ", key, call. = FALSE)
    }
    vers[[key]]
}

BIOC <- need("BIOC_VERSION")
R_WANT <- need("R_VERSION_FULL")
# ROXYGEN_VERSION is a later addition to the resolver; a branch predating it
# resolves fine and only skips the roxygen2 pin below.
ROXYGEN <- if ("ROXYGEN_VERSION" %in% names(vers)) vers[["ROXYGEN_VERSION"]] else ""
say("resolved: R", R_WANT, "| Bioconductor", BIOC, "| roxygen2",
    if (nzchar(ROXYGEN)) ROXYGEN else "(not pinned on this branch)")

# --- 1. right R? ---------------------------------------------------------
have <- paste(R.version$major, R.version$minor, sep = ".")
say("R", have, "at", R.home())
if (have != R_WANT) {
    stop("expected R ", R_WANT, " but this is ", have,
         ". Run `rig default ", R_WANT, "` first.", call. = FALSE)
}

# --- 2. a writable user library, chosen rather than prompted for ---------
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

# --- 3. binary repos: the same CRAN mirror CI uses -----------------------
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

# --- 4. mesa's dependencies ---------------------------------------------
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

# --- 5. the check tooling CI uses ---------------------------------------
# rcmdcheck and BiocCheck are what CI itself runs; devtools and covr complete
# the ladder locally (devtools::test, covr::package_coverage).
for (pkg in c("rcmdcheck", "BiocCheck", "remotes", "devtools", "covr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg, ask = FALSE, update = FALSE, lib = lib)
    }
}

# roxygen2 is pinned, not merely present: DESCRIPTION records the version that
# generated the committed man/ and NAMESPACE, and regenerating with a different
# one rewrites every page. Two details matter here. It comes after devtools,
# which depends on roxygen2 and would otherwise pull the current release over
# the pin; and presence is not enough, so the version is compared — the same
# reasoning as the roxygen2 block in .devcontainer/install.R.
if (!nzchar(ROXYGEN)) {
    say("no roxygen2 version resolved — skipping the pin.",
        "Do not regenerate man/ from this machine.")
} else if (!requireNamespace("roxygen2", quietly = TRUE) ||
           !identical(as.character(packageVersion("roxygen2")), ROXYGEN)) {
    say("installing pinned roxygen2", ROXYGEN)
    remotes::install_version("roxygen2", version = ROXYGEN,
                             upgrade = "never", lib = lib)
}

# --- 6. report -----------------------------------------------------------
say("BiocManager::valid():")
print(BiocManager::valid())
if (length(failed)) {
    say("COULD NOT INSTALL:", paste(failed, collapse = ", "))
    say("These are declared dependencies (Depends/Imports/Suggests); the full",
        "check ladder cannot run until they install.")
} else {
    say("every dependency installed")
}
rox_have <- if (requireNamespace("roxygen2", quietly = TRUE)) {
    as.character(packageVersion("roxygen2"))
} else "absent"
say("toolchain: R", have, "| Bioconductor",
    as.character(BiocManager::version()), "| roxygen2", rox_have,
    if (nzchar(ROXYGEN) && !identical(rox_have, ROXYGEN)) {
        sprintf("(MISMATCH — man/ must be regenerated with %s only)", ROXYGEN)
    } else "")
say("next: R CMD build . && R CMD check --no-manual mesa_*.tar.gz")
# Fail last, after the toolchain line, so the summary is still printed.
if (length(failed)) quit(status = 1L)
