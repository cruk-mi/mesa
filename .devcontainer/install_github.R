# ============================================================
# install_github.R
# Purpose: Install the GitHub-only packages into the devcontainer
#          image, after install.R has installed the main stack.
#
# Every package is pinned to a commit and fetched as a plain
# github.com archive, never through the GitHub API. Anonymous API
# calls are capped at 60/hour per IP, which shared CI runners
# exhaust; archive downloads are not API calls, so the build needs
# no token (and so exposes none to the packages it installs).
#
# immunedeconv declares its GitHub dependencies in `Remotes:`, which
# remotes would resolve through the API, at HEAD. So they are pinned
# here too and installed first, and each package is installed with
# `dependencies = FALSE` once its CRAN/Bioc dependencies (read from
# its DESCRIPTION) are in place -- `Remotes:` is never consulted.
# When bumping immunedeconv, re-check its `Remotes:` against the pins.
# ============================================================

# Dependencies come from the image's CRAN mirror plus the Bioc repos
# for the pinned Bioc version.
options(repos = BiocManager::repositories())

# Installed in this order: immunedeconv's `Remotes:` before immunedeconv.
pinned <- c(
  # ggtree dev version (needs ggplot2 >= 4.0.0); no formal releases.
  "YuLab-SMU/ggtree"        = "9f645a2b89e4150d9748547b3ea1b03906275c27",
  # immunedeconv's `Remotes:`, pinned to the commits current at v2.1.4.
  "dviraran/xCell"          = "20e2919eefd37e15af35f29f4944e30697098a28",
  "GfellerLab/EPIC"         = "50a4f404f96c2842b2891b517b4e3bfaa6c64b8f",
  "grst/MCPcounter"         = "7ea6e040af68d8e0641d80ddf6b0d06226752b87",
  "grst/mMCPcounter"        = "cec69714e0028a0b46f19b3192c48d26f8eaa177",
  "cansysbio/ConsensusTME"  = "6e14ba3d09e39d48b5b9a6b5fa9cd10296cce06a",
  # immunedeconv is only available from GitHub; v2.1.4.
  "omnideconv/immunedeconv" = "e625e6c28ed14a30f9f40f159925cc9f0df4fa49"
)
archive_url <- function(repo) {
  sprintf("https://github.com/%s/archive/%s.tar.gz", repo, pinned[[repo]])
}

for (repo in names(pinned)) {
  # CRAN/Bioc dependencies first, read from the pinned DESCRIPTION.
  desc_dir <- tempfile()
  dir.create(desc_dir)
  download.file(sprintf("https://raw.githubusercontent.com/%s/%s/DESCRIPTION",
                        repo, pinned[[repo]]),
                file.path(desc_dir, "DESCRIPTION"), quiet = TRUE)
  deps <- remotes::local_package_deps(desc_dir, dependencies = NA)
  missing <- setdiff(deps, c("R", rownames(installed.packages())))
  if (length(missing)) {
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  }
  remotes::install_url(archive_url(repo), dependencies = FALSE,
                       upgrade = "never")
}

# remotes only warns when an install fails (R_REMOTES_NO_ERRORS_FROM_WARNINGS
# is set in the Dockerfile), so a failed install -- e.g. a dependency download
# timing out -- would still pass the build and ship an image without the
# package. Fail the build unless each package loads from its pinned archive.
ok <- vapply(names(pinned), function(repo) {
  pkg <- basename(repo)
  requireNamespace(pkg, quietly = TRUE) &&
    identical(packageDescription(pkg)$RemoteUrl, archive_url(repo))
}, logical(1))
if (!all(ok)) {
  stop("install_github.R: not installed at the pinned SHA: ",
       paste(names(pinned)[!ok], collapse = ", "),
       ". See the install log above for the cause.")
}

message("✅ GitHub-only packages ready")
