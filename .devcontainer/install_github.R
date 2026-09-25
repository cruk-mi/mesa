# ============================================================
# install_github.R
# Purpose: Install the GitHub-only packages into the devcontainer
#          image, after install.R has installed the main stack.
#
# Kept separate from install.R so that the GitHub API token is
# only ever mounted into this build step (see the Dockerfile).
# Anonymous GitHub API calls are capped at 60/hour per IP, which
# shared CI runners exhaust; immunedeconv alone resolves five more
# GitHub repos through its `Remotes:` field.
#
# The token is read from the BuildKit secret and set as GITHUB_PAT
# for this process only; remotes reads that variable for the
# `Remotes:` dependencies, which ignore an explicit auth_token. It
# never lands in an image layer and is never printed. Only this
# step can see it: the hundreds of CRAN/Bioc packages are installed
# by install.R in an earlier step with no secret mounted. Without
# the secret (e.g. a local `docker build`) the installs fall back
# to anonymous requests.
# ============================================================

gh_pat_file <- "/run/secrets/github_pat"
gh_token <- if (file.exists(gh_pat_file)) {
  trimws(readLines(gh_pat_file, n = 1L, warn = FALSE))
} else {
  ""
}
if (nzchar(gh_token)) Sys.setenv(GITHUB_PAT = gh_token)
message(sprintf("── GitHub API: %s ──",
                if (nzchar(gh_token)) "authenticated" else "anonymous"))
rm(gh_token)

# Any dependency not already installed by install.R comes from the
# image's CRAN mirror plus the Bioc repos for the pinned Bioc version.
options(repos = BiocManager::repositories())

# ggtree dev version (needs ggplot2 >= 4.0.0); no formal releases on GitHub,
# pinned to a specific SHA for reproducibility.
# immunedeconv is only available from GitHub; pinned to the SHA for v2.1.4.
pinned <- c(
  "YuLab-SMU/ggtree"        = "9f645a2b89e4150d9748547b3ea1b03906275c27",
  "omnideconv/immunedeconv" = "e625e6c28ed14a30f9f40f159925cc9f0df4fa49"
)
for (repo in names(pinned)) {
  remotes::install_github(repo, ref = pinned[[repo]], upgrade = "never")
}

# remotes only warns when an install fails (R_REMOTES_NO_ERRORS_FROM_WARNINGS
# is set in the Dockerfile), so a failed install -- e.g. a dependency download
# timing out -- would still pass the build and ship an image without the
# package. Fail the build unless each package loads at its pinned SHA.
ok <- vapply(names(pinned), function(repo) {
  pkg <- basename(repo)
  requireNamespace(pkg, quietly = TRUE) &&
    identical(packageDescription(pkg)$RemoteSha, pinned[[repo]])
}, logical(1))
if (!all(ok)) {
  stop("install_github.R: not installed at the pinned SHA: ",
       paste(names(pinned)[!ok], collapse = ", "),
       ". See the install log above for the cause.")
}

message("✅ GitHub-only packages ready")
