# check=skip=InvalidDefaultArgInFrom
# ============================================================
# Devcontainer image for mesa — slim and full variants.
#
# Based on the official Bioconductor Docker image, the SAME image
# family CI uses (.github/workflows/check-bioc.yml). R, the core
# Bioconductor stack, pandoc and every system library for Bioc
# packages come from the base image — no conda.
#
# The image tag is derived from versions, never hardcoded:
#   BIOC_RELEASE = RELEASE_<BIOC_VERSION with . -> _>
# is computed by .devcontainer/resolve_versions.sh and passed in as
# a build-arg by .github/workflows/build-image.yml.
# ============================================================

# ── Single source of truth (passed in by the resolver / workflow) ──────────────
# ARG before FROM is only usable in FROM; re-declare after FROM to reuse below.
# No default: any hardcoded value would need updating on each Bioc upgrade,
# breaking the single-file (DESCRIPTION-only) upgrade path. The workflow
# always supplies BIOC_RELEASE via --build-arg.
ARG BIOC_RELEASE
FROM bioconductor/bioconductor_docker:${BIOC_RELEASE}

# Variant selector: "slim" (default) or "full"
ARG VARIANT=slim

# Versions, for use inside the image (install.R reads these).
ARG R_VERSION
ARG BIOC_VERSION
ENV R_VERSION=${R_VERSION}
ENV BIOC_VERSION=${BIOC_VERSION}
ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

# Install the dependency stack into the image. Only DESCRIPTION (the source of
# truth for *which* packages) and install.R are needed — copying just these
# keeps this layer cached unless they change. mesa itself is installed later,
# from the live workspace, by devcontainer.json's postCreateCommand.
COPY DESCRIPTION /tmp/mesa/DESCRIPTION
COPY .devcontainer/install.R /tmp/install.R
RUN MESA_VARIANT=${VARIANT} MESA_PKG_DIR=/tmp/mesa Rscript /tmp/install.R \
 && rm -rf /tmp/downloaded_packages /tmp/Rtmp*

# Non-root user for Codespaces. Builds above ran as root and installed into
# /usr/local/lib/R/site-library; hand that library to vscode so the
# postCreateCommand (run as vscode) can install mesa into it.
RUN useradd -m -s /bin/bash vscode \
 && echo 'vscode ALL=(ALL) NOPASSWD:ALL' >/etc/sudoers.d/99-vscode \
 && chmod 0440 /etc/sudoers.d/99-vscode \
 && chown -R vscode /usr/local/lib/R/site-library

USER vscode
WORKDIR /workspaces/mesa
CMD ["R"]
