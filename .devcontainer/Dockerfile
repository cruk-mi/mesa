FROM continuumio/miniconda3
ARG DEBIAN_FRONTEND=noninteractive

# Variant selector: "slim" (default) or "full"
ARG VARIANT=slim

# ── Single source of truth ───────────────────────────────────────────────────
# Passed in from versions.env via build-args in build-image.yml and
# devcontainer.json. Change versions.env only — never edit these defaults here.
ARG R_VERSION=4.5
ARG BIOC_VERSION=3.22
# ─────────────────────────────────────────────────────────────────────────────

# System dependencies (R-version agnostic — rarely need changing)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential ca-certificates curl git sudo \
    libbz2-dev libfontconfig1-dev libfreetype6-dev libfribidi-dev \
    libharfbuzz-dev libicu-dev libjpeg-dev liblzma-dev \
    libncurses5-dev libncursesw5-dev libpng-dev libssl-dev \
    libtiff-dev libtiff5-dev libwebp-dev libxml2-dev \
    libcurl4-openssl-dev zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

# Conda + mamba; conda-forge/bioconda only
RUN conda update -n base -c conda-forge conda -y \
 && conda install -n base -c conda-forge mamba -y \
 && conda config --system --remove-key channels || true \
 && conda config --system --remove-key default_channels || true \
 && conda config --system --add channels bioconda \
 && conda config --system --add channels conda-forge \
 && conda config --system --add default_channels conda-forge \
 && conda config --system --set channel_priority strict \
 && conda config --system --set add_pip_as_python_dependency false

# Copy both package lists
COPY .devcontainer/conda/package_list.slim.txt /opt/package_list.slim.txt
COPY .devcontainer/conda/package_list.full.txt /opt/package_list.full.txt
RUN cp /opt/package_list.${VARIANT}.txt /opt/package_list.txt \
 && echo "Using package list: /opt/package_list.${VARIANT}.txt" && head -n 40 /opt/package_list.txt

# Create conda env — r-base version injected from ARG
# Note: r-base is NOT in the package list files (see comment there)
RUN cp /opt/package_list.${VARIANT}.txt /opt/package_list.txt \
 && echo "Building ${VARIANT} image with R=${R_VERSION}, Bioc=${BIOC_VERSION}" \
 && mamba create -vv -y -n mesa_env \
      --override-channels -c conda-forge -c bioconda \
      --file /opt/package_list.txt \
      r-base=${R_VERSION} \
 && conda clean -afy && mamba clean -afy

# Write R and Bioc versions into the environment activation script
# so they're always available as env vars inside the container
RUN mkdir -p /opt/conda/envs/mesa_env/etc/conda/activate.d \
 && printf 'export R_VERSION=%s\nexport BIOC_VERSION=%s\n' \
      "${R_VERSION}" "${BIOC_VERSION}" \
      > /opt/conda/envs/mesa_env/etc/conda/activate.d/mesa_env_vars.sh
      
# Dev tools via conda (both variants)
RUN mamba install -y -n mesa_env -c conda-forge \
      r-devtools r-roxygen2 r-rcmdcheck r-pkgload r-usethis \
 && conda clean -afy && mamba clean -afy

# Put env on PATH
ENV PATH="/opt/conda/envs/mesa_env/bin:${PATH}"
ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

# Extra R bits not reliably available via conda
RUN R -e "install.packages(c('remotes','languageserver','imsig'), \
            repos='https://cloud.r-project.org')" \
 && R -e "remotes::install_github('omnideconv/immunedeconv')"

# Non-root user for Codespaces
RUN useradd -m -s /bin/bash vscode \
 && echo 'vscode ALL=(ALL) NOPASSWD:ALL' >/etc/sudoers.d/99-vscode \
 && chmod 0440 /etc/sudoers.d/99-vscode

USER vscode
WORKDIR /workspaces/mesa
CMD ["R"]

