FROM continuumio/miniconda3

ARG DEBIAN_FRONTEND=noninteractive
# Variant selector: "slim" (default) or "full"
ARG VARIANT=slim

# System dependencies (trim as needed)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential ca-certificates curl git sudo \
    libbz2-dev libfontconfig1-dev libfreetype6-dev libfribidi-dev \
    libharfbuzz-dev libicu-dev libjpeg-dev liblzma-dev \
    libncurses5-dev libncursesw5-dev libpng-dev libssl-dev \
    libtiff-dev libtiff5-dev libwebp-dev libxml2-dev \
    libcurl4-openssl-dev zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

# Conda + mamba; conda-forge/bioconda only
RUN conda update -n base -c defaults conda -y \
 && conda install -n base -c conda-forge mamba -y \
 && conda config --system --remove-key channels || true \
 && conda config --system --remove-key default_channels || true \
 && conda config --system --add channels conda-forge \
 && conda config --system --add channels bioconda \
 && conda config --system --set channel_priority strict \
 && conda config --system --set add_pip_as_python_dependency false

# Bring in both lists; select based on VARIANT
COPY .devcontainer/conda/package_list.slim.txt /opt/package_list.slim.txt
COPY .devcontainer/conda/package_list.full.txt /opt/package_list.full.txt
RUN cp /opt/package_list.${VARIANT}.txt /opt/package_list.txt \
 && echo "Using package list: /opt/package_list.${VARIANT}.txt" && head -n 40 /opt/package_list.txt

# Create env
RUN mamba create -vv -y -n mesa_env \
      --override-channels -c conda-forge -c bioconda \
      --file /opt/package_list.txt \
 && conda clean -afy && mamba clean -afy

# Dev tools (both variants): devtools + helpers via conda
RUN mamba install -y -n mesa_env -c conda-forge \
      r-devtools r-roxygen2 r-rcmdcheck r-pkgload r-usethis \
 && conda clean -afy && mamba clean -afy

# Put env on PATH
ENV PATH="/opt/conda/envs/mesa_env/bin:${PATH}"
ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

# Extra R bits not always great via conda
RUN R -e "install.packages(c('remotes','languageserver','imsig'), repos='https://cloud.r-project.org')" \
 && R -e "remotes::install_github('omnideconv/immunedeconv')"

# Non-root user for Codespaces
RUN useradd -m -s /bin/bash vscode \
 && echo 'vscode ALL=(ALL) NOPASSWD:ALL' >/etc/sudoers.d/99-vscode \
 && chmod 0440 /etc/sudoers.d/99-vscode
USER vscode
WORKDIR /workspaces/mesa

CMD ["R"]

