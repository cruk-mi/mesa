FROM continuumio/miniconda3

ARG DEBIAN_FRONTEND=noninteractive

# System deps (trim or extend as you need)
# Install system dependencies
RUN apt-get update && apt-get install -y -q \
    libnss-ldapd \
    acl \
    libnetcdf-dev \
    gcc-multilib \
    build-essential \
    ca-certificates \
    curl \
    git \
    libbz2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libicu-dev \
    libjpeg-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libpng-dev \
    libssl-dev \
    libtiff-dev \
    libtiff5-dev \
    libwebp-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    sudo \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


# --- Conda + mamba; only conda-forge + bioconda (no defaults)
RUN conda update -n base -c defaults conda -y \
 && conda install -n base -c conda-forge mamba -y \
 && conda config --system --remove-key channels || true \
 && conda config --system --remove-key default_channels || true \
 && conda config --system --add channels conda-forge \
 && conda config --system --add channels bioconda \
 && conda config --system --set channel_priority strict \
 && conda config --system --set add_pip_as_python_dependency false \
 && conda config --show-sources && conda config --show channels

# --- Bring in the conda package list generated from DESCRIPTION
COPY package_list.txt /opt/package_list.txt

# --- Create the env using ONLY CF+Bioconda; clean caches
RUN ls -l /opt && head -n 40 /opt/package_list.txt \
 && mamba create -vv -y -n shiny_env \
      --override-channels -c conda-forge -c bioconda \
      --file /opt/package_list.txt \
 && conda clean -afy && mamba clean -afy

# --- Put the env on PATH so R is available without activation
ENV PATH="/opt/conda/envs/shiny_env/bin:${PATH}"
ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

# --- Install extra R packages (inside shiny_env)
#     - CRAN 'imsig'
#     - 'remotes' (for postCreate install_local)
#     - 'languageserver' for IDE features
#     - GitHub 'omnideconv/immunedeconv'
RUN R -e "install.packages(c('imsig','remotes','languageserver'), repos='https://cloud.r-project.org')" \
 && R -e "remotes::install_github('omnideconv/immunedeconv')"

# --- Create an unprivileged user 'vscode' (Codespaces will use this)
RUN useradd -m -s /bin/bash vscode \
 && echo 'vscode ALL=(ALL) NOPASSWD:ALL' >/etc/sudoers.d/99-vscode \
 && chmod 0440 /etc/sudoers.d/99-vscode

# (Optional) a friendly workdir; Codespaces will mount the repo here
WORKDIR /workspaces/mesa

CMD ["R"]
