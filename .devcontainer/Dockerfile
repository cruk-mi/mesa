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

# Conda + mamba + strict channel priority
RUN conda update -n base -c defaults conda -y \
 && conda install -n base -c conda-forge mamba -y \
 && conda config --remove-key channels || true \
 && conda config --add channels conda-forge \
 && conda config --add channels bioconda \
 && conda config --set channel_priority strict

# Copy package list (generated from DESCRIPTION)
COPY package_list.txt /opt/package_list.txt

# Sanity print, then create env (verbose to surface any missing names)
RUN ls -l /opt && head -n 40 /opt/package_list.txt \
 && mamba create -vv -y -n shiny_env --file /opt/package_list.txt \
 && conda clean -afy && mamba clean -afy

# From here on, run inside the env
SHELL ["conda", "run", "-n", "shiny_env", "/bin/bash", "-c"]
ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

# Extras not on conda (keep only what you need)
RUN R -e "install.packages('imsig', repos='https://cloud.r-project.org')" \
 && R -e "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='https://cloud.r-project.org')" \
 && R -e "remotes::install_github('omnideconv/immunedeconv')"

CMD ["R"]

