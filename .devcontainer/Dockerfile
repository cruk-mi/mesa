# Start from official Rocker R image
FROM rocker/r-ver:4.3.3

# Install system dependencies for R packages and Git
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    curl \
    ca-certificates \
    sudo \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff-dev \
    libicu-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfontconfig1-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create vscode user (UID 1000) with sudo privileges
RUN useradd -ms /bin/bash vscode \
    && echo "vscode ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers.d/90-vscode

# Preinstall bootstrap packages globally (needed before user lib init)
RUN R -e "install.packages(c('jsonlite','BiocManager','languageserver','remotes'), repos='https://cloud.r-project.org')"

# Switch to vscode user
USER vscode
WORKDIR /workspaces/mesa
