# Start from official Rocker R image
FROM rocker/r-ver:4.3.3

# Install system dependencies for R packages and Git
RUN apt-get update && apt-get install -y --no-install-recommends \
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


# Create vscode user (UID 1000) with sudo privileges
RUN useradd -ms /bin/bash vscode \
    && echo "vscode ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers.d/90-vscode

# Preinstall bootstrap packages globally (needed before user lib init)
RUN R -e "install.packages(c('jsonlite','BiocManager','languageserver','remotes'), repos='https://cloud.r-project.org')"

# Switch to vscode user
USER vscode
WORKDIR /workspaces/mesa
