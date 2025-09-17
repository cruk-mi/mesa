# Start from official Rocker R image
FROM rocker/r-ver:4.3.3

# Install system dependencies for R packages and Git
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    curl \
    ca-certificates \
    sudo \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfontconfig1-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Preinstall bootstrap R packages (for smoother Codespace startup)
RUN R -e "install.packages(c('jsonlite','BiocManager'), repos='https://cloud.r-project.org')"
