FROM continuumio/miniconda3

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
    
# Update conda
RUN conda update -n base -c defaults conda -y

# Install mamba
RUN conda install -n base -c conda-forge mamba -y

# Copy package list
COPY package_list.txt /opt/package_list.txt

# Create conda environment with all packages (Imports + Depends + Suggests)
RUN mamba create -y -n shiny_env -c conda-forge -c bioconda --file /opt/package_list.txt

# Activate env and install extra CRAN/GitHub-only packages
RUN /bin/bash -c "source activate shiny_env && R -e \"install.packages('imsig', repos='https://cloud.r-project.org')\""
RUN /bin/bash -c "source activate shiny_env && R -e \"remotes::install_github('omnideconv/immunedeconv')\""

# Set default environment for subsequent RUN/CMD
SHELL ["conda", "run", "-n", "shiny_env", "/bin/bash", "-c"]

CMD ["R"]
