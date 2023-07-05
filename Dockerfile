FROM rocker/ml-verse

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

#install BiocManager and other packages required for the classifier pipeline
RUN install2.r --skipinstalled BiocManager tidymodels xgboost optparse mlbench vip cpp11 butcher && rm -rf /tmp/downloaded_packages

#install BiocParallel first because it was causing failures
RUN Rscript -e "BiocManager::install('BiocParallel')"

#install the package into the container
RUN Rscript -e "devtools::install('.', dependencies=TRUE, build_vignettes=TRUE, repos = BiocManager::repositories())" && rm -rf /tmp/downloaded_packages
