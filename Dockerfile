FROM rocker/ml-verse

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

RUN install2.r --skipinstalled BiocManager tidymodels xgboost optparse mlbench vip cpp11 && rm -rf /tmp/downloaded_packages

RUN Rscript -e "BiocManager::install('BiocParallel')"

RUN Rscript -e "devtools::install('.', dependencies=TRUE, build_vignettes=TRUE, repos = BiocManager::repositories())"

#RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); install.packages('tidymodels', 'xgboost', 'optparse', 'mlbench', 'vip')"
