FROM crukmi/mesa_dependencies

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

#RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); BiocManager::install(ask=FALSE)"

RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); devtools::install('.', dependencies=TRUE, build_vignettes=TRUE, repos = BiocManager::repositories())"

RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); install.packages('tidymodels', 'xgboost', 'optparse', 'mlbench', 'vip')"
