FROM rocker/ml-verse

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

#install BiocManager and other packages required for the classifier pipeline
RUN install2.r --skipinstalled BiocManager tidymodels xgboost optparse mlbench vip cpp11 butcher qs ggpubr ggrepel && rm -rf /tmp/downloaded_packages

#install BiocParallel first because it was causing failures
RUN Rscript -e "BiocManager::install('BiocParallel')"

#install other BSgenomes
RUN Rscript -e "BiocManager::install('BSgenome.Mmusculus.UCSC.mm39')"
RUN Rscript -e "BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')"
RUN Rscript -e "BiocManager::install('BSgenome.Rnorvegicus.UCSC.rn7')"

#install sesame
RUN Rscript -e "BiocManager::install('sesame')"
RUN Rscript -e "sesame::sesameDataCache()"
