## code to prepare `hg38CpGIslands` dataset goes here

#Generated via:
# ml apps/bedtools2/2.27.1-7
# wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz \
# | gunzip -c \
# | awk 'BEGIN{ OFS="\t"; }{ print $2, $3, $4, $5$6, $7, $8, $9, $10, $11, $12 ; }' \
# | sed 's/chr//g' \
# | bedtools sort -i -  > cpgIslandExt.hg38.noChr.noHead.bed
# echo $'chrom\tstart\tend\tname\tlength\tcpgNum\tgcNum\tperCpg\tperGc\tobsExp' | cat - cpgIslandExt.hg38.noChr.noHead.bed > cpgIslandExt.hg38.noChr.bed

hg38CpGIslands <- regioneR::toGRanges(read.table("/data/cep/Methylation/refData/cpgIslandExt.hg38.noChr.bed",
                                                 sep = "\t", head = TRUE)) %>%
  dplyr::select(-name)

usethis::use_data(hg38CpGIslands, overwrite = TRUE, internal = TRUE)
