## code to prepare `exampleTumourNormal` dataset goes here

## Example data for showing the package's functions, not a dataset for
## analysis. Tumour and normal colon and lung tissue bought from OriGene and
## processed at CRUK-MI CBC (project M022) into a qseaSet on GRCh38 with
## makeQset(): insert size 50-1000 bp, 300 bp windows, MAPQ 10. The subset is
## unpublished and the full qseaSet is internal, so this script cannot be
## re-run outside CRUK-MI.
## Originally built from the CRUK internal file
## /data/cep/Methylation/pipelineOutput/M022/qseaSets/
## M022_All_min50_max1000_w300_q10.rds.
## Subset: chr7 25-28 Mb, windows with CpG density >= 5, and 5 tumour/normal
## pairs renamed Colon1-2 and Lung1-3. The input counts are dropped.

qsetFile <- "M022_All_min50_max1000_w300_q10.rds"

exampleTumourNormal <- read_rds(qsetFile) %>%
  filterWindows(seqnames == 7, start >= 25000000, end <= 28000000, CpG_density >= 5) %>%
  filter(sample_name %in% c("MR006_T","MR006_N","MR007_T","MR007_N","MR059_T","MR059_N","MR055_N",
                                "MR060_N","MR060_T","MR062_T","MR062_N")) %>%
  selectQset(sample_name, group, patient, type, tumour, tissue, age, gender, stage) %>%
  renameQsetNames("MR006_T","Colon1_T") %>%
  renameQsetNames("MR007_T","Colon2_T") %>%
  renameQsetNames("MR006_N","Colon1_N") %>%
  renameQsetNames("MR007_N","Colon2_N") %>%
  renameQsetNames("MR059_T","Lung1_T") %>%
  renameQsetNames("MR059_N","Lung1_N") %>%
  renameQsetNames("MR060_T","Lung2_T") %>%
  renameQsetNames("MR060_N","Lung2_N") %>%
  renameQsetNames("MR062_T","Lung3_T") %>%
  renameQsetNames("MR062_N","Lung3_N") %>%
  mutate(patient = str_remove(sample_name,"_[NT]"),
             group = sample_name)

exampleTumourNormal@libraries$input_counts <- NULL

usethis::use_data(exampleTumourNormal, overwrite = TRUE, compress = "xz")
