## code to prepare `examplePairedTumourQset` dataset goes here

examplePairedTumourQset <- read_rds("/data/cep/Methylation/pipelineOutput/M022/qseaSets/M022_All_min50_max1000_w300_q10.rds") %>%
  dropPooledControl() %>%
  filterRegions(seqnames == 7, start >= 25000000, end <= 28000000, CpG_density >= 5) %>%
  filterQset(sample_name %in% c("MR006_T","MR006_N","MR007_T","MR007_N","MR059_T","MR059_N","MR055_N",
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
  mutateQset(patient = str_remove(sample_name,"_[NT]"),
             group = sample_name)

examplePairedTumourQset@libraries$input_counts <- NULL

usethis::use_data(examplePairedTumourQset, overwrite = TRUE)
