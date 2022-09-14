#' Combine two qseaSets
#'
#' This function takes a list of qseaSets (or strings to read rds files from) and merges them together
#'
#' @param qseaSets A list of qseaSet objects or strings to rds files with them in
#' @param firstQset An first qseaSet object to combine the rest into (optional)
#' @param dropDuplicates Whether to drop samples with the same name. Else renames them by adding "_Dup" to the end of the name.
#' @param checkParams Whether to verify that all the parameters are identical in both objects
#' @param regionsToKeep A GRanges object (or table coercible to one) to use to subset the samples. Helps keep RAM use down.
#' @return A qseaSet object, containing all the samples from both qseaSet objects.
#' @export
combineQsetsList <- function(qseaSets, firstQset = NULL, dropDuplicates = TRUE, checkParams = TRUE, regionsToKeep = NULL) {
  if (is.character(firstQset)) {
    if(length(firstQset) == 1 & tools::file_ext(firstQset) == "rds"){#has two be afterwards, else errors if it is a qseaSet
      message(glue::glue("Character string given as firstQset, loading {firstQset}"))
      firstQset <- readr::read_rds(firstQset)
    }
  }

  #TODO catch errors better
  if(is.null(firstQset) & length(qseaSets) >= 2){
    message(glue::glue("No initial qseaSet given, using first element as initial qseaSet"))
    firstQset <- qseaSets[[1]]
    qseaSets <- utils::tail(qseaSets, n = -1)
  }

  combinedQset <- firstQset

  for(i in 1:length(qseaSets)){
    combinedQset <- combineQsets(combinedQset, qseaSets[[i]],
                                 checkParams = checkParams,
                                 regionsToKeep = regionsToKeep,
                                 dropDuplicates = dropDuplicates)
  }
  return(combinedQset)
}

#' Combine two qseaSets
#'
#' This function takes two qseaSets and combines them into one.
#'
#' @param qseaSet1 The first qseaSet object.
#' @param qseaSet2 The second qseaSet object.
#' @param checkParams Whether to verify that all the parameters are identical in both objects
#' @param dropDuplicates Whether to drop samples with the same name. Else renames them by adding "_Dup" to the end of the name.
#' @param regionsToKeep A GRanges object (or table coercible to one) to use to subset the samples. Helps keep RAM use down.
#' @return A qseaSet object, containing all the samples from both qseaSet objects.
#' @export
combineQsets <- function(qseaSet1, qseaSet2, checkParams = FALSE, regionsToKeep = NULL, dropDuplicates = FALSE) {

  if (is.character(qseaSet1)) {
    if(length(qseaSet1) == 1 & tools::file_ext(qseaSet1) == "rds"){#has to be after first if, else errors if it is a qseaSet
      message(glue::glue("Character string given, loading {qseaSet1}"))
      qseaSet1 <- readr::read_rds(qseaSet1)
    }
  }

  if (!is.qseaSet(qseaSet1)) {
    stop("Please supply a qseaSet object in the first position.")
  }

  if(!is.null(regionsToKeep)){
    qseaSet1 <- qseaSet1 %>% filterByOverlaps(regionsToKeep)
  }

  if (is.qseaSet(qseaSet2)) {
    if(length(qseaSet2) == 1 & tools::file_ext(qseaSet2) == "rds") {#has to be after first if, else errors if it is a qseaSet
      message(glue::glue("Character string given, loading {qseaSet2}"))
      qseaSet2 <- readr::read_rds(qseaSet2)
    }
  }

  if (!is.qseaSet(qseaSet2)) {
    stop("Please supply a qseaSet object in the second position.")
  }

  if(!is.null(regionsToKeep)){
    qseaSet2 <- qseaSet2 %>% filterByOverlaps(regionsToKeep)
  }

  sampleNames1 <- qsea::getSampleNames(qseaSet1)
  sampleNames2 <- qsea::getSampleNames(qseaSet2)

  commonNames <- intersect(sampleNames1,sampleNames2)

  if (length(commonNames) > 0 & !dropDuplicates) {
    message(glue::glue("Samples exist with the same name, adding _Dup to their name"))
    message(glue::glue("{commonNames} "))

    for (i in commonNames) {
      qseaSet2 <- renameQsetNames(qseaSet2, paste0("^",i,"$"), paste0(i,"_Dup"))
      sampleNames2 <- qsea::getSampleNames(qseaSet2)
      commonNames <- intersect(sampleNames1,sampleNames2)
    }
  }

  if (all(sampleNames2 %in% sampleNames1 )) {

    message("All names in common, returning first qseaSet")
    return(qseaSet1)
  }

  if (all(sampleNames1 %in% sampleNames2 )) {

    message("All names in common, returning second qseaSet")
    return(qseaSet2)
  }




  if (length(commonNames) > 0) {
    if (length(commonNamesNoControl) > 0 & dropDuplicates) {
      message(
        glue::glue(
          "Dropping {length(commonNamesNoControl)} common names (use dropDuplicates=FALSE to stop this)"
        )
      )
    }
    qseaSet2 <-
      subsetQset(qseaSet2, samplesToDrop = commonNames)
  }

  slots1 <- methods::slotNames(qseaSet1)
  slots2 <- methods::slotNames(qseaSet2)

  if (!identical(slots1,slots2)) {stop("Objects have different slots!")}
  #if (!identical(qseaSet1@parameters, qseaSet2@parameters)) {stop("Parameters are not the same, use checkParams = FALSE to combine anyway.")}

  parameters1 <- tibble::as_tibble(qseaSet1@parameters) %>% dplyr::select(order(colnames(.)))
  parameters2 <- tibble::as_tibble(qseaSet2@parameters) %>% dplyr::select(order(colnames(.)))


  if (!all(identical(GenomeInfoDb::seqnames(qseaSet1@regions),
                     GenomeInfoDb::seqnames(qseaSet2@regions)),
           identical(IRanges::ranges(qseaSet1@regions),
                     IRanges::ranges(qseaSet2@regions)))) {

    regions1 <- qseaSet1 %>% qsea::getRegions() %>% tibble::as_tibble() %>% plyranges::as_granges()
    regions2 <- qseaSet2 %>% qsea::getRegions() %>% tibble::as_tibble() %>% plyranges::as_granges() #stop complaint if genome name not identical

    qseaSet1 <- qseaSet1 %>% filterByOverlaps(regions1 %>% plyranges::filter_by_overlaps(regions2))
    qseaSet2 <- qseaSet2 %>% filterByOverlaps(regions1 %>% plyranges::filter_by_overlaps(regions2))

    ##TODO Check intersection properly to make sure they are exactly start/end together
    message(glue::glue("Regions are not identical: {length(regions1)} and {length(regions2)} regions.
                           Taking intersection of {length(regions1 %>% plyranges::filter_by_overlaps(regions2))} regions."))
  }

  if (checkParams) {

    if (!all(parameters1 == parameters2)) {stop("Parameters are not the same!")}


    if (!identical(names(qseaSet1@enrichment),names(qseaSet2@enrichment))) {stop("Enrichment entries are not the same")}
    if (qseaSet1@enrichment$pattern_name != qseaSet2@enrichment$pattern_name) {stop("Pattern names are different")}
    if (!identical(qseaSet1@enrichment$density,qseaSet2@enrichment$density)) {stop("Enrichment density is different.")}
    if (!identical(qseaSet1@enrichment$n,qseaSet2@enrichment$n)) {stop("Enrichment n is different")}

  }

  newQSet <- qseaSet1

  newQSet@regions <- qseaSet1@regions %>%
    tibble::as_tibble() %>%
    dplyr::left_join(tibble::as_tibble(qseaSet2@regions)) %>%
    dplyr::select(CpG_density, tidyselect::everything()) %>%
    plyranges::as_granges()

  GenomeInfoDb::seqinfo(newQSet@regions) <- GenomeInfoDb::seqinfo(qseaSet1@regions)

  newQSet@zygosity <- rbind(qseaSet1@zygosity,qseaSet2@zygosity)
  newQSet@sampleTable <- dplyr::bind_rows(qseaSet1@sampleTable,qseaSet2@sampleTable)

  df1 <- qseaSet1@libraries$file_name
  df2 <- qseaSet2@libraries$file_name

  newQSet@libraries$file_name <- dplyr::bind_rows(as.data.frame(df1), as.data.frame(df2))

  if ("input_file" %in% names(qseaSet1@libraries)) {

    df1 <- qseaSet1@libraries$input_file
    df2 <- qseaSet2@libraries$input_file

    newQSet@libraries$input_file <- dplyr::bind_rows(as.data.frame(df1), as.data.frame(df2))
  }

  newQSet@cnv <- qseaSet1@cnv %>%
    data.frame(check.names = FALSE) %>%
    dplyr::full_join(data.frame(qseaSet2@cnv, check.names = FALSE ), by = c("seqnames", "start", "end", "width", "strand")) %>%
    plyranges::as_granges()

  newQSet@count_matrix <- cbind(qseaSet1@count_matrix, qseaSet2@count_matrix)

  newQSet@enrichment$parameters <- rbind(qseaSet1@enrichment$parameters, qseaSet2@enrichment$parameters)
  newQSet@enrichment$factors <- cbind(qseaSet1@enrichment$factors, qseaSet2@enrichment$factors)

  return(newQSet)
}
