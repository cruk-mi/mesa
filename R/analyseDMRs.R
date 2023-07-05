#' This function transforms a DMR data frame from the wide format to a long format
#' @param DMRtable A data frame with multiple comparisons inside
#' @param FDRthres FDR threshold to apply to each comparison
#' @param makePositive Whether to reverse the contrast when the window is hypomethylated in the contrast
#' @export
#'
pivotDMRsLonger <- function(DMRtable, FDRthres = 0.05, makePositive = FALSE){
  pivotedDMRs <- DMRtable %>%
    dplyr::rename_with(~ stringr::str_replace(.x, "_adjPval", ":adjPval")) %>%
    dplyr::rename_with(~ stringr::str_replace(.x, "_log2FC", ":log2FC")) %>%
    dplyr::rename_with(~ stringr::str_replace(.x, "_deltaBeta", ":deltaBeta")) %>%
    tidyr::pivot_longer(tidyselect::matches("_vs_"),
                        names_to = c(".comparison", ".ext"),
                        values_to = ".value",
                        names_sep = ":") %>%
    tidyr::pivot_wider(names_from = .ext,
                       values_from = .value) %>%
    dplyr::filter(adjPval <= FDRthres) %>%
    tidyr::separate(.comparison, into = c("group1","group2"), sep = "_vs_")

  if(makePositive){
    pivotedDMRs <- pivotedDMRs %>%
      dplyr::mutate(group1new = ifelse(log2FC < 0, group2, group1),
                    group2new = ifelse(log2FC < 0, group1, group2),
                    deltaBeta = sign(log2FC)*deltaBeta,
                    log2FC = sign(log2FC)*log2FC,
      ) %>%
      dplyr::select(-group1, -group2) %>%
      dplyr::rename(group1 = group1new, group2 = group2new)
  }

  pivotedDMRs %>%
    return()
}

#' This function summarises the number of windows called as up/down in each comparison in a .
#' Internally transforms the data to a long format. Contrasts with no significant windows will be dropped!
#' @param DMRtable A data frame with multiple comparisons inside
#' @param FDRthres FDR threshold to apply to each comparison
#' @param log2FCthres A log2FC threshold to apply to each comparison (absolute)
#' @param deltaBetaThres A deltaBeta (change in average beta values) threshold to apply to each comparison (absolute)
#' @export
#'
summariseDMRsByContrast <- function(DMRtable, FDRthres = 0.05, log2FCthres = 0, deltaBetaThres = 0){

  if (!("adjPval" %in% colnames(DMRtable))) {

    contrastsDF <- DMRtable %>%
      colnames() %>%
      stringr::str_subset(".*_vs_.*_adjPval$") %>%
      tibble::enframe(name = NULL) %>%
      tidyr::separate(value, into = c("group1","group2",NA), sep = "_vs_|_adjPval") %>%
      dplyr::arrange(group1, group2)

    DMRtable <- DMRtable %>%
      dplyr::select(-tidyselect::matches("_nrpm$|_beta$|_means$")) %>%
      pivotDMRsLonger(FDRthres = FDRthres)

  } else {
    contrastsDF <- DMRtable %>%
      dplyr::distinct(group1, group2) %>%
      dplyr::arrange(group1, group2)
  }

  DMRsummary <- DMRtable %>%
    tibble::as_tibble() %>%
    dplyr::filter(adjPval <= FDRthres,
                  abs(log2FC) >= log2FCthres,
                  abs(deltaBeta) >= deltaBetaThres,
                  ) %>%
    dplyr::mutate(.up = ifelse(log2FC > 0, "nUp","nDown")) %>%
    dplyr::group_by(group1, group2, .up) %>%
    dplyr::tally() %>%
    tidyr::pivot_wider(names_from = .up, values_from = n, values_fill = 0) %>%
    dplyr::ungroup()

  contrastsDF %>%
    dplyr::left_join(DMRsummary) %>%
    dplyr::mutate(nDown = tidyr::replace_na(nDown, 0),
           nUp = tidyr::replace_na(nUp, 0)) %>%
    dplyr::select(group1, group2, nUp, nDown) %>%
    return()
}

#' This function summarises a data frame by gene (using annotateWindows)
#' @param DMRtable A data frame or GRanges object, may already be annotated
#' @export
#'
summariseDMRsByGene <- function(DMRtable){

  if(!("ENSEMBL" %in% colnames(DMRtable))){
    DMRtable <- DMRtable %>% annotateWindows()
  }

  DMRtable %>%
    tibble::as_tibble() %>%
    dplyr::group_by(ENSEMBL, SYMBOL, GENENAME, geneChr, geneStart, geneEnd, geneLength) %>%
    dplyr::summarise(n = length(width)) %>%
    dplyr::arrange(dplyr::desc(n)) %>%
    return()

}


#' Write the qsea result to an Excel file
#'
#' This writes an Excel file with the results, with one sheet per contrast.
#'
#' @param dataTable An (optionally annotated) data table with contrasts
#' @param path The path to write the excel file to.
#' @param fdrThres FDR rate threshold to apply to each contrast
#' @return Returns the original table invisibly, so can be used in a pipe.
#' @export
writeDMRsToExcel <- function(dataTable, path, fdrThres = 0.05) {

  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop(
      "Package \"openxlsx\" must be installed to use this function.",
      call. = FALSE
    )
  }

  wb_DMR <- openxlsx::createWorkbook()

  contrasts <- dataTable %>% colnames() %>% stringr::str_subset("adjPval$")

  purrr::walk(contrasts, function(x){

    ##TODO Make this cleaner, particularly for multiple tabs to stop them being the same name
    ## Probably by using walk2 and the name.
    sheetName <- stringr::str_trunc(x,31,"right", ellipsis = "")

    openxlsx::addWorksheet(wb_DMR, sheetName)
    dataTable %>%
      dplyr::filter(!!dplyr::sym(x) <= fdrThres) %>%
      openxlsx::writeData(wb_DMR, sheetName, .)
    openxlsx::addFilter(wb_DMR, sheetName, row = 1, cols = 1:ncol(dataTable))
    openxlsx::freezePane(wb_DMR, sheetName,  firstRow = TRUE)

  }
  )

  openxlsx::saveWorkbook(wb_DMR, file = path, overwrite = TRUE)
  invisible(dataTable)
}

#' Write the qsea result to a set of bed file
#'
#' This writes a set of bed files with the results, with one bed file per contrast.
#'
#' @param dataTable  An (annotated) data table with contrasts
#' @param folder The folder to place the bed files.
#' @param fdrThres FDR rate threshold to apply to each contrast
#' @return Returns the original table invisibly, so can be used in a pipe.
#' @export
writeDMRsToBed <- function(dataTable, folder, fdrThres = 0.05) {

  dir.create(folder, showWarnings = TRUE, recursive = TRUE)

  contrasts <- dataTable %>% colnames() %>% stringr::str_subset("adjPval$")

  purrr::walk(contrasts, function(x){

    contrastName <- stringr::str_remove(x,"_adjPval$")

    dataTable %>%
      dplyr::filter(!!rlang::sym(x) <= fdrThres) %>%
      plyranges::as_granges() %>%
      rtracklayer::export.bed(file = file.path(folder, paste0(contrastName,".bed")))

  }
  )

  invisible(dataTable)
}
