#' Transform DMR results to long format
#'
#' Convert a wide-format DMR results table into a long-format table
#' with explicit columns for contrasts, effect sizes, and significance values.
#'
#' @param DMRtable A `data.frame` of DMR results in wide format, typically the
#'   output of [calculateDMRs()].
#' @param FDRthres `numeric(1)` False discovery rate threshold to filter significant
#'   windows.
#' @param makePositive `logical(1)` If `TRUE`, reverse the direction of contrasts
#'   such that all retained windows are positively associated with `group1`.
#'
#' @return A tibble in long format with columns including `group1`, `group2`,
#'   `deltaBeta`, `log2FC`, and `adjPval`.
#'
#' @seealso [tidyr::pivot_longer()], [summariseDMRsByContrast()]
#' @family DMR-helpers
#' 
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' dmr <- calculateDMRs(exampleTumourNormal, variable = "type",
#'                      contrasts = "LUAD_vs_NormalLung", fdrThres = 0.1)
#' head(pivotDMRsLonger(dmr, FDRthres = 0.1))
#'
#' @export
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


#' Summarise DMRs by contrast
#'
#' Count the number of up- and down-regulated windows per contrast.
#' Internally transforms the DMR results to long format.
#'
#' @param DMRtable A `data.frame` of DMR results, typically from [calculateDMRs()].
#'   If in wide format, it will be converted with [pivotDMRsLonger()].
#' @param FDRthres `numeric(1)` False discovery rate threshold.
#' @param log2FCthres `numeric(1)` Absolute log2 fold-change threshold.
#' @param deltaBetaThres `numeric(1)` Absolute delta-beta threshold.
#'
#' @return A tibble with one row per contrast and columns `nUp` and `nDown`.
#'
#' @seealso [pivotDMRsLonger()], [summariseDMRsByGene()]
#' @family DMR-helpers
#' 
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' dmr <- calculateDMRs(exampleTumourNormal, variable = "type",
#'                      contrasts = "LUAD_vs_NormalLung", fdrThres = 0.1)
#' summariseDMRsByContrast(dmr, FDRthres = 0.1)
#'
#' @export
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


#' Summarise DMRs by gene
#'
#' Aggregate DMR windows by gene annotation.
#'
#' @param DMRtable A `data.frame` or `GRanges` of annotated DMRs.
#'   Must include columns `ENSEMBL`, `SYMBOL`, and `GENENAME`
#'   (typically added with [annotateWindows()]).
#'
#' @return A tibble with one row per gene and a count of associated DMRs.
#'
#' @seealso [annotateWindows()], [summariseDMRsByContrast()]
#' @family DMR-helpers
#' 
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#' dmr <- calculateDMRs(exampleTumourNormal, variable = "type",
#'                      contrasts = "LUAD_vs_NormalLung", fdrThres = 0.1)
#' dmr_annot <- annotateWindows(dmr)  # requires TxDb/annotation
#' summariseDMRsByGene(dmr_annot)
#'
#' @export
summariseDMRsByGene <- function(DMRtable){

  if(!("ENSEMBL" %in% colnames(DMRtable))){
    stop("Please call annotateWindows first on the DMRs prior to summariseDMRsByGene.")
  }

  summary <- DMRtable %>%
    tibble::as_tibble() %>%
    dplyr::group_by(ENSEMBL, SYMBOL, GENENAME) %>%
    dplyr::tally() %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(n))
  
  return(summary)

}


#' Write DMR results to Excel
#'
#' Export DMR results to an Excel workbook with one sheet per contrast.
#'
#' @param dataTable A `data.frame` of DMR results, typically from [calculateDMRs()].
#' @param path `character(1)` File path to write the Excel workbook.
#' @param fdrThres `numeric(1)` False discovery rate threshold for filtering.
#'
#' @return Invisibly returns the input `dataTable`, enabling use in a pipe.
#'
#' @seealso [writeDMRsToBed()]
#' @family DMR-helpers
#' 
#' @examples
#' \dontrun{
#' data(exampleTumourNormal, package = "mesa")
#' dmr <- calculateDMRs(exampleTumourNormal, variable = "type",
#'                      contrasts = "LUAD_vs_NormalLung")
#' writeDMRsToExcel(dmr, path = "dmr_results.xlsx")
#' }
#'
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


#' Write DMR results to BED files
#'
#' Export DMR results to a set of BED files, one file per contrast.
#'
#' @param dataTable A `data.frame` of DMR results, typically from [calculateDMRs()].
#' @param folder `character(1)` Directory in which to write BED files.
#' @param fdrThres `numeric(1)` False discovery rate threshold for filtering.
#'
#' @return Invisibly returns the input `dataTable`, enabling use in a pipe.
#'
#' @seealso [writeDMRsToExcel()]
#' @family DMR-helpers
#' 
#' @examples
#' \dontrun{
#' data(exampleTumourNormal, package = "mesa")
#' dmr <- calculateDMRs(exampleTumourNormal, variable = "type",
#'                      contrasts = "LUAD_vs_NormalLung")
#' writeDMRsToBed(dmr, folder = "bed_results")
#' }
#'
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
