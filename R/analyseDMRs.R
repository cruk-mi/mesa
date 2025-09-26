#' Transform DMR results to long format
#'
#' Convert wide-format DMR results into a long-format table with explicit
#' columns for contrasts, effect sizes, and significance values.
#'
#' @param DMRtable `data.frame`  
#'   DMR results in wide format, typically produced by [calculateDMRs()].
#'
#' @param FDRthres `numeric(1)`  
#'   False discovery rate threshold for filtering significant windows.  
#'   **Default:** `0.05`.
#'
#' @param makePositive `logical(1)`  
#'   If `TRUE`, reverse the direction of contrasts so that all retained windows
#'   are positively associated with `group1`.  
#'   **Default:** `FALSE`.
#'
#' @return A `tibble`:
#'
#' * **pivotDMRsLonger()**: returns a long-format table with one row per
#'   DMR–contrast combination, including columns `group1`, `group2`,
#'   `deltaBeta`, `log2FC`, and `adjPval`.
#'
#' @seealso
#'   [tidyr::pivot_longer()],  
#'   [summariseDMRsByContrast()]
#'
#' @family DMR-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Convert results to long format and filter at custom FDR threshold
#' exampleTumourNormal %>%
#'   calculateDMRs(variable = "type",
#'                 contrasts = "LUAD_vs_NormalLung",
#'                 FDRthres = 0.1) %>%
#'   pivotDMRsLonger(FDRthres = 0.1)
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
#' Count the number of up- and down-regulated DMR windows per contrast.
#' Internally, results are reshaped to long format before summarisation.
#'
#' @param DMRtable `data.frame` or `GRanges`  
#'   DMR results, typically returned by [calculateDMRs()]. If in wide format,
#'   it will be converted with [pivotDMRsLonger()].
#'
#' @param FDRthres `numeric(1)`  
#'   False discovery rate threshold for significance.  
#'   **Default:** `0.05`.
#'
#' @param log2FCthres `numeric(1)`  
#'   Absolute log2 fold-change threshold for calling up/down regulation.  
#'   **Default:** `0`.
#'
#' @param deltaBetaThres `numeric(1)`  
#'   Absolute delta-beta threshold for calling up/down regulation.  
#'   **Default:** `0`.
#'
#' @return A `tibble`:
#'
#' * **summariseDMRsByContrast()**: returns one row per contrast, with
#'   counts of up-regulated (`nUp`) and down-regulated (`nDown`) DMRs that
#'   pass the specified thresholds.
#'
#' @seealso
#'   [pivotDMRsLonger()],  
#'   [summariseDMRsByGene()]
#'
#' @family DMR-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Summarise DMRs for a single contrast, using custom FDR threshold
#' exampleTumourNormal %>%
#'   calculateDMRs(variable = "type",
#'                 contrasts = "LUAD_vs_NormalLung",
#'                 FDRthres = 0.1) %>%
#'   summariseDMRsByContrast(FDRthres = 0.1)
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
    dplyr::left_join(DMRsummary, by = dplyr::join_by(group1, group2)) %>%
    dplyr::mutate(nDown = tidyr::replace_na(nDown, 0),
           nUp = tidyr::replace_na(nUp, 0)) %>%
    dplyr::select(group1, group2, nUp, nDown) %>%
    return()
}


#' Summarise DMRs by gene
#'
#' Aggregate differentially methylated region (DMR) windows by gene annotation.
#'
#' @param DMRtable `data.frame` or `GRanges`  
#'   Annotated DMRs. Must include the columns `ENSEMBL`, `SYMBOL`, and `GENENAME`
#'   (typically added with [annotateWindows()]).
#'
#' @return A `tibble`:
#'
#' * **summariseDMRsByGene()**: returns one row per gene, with counts of
#'   associated DMRs and any aggregated metadata.
#'
#' @seealso
#'   [annotateWindows()],  
#'   [summariseDMRsByContrast()]
#'
#' @family DMR-helpers
#'
#' @examples
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Summarise DMRs with explicit annotation databases
#' exampleTumourNormal %>%
#'   calculateDMRs(variable = "tumour", contrasts = "all") %>%
#'   annotateWindows(TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
#'                   annoDb = "org.Hs.eg.db") %>%
#'   summariseDMRsByGene()
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


#' Write DMR results to an Excel workbook
#'
#' Export differentially methylated region (DMR) results to an Excel workbook,
#' creating one worksheet per contrast.
#'
#' @param dataTable `data.frame` or `GRanges`  
#'   DMR results, typically produced by [calculateDMRs()]. If a `GRanges`, it
#'   will be coerced to a data frame internally.
#'
#' @param path `character(1)`  
#'   File path of the Excel workbook to write (e.g., `"dmr_results.xlsx"`).
#'
#' @param FDRthres `numeric(1)`  
#'   False discovery rate threshold used to filter DMRs before writing. 
#'   **Default:** `0.05`.
#'
#' @return (Invisibly) returns the input object:
#'
#' * **writeDMRsToExcel()**: returns `dataTable` (invisibly), enabling use in
#'   a pipeline.
#'
#' @details
#' One worksheet is created per contrast found in `dataTable`. If no contrast
#' column is present, all rows are written to a single worksheet. If file I/O is
#' restricted on the system (e.g., some HPC build environments), consider
#' writing to `tempfile(fileext = ".xlsx")` in examples or tests.
#'
#' @seealso
#'   [writeDMRsToBed()]
#'
#' @family DMR-helpers
#'
#' @examples
#' \donttest{
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Minimal example: write results for a single contrast
#' exampleTumourNormal %>%
#'   calculateDMRs(variable = "type", contrasts = "LUAD_vs_NormalLung") %>%
#'   writeDMRsToExcel(path = file.path(tempdir(),"test.xlsx"))
#'   
#' # With multiple contrasts and annotation 
#' exampleTumourNormal %>%
#'   calculateDMRs(variable = "tumour", contrasts = "all") %>%
#'   annotateWindows(TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene", annoDb = "org.Hs.eg.db") %>%
#'   writeDMRsToExcel(path = file.path(tempdir(),"test.xlsx"))
#' }
#'
#' @export
writeDMRsToExcel <- function(dataTable, path, FDRthres = 0.05) {

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
      dplyr::filter(!!dplyr::sym(x) <= FDRthres) %>%
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
#' Export differentially methylated region (DMR) results to a set of BED files,
#' with one file generated per contrast.
#'
#' @param dataTable `data.frame`  
#'   DMR results, typically returned by [calculateDMRs()]. If a `GRanges`, it
#'   will be coerced to a data frame internally.
#'
#' @param folder `character(1)`  
#'   Directory in which to write BED files. Must exist and be writable.
#'
#' @param FDRthres `numeric(1)`  
#'   False discovery rate threshold used to filter DMRs before writing.  
#'   **Default:** `0.05`.
#'
#' @return (Invisibly) returns the input object:
#'
#' * **writeDMRsToBed()**: returns `dataTable` (invisibly), enabling use in a
#'   pipeline.
#'
#' @seealso
#'   [writeDMRsToExcel()]
#'
#' @family DMR-helpers
#'
#' @examples
#' \donttest{
#' data(exampleTumourNormal, package = "mesa")
#'
#' # Export DMRs for a single contrast to BED files in a temp directory
#' exampleTumourNormal %>%
#'   calculateDMRs(variable = "type",
#'                 contrasts = "LUAD_vs_NormalLung") %>%
#'   writeDMRsToBed(folder = tempdir())
#'
#' # Export DMRs with explicit annotation (requires TxDb/annotation packages)
#' exampleTumourNormal %>%
#'   calculateDMRs(variable = "tumour", contrasts = "all") %>%
#'   annotateWindows(TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
#'                   annoDb = "org.Hs.eg.db") %>%
#'   writeDMRsToBed(folder = tempdir(), FDRthres = 0.1)
#' }
#'
#' @export
writeDMRsToBed <- function(dataTable, folder, FDRthres = 0.05) {

  dir.create(folder, showWarnings = FALSE, recursive = TRUE)

  contrasts <- dataTable %>% colnames() %>% stringr::str_subset("adjPval$")

  if(length(contrasts) == 0) {
    stop("No columns ending in adjPval found! Is this a set of DMRs generated by calculateDMRs?")
  }
  
  purrr::walk(contrasts, function(x){

    contrastName <- stringr::str_remove(x,"_adjPval$")

    dataTable %>%
      dplyr::filter(!!rlang::sym(x) <= !!FDRthres) %>%
      plyranges::as_granges() %>%
      rtracklayer::export.bed(file.path(folder, paste0(contrastName,".bed")))

    }
  )

  invisible(dataTable)
}

#' Take the top most DMRs per contrast, based on those with the largest value of the selected metric
#'
#' @param DMRs  A data frame containing the output of calculateDMRs, potentially with multiple contrasts
#' @param n How many DMRs to take of each contrast (if that many exists)
#' @param FDRthres Threshold on the adjusted p values
#' @param metric Which metric to use to select the top DMRs. Options are deltaBeta, log2FC, adjPval, CpG_density, position (using seqnames and start columns) or any other column in the data frame. 
#' If `adjPval` or `position` are used, then the window with the smallest value will be chosen, otherwise the largest value will be used.
#' @param makePositive Whether to reverse the contrast when the window is hypomethylated in the contrast.
#' @return A data frame with the DMRs with the largest value of the selected metrics
#' @examples
#' # calculate some DMRs
#' DMRs <- exampleTumourNormal %>% calculateDMRs(variable = "type", contrasts = "all", keepContrastMeans = FALSE)
#' # Find the DMRs with the largest deltaBeta between each comparison:
#' DMRs %>% sliceDMRs(n = 1)
#' # Or the windows with the largest log2FC: 
#' DMRs %>% sliceDMRs(n = 1, metric = log2FC)
#' # Or the windows with the largest CpG_density: 
#' DMRs %>% sliceDMRs(n = 1, metric = CpG_density)
#' # If adjPval is used, then the smallest value is chosen instead:
#' DMRs %>% sliceDMRs(n = 1, metric = CpG_density)
#' # If position is used, then the windows are sorted by genomic position:
#' DMRs %>% sliceDMRs(n = 1, metric = position)
#' @export
sliceDMRs <- function(DMRs, n = 1, metric = deltaBeta, makePositive = TRUE, FDRthres = 0.05) {
  
  positiveDMRs <- DMRs %>% 
    pivotDMRsLonger(makePositive = makePositive, FDRthres = FDRthres) %>% 
    dplyr::relocate(tidyselect::matches("means$"), .after = tidyselect::last_col()) %>%
    dplyr::group_by(group1, group2)

  if(deparse(substitute(metric)) %in% c("position")){
    out <- positiveDMRs %>% 
      dplyr::arrange(seqnames, start) %>% 
      dplyr::slice(1:(!!n)) %>%
      dplyr::ungroup()
  } else if (deparse(substitute(metric)) %in% c("adjPval")){
    message(glue::glue("Choosing the windows with the smallest value of {rlang::ensym(metric)}"))
    out <- positiveDMRs %>% 
      dplyr::arrange({{metric}}) %>% 
      dplyr::slice(1:(!!n)) %>%
      dplyr::ungroup()
  } else {
    out <- positiveDMRs %>% 
      dplyr::arrange(dplyr::desc({{metric}})) %>% 
      dplyr::slice(1:(!!n)) %>%
      dplyr::ungroup()
  }

  return(out)
}
