#' Construct an initial qseaSet from BAMs and metadata
#'
#' Build a `qseaSet` from a `sampleTable` that lists BAM files and sample
#' metadata, tile the genome into fixed windows, optionally remove blacklisted
#' regions, and load coverage and (optionally) CNV information according to the
#' selected methods. This is typically the first step in a mesa/qsea workflow.
#'
#' @details
#' **Inputs & requirements**
#' - `sampleTable` must include at least:
#'   - `sample_name` (unique per sample),
#'   - `file_name` (path to MeDIP/capture BAM),
#'   - `group` (string label).  
#'   If the selected `CNVmethod` uses inputs (`"HMMdefault"`, `"qseaInput"`), 
#'   also include `input_file` (path to matched Input BAM).
#' - BAMs should be coordinate-sorted. For paired data, running `samtools fixmate`
#'   is recommended so MAPQ tags are available (see `minMapQual` note below).
#' - `BSgenome` must be an installed genome package string (e.g.
#'   `"BSgenome.Hsapiens.NCBI.GRCh38"`), which provides chromosome lengths and
#'   sequences for tiling windows.
#'
#' **Windows & blacklist**
#' - Chromosomes in `chrSelect` are tiled into non-overlapping windows of size
#'   `windowSize` (bp). Windows overlapping `badRegions` (if provided) are removed.
#'
#' **Coverage loading**
#' - `coverageMethod = "PairedAndR1s"` loads proper pairs and, optionally,
#'   high-MAPQ R1 singletons; governed by `minMapQual`, `minInsertSize`,
#'   `maxInsertSize`, `properPairsOnly`, `minReferenceLength`.
#' - `coverageMethod = "qseaPaired"` delegates to [qsea::addCoverage()] with 
#'   paired settings.
#' - `fragmentType` can set defaults for `fragmentLength`/`fragmentSD` (`"Sheared"`
#'   or `"cfDNA"`); otherwise supply both explicitly.
#'
#' **CNV options**
#' - `"qseaInput"`: run [qsea::addCNV()] using `input_file` BAMs.
#' - `"HMMdefault"`: use mesa's [addHMMcopyCNV()] with supplied defaults for GRCh38
#'   GC/mappability tracks (objects like `gc_hg38_1000kb`, `map_hg38_1000kb`,
#'   chosen by `CNVwindowSize`).
#' - `"MeCap"`: estimate CNV from MeCap BAMs (`file_name`).
#' - `"None"`: do not compute CNV.
#'
#' **Parallelisation**
#' - If `parallel = TRUE`, computation uses the currently registered
#'   BiocParallel backend (see [BiocParallel::register()]); otherwise runs serially.
#'   Use `setMesaParallel()` to toggle package-wide.
#'
#' **Additional processing**
#' - Adds CpG density via `qsea::addPatternDensity(pattern = "CG")`, and 
#'   enrichment parameters, allowing for estimation of beta values.
#'
#' @param sampleTable `data.frame`.  
#'   Must contain `sample_name`, `file_name`, and `group`. If `CNVmethod`
#'   requires inputs, also `input_file`. Additional columns are preserved as
#'   sample metadata.  
#'   **Default:** none (must be supplied).
#'
#' @param BSgenome `character(1)` or `NULL`.  
#'   BSgenome package string (e.g., `"BSgenome.Hsapiens.NCBI.GRCh38"`).
#'   See [BSgenome::available.genomes()].  
#'   **Default:** `NULL` (must be supplied).
#'
#' @param chrSelect `integer()` or `character()`.  
#'   Chromosomes to include.
#'   **Default:** `1:22`.
#'
#' @param windowSize `integer(1)`.  
#'   Window size in bp for tiling the genome.  
#'   **Default:** `300`.
#'
#' @param CNVwindowSize `integer(1)`.  
#'   Copy number variation (CNV) window size (bp).
#'   **Default:** `1e6`.
#'
#' @param fragmentType `character(1)` or `NULL`.  
#'   If `"Sheared"` or `"cfDNA"`, sets defaults for `fragmentLength`/`fragmentSD`;
#'   otherwise both must be supplied explicitly.  
#'   **Default:** `NULL`.
#'
#' @param fragmentLength `numeric(1)` or `NULL`.  
#'   Average fragment length (bp). This is used to extend single-end reads if 
#'   used, but also strongly affects the CpG density calculation, which will 
#'   strongly impact the beta value estimation.
#'   **Default:** `NULL` (inferred from `fragmentType` or must be provided).
#'
#' @param fragmentSD `numeric(1)` or `NULL`.  
#'   Standard deviation of fragment length (bp). Used to calculate the window 
#'   based CpG density calculation, which will strongly impact the beta value 
#'   estimation.
#'   **Default:** `NULL` (inferred from `fragmentType` or must be provided).
#'
#' @param CNVmethod `character(1)`.  
#'   One of `"HMMdefault"`, `"qseaInput"`, `"MeCap"`, `"None"`.  
#'   **Default:** `"HMMdefault"`.
#'
#' @param coverageMethod `character(1)`.  
#'   One of `"PairedAndR1s"` (mesa) or `"qseaPaired"` (qsea).  
#'   **Default:** `"PairedAndR1s"`.
#'
#' @param minMapQual `integer(1)`.  
#'   Minimum MAPQ to retain a read. For `"PairedAndR1s"`, a pair is kept if
#'   either end passes when properly paired and MAPQ tags are set.  
#'   **Default:** `10`.
#'
#' @param minInsertSize `integer(1)`.  
#'   Minimum absolute insert size for proper pairs (bp).
#'   Only applies to coverageMethod = "PairedAndR1s", and applies to Input 
#'   samples as well as MeCap.
#'   **Default:** `70`.
#'
#' @param maxInsertSize `integer(1)`.  
#'   Maximum absolute insert size for proper pairs (bp). 
#'   Only applies to coverageMethod = "PairedAndR1s", and applies to Input 
#'   samples as well as MeCap. 
#'   **Default:** `1000`.
#'
#' @param minReferenceLength `integer(1)`.  
#'   A minimum mapped distance on the genome to keep an read. This is to avoid
#'   very short reads that may be due to mapping artefacts. Note that with 
#'   defaultparameters bwa will allow reads where only 19bps have mapped to the
#'   genome.
#'   Only applies to coverageMethod = "PairedAndR1s", and applies to Input 
#'   samples as well as MeCap. 
#'   **Default:** `30`.
#'
#' @param badRegions `GRanges` or `NULL`.  
#'   Genomic regions to exclude (blacklist). These regions will be removed from 
#'   the output. Note a set of GRCh38 windows identified by ENCODE are provided 
#'   in the [ENCODEbadRegions] data object.
#'   **Default:** `NULL`.
#'
#' @param properPairsOnly `logical(1)`.  
#'   If `TRUE`, use only proper pairs (stricter size selection).  
#'   **Default:** `FALSE`.
#'
#' @param hmmCopyGC `data.frame` or `NULL`.  
#'   GC content per CNV bin (size = `CNVwindowSize`) for HMMcopy.
#'   Required if using the `CNVmethod = "HMMdefault"` option *unless* you are 
#'   using hg38/GRCh38 with provided default bin sizes (50kb, 500kb, 1Mb).
#'   Must be a dataframe with columns `chr`, `start`, `end`, and `gc`.
#'   **Default:** `NULL`.
#'
#' @param hmmCopyMap `data.frame` or `NULL`.  
#'   Mappability per CNV bin (size = `CNVwindowSize`) for HMMcopy.  
#'   Required if using the `CNVmethod = "HMMdefault"` option *unless* you are 
#'   using hg38/GRCh38 with provided default bin sizes (50kb, 500kb, 1Mb).
#'   Must be a dataframe with columns `chr`, `start`, `end`, and `map`.
#'   **Default:** `NULL`.
#'
#' @param maxPatternDensity `numeric(1)`.
#'  Maximum pattern density (e.g., CpG) to consider a window for background 
#'  offset calculation. Regions above this are excluded from the calculation.
#'  This may need to be set higher if you are not considering many windows or 
#'  have a high average genomic CG content.
#'  **Default:** `0.05`.
#'  
#' @param enrichmentMethod `character(1)`.
#'  Method for calculating enrichment for target pattern (e.g. CpG). 
#'  See [addNormalisation()] for details.
#'  **Default:** `"blind1-15"`.
#'
#' @param parallel `logical(1)`.  
#'   Use the registered BiocParallel backend (`TRUE`) or run serially (`FALSE`).  
#'   See [BiocParallel::register()].  
#'   **Default:** `getMesaParallel()`.
#'
#' @return A `qseaSet` containing:  
#' * tiled (and optionally blacklisted) genomic windows;  
#' * the input `sampleTable` as sample metadata;  
#' * loaded coverage according to `coverageMethod`;  
#' * CNV results depending on `CNVmethod`;  
#' * CpG density and enrichment parameters;  
#' * recorded processing parameters in `@parameters`.
#'
#' @seealso
#' [qsea::createQseaSet()], [qsea::addCoverage()], [qsea::addCNV()],  
#' [addHMMcopyCNV()], [setMesaParallel()], [BSgenome::available.genomes()]
#'
#' @examples
#' # Minimal runnable sketch if MEDIPSData is installed
#' if (requireNamespace("MEDIPSData", quietly = TRUE)) {
#'   sampleTable <- data.frame(
#'     sample_name = c("Normal1", "Tumour1"),
#'     group       = c("Normal",  "Tumour"),
#'     file_name   = c(
#'       system.file("extdata", "NSCLC_MeDIP_1N_fst_chr_20_21_22.bam",
#'                   package = "MEDIPSData", mustWork = TRUE),
#'       system.file("extdata", "NSCLC_MeDIP_1T_fst_chr_20_21_22.bam",
#'                   package = "MEDIPSData", mustWork = TRUE)
#'     )
#'   )
#'
#'   sampleTable %>%
#'     makeQset(
#'       BSgenome          = "BSgenome.Hsapiens.UCSC.hg19",
#'       chrSelect         = paste0("chr", 20:22),
#'       windowSize        = 300,
#'       fragmentLength    = 200,
#'       fragmentSD        = 50,
#'       CNVmethod         = "None",
#'       coverageMethod    = "PairedAndR1s",
#'       properPairsOnly   = FALSE,
#'       minMapQual        = 10
#'     )
#' }
#' @export
makeQset <- function(sampleTable,
                     BSgenome,
                     chrSelect,
                     windowSize = 300,
                     CNVwindowSize = 1000000,
                     fragmentType = NULL,
                     fragmentLength = NULL,
                     fragmentSD = NULL,
                     CNVmethod = "HMMdefault",
                     coverageMethod = "PairedAndR1s",
                     minMapQual = 10,
                     minInsertSize = 70,
                     maxInsertSize = 1000,
                     minReferenceLength = 30,
                     badRegions = NULL,
                     properPairsOnly = FALSE,
                     hmmCopyGC = NULL,
                     hmmCopyMap = NULL,
                     maxPatternDensity = 0.05,
                     enrichmentMethod = "blind1-15",
                     parallel = getMesaParallel()) {
  
  if(parallel) {
    if(BiocParallel::bpworkers() == 1){
      message("No configured parallelisation, use e.g. register(MulticoreParam(workers = 4)) to process multiple files at once.")
      parallel <- FALSE
    } else {
      message(glue::glue("Detected parallel setup with {BiocParallel::bpworkers()} workers."))
    }

  }


  if (!is.null(fragmentType)) {
    
    if(!is.null(fragmentLength)) {
      stop("Please specify one or the other of fragmentType and fragmentLength.")
    }
    
    if (fragmentType %in% c("Sheared","sheared") ) {
      fragmentLength <- 213
      fragmentSD <- 60
    } else if (fragmentType == "cfDNA") {
      fragmentLength <- 167
      fragmentSD <- 38
    } else {
      stop("fragmentType should be either Sheared or cfDNA.")}
  }

  if (is.null(fragmentLength) | is.null(fragmentSD)) {
    stop("fragmentLength and fragmentSD must be specified, or fragmentType can be specified for some defaults.")
  }

  # convert sampleTable to data.frame as qseaSet doesn't like tibbles.
  sampleTable <- as.data.frame(sampleTable)
  
  if (!("sample_name" %in% colnames(sampleTable)))  {
    stop("Required column sample_name not included in the sampleTable.")
  }

  if (!("file_name" %in% colnames(sampleTable))) {
    stop("Required column file_name not included in the sampleTable.")
  }

  if (!("group" %in% colnames(sampleTable)))  {
    stop("Required column group included in the sampleTable.")
  }

  asValidNames <- base::make.names(sampleTable$sample_name)
  
  if(any(asValidNames != sampleTable$sample_name)){
    stop(glue::glue("sample_name column must be valid names for columns in R without quoting.
  See the help for base::make.names, but generally use only letters, numbers, 
  underscores and dots, and names can't start with a number. 
  Issues were found with: 
    {paste(sampleTable$sample_name[sampleTable$sample_name != asValidNames], collapse = '\n    ')}
   "))
  }
  
  if(!all(file.exists(sampleTable$file_name))){
    stop(glue::glue(" MeCap file not found for: {dplyr::filter(sampleTable,!file.exists(file_name)) %>% dplyr::pull(sample_name)}. "))
  }

  if ((!("input_file" %in% colnames(sampleTable))) & CNVmethod == "HMMdefault")  {
    stop("Input_file column needed in the sampleTable for calculating Copy Number Variation (CNV), set CNVmethod = none to not call CNV.")
  }

  if (("input_file" %in% colnames(sampleTable)) & CNVmethod != "none")  {
    if(!all(file.exists(sampleTable$input_file))){
      stop(glue::glue("Input file not found for: {dplyr::filter(sampleTable,!file.exists(input_file)) %>% dplyr::pull(input_file)}. "))
    }
  }

  if (is.null(badRegions)) {
    badRegions <- GenomicRanges::GRanges()
  }

  # Get data from the BSgenome
  refGenome <- BSgenome::getBSgenome(BSgenome)
  
  # check that the chromosome names match the BS genome
  bsNames <- GenomeInfoDb::seqnames(refGenome)
  unknownChr <- setdiff(chrSelect, bsNames)
  if (length(unknownChr) > 0) {
    
    if(stringr::str_detect(BSgenome,"UCSC") && all(stringr::str_detect(chrSelect,"chr", negate = TRUE))) {
      stop(glue::glue("UCSC genome requires 'chr' prefixes which were not found. \\
            Add these or consider swapping to a BSgenome without 'chr', \\
            e.g. BSgenome.Hsapiens.NCBI.GRCh38"))
    } 
    
    if(stringr::str_detect(BSgenome,"NCBI") && all(stringr::str_detect(chrSelect,"chr"))) {
      stop(glue::glue("NCBI genome does not use 'chr' prefixes. \\
            Remove these or consider swapping to e.g. BSgenome.Hsapiens.UCSC.hg38"))
    }
      
    stop(glue::glue(
      "Chromosomes provided not found in the given BSgenome!\n",
      "Showing first {length(unknownChr)} not found:\n",
      "    {glue::glue_collapse(head(unknownChr), sep = '\n    ')}"
    ))
  }
  
  # chromosome lengths, and then a Seqinfo object with that
  chrLength <- GenomeInfoDb::seqlengths(refGenome)[chrSelect]
  seqinfo <- GenomeInfoDb::Seqinfo(as.character(chrSelect),chrLength, NA, BSgenome)

  # number of windows of size windowSize on each chromosome
  numWindows <- floor(chrLength/windowSize)

  # starting point of each window on each chromosome, then make a GRanges object with all those windows
  windowStart <- unlist(lapply(FUN = seq, X = chrLength - windowSize + 1,
                            from = 1, by = windowSize), FALSE, FALSE)

  windowsGRanges <- GenomicRanges::GRanges(seqnames = rep(factor(chrSelect), numWindows),
                                           ranges = IRanges::IRanges(start = windowStart,
                                                                     width = windowSize),
                                           seqinfo = seqinfo)

  # remove both sets of blacklisted windows from the full set of windows
  windowsWithoutBlacklist <- windowsGRanges %>%
    plyranges::filter_by_non_overlaps(badRegions)

  message(
    "Considering ", length(windowsWithoutBlacklist),
    " regions with total size ",
    sum(BiocGenerics::width(windowsWithoutBlacklist))
  )

  #make the initial Qsea object, with the reduced set of windows, using the sampleTable
  qseaSet <- qsea::createQseaSet(sampleTable = sampleTable,
                                 BSgenome = BSgenome,
                                 chr.select = chrSelect,
                                 Regions = windowsWithoutBlacklist,
                                 window_size = windowSize)

  if(any(GenomeInfoDb::genome(qsea::getRegions(qseaSet)) != BSgenome)) {
    regions <- qseaSet %>% qsea::getRegions() 
    GenomeInfoDb::genome(regions) <- BSgenome
    qseaSet@regions <- regions
  }

  # store the extra blacklist file location
  #qseaSet@parameters$badRegions2 <- blacklistBed

  #TODO add a single-end coverage method.

  if (coverageMethod == "qseaPaired") {

    #load the coverage from each bam file, using qsea default method
    qseaSet <- qsea::addCoverage(qseaSet,
                                 uniquePos = FALSE,
                                 paired = TRUE,
                                 parallel = parallel,
                                 minMapQual = minMapQual
    )

    qseaSet <- addMedipsEnrichmentFactors(qseaSet, nCores = ifelse(parallel, BiocParallel::bpworkers(), 1), nonEnrich = FALSE)

  } else if (coverageMethod == "PairedAndR1s") {

    # load the coverage from each bam file, including using R1s from high MAPQ reads that aren't in perfect pairs.
    qseaSet <- addBamCoveragePairedAndUnpaired(qseaSet,
                                               fragmentLength = fragmentLength,
                                               parallel = parallel,
                                               minMapQual = minMapQual,
                                               minReferenceLength = minReferenceLength,
                                               maxInsertSize = maxInsertSize,
                                               minInsertSize = minInsertSize,
                                               properPairsOnly = properPairsOnly
    )
  } else {stop(glue::glue("Unknown coverageMethod {coverageMethod}. Use PairedAndR1s or qseaPaired."))}

  qseaSet@parameters$minInsertSize <- minInsertSize
  qseaSet@parameters$maxInsertSize <- maxInsertSize
  qseaSet@parameters$properPairsOnly <- properPairsOnly
  qseaSet@parameters$minReferenceLength <- minReferenceLength
  qseaSet@parameters$minMapQual <- minMapQual

  numEmpty <- qseaSet@libraries$file_name %>%
    as.data.frame() %>%
    dplyr::select(total_fragments) %>%
    dplyr::filter(total_fragments == 0)

  if (nrow(numEmpty) > 0) {stop(glue::glue("Empty sample: {rownames(numEmpty)}"))}

  if (CNVmethod == "qseaInput") {

    # calculate the CNV, using the input files included in the sampleTable CSV.
    # uses HMMCopy behind the scenes
    # note that apparently if you give any samples with normal or control in the name it will try and use those to normalise!
    qseaSet <- qsea::addCNV(qseaSet,
                            file_name = "input_file",
                            window_size = CNVwindowSize,
                            fragment_length = fragmentLength,
                            paired = TRUE,
                            parallel = TRUE,
                            MeDIP = FALSE
    )

    #this is included in the addHMMcopyCNV method more efficiently, don't need to call it there.
    qseaSet <- addMedipsEnrichmentFactors(qseaSet, nCores = ifelse(parallel, BiocParallel::bpworkers(), 1), nonEnrich = TRUE)

  } else if (CNVmethod == "HMMdefault") {

    if(is.null(hmmCopyGC) && BSgenome == "BSgenome.Hsapiens.NCBI.GRCh38") {

      if(CNVwindowSize == 1000000) {
        hmmCopyGC <- gc_hg38_1000kb
      } else if (CNVwindowSize == 500000) {
        hmmCopyGC <- gc_hg38_500kb
      } else if (CNVwindowSize == 50000) {
        hmmCopyGC <- gc_hg38_50kb
      } else {
        stop("Please supply gc data for this CNVwindowSize via the hmmCopyGC argument")
      }

    }

    if(is.null(hmmCopyGC) && BSgenome == "BSgenome.Hsapiens.UCSC.hg38") {
      
      if(CNVwindowSize == 1000000) {
        hmmCopyGC <- gc_hg38_1000kb %>%
          dplyr::mutate(chr = paste0("chr",chr))
      } else if (CNVwindowSize == 500000) {
        hmmCopyGC <- gc_hg38_500kb %>%
          dplyr::mutate(chr = paste0("chr",chr))
      } else if (CNVwindowSize == 50000) {
        hmmCopyGC <- gc_hg38_50kb %>%
          dplyr::mutate(chr = paste0("chr",chr))
      } else {
        stop("Please supply gc data for this CNVwindowSize via the hmmCopyGC argument")
      }
      
    }    
    
    if(is.null(hmmCopyMap) && BSgenome == "BSgenome.Hsapiens.NCBI.GRCh38") {
      if(CNVwindowSize == 1000000) {
        hmmCopyMap <- map_hg38_1000kb
      } else if (CNVwindowSize == 500000) {
        hmmCopyMap <- map_hg38_500kb
      } else if (CNVwindowSize == 50000) {
        hmmCopyMap <- map_hg38_50kb
      } else {
        stop("Please supply mapability data for this CNVwindowSize via the hmmCopyGC argument")
      }

    }

    if(is.null(hmmCopyMap) && BSgenome == "BSgenome.Hsapiens.UCSC.hg38") {
      if(CNVwindowSize == 1000000) {
        hmmCopyMap <- map_hg38_1000kb %>%
          dplyr::mutate(chr = paste0("chr",chr))
      } else if (CNVwindowSize == 500000) {
        hmmCopyMap <- map_hg38_500kb %>%
          dplyr::mutate(chr = paste0("chr",chr))
      } else if (CNVwindowSize == 50000) {
        hmmCopyMap <- map_hg38_50kb %>%
          dplyr::mutate(chr = paste0("chr",chr))
      } else {
        stop("Please supply mapability data for this CNVwindowSize via the hmmCopyGC argument")
      }
      
    }    
    
      if (is.null(hmmCopyGC)) {
      stop("No hmmCopy GC content object provided!")
      } else {
        requiredColumns <- c("chr","start","end","gc")
        columnDiff <- setdiff(colnames(hmmCopyGC), requiredColumns)
      if(length(columnDiff) > 0) {
        stop(glue::glue("hmmCopyGC object missing required columns: {paste(columnDiff, collapse = ' ')}"))
      } 
      objectWindowSize <- hmmCopyGC %>% mutate(size = end - start + 1) %>% pull(size) %>% unique()
      
      if(length(objectWindowSize) != 1) {
        stop(glue::glue("hmmCopyGC object should have constant window size for all windows."))
      } else if (objectWindowSize != CNVwindowSize) {
        stop(glue::glue("hmmCopyGC object window size should be the same as the CNVwindowSize."))
      }
    }

    if (is.null(hmmCopyMap)) {
      stop("No hmmCopy Mapability file provided!")
    } else {
      requiredColumns <- c("chr","start","end","map")
      columnDiff <- setdiff(colnames(hmmCopyMap), requiredColumns)
      if(length(columnDiff) > 0) {
        stop(glue::glue("hmmCopyMap object missing required columns: {paste(columnDiff, collapse = ' ')}"))
      } 
      objectWindowSize <- hmmCopyMap %>% mutate(size = end - start + 1) %>% pull(size) %>% unique()
      
      if(length(objectWindowSize) != 1) {
        stop(glue::glue("hmmCopyMap object should have constant window size for all windows."))
      } else if (objectWindowSize != CNVwindowSize) {
        stop(glue::glue("hmmCopyMap object window size should be the same as the CNVwindowSize."))
      }
    }

    # use HMMCopy directly, with default settings
    qseaSet <- addHMMcopyCNV(qseaSet,
                             inputColumn = "input_file",
                             windowSize = CNVwindowSize,
                             fragmentLength = fragmentLength,
                             minMapQual = minMapQual,
                             minReferenceLength = minReferenceLength,
                             maxInsertSize = maxInsertSize,
                             minInsertSize = minInsertSize,
                             properPairsOnly = properPairsOnly,
                             hmmCopyGC = hmmCopyGC,
                             hmmCopyMap = hmmCopyMap,
                             parallel = parallel)
  }  else if (CNVmethod == "MeCap") {
    message("Adding CNV using MeCap samples")

    qseaSet <- qsea::addCNV(qseaSet,
                            file_name = "file_name",
                            window_size = CNVwindowSize,
                            fragment_length = fragmentLength,
                            paired = TRUE,
                            parallel = parallel,
                            MeDIP = TRUE
    )

    #qseaSet <- addMedipsEnrichmentFactors(qseaSet, nCores = ifelse(parallel, BiocParallel::bpworkers(), 1), nonEnrich = TRUE)


  } else if (CNVmethod == "None") {
    message("No CNV being calculated")
  }

  # calculate the average CG density across the genome. Does not use the reads at all.
  qseaSet <- qsea::addPatternDensity(qseaSet,
                                     pattern = "CG",
                                     name = "CpG",
                                     fragment_length = fragmentLength,
                                     fragment_sd = fragmentSD,
                                     masks = NULL)

  # store the values of the fragment length and SD used into the object
  qseaSet@parameters$fragmentLength <- fragmentLength
  qseaSet@parameters$fragmentSD <- fragmentSD

  #do not set library factors via TMM, just set to be 1. Makes no difference to beta values as gets normalised out anyway, only nrpms are affected.

  qseaSet <- addNormalisation(qseaSet, enrichmentMethod = enrichmentMethod, maxPatternDensity = maxPatternDensity)
  
  qseaSet@parameters$coverageMethod <- coverageMethod
  qseaSet@parameters$cnvMethod <- CNVmethod

  message("qseaSet object generated successfully")

  return(qseaSet)
}
