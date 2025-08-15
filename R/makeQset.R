#' Construct an initial qseaSet from BAMs and metadata
#'
#' Builds a `qseaSet` from a `sampleTable` that lists BAM files and sample
#' metadata, lays out the genome into fixed windows, removes blacklisted regions,
#' and then loads coverage and (optionally) CNV information according to the
#' selected methods. Typical use is the first step in a mesa/qsea workflow.
#'
#' @details
#' **Inputs & requirements**
#' - `sampleTable` must include at least the columns:
#'   - `sample_name` (unique per sample),
#'   - `file_name` (path to MeDIP/capture BAM),
#'   - `group` (string label).  
#'   If `CNVmethod` uses inputs (`"HMMdefault"`, `"qseaInput"`), the table must
#'   also contain `input_file` (path to matched input BAM).
#' - BAMs should be coordinate-sorted; for paired data it is recommended to run
#'   `samtools fixmate` so MAPQ tags are available (see `minMapQual` note below).
#' - `BSgenome` must be a valid genome package string (e.g.
#'   `"BSgenome.Hsapiens.NCBI.GRCh38"`). The genome provides chromosome lengths
#'   and sequence info used to tile windows.
#'
#' **Windows & blacklist**
#' - Chromosomes in `chrSelect` are tiled into non-overlapping windows of size
#'   `windowSize` (bp). Windows overlapping `badRegions` (if provided) are removed.
#'
#' **Coverage loading**
#' - `coverageMethod = "PairedAndR1s"` loads paired reads and allows high-MAPQ
#'   singleton R1s from properly paired fragments (controlled by `minMapQual`,
#'   `minInsertSize`, `maxInsertSize`, `properPairsOnly`, `minReferenceLength`).
#' - `coverageMethod = "qseaPaired"` delegates to `qsea::addCoverage()` with paired
#'   settings.
#' - `fragmentType` can set default `fragmentLength`/`fragmentSD` for `"Sheared"`
#'   or `"cfDNA"`; otherwise both must be supplied explicitly.
#'
#' **CNV options**
#' - `"qseaInput"`: runs `qsea::addCNV()` using `input_file` BAMs.
#' - `"HMMdefault"`: uses mesa's `addHMMcopyCNV()` with supplied or default GC/
#'   mappability tracks for GRCh38 (objects like `gc_hg38_1000kb`, `map_hg38_1000kb`,
#'   chosen by `CNVwindowSize`).  
#' - `"MeCap"`: estimates CNV from MeCap BAMs (`file_name`).  
#' - `"None"`: no CNV is computed.
#'
#' **Parallelisation**
#' - If `parallel = TRUE`, computation uses the currently registered
#'   BiocParallel backend (set with `BiocParallel::register()`); otherwise runs
#'   serially. `setMesaParallel()` can help toggle this package-wide.
#'
#' **Additional processing**
#* - Adds CpG density (`qsea::addPatternDensity()` with pattern `"CG"`),
#'   enrichment parameters, sets library factors to 1 (so betas are unchanged by
#'   library scaling), and stores key parameters in `qseaSet@parameters`.
#'
#' **Errors & validations**
#' - The function checks required columns, existence of listed BAMs, and presence
#'   of required GC/mappability tracks for `"HMMdefault"`; it stops with a
#'   descriptive message if any prerequisite is missing.
#'
#' @param sampleTable A data frame with at least `sample_name`, `file_name`,
#'   and `group`; if `CNVmethod` requires inputs, also `input_file`. Additional
#'   columns are kept as sample metadata.
#' @param BSgenome Character. BSgenome package string (e.g.
#'   `"BSgenome.Hsapiens.NCBI.GRCh38"`). See [BSgenome::available.genomes()].
#' @param chrSelect Integer or character vector of chromosomes to include (default `1:22`).
#' @param windowSize Integer window size in bp for tiling the genome (default `300`).
#' @param fragmentType Character. Either `"Sheared"` or `"cfDNA"` to set default
#'   `fragmentLength`/`fragmentSD`. If `NULL`, both must be supplied manually.
#' @param fragmentLength Numeric. Average DNA fragment length (bp).
#' @param fragmentSD Numeric. Standard deviation of fragment length (bp).
#' @param CNVwindowSize Integer CNV bin size (bp), default `1e6`.
#' @param CNVmethod Character. One of `"HMMdefault"`, `"qseaInput"`, `"MeCap"`, or `"None"`.
#' @param coverageMethod Character. One of `"PairedAndR1s"` (mesa) or `"qseaPaired"` (qsea).
#' @param minMapQual Integer. Minimum MAPQ to keep a read. For `PairedAndR1s`, a pair is
#'   kept if either R1 or R2 passes when properly paired and MAPQ tags are set.
#' @param minInsertSize,maxInsertSize Integer bounds for paired insert sizes
#'   (bp). Useful for cfDNA size selection.
#' @param properPairsOnly Logical; if `TRUE`, keep only properly paired reads
#'   (recommended for strict size selection).
#' @param minReferenceLength Integer. Minimum aligned length (bp) to keep a read.
#' @param badRegions Optional `GRanges` with regions to exclude (blacklist).
#' @param hmmCopyGC data.frame or `NULL`. GC content per CNV bin (size = `CNVwindowSize`) for HMMcopy.
#' @param hmmCopyMap data.frame or `NULL`. Mappability per CNV bin (size = `CNVwindowSize`) for HMMcopy.
#' @param parallel Logical; use the registered BiocParallel backend (`TRUE`) or
#'   run serially (`FALSE`). See [BiocParallel::register()].
#'
#' @return
#' A `qseaSet` containing:
#' - tiled and blacklisted genomic windows,
#' - the input `sampleTable` (as sample metadata),
#' - loaded coverage according to `coverageMethod`,
#' - CNV results depending on `CNVmethod`,
#' - CpG density and enrichment parameters, and
#' - recorded processing parameters in `@parameters`.
#'
#' @seealso
#' [qsea::createQseaSet()], [qsea::addCoverage()], [qsea::addCNV()],
#' [addHMMcopyCNV()], [setMesaParallel()], [BSgenome::available.genomes()]
#'
#' @examples
#' \donttest{
#' ## Skeleton workflow (requires existing BAMs; paths below are examples)
#' # sampleTable <- data.frame(
#' #   sample_name = c("S1_T", "S1_N"),
#' #   file_name   = c("/path/to/S1_T.meCap.bam", "/path/to/S1_N.meCap.bam"),
#' #   input_file  = c("/path/to/S1_T.input.bam", "/path/to/S1_N.input.bam"),
#' #   group       = c("Tumour", "Normal"),
#' #   stringsAsFactors = FALSE
#' # )
#' # qs <- makeQset(
#' #   sampleTable   = sampleTable,
#' #   BSgenome      = "BSgenome.Hsapiens.NCBI.GRCh38",
#' #   chrSelect     = 1:22,
#' #   windowSize    = 300,
#' #   fragmentType  = "cfDNA",   # sets fragmentLength/fragmentSD defaults
#' #   CNVmethod     = "HMMdefault",
#' #   coverageMethod= "PairedAndR1s",
#' #   parallel      = FALSE
#' # )
#' }
#'
#' @export
makeQset <- function(sampleTable,
                     BSgenome = NULL,
                     chrSelect = 1:22,
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
                     parallel = getMesaParallel()) {

  if(parallel) {
    if(BiocParallel::bpworkers() == 1){
      message("No configured parallelisation, use e.g. register(MulticoreParam(workers = 4)) to process multiple files at once.")
      parallel = FALSE
    } else {
      message(glue::glue("Detected parallel setup with {BiocParallel::bpworkers()} workers."))
    }

  }

  if (!is.null(fragmentType)) {
    if (fragmentType %in% c("Sheared","sheared") ) {
      fragmentLength = 213
      fragmentSD = 60
    } else if (fragmentType == "cfDNA") {
      fragmentLength = 167
      fragmentSD = 38
    } else {
      stop("fragmentType should be either Sheared or cfDNA.")}

  }

  if (is.null(fragmentLength) | is.null(fragmentSD)) {
    stop("fragmentLength and fragmentSD must be specified, or fragmentType can be specified for some defaults.")
  }

  if (!("sample_name" %in% colnames(sampleTable)))  {
    stop("Required column sample_name not included in the sampleTable.")
  }

  if (!("file_name" %in% colnames(sampleTable))) {
    stop("Required column file_name not included in the sampleTable.")
  }

  if (!("group" %in% colnames(sampleTable)))  {
    stop("Required column group included in the sampleTable.")
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
  # chromosome lengths, and then a Seqinfo object with that
  chr_length <- GenomeInfoDb::seqlengths(refGenome)[chrSelect]
  seqinfo <- GenomeInfoDb::Seqinfo(as.character(chrSelect),chr_length, NA, BSgenome)

  # number of windows of size windowSize on each chromosome
  nr_wd <- floor(chr_length/windowSize)

  # starting point of each window on each chromosome, then make a GRanges object with all those windows
  wd_start <- unlist(lapply(FUN = seq, X = chr_length - windowSize + 1,
                            from = 1, by = windowSize), FALSE, FALSE)

  windowsGRanges <- GenomicRanges::GRanges(seqnames = rep(factor(chrSelect), nr_wd),
                                           ranges = IRanges::IRanges(start = wd_start,
                                                                     width = windowSize),
                                           seqinfo = seqinfo)

  # remove both sets of blacklisted windows from the full set of windows
  windowsWithoutBlacklist <- windowsGRanges %>%
    plyranges::filter_by_non_overlaps(badRegions)

  print(paste0("Considering ", length(windowsWithoutBlacklist), " regions with total size ",
               sum(BiocGenerics::width(windowsWithoutBlacklist))))

  #make the initial Qsea object, with the reduced set of windows, using the sampleTable
  qseaSet <- qsea::createQseaSet(sampleTable = sampleTable,
                                 BSgenome = BSgenome,
                                 chr.select = chrSelect,
                                 Regions = windowsWithoutBlacklist,
                                 window_size = windowSize)

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

      if (is.null(hmmCopyGC)) {
      stop("No hmmCopy GC content object provided!")
    }

    if (is.null(hmmCopyMap)) {
      stop("No hmmCopy Mapability file provided!")
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
                            paired = FALSE,
                            parallel = parallel,
                            MeDIP = TRUE
    )

    qseaSet <- addMedipsEnrichmentFactors(qseaSet, nCores = ifelse(parallel, BiocParallel::bpworkers(), 1), nonEnrich = TRUE)


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
  #if you do use TMM, then it depends on the samples in the qseaSet (with one being the reference to compare with)
  qseaSet <- qsea::addLibraryFactors(qseaSet, factors = 1)

  qseaSet <- qsea::addOffset(qseaSet, enrichmentPattern = "CpG", maxPatternDensity = 0.05)

  wd <- which(qsea::getRegions(qseaSet)$CpG_density >= 1 &
                qsea::getRegions(qseaSet)$CpG_density <= 15)
  signal <- (15 - qsea::getRegions(qseaSet)$CpG_density[wd])*.55/15 + .25

  qseaSet <- qsea::addEnrichmentParameters(qseaSet,
                                           enrichmentPattern = "CpG",
                                           windowIdx = wd,
                                           signal = signal,
                                           min_wd = 5,
                                           bins = seq(.5,40.5,0.5)
  )

  # store the enrichment method used
  qseaSet@parameters$enrichmentMethod <- "blind1-15"
  qseaSet@parameters$coverageMethod <- coverageMethod
  qseaSet@parameters$cnvMethod <- CNVmethod

  print("qseaSet object generated successfully")

  return(qseaSet)
}
