#' This function takes a bam file and gets the coverage over a set of regions. It keeps pairs or high quality R1s.
#' @param fileName Path to the bam file
#' @param BSgenome String denoting the BSgenome to be used.
#' @param regions GRanges object with the regions to calculate the coverage for.
#' @param fragmentLength Length to extend unpaired reads to.
#' @param minMapQual Minimum mapping quality to include a read (on either end if proper pair)
#' @param maxInsertSize Maximum insert size for a proper paired read
#' @param minInsertSize Minimum insert size for a proper paired read
#' @param properPairsOnly Whether to only keep proper pairs, else high quality R1s are kept. Also rescues pairs that are almost proper pairs.
#' @param minReferenceLength Minimum length on the reference coordinates for a non-paired read to be included
#' @return A list with two elements, the first containing the genomic regions with read numbers, the second the distribution related parameters.
#'
getBamCoveragePairedAndUnpairedR1 <- function(fileName = NULL, BSgenome = NULL, regions = NULL, fragmentLength = NULL,
                                              minMapQual = 30, maxInsertSize = 1000, minInsertSize = 50,
                                              minReferenceLength = 30, properPairsOnly = FALSE){

  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    stop(
      "Package \"Rsamtools\" must be installed to use this function.",
      call. = FALSE
    )
  }

  #sampleName <- stringr::str_remove(basename(fileName),"_MeCap.bam")

  message(glue::glue("Reading {basename(fileName)}, with minimum mapping quality {minMapQual}."))

  #if (is.null(fragmentLength) ){
  #  stop("Please provide an estimate for the fragment length")
  #}

  if (!tools::file_ext(fileName) %in% c("bam","BAM")) {
    stop(glue::glue("Please supply a bam file, not {file_name}"))
  }

  if(properPairsOnly){

    myParam <- Rsamtools::ScanBamParam(what = c("qname","flag","strand","rname","pos",
                                              "mapq","cigar","isize"), tag = "MQ")
  } else {

    myParam <- Rsamtools::ScanBamParam(what = c("flag","strand","rname","pos",
                                                "mapq","cigar","isize"), tag = "MQ")

  }

  readBamData <- Rsamtools::scanBam(file = fileName, param = myParam)

  if (is.null(readBamData[[1]]$tag$MQ)) {
    #MQ is the mate quality tag, coming from samtools fixmate. If not present, then add as NA.
    message(glue::glue("No samtools fixmate MQ (mate quality) tags on the bam file, using R1 MAPQ only."))
    readBamData[[1]]$tag$MQ <- rep(NA, length(readBamData[[1]]$flag))
  }

  readDF <- readBamData %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::bind_cols(., tibble::as_tibble(Rsamtools::bamFlagAsBitMatrix(.$flag)))

  #remove the chr if present in the bam file but not in the regions
  if(!any(stringr::str_detect(regions %>% GenomeInfoDb::seqlevels(),"chr"))) {
    readDF <- readDF  %>%
      dplyr::mutate(rname = stringr::str_replace(rname,"chrM","chrMT"),
                    rname = stringr::str_remove(rname,"chr")
      )

  }

  numInitialReads <- nrow(readDF)

  rm(readBamData)

  # Finds "proper pairs" with at least one read over the quality threshold.
  properPairsGRanges <- readDF %>%
    dplyr::filter(isProperPair == 1, isUnmappedQuery == 0, isSupplementaryAlignment == 0,
                  isSecondaryAlignment == 0, isDuplicate == 0, hasUnmappedMate == 0,
                  isNotPassingQualityControls == 0) %>% # get proper paired reads
    dplyr::filter(strand == "+") %>% # keep only the part of the pair that is on the positive strand
    dplyr::mutate(seqnames = rname, start = pos, end = pos + isize - 1) %>%
    plyranges::as_granges()

  if(properPairsOnly){

  shouldBeProperPairsGRanges <- readDF %>%
    dplyr::filter(isProperPair == 0, isUnmappedQuery == 0, isSupplementaryAlignment == 0,
           isSecondaryAlignment == 0, isDuplicate == 0, hasUnmappedMate == 0,
           isNotPassingQualityControls == 0,
           abs(isize) <= maxInsertSize, abs(isize) >= minInsertSize) %>%
    dplyr::group_by(qname) %>% #this grouping is a expensive operation, so filter as much as possible beforehand.
    dplyr::mutate(nStrands = length(unique(strand)),
           nSign = length(sign(isize)),
           inOut = sum(ifelse(strand == "+", 1, -1) * sign(isize))) %>%
    dplyr::filter(inOut == 2, nSign == 2, nStrands == 2) %>% #filter to just reads that meet the proper pair criteria
    dplyr::ungroup() %>%
    dplyr::filter(strand == "+") %>% # keep only the part of the pair that is on the positive strand
    dplyr::mutate(seqnames = rname, start = pos, end = pos + isize - 1) %>%
    plyranges::as_granges() %>%
    dplyr::filter(MQ >= minMapQual | mapq >= minMapQual)

    message(glue::glue("Rescued {length(shouldBeProperPairsGRanges)} paired reads that met all but the length criteria of being proper paired reads (but 70-1000bp)"))

  properPairsGRanges <- properPairsGRanges %>%
    plyranges::bind_ranges(shouldBeProperPairsGRanges)

  }

  numPairsInit <- length(properPairsGRanges)

  properPairsGRanges <- properPairsGRanges %>%
    dplyr::filter(MQ >= minMapQual | mapq >= minMapQual)  #keep proper pair fragments if either mapq passes the quality

  numPairsAfterMAPQ <- length(properPairsGRanges)

  if (numPairsAfterMAPQ < numPairsInit) {
    message(glue::glue("Dropping {numPairsInit - numPairsAfterMAPQ} proper pairs that have neither read with mapq >= {minMapQual}"))
  }

  #Filters out ones that are too big.
  properPairsGRanges <- properPairsGRanges %>%
    dplyr::filter(abs(isize) <= maxInsertSize,
                  abs(isize) >= minInsertSize)

  numPairsWithinSize <- length(properPairsGRanges)

  if (numPairsWithinSize < numPairsInit) {
    message(glue::glue("Dropping {numPairsAfterMAPQ - numPairsWithinSize} proper pairs that do not lie within [{minInsertSize}, {maxInsertSize}]"))
  }

  R1GRanges <- readDF %>%
    dplyr::filter(isProperPair == 0, isUnmappedQuery == 0, isFirstMateRead == 1 | isPaired == 0, isSupplementaryAlignment == 0,
                  isSecondaryAlignment == 0, isDuplicate == 0, isNotPassingQualityControls == 0) %>% #get R1s that aren't proper pairs
    dplyr::mutate(seqnames = rname, start = pos, end = pos) %>%
    plyranges::as_granges()

  numR1sInit <- length(R1GRanges)

  R1GRanges <- R1GRanges %>%
    dplyr::filter(mapq >= minMapQual,  #require good mapq and at least minReferenceLength distance mapped onto the genome
                  GenomicAlignments::cigarWidthAlongReferenceSpace(cigar) >= minReferenceLength)

  numR1sMAPQAndLength <- length(R1GRanges)

  rm(readDF)

  if (!properPairsOnly) {

    message("Using high quality R1 unpaired reads as well as pairs")

    if (is.null(fragmentLength) | length(properPairsGRanges) >= length(R1GRanges) ) {
      fragmentLength = stats::median(BiocGenerics::width(properPairsGRanges))
      message(glue::glue("Estimating fragment length for {length(R1GRanges)} unpaired R1 reads as {fragmentLength}, based on {length(properPairsGRanges)} proper pairs."))
    } else{
      message(glue::glue("Setting unpaired fragment lengths to {fragmentLength}."))
    }

    R1GRanges <- R1GRanges %>%
      plyranges::anchor_5p() %>%
      plyranges::stretch(extend = fragmentLength)

    readGRanges <-  c(properPairsGRanges, R1GRanges)

  } else{

    message("Ignoring R1 unpaired reads, using proper pairs only")
    readGRanges <- properPairsGRanges
  }



  regions <- regions %>%
    tibble::as_tibble() %>%
    tibble::rowid_to_column("id") %>%
    plyranges::as_granges()

  # regionAvgs <- regions %>%
  #   plyranges::join_overlap_left(readGRanges, .) %>%
  #   tibble::as_tibble() %>%
  #   dplyr::group_by(id) #%>%
  #dplyr::summarise(avgFragmentLength = median(width))

  # Count the number of reads within each region. Use the midpoints to be consistent with qsea.
  midPointCounts <- readGRanges %>%
    dplyr::mutate(start = (start + end)/2, end = start ) %>%
    plyranges::count_overlaps(regions,.)

  regionsOut <- regions %>%
    tibble::as_tibble() %>%
    dplyr::mutate(counts = midPointCounts) %>%
    #dplyr::left_join(regionAvgs, by = "id") %>%
    #dplyr::mutate(nProperPairs = tidyr::replace_na(nProperPairs, 0) ) %>%
    plyranges::as_granges()

  enrichment = calculateCGEnrichmentGRanges(readGRanges,
                                             BSgenome,
                                             chr.select = regions %>% GenomeInfoDb::seqinfo() %>% GenomeInfoDb::seqnames())

  distributionParameters <- data.frame(
    qsea_initial_reads = numInitialReads,
    qsea_initial_pairs = numPairsInit,
    qsea_initial_r1s = numR1sInit,
    qsea_mapq_filtered_pairs = numPairsAfterMAPQ,
    qsea_size_filtered_pairs = numPairsWithinSize,
    qsea_filtered_r1s = numR1sMAPQAndLength,
    total_fragments = length(readGRanges), # total number of valid reads
    valid_fragments = sum(midPointCounts), #  number of reads in the regions under consideration
    library_factor = NA, # these two are filled in later by qsea
    offset = NA,
    fragment_median = readGRanges %>% width() %>% stats::median(),
    fragment_mean = readGRanges %>% width() %>% mean(),
    fragment_sd = readGRanges %>% width() %>% stats::sd(),
    relH = enrichment$relH,
    GoGe = enrichment$GoGe,
    fragments_without_pattern = enrichment$nReadsWithoutPattern,
    prop_without_pattern = round(enrichment$nReadsWithoutPattern/length(readGRanges), 3)
  )

  if (properPairsOnly) {
    distributionParameters <- distributionParameters %>%
      dplyr::select(-tidyselect::matches("r1s"))
  }

  return(list(regionsGR = regionsOut, library = distributionParameters))

}



#' This function adds coverage to the qseaSet, using both paired reads and the high-quality R1s.
#' @param qs qseaSet object to add coverage to.
#' @param fragmentLength Length to extend unpaired reads to
#' @param minMapQual Minimum mapping quality to include a read (on either end if proper pair)
#' @param maxInsertSize Maximum insert size for a proper paired read
#' @param minInsertSize Minimum insert size for a proper paired read
#' @param properPairsOnly Whether to only keep proper pairs, else high quality unpaired R1s are also kept.
#' @param minReferenceLength Minimum length on the reference coordinates for a non-paired read to be included
#' @param parallel Whether to use parallelisation
#' @return A qseaSet with the coverage added.
#' @export
#'
#'
addBamCoveragePairedAndUnpaired <- function(qs,
                                            fragmentLength = NULL,
                                            minMapQual = 30,
                                            maxInsertSize = 1000,
                                            minInsertSize = 100,
                                            minReferenceLength = 30,
                                            parallel = TRUE,
                                            properPairsOnly = FALSE) {
  sampleTable = qsea::getSampleTable(qs)
  Regions = qsea::getRegions(qs)

  if (parallel & (BiocParallel::bpworkers() > 1)) {
    BPPARAM = BiocParallel::bpparam()
    message("Scanning up to ", BiocParallel::bpnworkers(BPPARAM), " files in parallel.")
  } else {
    BPPARAM = BiocParallel::SerialParam()
    message("Scanning one file at a time. Use e.g. register(MulticoreParam(workers = 4)) to parallelise this process.")
    }

  getCGPositions(qsea:::getGenome(qs), qs %>% qsea::getRegions()  %>%
                   GenomeInfoDb::seqnames() %>%
                   levels()) #memoise this once, for using in all the getBamCoveragePairedAndUnpairedR1 function.

  bamOutList <- BiocParallel::bplapply(X = sampleTable[, "file_name"],
                                       FUN = getBamCoveragePairedAndUnpairedR1,
                                       BSgenome = qsea:::getGenome(qs),
                                       regions = Regions,
                                       fragmentLength = fragmentLength, maxInsertSize = maxInsertSize,
                                       minInsertSize = minInsertSize, minReferenceLength = minReferenceLength,
                                       minMapQual = minMapQual, BPPARAM = BPPARAM, properPairsOnly = properPairsOnly)

  names(bamOutList) <- sampleTable$sample_name

  message("Generated objects")

  coverage =  bamOutList %>%
    purrr::map_dfc(function(x) {
      purrr::pluck(x,"regionsGR") %>%
        tibble::as_tibble() %>%
        dplyr::pull(counts)  }
    ) %>% as.matrix()
  colnames(coverage) <- sampleTable$sample_name

  message("Calculated coverage")

  # regionAvgs =  bamOutList %>%
  #   purrr::map_dfc(function(x) {
  #     purrr::pluck(x,"regionsGR") %>%
  #       tibble::as_tibble() %>%
  #       dplyr::select(avgFragmentLength, avgFragmentMAPQ)
  #     }
  #   ) %>% as.matrix()
  #
  # colnames(regionAvgs) <- paste0(rep(sampleTable$sample_name,rep(2,nrow(sampleTable))), "_",
  #                                c("avgFragmentLength", "avgFragmentMAPQ"))
  #
  #   message("Calculated read coverage")

  libraries = bamOutList %>% purrr::map_dfr(~ purrr::pluck(.,"library")) %>% as.data.frame()
  rownames(libraries) <- sampleTable$sample_name

  param = list(minMapQual = minMapQual, minReferenceLength = minReferenceLength,
               maxInsertSize = maxInsertSize, minInsertSize = minInsertSize,
               properPairsOnly = properPairsOnly)
  qs = qsea:::addParameters(qs, param)
  qs = qsea:::setCounts(qs, count_matrix = coverage)
  qs = qsea:::setLibrary(qs, "file_name", libraries)

  #GenomicRanges::mcols(qs@regions) <- cbind(GenomicRanges::mcols(qs@regions), regionAvgs)


  return(qs)
}




#' This function takes a qseaSet and adds several normalisation steps from qsea, with default values
#' @param qseaSet The qseaSet object.
#' @param enrichmentMethod What method to use to calculate the enrichment step
#' @param maxPatternDensity Maximum pattern density in a window to consider it for the background calculation
#' @return A qseaSet object with the
#' @examples
#' getExampleQseaSet(expSamplingDepth = 100000) %>% addNormalisation(maxPatternDensity = 0.5)
#' @export
#'
addNormalisation <- function(qseaSet, enrichmentMethod = "blind1-15", maxPatternDensity = 0.05){

  # do not calculate the TMM normalisation, set all to be 0.
  qseaSet <- qsea::addLibraryFactors(qseaSet, factors = 1)

  qseaSet <- qsea::addOffset(qseaSet, enrichmentPattern = getPattern(qseaSet), maxPatternDensity = maxPatternDensity)

  wd <- which(qsea::getRegions(qseaSet)$CpG_density >= 1 &
                qsea::getRegions(qseaSet)$CpG_density <= 15)
  signal <- (15 - qsea::getRegions(qseaSet)$CpG_density[wd])*.55/15 + .25

  qseaSet <- qsea::addEnrichmentParameters(qseaSet,
                                           enrichmentPattern = getPattern(qseaSet),
                                           windowIdx = wd,
                                           signal = signal,
                                           min_wd = 5,
                                           bins = seq(.5,40.5,0.5)
  )

  # store the enrichment method used
  qseaSet@parameters$enrichmentMethod <- "blind1-15"

  return(qseaSet)
}
