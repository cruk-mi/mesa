#' Make the initial qseaSet object from a set of samples
#'
#' This function makes the initial qseaSet from a sampleTable that includes the file paths for the bam files and their metadata
#'
#' @param sampleTable A data frame with the sample information to be passed to qsea.
#' @param chrSelect Which chromosomes to use (default 1:22).
#' @param windowSize What window size (in bp) to use in the genome (default 300).
#' @param fragmentType What type of procedure generated the library, Sheared or cfDNA (used to set the average fragment length/SD to a default for CBC)
#' @param fragmentLength Average DNA fragment length. Can be set by fragmentType.
#' @param fragmentSD Standard deviation of the DNA fragment lengths. Can be set by fragmentType.
#' @param CNVwindowSize What window size (in bp) to use in the calculation of CNV (default 1000000).
#' @param CNVmethod Which method to use for calculation of CNV. Options include "HMMdefault" (hmmcopy with default parameters) and "None"
#' @param coverageMethod Whether to use custom methtools method for reading coverage in (set to PairedAndR1s), rather than qsea's (qseaPaired).
#' @param minMapQual Minimum MAPQ score for a read to be kept. For PairedAndR1s, will keep if either R1 or R2 meet this cutoff in a properly paired read (if MQ tags are set in the bam file with samtools fixmate).
#' @param minInsertSize For paired reads, only keep them if they are above a minimum length. Can be used for cfDNA size selection. Applies to Input samples as well as MeCap.
#' @param maxInsertSize For paired reads, only keep them if they are below a maximum length. Can be used for cfDNA size selection. Applies to Input samples as well as MeCap.
#' @param properPairsOnly Whether to only keep properly paired reads, or to keep high-quality (MAPQ 30+) unpaired R1s as well. Set to TRUE for size selection.
#' @param minReferenceLength A minimum distance on the genome to keep the read. bwa by default gives 19bp as minimum for a read, which is quite short.
#' @param badRegions A GRanges object containing regions to filter out from the result. Defaults to the ENCODE blocklist for hg38.
#' @param badRegions2 A GRanges object containing additional regions to filter out from the result.
#' @param nCores How many cores to use to read the data in. Set as high as possible on your machine for speedy results.
#' @return A qseaSet object, containing all the information required.
#' @export
makeQset <- function(sampleTable,
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
                     badRegions = encodeBlacklist,
                     badRegions2 = NULL,
                     properPairsOnly = FALSE,
                     nCores = 1) {

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
    stop("fragmentLength and fragmentSD must be specified, or fragmentType")
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
    stop(glue::glue(" MeCap file not found for: {filter(sampleTable,!file.exists(file_name)) %>% pull(sample_name)}. "))
  }

  if(!all(file.exists(sampleTable$input_file))){
    stop(glue::glue(" Input file not found for: {filter(sampleTable,!file.exists(input_file)) %>% pull(input_file)}. "))
  }


  # Additional file of bad mapping regions to remove
  ## TODO allow to be a GRanges object or bed

  if (!is.null(badRegions2)) {
    extraBlacklistGRanges <- plyranges::read_bed(badRegions2)

  } else {
    extraBlacklistGRanges <- GenomicRanges::GRanges()
  }

  # Get data for hg38
  refGenome <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
  # chromosome lengths, and then a Seqinfo object with that
  chr_length <- GenomeInfoDb::seqlengths(refGenome)[chrSelect]
  seqinfo <- GenomeInfoDb::Seqinfo(as.character(chrSelect),chr_length, NA, "BSgenome.Hsapiens.NCBI.GRCh38")

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
    plyranges::filter_by_non_overlaps(badRegions) %>%
    plyranges::filter_by_non_overlaps(extraBlacklistGRanges)

  print(paste0("Considering ", length(windowsWithoutBlacklist), " regions with total size ",
               sum(BiocGenerics::width(windowsWithoutBlacklist))))

  #make the initial Qsea object, with the reduced set of windows, using the sampleTable
  qseaSet <- qsea::createQseaSet(sampleTable = sampleTable,
                                 BSgenome = "BSgenome.Hsapiens.NCBI.GRCh38",
                                 chr.select = chrSelect,
                                 Regions = windowsWithoutBlacklist,
                                 window_size = windowSize)

  # store the extra blacklist file location
  #qseaSet@parameters$badRegions2 <- blacklistBed

  if (coverageMethod == "qseaPaired") {

    #load the coverage from each bam file, using qsea default method
    qseaSet <- qsea::addCoverage(qseaSet,
                                 uniquePos = FALSE,
                                 paired = TRUE,
                                 parallel = TRUE,
                                 minMapQual = minMapQual
    )

    #this is included in the PairedAndR1s method more efficiently, don't need to call it there.
    qseaSet <- addMedipsEnrichmentFactors(qseaSet, nCores = nCores)

  } else if (coverageMethod == "PairedAndR1s") {

    # load the coverage from each bam file, including using R1s from high MAPQ reads that aren't in perfect pairs.
    qseaSet <- addBamCoveragePairedAndUnpaired(qseaSet,
                                               fragmentLength = fragmentLength,
                                               parallel = TRUE,
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
    dplyr::select(total_fragments) %>%
    dplyr::filter(total_fragments == 0)

  if (nrow(numEmpty) > 0) {stop(glue::glue("Empty sample: {rownames(numEmpty)}"))}

  if (CNVmethod == "qseaPooled") {

    # calculate the CNV, using the input files included in the sampleTable CSV.
    # uses HMMCopy behind the scenes, but with the PooledControl sample as a normalisation
    # note that apparently if you give any samples with normal or control in the name it will try and use those to normalise!
    qseaSet <- qsea::addCNV(qseaSet,
                            file_name = "input_file",
                            window_size = CNVwindowSize,
                            fragment_length = fragmentLength,
                            paired = TRUE,
                            parallel = TRUE,
                            MeDIP = FALSE
    )
  } else if (CNVmethod == "HMMdefault") {

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
                             parallel = TRUE)
  }  else if (CNVmethod == "MeCap") {
    message("Adding CNV using MeCap samples")

    qseaSet <- qsea::addCNV(qseaSet,
                            file_name = "file_name",
                            window_size = CNVwindowSize,
                            paired = FALSE,
                            parallel = TRUE,
                            MeDIP = TRUE
    )

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

  #do not set library factors via TMM, just set to be 1. Makes no difference to beta values, only nrpms are affected.
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
