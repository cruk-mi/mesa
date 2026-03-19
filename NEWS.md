---
editor_options: 
  markdown: 
    wrap: sentence
---

# dev

### ADDED
* Many more examples for individual functions on their help pages, as well as a set of vignettes. [!14](https://github.com/cruk-mi/mesa/pull/14)
* Added validity checks for mesa classes (mesaDimRed, mesaPCA, mesaUMAP). [!43](https://github.com/cruk-mi/mesa/pull/43)
* Objects with invalid slots will now throw informative errors. [!43](https://github.com/cruk-mi/mesa/pull/43)
* Added a function `sliceDMRs` to take the 'top' DMRs in each contrast, based on a specified ranking column. [!24](https://github.com/cruk-mi/mesa/pull/24)
* Added `getMesaGenome()`, `getMesaTxDb()`, and `getMesaAnnoDb()` functions to complete mesa's genome management system. Supports hg38, hg19, and mm10 with extensible architecture. [!68](https://github.com/cruk-mi/mesa/pull/68)
* `plotGenomicFeatureDistribution()` now supports multiple genomes via `genome`, `TxDb`, and `annoDb` parameters. Function now works with mesa's genome system via `setMesaGenome()` or custom genome builds. No longer limited to hg38. [!68](https://github.com/cruk-mi/mesa/pull/68)

### CHANGES
* Converted `plotPCA` into a submethod for the `qsea` defined method. [!11](https://github.com/cruk-mi/mesa/pull/11)
* `plotPCA` gains a `verbose` option to turn off most of the messages produced. [!11](https://github.com/cruk-mi/mesa/pull/11)
* `plotPCA` and `plotUMAP` now default to filled shapes rather than empty shapes if less than 6 shapes are required. [!36](https://github.com/cruk-mi/mesa/pull/36)
* `plotPCA` and `plotUMAP` can now specify the shapes used via `shapePalette`, even if no shape annotation is being given. [!36](https://github.com/cruk-mi/mesa/pull/36)
* `getSampleTable` is now defined for PCA/UMAP objects. [!11](https://github.com/cruk-mi/mesa/pull/11)
* `plotGeneHeatmap` now automatically retries if it fails to connect to biomaRt, and fails with a clear error message if it cannot connect. [!44](https://github.com/cruk-mi/mesa/pull/44)
* Updated function documentation with examples, default parameter values, and more detailed descriptions. [!43](https://github.com/cruk-mi/mesa/pull/43)
* `fdrThres` changed to `FDRthres` in `calculateDMRs` and `subsetWindowsOverBackground`. [!48](https://github.com/cruk-mi/mesa/pull/48)
* `makeQset` now checks that the chromosomes provided match with those present in the BSgenome. [!32](https://github.com/cruk-mi/mesa/pull/32)
* The `"PairedAndR1s"` coverage method for `makeQset` will now only process reads that are present in the regions, this should reduce memory requirements and fix issues when `_alt` chromosomes exist in the bam files. This means the fragment size measurements are now only calculated over the selected regions. [!32](https://github.com/cruk-mi/mesa/pull/32)
* `plotPCA` and `plotUMAP` no longer take the `qseaSet` as an input. The PCA/UMAP object has a copy of the `sampleTable` which may be modified instead. [!58](https://github.com/cruk-mi/mesa/pull/58)
* Swapped to use of `seq_along` rather than `1:n` throughout. [!67](https://github.com/cruk-mi/mesa/pull/67)
* Be more specific in the use of message suppression. [!66](https://github.com/cruk-mi/mesa/pull/66)

### REMOVED
* Made `plotGenomicFeatureDistribution` and `getGenomicFeatureDistribution` internal as they currently only work for hg38. [!14](https://github.com/cruk-mi/mesa/pull/14)
* Made `calculateFractionReadsInGRanges` internal as it seems to be returning the fraction of windows that overlap not reads. [!14](https://github.com/cruk-mi/mesa/pull/14)
* Made `countWindowsAboveCutoff` internal as it needs the arguments renaming and better documentation.  [!14](https://github.com/cruk-mi/mesa/pull/14)
* Removed internal functions`getAnnotationDataFrame` and `getAnnotationDataFrameIndividual` as they are superseded by `getAnnotation` and the shift to tidy evaluation via `sampleAnnotation` in the plotting functions. [!14](https://github.com/cruk-mi/mesa/pull/14)
* Removed `colnames` function definion on a qseaSet, which was not working anyway. [!14](https://github.com/cruk-mi/mesa/pull/14)
* Removed `dropAvgFragDetails` as no longer required. [!63](https://github.com/cruk-mi/mesa/pull/63)

### BUG FIXES
* `makeTransposedTable` no longer adds `chr` to the window names even if they already had a `chr` prefix. [!14](https://github.com/cruk-mi/mesa/pull/14)
* Correctly pass the `...` inside `plotGeneHeatmap` and `plotRegionsHeatmap`. [!14](https://github.com/cruk-mi/mesa/pull/14)
* `writeDMRsToBed` should now correctly export the files. [!14](https://github.com/cruk-mi/mesa/pull/14)
* Fixed error when `plotRegionsHeatmap` was given more than one region that overlapped one window. [!21](https://github.com/cruk-mi/mesa/pull/21)
* Correct the message produced by `addMedipsEnrichmentFactors` (thanks @daonslog for reporting). [!14](https://github.com/cruk-mi/mesa/pull/14)
* `makeQset`, `renameSamples` and `renameQsetNames` will no longer accept sample names that are not valid column names in R without quotation. [!31](https://github.com/cruk-mi/mesa/pull/31)
* Correctly pass the `fragmentLength` when calling `makeQset` with the `CNVmethod = "MeCap"` option, and fix an issue with hg19 GRanges. [!21](https://github.com/cruk-mi/mesa/pull/21)
* `plotPCA` now plots a shape column that contains NA values without needing to also specifying `NAshape`. [!36](https://github.com/cruk-mi/mesa/pull/36)
* When no colour or shape annotation is provided, `plotPCA` and `plotUMAP` no longer print a NULLcol or NULLshape column in the legend when using filled shapes. [!36](https://github.com/cruk-mi/mesa/pull/36)
* `plotPCA` and `plotUMAP` now show the fill colour as opposed to black points in the legend when using filled shapes.  [!36](https://github.com/cruk-mi/mesa/pull/36)
* Correctly pass the `fragmentLength` when calling `makeQset` with the `CNVmethod = "MeCap"` option, and fix an issue with hg19 GRanges.   [!21](https://github.com/cruk-mi/mesa/pull/21)
* Prevent exponentially increasing numbers of rows in CNV object when incorrect hmmCopy objects are provided, [fixes issue #26](https://github.com/cruk-mi/mesa/issues/26) reported by @lbeltrame. [!29](https://github.com/cruk-mi/mesa/pull/29)
* Fixed an issue where `makeQset` printed the wrong number of paired reads being filtered out due to having an insert size outside of the selected size range. [!32](https://github.com/cruk-mi/mesa/pull/32)
* Fixed plotting a shape inside `plotPCA` when using ggplot2 4.0.0. [!55](https://github.com/cruk-mi/mesa/pull/55)
* Fixed GRanges conversion error in `plotGenomicFeatureDistribution` that occurred with Bioconductor 3.21 when multiple chromosome columns existed after ChIPseeker annotation. [!68](https://github.com/cruk-mi/mesa/pull/68)
* Fixed makeQset validation tests to handle updated annotation database versions and focus on parameter validation rather than computational integration tests. [!68](https://github.com/cruk-mi/mesa/pull/68)
* Updated network error patterns in testPlotGeneHeatmap, preventing biomart HTP 503 error in test-makeQset.R:83:3. [!68](https://github.com/cruk-mi/mesa/pull/68)

# mesa 0.5.1

This is the first released version on github, following a lengthy period of internal development. Many things have changed in a major overhaul of the package.

### BUG FIXES
* `summariseDMRsByGene` now correctly summarises genes when they have different positions or lengths in different windows.  [!5](https://github.com/cruk-mi/mesa/pull/5)
* `plotCorrelationMatrix` no longer crashes when no annotation is being used.  [!5](https://github.com/cruk-mi/mesa/pull/5)

# mesa 0.5.0 [!6](https://github.com/cruk-mi/mesa/pull/6)

### ADDED
* Examples added for some functions, and documentation for some related functions merged.
* Added more tests following inspection of `covr` coverage.
* `getPCA()` and `getUMAP()` now return an object of class `mesaDimRed`, which wraps the previous list based output. This object now contains a copy of the sampleTable, which may be edited using `mutate` and `left_join`.
* A new example qseaSet for a small portion of the mouse genome.
* Added the ability to specifying a biomaRt object directly on a qseaSet for used for finding gene annotation in `plotGeneHeatmap`.

### CHANGES
* `plotPCA()` and `plotUMAP()` no longer require passing the qseaSet, as the sampleTable is stored in the object.
* Added functions `setMesaGenome`, `setMesaTxDb` and `setMesaAnnoDb` to set global defaults for the annotation packages required by `annotateWindows`. This has the effect that `annotateWindows` will no longer assume hg38 by default with no arguments.
* `summariseDMRsByGene` function now requires `annotateWindows` to have already been called on the DMRs (previously it called this internally if necessary).
* Reduced the number of messages produced when generating tables of data.

### BUG FIXES
* `poolSamples` now returns data frames instead of matrices in the libraries slot.
* `select` now works when dropping a column. [!6](https://github.com/cruk-mi/mesa/pull/6)
*  PCA/UMAP functions now give a more informative error if there are insufficient regions or samples.
* `downSample` now works correctly again following changes to how `table` and `enframe` interact.
* `calculateDMRs` now correctly returns an empty data frame if no DMRs are found.

# mesa 0.4.1

### BUG FIXES
* `plotPCA` now works when colouring by a factor variable.
* `mutate` now cannot be used to change the sample_name in the sample table of a qseaSet which breaks the object. `renameQsetNames` or `renameSamples` must be used for this.

### CHANGES
* HMMCopy related functions now take explicit input of the GC and mappability tracks, rather than having hardcoded internal hg38 objects, allowing for different genomes.
* Removed the data objects `gc_hg38_10kb` and `map_hg38_10kb` (GC and mappability tracks for hg38 at 10kb resolution) for package size reasons.

# mesa 0.4.0

### BUG FIXES

-   Fixed how the relative enrichment (relH/GoGe) calculation is performed (`addMedipsEnrichmentFactors`) when using the qsea default method for reads (#1).

-   Filtering a `qseaSet` now works correctly when a `qseaSet` has no `cnv` data.

-   Fixed `pull` to work directly on a `qseaSet` to return a column of the sampleTable.

### CHANGES

-   Removed `mutateQset`, `filterQset`, `arrangeQset`, `sortQset`, `pullQset`, `leftJoinQset` as these just work on the `dplyr` verbs directly.

-   Combined `getBetaMeans` and `getNormalisedReadSum` into new function called `summariseAcrossWindows`.
    This allows an arbitrary function to be (e.g. mean or sd) to be applied over all the windows in each sample, and returns a data frame with one result per sample.

-   Added a new function `addSummaryAcrossWindows` that performs `summariseAcrossWindows` but returns the `qseaSet` with the results added to the `sampleTable`.

-   Added a function `asValidGranges` that checks if an object can be coerced to a GRanges object and does so if possible.

-   Added `setMesaGenome` function to globally set which genome to be used for plotting functions.

-   Standardised on using "Windows" instead of "Regions".
    For instance, `filterRegions` is now `filterWindows`, and `plotGRangesHeatmap` is now `plotRegionsHeatmap`.
    Added a function `getWindows` to mimic `getRegions` from qsea.

-   Complete overhaul of the PCA functionality.

    -   This is now split into two functions, `getPCA` and `plotPCA`, and the `getPCAwithBatch` function has been removed.

-   UMAP functionality added, using the same framework as PCAs, with `getUMAP` and `getPCA`.
    This is still experimental.

-   Changed `pivotDMRsLonger` to return columns "`group1`" and "`group2`" instead of "`sample1`" and "`sample2`".

-   In DMR results, moved the `betaDelta` columns to be after the equivalent `adjPval` columns instead of all being at the end.

-   Renamed the `contrastsToDo` argument of `makeDMRs` to just `contrasts`.

-   `mixSamples` (formerly `mixQsetSamples`) no longer sets `type` and `tumour` in the `sampleTable` of the returned `qseaSet`.

-   Changes to heatmap plotting functions (`plotRegionsHeatmap` and `plotGeneHeatmap`, `plotCNVHeatmap`):

    -   `plotGRangesHeatmap` is now renamed as `plotRegionsHeatmap`

    -   Both support a `useGroups` argument, to determine whether to average over the "group" column of the sampleTable or not.

    -   Argument `signatureGR` of `plotRegionsHeatmap` changed to `regionsToOverlap`.

    -   Changed from `pheatmap` package to `ComplexHeatmap` package for the plotting.

    -   The method of adding sample annotations has been overhauled.
        This is now specified using `sampleAnnotation` (instead of `annotationCol`), and does not support the use of the `getAnnotationDataFrame` function.
        These options now directly take the name of the features to be plotted, supporting unquoted tidy evaluation, i.e. all of these work:

        -   `sampleAnnotation = group`

        -   `sampleAnnotation = c(group, type)`

        -   `sampleAnnotation = c("group","type")`

    -   Heatmaps can now have windows annotated with the use of `windowAnnotation`.
        This supports unquoted tidy evaluation as well, and can use anything on the qseaSet regions slot or on the `regionsToOverlap` argument.

    -   The `getAnnotationDataFrame` and `getAnnotationDataFrameIndividual` functions have been removed.

    -   The colours for the annotations have been overhauled.
        These now try to choose appropriate colour scales for continuous data based on whether they are strictly positive or not.
        For discrete annotations, the colours are set globally for the entire annotation set, to try and prevent the occurrance of different annotations using very similar colours.

-   `getDataTable` function can now keep the suffix (e.g. `_nrpm` or `_beta`) if required.

-   `removeWindowsOverCutoff` and `keepWindowsOverCutoff` have been replaced by a more general function called `subsetWindowsBySignal`.
    This allows for a general function to be specified (for instance mean or max) as well as a threshold, enabling for the filtering of windows based on that function.

-   Parallelisation of functions is now controlled explicitly by the `setMesaParallel` function.
    This must be used to enable multicore evaluation.

-   When calculating DMRs, by default the `qseaSet` will now be filtered to only include samples which are present in the contrasts for the calculating of the dispersion estimates in the generalised linear model.
    This can be controlled by the `calcDispersionAll` option.
    The pre-filtering of windows based on expression is unchanged - this only uses samples included in the contrasts.

-   The difference in mean beta values in the DMR output has been renamed from `betaDelta` to `deltaBeta` to reflect the mathematical usage \tex{\Delta \beta}.

-   Updated documentation in many functions.

-   Authors updated to include Paddy Harker, Kevin Brennan and Katarzyna Kamieniecka.

-   Removed exporting of internal functions: `dropAvgFragDetails`, `getBamCoveragePairedAndUnpairedR1`, `mixThreeQsetSamples`, `getCGPositions`, `fitQseaGLM`, `getDMRsData`, `subsetWindowsOverBackground`.

-   Renamed many functions throughout the package:

    | Previous Function Name          | New Function Name              |
    |---------------------------------|--------------------------------|
    | `addQseaNormalisationSteps`     | `addNormalisation`             |
    | `annotateData`                  | `annotateWindows`              |
    | `calculateArrayBetas`           | `convertToArrayBetaTable`      |
    | `calculateCpGEnrichment`        | `calculateCGEnrichment`        |
    | `calculateCpGEnrichmentGRanges` | `calculateCGEnrichmentGRanges` |
    | `countWindowsOverCutoff`        | `countWindowsAboveCutoff`      |
    | `downsampleQsea`                | `downSample`                   |
    | `filterRegions`                 | `filterWindows`                |
    | `liftOverhg19`                  | `liftOverHg19`                 |
    | `mixQsetSamples`                | `mixSamples`                   |
    | `plotDMRUpSet`                  | `plotDMRUpset`                 |
    | `qseaWriteDMRsBeds`             | `writeDMRsToBed`               |
    | `qseaWriteDMRsExcel`            | `writeDMRsToExcel`             |
    | `relabelQset`                   | `renameSamples`                |
    | `summariseByGene`               | `summariseDMRsByGene`          |

-   Six functions have been replaced by different functions (detailed above).

    +---------------------------+--------------------------+
    | Previous Function         | Replacement Function(s)  |
    +===========================+==========================+
    | `getQseaPCA`              | `getPCA`/`plotPCA`       |
    |                           |                          |
    | `getPCAwithBatch`         |                          |
    +---------------------------+--------------------------+
    | `getBetaMeans`            | `summariseAcrossWindows` |
    |                           |                          |
    | `getNormalisedReadSum`    |                          |
    +---------------------------+--------------------------+
    | `removeWindowsOverCutoff` | `subsetWindowsBySignal`  |
    |                           |                          |
    | `keepWindowsOverCutoff`   |                          |
    +---------------------------+--------------------------+

-   Renamed data objects:

    | Previous data name        | New data name         |
    |---------------------------|-----------------------|
    | `encodeBlacklist`         | `ENCODEbadRegions`    |
    | `examplePairedTumourQset` | `exampleTumourNormal` |

# mesa 0.2.3

### CHANGES

-   `plotGeneHeatmap` now uses `bioMaRt` to find the gene details, based on the genome or mart supplied.

-   `filterByOverlaps` now works when the \`qseaSet\` regions contains "chr" but the regions to filter by do not or vice versa.

-   Removed all references to `PooledControl` sample (previously suggested for normalisation).

-   Various R CMD check warning issues fixed.[!6](https://github.com/cruk-mi/mesa/pull/6)

-   `filterByOverlaps` no longer allows the option to take a list of window ID numbers.

-   Fixed `addMedipsEnrichmentFactors` function inside `makeQset` when using the `"qseaPaired"` coverageMethod.

-   Fixed tests for `addMedipsEnrichmentFactors`.

# mesa 0.2.2

### BUG FIXES

-   Update mixQsetSamples to match the mixArrayWithQset code (and now runs with latest version of R).
-   Less stringent checking on the p values being zero bug in qsea, to minimise false positives.

### KNOWN ISSUES

\* `pull` is not working directly, `pullQset` must be used.

# mesa 0.2.1

Version set for initial release in v2.0.1 of the internal CBC Nextflow pipeline.
