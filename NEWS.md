---
editor_options: 
  markdown: 
    wrap: sentence
---
# dev

### Added
* Examples added for some functions, and documentation for some related functions merged.
* Added more tests following inspection of `covr` coverage.
* A new example qseaSet for a small portion of the mouse genome. 
* Added the ability to specifying a biomaRt object directly on a qseaSet for used for finding gene annotation in `plotGeneHeatmap`. 

### Changes
* Added functions `setMesaGenome`, `setMesaTxDb` and `setMesaAnnoDb` to set global defaults for the annotation packages required by `annotateWindows`. This has the effect that `annotateWindows` will no longer assume hg38 by default with no arguments.
* `summariseDMRsByGene` function now requires `annotateWindows` to have already been called on the DMRs (previously it called this internally if necessary). 

### BUG FIXES
* `poolSamples` now returns data frames instead of matrices in the libraries slot. 
* `select` now works when dropping a column.
*  PCA/UMAP functions now give a more informative error if there are insufficient regions or samples.
* `downSample` now works correctly again following changes to how `table` and `enframe` interact.

# mesa 0.4.1

### BUG FIXES
* `plotPCA` now works when colouring by a factor variable
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

-   Removed exporting of internal functions: `dropAvgFragDetails`, `getBamCoveragePairedAndUnpairedR1`, `mixThreeQsetSamples`, `getCGPositions`, `fitQseaGLM`, `getDMRsData`, `subsetWindowsOverBackground`

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

-   Six functions have been replaced by different functions (detailed above)

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

-   Various R CMD check warning issues fixed.

-   `filterByOverlaps` no longer allows the option to take a list of window ID numbers.

-   Fixed `addMedipsEnrichmentFactors` function inside `makeQset` when using the `"qseaPaired"` coverageMethod.

-   Fixed tests for `addMedipsEnrichmentFactors`

# mesa 0.2.2

### BUG FIXES

-   Update mixQsetSamples to match the mixArrayWithQset code (and now runs with latest version of R).
-   Less stringent checking on the p values being zero bug in qsea, to minimise false positives.

### KNOWN ISSUES

\* `pull` is not working directly, `pullQset` must be used.

# mesa 0.2.1

Version set for initial release in v2.0.1 of the internal CBC Nextflow pipeline.
